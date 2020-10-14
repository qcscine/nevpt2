!
!  Module to handle reading of higher-order RDMs from (QCMaquis) HDF5
!

module rdm_utils

    use hdf5_utils
    use swap_util, only: swap
    use indices_norm, only: norm
#ifdef DMRG_NEVPT
    use qcmaquis_interface
    use qcmaquis_interface_utility_routines, only: str
#endif
    implicit none

    public hdf5_read_rdm
    public hdf5_read_rdm_from_resfile
    private hdf5_qcmaquis_result_from_checkpointfile
contains

    !! Get a HDF5 QCMaquis result file name from a checkpoint file
    function hdf5_qcmaquis_result_from_checkpointfile(chkpname) result (resultname)
      character(len=*),intent(in)   :: chkpname
      character(len=255) :: resname_tmp
      character(len=:),allocatable  :: resultname

      character(len=:),allocatable :: chkp_props_name
      integer(HID_T) :: file_id
      logical file_exists

      resname_tmp = ''

      chkp_props_name = trim(chkpname)//"/props.h5"
      inquire(file=chkp_props_name, exist=file_exists)
      if (.not.file_exists) then
        write(*,*) "Cannot find file "//chkp_props_name
        stop
      end if
      call hdf5_open(chkp_props_name, file_id)
      ! Get the result filename and close the checkpoint name
      ! Use a temporary fixed-length string resname_tmp because the HDF5 call does not support
      ! variable-length strings
      call hdf5_get_data_str(file_id, "/parameters/resultfile", [1], resname_tmp) ! error handling?
      call hdf5_close(file_id)
      resultname = trim(resname_tmp)

    end function hdf5_qcmaquis_result_from_checkpointfile

    !! Read 3- and 4-RDM by specifying the corresponding QCMaquis checkpoint file
    !! The subroutine finds out the corresponding result file and reads the RDMs from there
    subroutine hdf5_read_rdm(chkpname, rdm3, rdm4, read_num_elements, trans_rdm)
      character(len=*), intent(in)               :: chkpname
      real*8, dimension(:,:,:,:,:,:), intent(inout), optional :: rdm3
      real*8, dimension(:,:,:,:,:),   intent(inout), optional :: rdm4 ! collapsed dimensions for 4-RDM
      logical, intent(in), optional :: trans_rdm

      integer, intent(out),optional :: read_num_elements ! number of read elements

      call hdf5_read_rdm_from_resfile( &
        hdf5_qcmaquis_result_from_checkpointfile(trim(chkpname)), rdm3, rdm4, read_num_elements, trans_rdm)
    end subroutine hdf5_read_rdm

    !! Routine to read (parts) of 3- and 4-RDM from QCMaquis HDF5 files
    !! This routine may be called multiple times with different HDF5 files
    !! to fill the rdm3 or rdm4 arrays only partially
    !! Optional parameter: read_num_elements: returns read
    !! trans_rdm: (optional, default=.false.) whether we read a transition 3-RDM. No effect for 4-RDM
    !! resultname: QCMaquis result file
    subroutine hdf5_read_rdm_from_resfile(result_name, rdm3, rdm4, read_num_elements, trans_rdm)
      character(len=*), intent(in)               :: result_name
      real*8, dimension(:,:,:,:,:,:), intent(inout), optional :: rdm3
      real*8, dimension(:,:,:,:,:),   intent(inout), optional :: rdm4 ! collapsed dimensions for 4-RDM


      logical, intent(in), optional :: trans_rdm
      integer, intent(out),optional :: read_num_elements ! number of read elements
      real*8, dimension(:), allocatable :: buffer ! buffer for values
      integer*4, dimension(:,:), allocatable :: idx_buffer ! buffer for indices

      integer(HID_T) :: file_id, dset_id, space_id
      integer(HERR_T) :: error
!       integer*8 :: L ! number of orbitals
      integer :: i ! loop variable
      integer(HSIZE_T) :: rank, nentries
      integer(HSIZE_T), dimension(1:2) :: ds_size
      character(len=:), allocatable :: ds_name_root, ds_name ! dataset name, variable-length string and QCMaquis result file

      character(len=:), allocatable :: chkp_props_name !name for checkpoint/props.h5 file
      logical :: file_exists

      logical :: trans_rdm_

      trans_rdm_ = .false.
      if (present(trans_rdm)) trans_rdm_ = trans_rdm

      if (present(rdm3).eqv.present(rdm4)) then
        write(*,*) "Error in hdf5_read_rdm: Only 3-RDM or 4-RDM variable may be specified at the same time."
        stop
      endif

      inquire(file=trim(result_name), exist=file_exists)
      if (.not.file_exists) then
        write(*,*) "Cannot find QCMaquis result file "//trim(result_name)//"."
        stop
      end if
      ! now open the result file
      call hdf5_open(trim(result_name), file_id)
!       ! get L
!       call hdf5_get_data(file_id, "/parameters/L",[1], L)
!       ! error handling?

      ! get the name of the labels dataset
      if (present(rdm4)) then
        ds_name_root = "/spectrum/results/fourptdm"
      else
        if (present(rdm3)) then
            if (trans_rdm_) then
              ds_name_root = "/spectrum/results/transition_threeptdm"
            else
              ds_name_root = "/spectrum/results/threeptdm"
            end if
        end if
      end if

      ds_name = trim(ds_name_root)//'/labels_num'
     ! If the 4-RDM is going to be too large, one should consider implement chunked reading, but
     ! for now we won't do it

      ! read values into buffer
      call h5dopen(file_id, ds_name, dset_id, error)

      if (error.ne.0) then
        write (*,*) "Could not open the RDM dataset. Have you calculated the RDMs with QCMaquis?"
        stop
      end if
      call h5dget_space(dset_id, space_id, error)

      ! obtain the dimensions of labels
      call h5sget_simple_extent_dims(space_id, ds_size, error)
      call h5sclose(space_id, error)
      rank = ds_size(1)
      nentries = ds_size(2)

      if (present(read_num_elements)) then
        read_num_elements = nentries
      end if
      if (present(rdm4)) then
        if (rank.ne.8) then
          write(*,*) "4-RDM should have 8 indices per entry stored in HDF5."
          stop
        end if
      end if
      if (present(rdm3)) then
        if (rank.ne.6) then
          write(*,*) "3-RDM should have 6 indices per entry stored in HDF5."
          stop
        end if
      end if
      ! read indices
      ! warning, indices are saved as integer*4
      allocate(idx_buffer(rank,nentries))
      call hdf5_get_data(file_id, ds_name, ds_size, idx_buffer)

      ! now read the values
      ! get the name of the value dataset

      ds_name = trim(ds_name_root)//"/mean/value"

      ! read values into buffer
      call h5dopen(file_id, ds_name, dset_id, error)

      if (error.ne.0) then
        write (*,*) "Could not open the RDM dataset. Have you calculated the RDMs with QCMaquis?"
        stop
      end if
      call h5dget_space(dset_id, space_id, error)

      ! obtain the dimensions of value dataset and make sure they match with dimensions of the labels
      call h5sget_simple_extent_dims(space_id, ds_size, error)
      call h5sclose(space_id, error)
      if (ds_size(1).ne.nentries) then
        write (*,*) "Dataset size mismatch between RDM labels and values. Your HDF5 file might be corrupt."
        stop
      end if

      allocate(buffer(nentries))
      call hdf5_get_data(file_id, ds_name, ds_size, buffer)
      call h5dclose(dset_id, error)
      call hdf5_close(file_id)

! calculate the size from L, might be interesting for a check but otherwise unneeded if we assume rdm arrays are
! properly allocated
      ! increase indices by 1 to make them Fortran-conformant

      idx_buffer = idx_buffer + 1
      ! fill RDM arrays
      do i=1, nentries
        if(present(rdm4)) then
          call permute_4rdm(rdm4, &
                            int(idx_buffer(1,i),8),& ! explicit conversion of integer*4 to integer*8 of indices,
                            int(idx_buffer(2,i),8),& ! otherwise the compiler will complain
                            int(idx_buffer(3,i),8),&
                            int(idx_buffer(4,i),8),&
                            int(idx_buffer(5,i),8),&
                            int(idx_buffer(6,i),8),&
                            int(idx_buffer(7,i),8),&
                            int(idx_buffer(8,i),8),&
                            buffer(i))
        else
         if(present(rdm3)) then
           call permute_3rdm(rdm3, &
                            int(idx_buffer(1,i),8),&
                            int(idx_buffer(2,i),8),&
                            int(idx_buffer(3,i),8),&
                            int(idx_buffer(4,i),8),&
                            int(idx_buffer(5,i),8),&
                            int(idx_buffer(6,i),8),&
                            buffer(i), &
                            trans_rdm_)
         end if
        end if
      end do

      ! clean up
      if (allocated(buffer)) deallocate(buffer)
      if (allocated(idx_buffer)) deallocate(idx_buffer)
    end subroutine hdf5_read_rdm_from_resfile

    ! Fills 4-RDM array with all symmetry permutations of i,j,k,l,m,n,o,p
    ! Essentially copy-pasted from QCMaquis fourptdm_perm.py by S. Knecht
    ! rdm4: 4-RDM array with first four indices collapsed
    ! i - p : indices (starting from 1)
    ! val: value
    ! TODO: this seems highly redundant, so it should be simplified!
#define dump_element(M, v, i, j, k, l, m, n, o, p) M(norm(i,j,k,l), m,n,o,p) = v
    subroutine permute_4rdm(rdm4,i,j,k,l,m,n,o,p,val)
      real*8, dimension(:,:,:,:,:), intent(inout) :: rdm4 ! first index is collapsed with norm(i,j,k,l)
      integer, intent(in) :: i,j,k,l,m,n,o,p
      real*8, intent(in) :: val
      real*8 :: valm05, valm2
      real*8, parameter :: rdm_thresh = 1.0d-12
      if (abs(val).lt.rdm_thresh) return

      valm05 = val*(-0.5d0)
      valm2  = val*(-2.0d0)


      if ((i.eq.j).and.(k.eq.l)) then ! case 1: i.eq.j k.eq.l
        if ((m.eq.n).and.(o.eq.p)) then
          !default
          dump_element(rdm4,val,i,j,k,l,m,n,o,o)
          dump_element(rdm4,val,i,j,k,l,o,o,m,m)
          dump_element(rdm4,valm05,i,j,k,l,m,o,o,m)
          dump_element(rdm4,valm05,i,j,k,l,m,o,m,o)
          dump_element(rdm4,valm05,i,j,k,l,o,m,m,o)
          dump_element(rdm4,valm05,i,j,k,l,o,m,o,m)
        else
          if ((m.eq.n).and.(o.ne.p)) then
              dump_element(rdm4,val,i,j,k,l,m,n,o,p)
              ! regular permutations
              dump_element(rdm4,val,i,j,k,l,m,m,p,o)
              dump_element(rdm4,val,i,j,k,l,o,p,m,m)
              dump_element(rdm4,val,i,j,k,l,p,o,m,m)
              ! remaining permutations
              dump_element(rdm4,valm05,i,j,k,l,o,m,m,p)
              dump_element(rdm4,valm05,i,j,k,l,o,m,p,m)
              dump_element(rdm4,valm05,i,j,k,l,p,m,m,o)
              dump_element(rdm4,valm05,i,j,k,l,p,m,o,m)
              dump_element(rdm4,valm05,i,j,k,l,m,o,p,m)
              dump_element(rdm4,valm05,i,j,k,l,m,p,o,m)
              dump_element(rdm4,valm05,i,j,k,l,m,p,m,o)
              dump_element(rdm4,valm05,i,j,k,l,m,o,m,p)
          else
            if ((m.ne.n).and.(o.eq.p)) then
                !default
                dump_element(rdm4,val,i,j,k,l,m,n,o,p)
                ! regular permutations
                dump_element(rdm4,val,i,j,k,l,n,m,o,o)
                dump_element(rdm4,val,i,j,k,l,o,o,m,n)
                dump_element(rdm4,val,i,j,k,l,o,o,n,m)
                ! remaining permutations
                dump_element(rdm4,valm05,i,j,k,l,m,o,o,n)
                dump_element(rdm4,valm05,i,j,k,l,n,o,o,m)
                dump_element(rdm4,valm05,i,j,k,l,o,m,o,n)
                dump_element(rdm4,valm05,i,j,k,l,o,n,o,m)
                dump_element(rdm4,valm05,i,j,k,l,o,m,n,o)
                dump_element(rdm4,valm05,i,j,k,l,o,n,m,o)
                dump_element(rdm4,valm05,i,j,k,l,m,o,n,o)
                dump_element(rdm4,valm05,i,j,k,l,n,o,m,o)
            else
              if ((m.ne.n).and.(o.ne.p)) then
                if (n.eq.o) then
                !default
                dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                dump_element(rdm4,valm2,i,j,k,l,n,n,m,p)
                dump_element(rdm4,valm2,i,j,k,l,n,n,p,m)
                dump_element(rdm4,valm2,i,j,k,l,m,p,n,n)
                dump_element(rdm4,valm2,i,j,k,l,p,m,n,n)

                dump_element(rdm4,val,i,j,k,l,p,n,n,m)
                dump_element(rdm4,val,i,j,k,l,p,n,m,n)
                dump_element(rdm4,val,i,j,k,l,m,n,p,n)
                dump_element(rdm4,val,i,j,k,l,n,p,m,n)
                dump_element(rdm4,val,i,j,k,l,n,m,p,n)
                dump_element(rdm4,val,i,j,k,l,n,m,n,p)
                dump_element(rdm4,val,i,j,k,l,n,p,n,m)
                else
                  if ((n.lt.o).and.(p.ne.m)) then
                  !default
                    dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                    dump_element(rdm4,val,i,j,k,l,m,n,p,o)
                    dump_element(rdm4,val,i,j,k,l,n,m,p,o)
                    dump_element(rdm4,val,i,j,k,l,n,m,o,p)
                    dump_element(rdm4,val,i,j,k,l,o,p,m,n)
                    dump_element(rdm4,val,i,j,k,l,o,p,n,m)
                    dump_element(rdm4,val,i,j,k,l,p,o,n,m)
                    dump_element(rdm4,val,i,j,k,l,p,o,m,n)
                  else
                    if (o.ne.p) then
                      dump_element(rdm4,val,i,j,k,l,m,n,o,p)
                      dump_element(rdm4,val,i,j,k,l,m,n,p,o)
                      dump_element(rdm4,val,i,j,k,l,n,m,p,o)
                      dump_element(rdm4,val,i,j,k,l,n,m,o,p)
                      dump_element(rdm4,val,i,j,k,l,o,p,m,n)
                      dump_element(rdm4,val,i,j,k,l,o,p,n,m)
                      dump_element(rdm4,val,i,j,k,l,p,o,n,m)
                      dump_element(rdm4,val,i,j,k,l,p,o,m,n)
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      else
        if ((i.eq.j).and.(k.ne.l)) then ! case 2: i.eq.j k.ne.l
          if ((m.eq.n).and.(o.eq.p)) then ! 2x2 equal
            !default
            dump_element(rdm4,val,i,j,k,l,m,n,o,p)

            dump_element(rdm4,val,i,j,k,l,o,o,m,m)

            dump_element(rdm4,valm05,i,j,k,l,m,o,o,m)
            dump_element(rdm4,valm05,i,j,k,l,o,m,m,o)
            dump_element(rdm4,valm05,i,j,k,l,o,m,o,m)
            dump_element(rdm4,valm05,i,j,k,l,m,o,m,o)
          endif
          if ((m.eq.n).and.(o.ne.p)) then ! 1 equal (first 2)
            !default
            dump_element(rdm4,val,i,j,k,l,m,n,o,p)

            dump_element(rdm4,valm05,i,j,k,l,m,o,m,p)
            dump_element(rdm4,valm05,i,j,k,l,m,p,o,m)
            dump_element(rdm4,valm05,i,j,k,l,o,m,m,p)
            dump_element(rdm4,valm05,i,j,k,l,p,m,o,m)
          endif
          if ((m.ne.n).and.(o.eq.p)) then ! 1 equal (latter 2)
            !default
            dump_element(rdm4,val,i,j,k,l,m,n,o,p)

            dump_element(rdm4,val,i,j,k,l,n,m,o,o)
          endif
          if ((m.ne.n).and.(n.ne.o).and.(o.ne.p).and.(m.ne.p)) then ! all different
            !default
            dump_element(rdm4,val,i,j,k,l,m,n,o,p)

            dump_element(rdm4,val,i,j,k,l,n,m,o,p)
          endif
        else
          if ((i.ne.j).and.(j.eq.k).and.(k.ne.l)) then ! case 3: i.ne.j j.eq.k k.ne.l
            if ((m.eq.n).and.(o.eq.p)) then ! 2x2 equal
                !default
                dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                dump_element(rdm4,val,i,j,k,l,o,o,m,m)
                dump_element(rdm4,val,i,j,k,l,o,m,o,m)
                dump_element(rdm4,val,i,j,k,l,m,o,m,o)

                dump_element(rdm4,valm2,i,j,k,l,m,o,o,m)
                dump_element(rdm4,valm2,i,j,k,l,o,m,m,o)
            endif
            if ((m.ne.n).and.(o.eq.p)) then ! 1 equal (latter 2)
                !default
                dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                dump_element(rdm4,val,i,j,k,l,o,o,m,n)
                dump_element(rdm4,val,i,j,k,l,o,m,o,n)
                dump_element(rdm4,val,i,j,k,l,m,o,n,o)

                dump_element(rdm4,valm2,i,j,k,l,m,o,o,n)
            endif
            if ((m.ne.n).and.(o.ne.p).and.(m.eq.p)) then ! 1 equal (first.and.last)
                !default
                dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                dump_element(rdm4,val,i,j,k,l,m,o,n,m)
            endif
            if ((m.ne.n).and.(n.ne.o).and.(o.ne.p).and.(m.ne.p)) then ! all different
                !default
                dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                dump_element(rdm4,val,i,j,k,l,m,o,n,p)
            endif
          else
            if ((i.ne.j).and.(k.eq.l)) then ! case 4: i.ne.j k.eq.l
              if ((m.eq.n).and.(o.eq.p)) then ! 2x2 equal
                  !default
                  dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                  dump_element(rdm4,val,i,j,k,l,o,o,m,m)

                  dump_element(rdm4,valm05,i,j,k,l,m,o,o,m)
                  dump_element(rdm4,valm05,i,j,k,l,o,m,m,o)
                  dump_element(rdm4,valm05,i,j,k,l,o,m,o,m)
                  dump_element(rdm4,valm05,i,j,k,l,m,o,m,o)
              endif
              if ((m.eq.n).and.(o.ne.p)) then ! 1 equal (first 2)
                  !default
                  dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                  dump_element(rdm4,val,i,j,k,l,m,m,p,o)
              endif
              if ((m.ne.n).and.(o.eq.p)) then ! 1 equal (latter 2)
                  !default
                  dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                  dump_element(rdm4,valm05,i,j,k,l,m,o,n,o)
                  dump_element(rdm4,valm05,i,j,k,l,m,o,o,n)
                  dump_element(rdm4,valm05,i,j,k,l,o,n,o,m)
                  dump_element(rdm4,valm05,i,j,k,l,o,n,m,o)
              endif
              if ((m.ne.n).and.(n.ne.o).and.(o.ne.p).and.(m.ne.p)) then ! all different
                  !default
                  dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                  dump_element(rdm4,val,i,j,k,l,m,n,p,o)
              endif
            else
              if ((i.ne.j).and.(k.ne.l).and.(j.ne.k)) then ! case 5: i.ne.j j.ne.k k.ne.l
                if (((m.eq.o).and.(n.eq.p)).or.((m.eq.p.and.n.eq.o))) then
                    !default
                    dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                    dump_element(rdm4,val,i,j,k,l,n,m,p,o)
                else
                  if ((m.eq.n).and.(o.eq.p)) then
                      !default
                      dump_element(rdm4,val,i,j,k,l,m,n,o,p)

                      dump_element(rdm4,val,i,j,k,l,o,o,m,m)
                  else
                      !default
                      dump_element(rdm4,val,i,j,k,l,m,n,o,p)
                  endif
                endif
              endif
            endif
          endif
        endif
      endif

    end subroutine permute_4rdm
#undef dump_element

    ! same for 3-RDMs, copy-paste from threeptdm.py
    ! trans_rdm: true if we have a transition density matrix, otherwise false
#define dump_element(M, v, i, j, k, l, m, n) M(i,j,k,l,m,n) = v
    subroutine permute_3rdm(rdm3,i,j,k,l,m,n, val, trans_rdm)
      real*8, dimension(:,:,:,:,:,:), intent(inout) :: rdm3
      integer, intent(in) :: i,j,k,l,m,n
      real*8, intent(in) :: val
      logical, intent(in), optional :: trans_rdm
      logical :: trans_rdm_

      real*8, parameter :: rdm_thresh = 1.0d-12
      real*8 :: tmp

      trans_rdm_ = .false.
      if (present(trans_rdm)) trans_rdm_ = trans_rdm

      ! invert the sign of the value if we don't have transition density matrix
      ! for whatever reason, but it was like this in the old QCMaquis interface
      ! (see threeptdm.py and transition_threeptdm.py)
      tmp = merge(val, -val, trans_rdm_)

      if (abs(val).lt.rdm_thresh) return
      if ((i == j).and.(i == k)) return
      if ((l == m).and.(l == n)) return

      ! 6 permutations for 3-RDM

      dump_element(rdm3,tmp,i,j,k,l,m,n)
      dump_element(rdm3,tmp,i,k,j,l,n,m)
      dump_element(rdm3,tmp,j,i,k,m,l,n)
      dump_element(rdm3,tmp,j,k,i,m,n,l)
      dump_element(rdm3,tmp,k,i,j,n,l,m)
      dump_element(rdm3,tmp,k,j,i,n,m,l)

      if(.not.trans_rdm_) then
        ! transpose elements
        dump_element(rdm3,tmp,l,m,n,i,j,k)
        dump_element(rdm3,tmp,l,n,m,i,k,j)
        dump_element(rdm3,tmp,m,l,n,j,i,k)
        dump_element(rdm3,tmp,m,n,l,j,k,i)
        dump_element(rdm3,tmp,n,l,m,k,i,j)
        dump_element(rdm3,tmp,n,m,l,k,j,i)
      end if

    end subroutine permute_3rdm

#undef dump_element

    ! Reads distributed n-RDM from a given path
    ! Traverses all subdirectories and searches for QCMaquis HDF5 file
    ! result_name, attempting to read n-RDM
    ! Specify one of
    ! rdm4: 4-RDM
    ! rdm3: (transition)-3RDM
    ! trans_rdm: optional, set to .true. for transition 3-RDMs
    subroutine read_distributed_rdm(rdm3, rdm4, path, result_name, trans_rdm)

      real*8, dimension(:,:,:,:,:,:), intent(inout), optional :: rdm3
      real*8, dimension(:,:,:,:,:),   intent(inout), optional :: rdm4 ! collapsed dimensions for 4-RDM
      character(*), intent(in) :: path, result_name
      character(len=:), allocatable :: tmp_string
      character(len=128), dimension(:), allocatable :: dir_list
      character(*), parameter :: tmpfile = '.fortran_dir'
      integer, parameter :: lu = 6366
      integer :: stat, num_dirs, acc_elements, num_elements, total_elements, i
      logical :: file_exists

      logical, intent(in), optional :: trans_rdm
      logical :: trans_rdm_


      trans_rdm_ = .false.
      if (present(trans_rdm)) trans_rdm_ = trans_rdm

      ! Only one out of 3- or 4-RDM may be specified as a parameter
      if (present(rdm3).eqv.present(rdm4)) then
        write(*,*) "Error in read_distributed_4rdm: Only 3-RDM or 4-RDM variable may be specified at the same time."
        stop
      endif

      ! check if path exists
      inquire(file=trim(path), exist=file_exists)
      if (.not.file_exists) then
        write(*,*) "Cannot find path "//trim(path)//"."
        stop
      end if

      ! the only way to list all directories in fortran seems
      ! to be the system call of ls -1 and reading its output
      call system("ls -1 "//trim(path)//" > "//tmpfile)
      open(lu, file=tmpfile, action="read", status="old", iostat=stat)
      if (stat.ne.0) stop "filesystem error"

      ! count the number of lines
      num_dirs = 0
      do
        read(lu,FMT='(a)',iostat=stat) tmp_string
        if (stat.ne.0) EXIT
        num_dirs = num_dirs+1
      end do
      rewind(lu)

      ! allocate the correct number of entries
      allocate(dir_list(num_dirs))

      do i=1, num_dirs
        read(lu,FMT='(a)',iostat=stat) dir_list(i)
        if (stat.ne.0) stop "filesystem error"
      end do
      close(lu, status="delete")
      write(*,*) "Reading distributed "//merge("4","3",present(rdm4))//"-RDM from directory "//trim(path)

      total_elements = merge(qcmaquis_interface_get_4rdm_elements(), &
                             qcmaquis_interface_get_3rdm_elements(bra_neq_ket=trans_rdm_), &
                             present(rdm4))
      write (*,*) "Total number of RDM elements: "//trim(str(total_elements))

      acc_elements = 0 ! counter to keep track how many 4-RDM elements we have read
      ! traverse all subdirectories and read the HDF5 results from each subdirectory
      do i=1, num_dirs
        tmp_string = trim(path)//"/"//trim(dir_list(i))//"/"//trim(result_name)
        inquire(file=tmp_string, exist=file_exists)
        ! check if file exists, otherwise do nothing
        if (file_exists) then
          call hdf5_read_rdm_from_resfile(result_name=tmp_string, &
              rdm3 = rdm3, rdm4 = rdm4, read_num_elements=num_elements, trans_rdm=trans_rdm)
          ! accumulate number of elements
          acc_elements = acc_elements + num_elements
          write (*,*) "Read file "//tmp_string//" ... "// &
            trim(str(acc_elements))//'/'//trim(str(total_elements))// &
            " elements"
        end if
        ! skip remaining files once all elements have been read.
        ! TODO: make a better check instead of checking only the # of elements?
        if (acc_elements.eq.total_elements) exit
      end do

      if (acc_elements.ne.total_elements) then
        if (acc_elements.eq.0) then
          write(*,*) "Did not find any QCMaquis HDF5 result files. Did you specify the correct path?"
        else
          write (*,*) "Mismatch in the number of RDM elements. Could read only "//trim(str(acc_elements))// &
              " out of "//trim(str(total_elements))//". Did all calculations succeed?"
        end if
        stop
      end if
      if (allocated(dir_list)) deallocate(dir_list)
    end subroutine

end module rdm_utils

