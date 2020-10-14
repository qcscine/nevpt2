module info_orbital_space

  use hdf5_utils

  implicit none

  public initialize_inforb
  public finalize_inforb
  public set_orbarray_pointers
  ! These routines should be called by the MOLCAS driver
  public initialize_inforb_molcas
  public finalize_inforb_molcas

  type type_inforb
  integer, allocatable   :: nish(:)
  integer, allocatable   :: nash(:)
  integer, allocatable   :: nssh(:)
  integer, allocatable   :: nosh(:)
  integer, allocatable   :: nfro(:)
  integer, allocatable   :: ndel(:)
  integer, allocatable   :: iosh(:)
  integer                :: nbast
  integer                :: nisht
  integer                :: nasht
  integer                :: nssht
  integer*8              :: norbtt_b
  integer*8              :: nsym_b
  integer*8, allocatable :: norb_b(:)
  integer*8, allocatable :: nfro_b(:)
  integer*8, allocatable :: ndel_b(:)
  real*8                 :: ecor_b
  end type type_inforb

  type type_inforb_short
  ! nfro/ndel is the number of
  ! frozen/deleted orbitals read from MOTRA
  ! norb -- number of employed orbitals
  ! in MOTRA -> check if it's the same as nbas?
  integer, allocatable   :: nfro(:)
  integer, allocatable   :: nish(:)
  integer, allocatable   :: nash(:)
  integer, allocatable   :: nssh(:)
  integer, allocatable   :: ndel(:)
  integer, allocatable   :: norb(:)
  integer, allocatable   :: nbas(:)
  integer                :: nbast
  integer                :: nbsqt
  integer                :: ncmo
  end type type_inforb_short

  type (type_inforb), save, public :: inforb

  ! Orbital info to be read by the MOLCAS interface
  ! We only care about nfro, nish, nash, nssh
  ! and some other stuff
  ! initialised by the MOLCAS interface!
  type (type_inforb_short), save, public :: inforb_molcas

  contains

  subroutine initialize_inforb(nsym,norb,from_molcas)

    integer, intent(out):: nsym
    integer, intent(out):: norb
    integer             :: i, j, k, l
    integer             :: ierr, ierr1, ierr2, ierr3
    integer             :: norb_tmp
    character(len=  8)  :: cflag2
    character(len= 12)  :: cflag3
    character(len=100)  :: guess
    real*8, allocatable :: readbuf(:,:)
    logical, intent(in) :: from_molcas

    !> get data from file
    datadim(1)    = 1; datadim_bound = 1
    call hdf5_get_data(file_id(2),"nsym  ",datadim,inforb%nsym_b)
    nsym = inforb%nsym_b

  ! Do not read nfro and ndel here if we're called from MOLCAS
  ! It's already filled in initialize_inforb_molcas()
    if (.not.from_molcas) then
      allocate(inforb%norb_b(inforb%nsym_b),inforb%nfro_b(inforb%nsym_b),inforb%ndel_b(inforb%nsym_b))
      inforb%norb_b = 0; inforb%nfro_b = 0; inforb%ndel_b = 0

      datadim(1)   = inforb%nsym_b
      allocate(readbuf(nsym,3)); readbuf = -1
      call hdf5_get_data(file_id(2),"norb  ",datadim,readbuf(1,1))
      call hdf5_get_data(file_id(2),"nfro  ",datadim,readbuf(1,2))
      call hdf5_get_data(file_id(2),"ndel  ",datadim,readbuf(1,3))
      do i = 1, nsym
        inforb%norb_b(i) = nint(readbuf(i,1))
        inforb%nfro_b(i) = nint(readbuf(i,2))
        inforb%ndel_b(i) = nint(readbuf(i,3))
      end do
      deallocate(readbuf)
      !print *,' norb from MOLCAS ',inforb%norb_b
      !print *,' nfro from MOLCAS ',inforb%nfro_b
      !print *,' ndel from MOLCAS ',inforb%ndel_b
    end if

    allocate(inforb%nosh(nsym),inforb%iosh(nsym))
    inforb%nosh = 0; inforb%iosh = 0
    allocate(inforb%nish(nsym),inforb%nash(nsym),inforb%nssh(nsym),inforb%nfro(nsym), &
             inforb%ndel(nsym))
    inforb%nish = 0; inforb%nash = 0; inforb%nssh = 0; inforb%nfro = 0; inforb%ndel = 0

    norb = 0




    if (from_molcas) then
    ! fill the values from MOLCAS variables
      do i = 1, nsym
        inforb%iosh(i) = norb
        norb           = norb           + inforb_molcas%norb(i)
        inforb%nosh(i) = inforb%nosh(i) + inforb_molcas%norb(i)
        inforb%nfro(i) = inforb%nfro(i) + inforb_molcas%nfro(i)
        inforb%ndel(i) = inforb%ndel(i) + inforb_molcas%ndel(i)
      end do

      inforb%nish(1:nSym)=inforb_molcas%nish(1:nSym)
      inforb%nash(1:nSym)=inforb_molcas%nash(1:nSym)
      inforb%nssh(1:nSym)=inforb_molcas%nssh(1:nSym)
    else

    ! fill the values from newly read ijkl.h5 file from MOTRA
      do i = 1, nsym
        inforb%iosh(i) = norb
        norb           = norb           + inforb%norb_b(i)
        inforb%nosh(i) = inforb%nosh(i) + inforb%norb_b(i)
        inforb%nfro(i) = inforb%nfro(i) + inforb%nfro_b(i)
        inforb%ndel(i) = inforb%ndel(i) + inforb%ndel_b(i)
      end do

    !> retrieve information about the active space (in RASSCF)
    ! Leon: this is error prone, especially when multiple RasOrb files are present in the scratch directory!
    ! If we use the MOLCAS interface, it reads the orbital space information from the JobIph file and places
    ! it into nevpt2_cfg arrays, so we should copy the info from them. Otherwise, use the old RasOrb reading
    ! routine to determine the orbital space occupancy

      call system("ls *.RasOrb > ORB_GSS_NEVPT")

      open(unit=101,file="ORB_GSS_NEVPT",iostat=ierr)
      if(ierr.ne.0)then
        write(*,*)" RASSCF orbital file missing in scratch!"
        stop
      end if
      read(101,*)guess
      close(101)
      open(unit=102,file=trim(guess))
      cflag2=""

      i=0
      do
  !     if(i.eq.nsym)exit
        read(102,"(A8)",iostat=ierr1)cflag2
        if(ierr1.ne.0)then
          exit
        end if
        if(cflag2(1:6).eq."#INDEX")then
          cflag3=""
          do
            read(102,"(A12)",iostat=ierr2)cflag3
            if(ierr2 /= 0) exit
            if(cflag3(3:12).eq."1234567890")then
              i=i+1
              cycle
            end if
            do j=3,len(trim(cflag3))
              !if(cflag3(j:j).eq."f")then
              !   inforb%nfro(i) = inforb%nfro(i) + 1
              !> assign inactive and frozen RASSCF orbitals to nish!
              if(cflag3(j:j).eq."i" .or. cflag3(j:j).eq."f")then
                inforb%nish(i) = inforb%nish(i) + 1
              else if(cflag3(j:j).eq."2")then
                inforb%nash(i) = inforb%nash(i) + 1
              !> assign secondary and deleted RASSCF orbitals to nssh!
              else if(cflag3(j:j).eq."s" .or. cflag3(j:j).eq."d")then
                inforb%nssh(i) = inforb%nssh(i) + 1
              end if
            end do
          end do
        end if
      end do
      close(102)
    end if

    !> possibly correct the inactive shell distribution for frozen orbitals (in MOTRA)
    !> possibly correct the secondary shell distribution for deleted orbitals (in MOTRA)
    !> norb and nosh are correct since they are derived from the MOLTRA
    !  vector...
    do i = 1, nsym
      inforb%nish(i) = inforb%nish(i) - inforb%nfro(i)
      inforb%nssh(i) = inforb%nssh(i) - inforb%ndel(i)
    end do

    !print *, 'nfro, nish, nash, nssh ... ', inforb%nfro, inforb%nish, inforb%nash, inforb%nssh
    !call flush(6)

  end subroutine initialize_inforb
! ----------------------------------------------------------------------
  subroutine initialize_inforb_molcas(nsym)
    integer, intent(in) :: nsym

    allocate(inforb_molcas%nfro(nsym),inforb_molcas%ndel(nsym),inforb_molcas%nish(nsym), &
      inforb_molcas%nash(nsym),inforb_molcas%nssh(nsym),inforb_molcas%nbas(nsym),inforb_molcas%norb(nsym))

    inforb_molcas%nfro = 0; inforb_molcas%ndel = 0; inforb_molcas%norb = 0;
    inforb_molcas%nish = 0; inforb_molcas%nash  = 0; inforb_molcas%nssh  = 0;
    inforb_molcas%nbas = 0; inforb_molcas%nbast = 0; inforb_molcas%ncmo  = 0; inforb_molcas%nbsqt = 0
  end subroutine initialize_inforb_molcas

! ----------------------------------------------------------------------

  subroutine finalize_inforb()
    if(allocated(inforb%norb_b)) deallocate(inforb%norb_b)
    if(allocated(inforb%nfro_b)) deallocate(inforb%nfro_b)
    if(allocated(inforb%ndel_b)) deallocate(inforb%ndel_b)

    deallocate(inforb%nosh,inforb%nish,inforb%nash,inforb%nssh,inforb%nfro,inforb%ndel,inforb%iosh)
  end subroutine finalize_inforb

  subroutine finalize_inforb_molcas()
    if (allocated(inforb_molcas%nfro)) deallocate(inforb_molcas%nfro)
    if (allocated(inforb_molcas%ndel)) deallocate(inforb_molcas%ndel)
    if (allocated(inforb_molcas%nish)) deallocate(inforb_molcas%nish)
    if (allocated(inforb_molcas%nash)) deallocate(inforb_molcas%nash)
    if (allocated(inforb_molcas%nssh)) deallocate(inforb_molcas%nssh)
    if (allocated(inforb_molcas%nbas)) deallocate(inforb_molcas%nbas)
    if (allocated(inforb_molcas%norb)) deallocate(inforb_molcas%norb)
  end subroutine finalize_inforb_molcas
! ----------------------------------------------------------------------

  subroutine set_orbarray_pointers(                              &
                                   nsym,                         &
                                   norb,                         &
                                   ncore,                        &
                                   nact,                         &
                                   nvirt,                        &
                                   itsym,                        &
                                   nord_mn_c,                    &
                                   nord_mn_a,                    &
                                   nord_mn_v,                    &
                                   nord_nev2mol,                 &
                                   nord_mol2nev,                 &
                                   nc,                           &
                                   na,                           &
                                   nv,                           &
                                   nbe,                          &
                                   nmo,                          &
                                   ncmax,                        &
                                   nvmax                         &
                                  )

    integer, intent(in) :: nsym
    integer, intent(in) :: norb
    integer, intent(in) :: ncore
    integer, intent(in) :: nact
    integer, intent(in) :: nvirt
    integer, intent(in) :: itsym(2*norb+1)
    integer, intent(out):: nord_mn_c(ncore,nsym)
    integer, intent(out):: nord_mn_a(nact,nsym)
    integer, intent(out):: nord_mn_v(nvirt,nsym)
    integer, intent(out):: nord_nev2mol(norb)
    integer, intent(out):: nord_mol2nev(norb)
    integer, intent(out):: nbe(nsym)
    integer, intent(out):: nmo(nsym)
    integer, intent(out):: nc(nsym)
    integer, intent(out):: na(nsym)
    integer, intent(out):: nv(nsym)
    integer, intent(out):: ncmax
    integer, intent(out):: nvmax

    integer             :: isym,i,idum

      ncmax=0
      do isym=1,nsym
        nc(isym)=0
        do i=1,ncore
          if (itsym(i).eq.isym) then
            nc(isym)=nc(isym)+1
            nord_mn_c(nc(isym),isym)=i
          endif
        enddo
        if (ncmax.lt.nc(isym)) ncmax=nc(isym)
      enddo

      do isym=1,nsym
        na(isym)=0
        do i=ncore+1,ncore+nact
          if (itsym(i).eq.isym) then
            na(isym)=na(isym)+1
            nord_mn_a(na(isym),isym)=i
          endif
        enddo
      enddo

      nvmax=0
      do isym=1,nsym
        nv(isym)=0
        do i=ncore+nact+1,norb
          if (itsym(i).eq.isym) then
            nv(isym)=nv(isym)+1
            nord_mn_v(nv(isym),isym)=i
          endif
        enddo
        if (nvmax.lt.nv(isym)) nvmax=nv(isym)
      enddo

      nbe(1)=0
      do isym=1,nsym
        nmo(isym)=nc(isym)+na(isym)+nv(isym)
        if (isym.gt.1) nbe(isym)=nbe(isym-1)+nmo(isym-1)
      enddo

      do isym=1,nsym
        idum=0
        do i=1,ncore
          if (itsym(i).eq.isym) then
            idum=idum+1
            nord_mol2nev(nbe(isym)+idum)=i
          endif
        enddo
      enddo

      do isym=1,nsym
        idum=0
        do i=ncore+1,ncore+nact
          if (itsym(i).eq.isym) then
            idum=idum+1
            nord_mol2nev(nbe(isym)+nc(isym)+idum)=i
          endif
        enddo
      enddo

      do isym=1,nsym
        idum=0
        do i=ncore+nact+1,norb
          if (itsym(i).eq.isym) then
            idum=idum+1
            nord_mol2nev(nbe(isym)+nc(isym)+na(isym)+idum)=i
          endif
        enddo
      enddo

      do i=1,norb
        nord_nev2mol(nord_mol2nev(i))=i
      enddo

  end subroutine set_orbarray_pointers

! ----------------------------------------------------------------------

end module info_orbital_space
