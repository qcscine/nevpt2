module koopro4QD

      use indices_norm
      use koop_matrices
      use ord_utils
      use hdf5_utils
      use info_orbital_space
      use info_symmetry
      use nevpt2_cfg ! input parameters
      use ijkl_utils
#include "version.h"
!> use implicit none requires more work to do (do it if time permits)
!implicit none
#ifdef DMRG_NEVPT
      use qcmaquis_interface
      use qcmaquis_interface_cfg, only: qcmaquis_param
      use qcmaquis_interface_utility_routines, only : str
      use rdm_utils
#endif

public koopro4QD_driver

contains

      subroutine koopro4QD_driver(rel_ham,t13,from_molcas)

!     notes for/by stknecht...

!     trou: array with length of ncf containing a list of orbitals considered as (singly/empty)
!     occupied in the simple-minded Fermi occupation limit counting
!     active orbitals by order not symmetry irreps

!     part: array with length of ncf containing a list of orbitals
!     to which we excite starting from the simple-minded Fermi occupation limit counting
!     active orbitals by order not symmetry irreps
!
!     notation: oocupation(symmetry irrep)
!     example for HF: CAS(4,3): core "2(1) 2(1) 2(1)" 2(2) 2(2) 0(3)

!     nd: offset array for each determinant to set the pointer in trou
!     and part. --> not needed for koopmatrix construction
!     ne: # of electrons moved from the reference determinant
!     (simple-minded Fermi occupation limit) to create a new determinant

!     iocc: occupation matrix (norb,ncf) of all orbitals for each determinant

!     icomp: array with entries always set to 1

!     itsym: array giving the point group irreps for each orbital

!     ispin: array of length 2*norb(+1): first 1,norb = 0 (alpha) norb+1,2*norb = 1 (beta)
!     iorb:  array of length 2*norb(+1): 1,2*norb = #orb number - not needed for koopmatrix construction

      implicit real*8 (a-h,o-y),logical*1 (z)
      logical, intent(in) :: rel_ham
      real*4, intent(out) :: t13
      logical, intent(in) :: from_molcas
      allocatable dal(:,:,:)
      allocatable daaa(:,:,:,:,:,:,:)
      allocatable taa(:,:,:,:,:)
      allocatable amat(:,:,:,:,:,:,:),bmat(:,:,:,:,:),cmat(:,:,:,:,:),dmat(:,:,:)
      allocatable atmat(:,:,:,:,:,:,:),btmat(:,:,:,:,:),ctmat(:,:,:,:,:),dtmat(:,:,:)
      allocatable ro4(:,:,:,:,:,:)
      allocatable atwo(:,:,:,:)
      INTEGER*2 NE,TROU,PART
      allocatable :: ne(:),ispin(:),iorb(:)
      integer, allocatable :: nd(:)
      real*8,  allocatable :: c(:,:)
      integer :: NORB,NOCA,NOCB,METAT,NCF,II,I,j,istate,nnrot,mmetat,nsym
      integer :: nele,ncoppie
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      allocatable trou(:),part(:)
      real*4 tarray(2),t,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12
      logical z
      allocatable f(:),aj(:)
      character (len=120) :: revision
      allocatable iocc(:,:)
      integer iocc, icomp
      allocatable icomp(:),zcapos(:)
      integer :: maxlong,idum,idum2,ncapos,inext,icapos
      !NAMELIST /LEGGI/ NCORE,NACT,THR,NSPIN,file04,file50,file26,file25,file32,ZORDER
      !NAMELIST /LEGGI/ NCORE,NACT,THR,NSPIN,file04,file50,ZORDER

!Cele QD_new
      integer, allocatable :: nord_nev2mol(:),nord_mol2nev(:)
      integer, allocatable :: nord_mn_c(:,:),nord_mn_v(:,:)
      integer, allocatable :: nmo(:),nc(:),na(:),nv(:),nbe(:)
      integer              :: ncmax, nvmax
      integer              :: nvirt
!Cele QD_new
      integer, allocatable :: nord_mn_a(:,:)
      real*8 :: escf

      call cpu_time(t1)

      !> enable HDF5 support and create the file nevpt.h5; open ijkl.h5 (written in Molcas)
      call hdf5_init()
      call hdf5_create(filename, file_id(1))
      call hdf5_open(ijklname, file_id(2))

      !> make sure the array is allocated for a minimal input
      if(.not.from_molcas)then
        if(.not.allocated(MultGroup%State))&
        allocate(MultGroup%State(nr_states))
        do i = 1, nr_states
          MultGroup%State(i) = i
        end do
      end if

      !> single-state calculation is always state-specific
      if(nr_states == 1) skip_effective_ham = .true.

      !> initialize
      norb = 0; nsym = 0
      call initialize_inforb(nsym,norb,from_molcas)

      !> set its and itsym
      call initialize_infsym(nsym,norb)

      ncore = inforb%nisht
      nact  = inforb%nasht
      nvirt = norb-ncore-nact

#ifdef DMRG_NEVPT
      if(do_dmrg_pt) qcmaquis_param%L = nact

#endif

      !> initialize indices and norm array
      call initialize_indices_norm(nact)

      if(.not. do_dmrg_pt)then

        if(file04.eq.' ')then
           print* ,' file04 is mandatory input'
           stop 'missing files'
        endif

        open(4,file=file04,form='UNFORMATTED',status='OLD')
        !open(50,file=file50,form='UNFORMATTED',status='OLD')
        !maxlong=200000000
        !call filesplit('OPEN',50,maxlong,idum,idum2)
        !open(25,file=file25,form='UNFORMATTED',status='OLD')
        !open(26,file=file26,form='UNFORMATTED',status='OLD')
        READ (4) NORB,NOCA,NOCB,METAT,NCF,II,I

        allocate(ne(ncf))                ! exclusively needed for determinant construction of RHOs
        allocate(nd(ncf))                ! exclusively needed for determinant construction of RHOs
        allocate(icomp(ncf))             ! exclusively needed for determinant construction of RHOs
        allocate(trou(ii))               ! exclusively needed for determinant construction of RHOs
        allocate(part(ii))               ! exclusively needed for determinant construction of RHOs
        allocate(c(ncf,metat))           ! exclusively needed for determinant construction of RHOs
        allocate(ispin(norb+norb+1))     ! exclusively needed for determinant construction of RHOs
        allocate(iorb(norb+norb+1))      ! exclusively needed for determinant construction of RHOs
        ncoppie=metat*(metat-1)/2
        print *, 'ncoppie is ...',ncoppie

        write (6,*) 'FILE04 informations:'
        write (6,*) 'NORB = ',norb
        write (6,*) 'NOCA = ',noca
        write (6,*) 'NOCB = ',nocb
        write (6,*) 'METAT= ',metat
        write (6,*) 'NCF  = ',ncf
        write (6,*) 'II   = ',ii
        write (6,*) 'I    = ',i


        !> read full information from file FILE04
        rewind 4
        READ (4) NORB,NOCA,NOCB,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),            &
             (TROU(J),PART(J),J=1,II),((C(J,istate),J=1,ncf),istate=1            &
             ,metat),(E,J=1,METAT),(EMP,J=1,METAT),(EPS,J=1,NORB),(EPS,J=1       &
             ,NORB),NNROT,(dummy,J=1,NNROT),MMETAT,(ZHEFF,J=1                    &
             ,MMETAT),zz,zz,FACTOR,escf

        !> debug print
        !print*,'determinanti appena letti'
        !call printtp(ne,nd,trou,part,ncf,norb,c)
        !print *, 'initial factor is ...',factor

        !print *, 'tr == ',trou(1:ii)
        !print *, 'pa == ',part(1:ii)
        !print *, 'nd == ',nd(1:ncf)
        !print *, 'ne == ',ne(1:ncf)

        !> calculating an ion/radical? --> odd number of electrons
        zion=.false.
        if (ne(1).ne.0.and.part(ne(1)).gt.2*norb) zion=.true.
        if(zion) print*,'Caso con numero dispari di elettroni'
        nele=(noca-ncore)*2
        if(zion)nele=nele-1

        !> ensure to have the same number of alpha and beta orbitals
        nocb=min(noca,nocb)
        noca=nocb

        do  i=1,norb
           ispin(i)=0
           ispin(i+norb)=1
           iorb(i)=i
           iorb(i+norb)=i
        enddo
        if(zion)then
           iorb(norb+norb+1)=norb+1
           ispin(norb+norb+1)=1
        endif

        call esclass(nd,ne,trou,part,icomp,ncf,nocb,norb,nspin)
        ncapos=0
        inext=1
        allocate(zcapos(ncf))
        do i=1,ncf
           if(i.eq.inext)then
              zcapos(i)=.true.
              if (zorder) then
               inext=i+icomp(i)
              else
               inext=i+1
!cele   16-03-2010
               icomp(i)=1
              endif
           else
              zcapos(i)=.false.
           endif
           if(zcapos(i))ncapos=ncapos+1
        enddo

        !print '(a,i6,a)','Ci sono',ncapos,' capostipiti'
        !print *, ' zorder ==',zorder
        !print *, ' zcapos ==',zcapos(1:ncf)
        !print *, ' icomp  ==',icomp(1:ncf)
        allocate(iocc(1:norb,1:ncapos))
        icapos=0
        do i=1,ncf
           if(zcapos(i))then
              icapos=icapos+1
              call giveocc(iocc(1,icapos),i,nd,ne,trou,part,nocb,norb)
           endif
        enddo
        deallocate(zcapos)

      else ! do_dmrg_pt
        metat   = nr_states
        ncoppie = metat*(metat-1)/2
        nele    = nr_active_electrons
        allocate(ne(1))
        allocate(nd(1))
        allocate(icomp(1))
        allocate(trou(1))
        allocate(part(1))
        allocate(c(1,1))
        allocate(ispin(1))
        allocate(iorb(1))
        allocate(iocc(1,1))
        ne = -1; nd = -1; icomp = -1; trou = -1; part = -1; c = -1; ispin = -1; iorb = -1; iocc = -1
      end if

      !> allocate memory for the 1-4 particle RDMs
      allocate(dal(metat,nact,nact),taa(metat,nact,nact,nact,nact),daaa(metat,nact,nact,nact,nact,nact,nact))
      if(.not.no_4rdm_terms) then
        allocate(ro4(nwords,nact,nact,nact,nact,metat))
        ro4  = 0;
      end if
      daaa = 0; taa  = 0; dal  = 0

      call cpu_time(t)

      call state_specific_densities_driver(                       &
                                           c,                     &
                                           ro4,                   &
                                           daaa,                  &
                                           taa,                   &
                                           dal,                   &
                                           ncore,                 &
                                           nact,                  &
                                           thr,                   &
                                           icomp,                 &
                                           iocc,                  &
                                           ncapos,                &
                                           trou,                  &
                                           part,                  &
                                           ne,                    &
                                           nd,                    &
                                           ispin,                 &
                                           iorb,                  &
                                           ncf,                   &
                                           metat,                 &
                                           MultGroup%State,       &
                                           MultGroup%h5_file_name,&
                                           nwords,                &
                                           norb,                  &
                                           nele,                  &
                                           do_dmrg_pt             &
                                          )

      call cpu_time(t2)
      print '(/a,f9.2,a )',' time needed for construction of RDM matrices: ',t2-t,' sec.'

      deallocate(icomp)

      !> put data to file
      datadim(1)   = metat; datadim(2:7) = nact; datadim_bound = 7
      call hdf5_put_data_dp(file_id(1), "3-RDM" , datadim, daaa)

      !> set reordering pointers (needed for integral reading)
      allocate(nord_mn_c(ncore,nsym),nord_mn_v(nvirt,nsym))
      allocate(nord_mn_a(nact,nsym))
      allocate(nord_nev2mol(norb),nord_mol2nev(norb))
      allocate(nv(nsym),nc(nsym),na(nsym),nbe(nsym),nmo(nsym))

      call set_orbarray_pointers(                              &
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

      !> allocate integral memory and read-in integrals
      call readint(norb,nsym,nord_mol2nev,escf)

      deallocate(nord_mn_c,nord_mn_a,nord_mn_v)
      deallocate(nord_nev2mol,nord_mol2nev)
      deallocate(nv,nc,na,nmo,nbe)

      !> build array (2e-integrals+fock matrix) in array atwo
      allocate(atwo(1:nact,1:nact,1:nact,1:nact),aj(norb*(norb+1)/2))
      atwo = 0; aj = 0
      call get_aj_atwo(atwo,nact,ncore,norb,aj)

      !> deallocate integral memory
      call finalize_ijkl()

      allocate(f(norb*(norb+1)/2))
      f = 0
!---calcolo Koopmans esteso doppie ionizzazioni---------------------
      call koop_matrices_ctl(metat,nact,ncore,atwo,daaa,taa,dal,f,aj)

!------------------------------------------------------------------------------------------------------------
      !> build and save the matrices A - D
!------------------------------------------------------------------------------------------------------------
      if(.not.no_4rdm_terms) then
        allocate(amat(metat,nact,nact,nact,nact,nact,nact))
        allocate(atmat(metat,nact,nact,nact,nact,nact,nact))

        amat  = 0; atmat = 0

        call cpu_time(t5)

        call bamat(amat,atmat,atwo,ro4,daaa,taa,dal,nact,ncore,f,aj,norb,metat,nsym)

        !> put data to file
        datadim(1)   = metat; datadim(2:7) = nact; datadim_bound = 7
        call hdf5_put_data_dp(file_id(1), "amatIP", datadim, amat)
        call hdf5_put_data_dp(file_id(1), "amatEA", datadim, atmat)

        deallocate(ro4, amat, atmat)

        call cpu_time(t6)
        print '(a,f10.2,a/)',' construction of AMAT took ',t6-t5,' seconds'



        allocate(bmat(metat,nact,nact,nact,nact),btmat(metat,nact,nact,nact,nact))

        bmat  = 0; btmat = 0

        call bbmat(atwo,bmat,btmat,daaa,taa,dal,nact,ncore,f,aj,metat)

        !> put data to file
        datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5
        call hdf5_put_data_dp(file_id(1), "bmatIP", datadim, bmat)
        call hdf5_put_data_dp(file_id(1), "bmatEA", datadim, btmat)

        call cpu_time(t7)
        print '(a,f10.2,a/)',' construction of BMAT took ',t7-t6,' seconds'

        allocate(cmat(metat,nact,nact,nact,nact),ctmat(metat,nact,nact,nact,nact))
        cmat  = 0; ctmat = 0

        call bcmat(atwo,cmat,ctmat,daaa,taa,dal,nact,ncore,f,aj,metat)

        !> put data to file
        datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5
        call hdf5_put_data_dp(file_id(1), "cmatIP", datadim, cmat)
        call hdf5_put_data_dp(file_id(1), "cmatEA", datadim, ctmat)

        call cpu_time(t8)
        print '(a,f10.2,a/)',' construction of CMAT took ',t8-t7,' seconds'

        call fdiff(bmat,cmat,btmat,ctmat,nact,metat)
        deallocate(bmat, btmat, cmat, ctmat, daaa)

        call cpu_time(t9)
        allocate(dmat(metat,nact,nact),dtmat(metat,nact,nact))
        dmat  = 0; dtmat = 0

        call bdmat(atwo,dmat,dtmat,taa,dal,nact,ncore,f,aj,metat)

        !> put data to file
        datadim(1)   = metat; datadim(2:3) = nact; datadim_bound = 3
        call hdf5_put_data_dp(file_id(1), "dmatIP", datadim, dmat)
        call hdf5_put_data_dp(file_id(1), "dmatEA", datadim, dtmat)

        deallocate(dmat, dtmat, taa, dal, aj, f)

        call cpu_time(t10)
        print '(a,f10.2,a/)',' construction of DMAT took ',t10-t9,' seconds'
      end if

!------------------------------------------------------------------------------------------------------------
!     end of state-specific part (if metat = 1)
!------------------------------------------------------------------------------------------------------------

      !> exit if we run a state-specific calculation
      if(ncoppie == 0) goto 999

      allocate(daaa(ncoppie,nact,nact,nact,nact,nact,nact),taa(ncoppie,nact,nact,nact,nact),dal(ncoppie,nact,nact))

      daaa = 0; taa = 0; dal = 0

      !> simple-minded version of multi-state state-specific calculations
      if(.not.skip_effective_ham)then
        call cpu_time(t)

        call transition_densities_driver(                       &
                                         c,                     &
                                         daaa,                  &
                                         taa,                   &
                                         dal,                   &
                                         ncore,                 &
                                         nact,                  &
                                         thr,                   &
                                         iocc,                  &
                                         ncapos,                &
                                         trou,                  &
                                         part,                  &
                                         ne,                    &
                                         nd,                    &
                                         ispin,                 &
                                         iorb,                  &
                                         ncf,                   &
                                         metat,                 &
                                         MultGroup%State,       &
                                         MultGroup%h5_file_name,&
                                         ncoppie,               &
                                         norb,                  &
                                         nele,                  &
                                         do_dmrg_pt             &
                                        )

        call cpu_time(t12)
        print '(/a,f9.2,a )',' time needed for construction of transition RDM matrices: ',t12-t,' sec.'
      end if

      !> put data to file
      datadim(1)   = ncoppie; datadim(2:7) = nact; datadim_bound = 7
      call hdf5_put_data_dp(file_id(1), "3-TRDM" , datadim, daaa)
      datadim_bound = 5
      call hdf5_put_data_dp(file_id(1), "2-TRDM" , datadim, taa)
      datadim_bound = 3
      call hdf5_put_data_dp(file_id(1), "1-TRDM" , datadim, dal)

      !> deallocate the transition RDM stuff
      deallocate(daaa,taa,dal)

      !> deallocate all stuff needed for RDM construction / read-in of integrals
 999  deallocate(ne,nd,trou,part,c,ispin,iorb)

      !> done...
      call finalize_indices_norm()
      call finalize_inforb()
      call finalize_infsym()

      !> close files...
      if(.not.do_dmrg_pt)then
        close(4,  status="keep")
        !close(50, status="keep")
      end if

      !> close the files nevpt.h5+ijkl.h5 and turn off HDF5 support.
      call hdf5_close(file_id(1)); call hdf5_close(file_id(2)); call hdf5_exit()

      call cpu_time(t13)

      end subroutine koopro4QD_driver
!-------------------------------------------------------------------

      subroutine state_specific_densities_driver(                &
                                                 c,              &
                                                 ro4,            &
                                                 daaa,           &
                                                 taa,            &
                                                 dal,            &
                                                 ncore,          &
                                                 nact,           &
                                                 thr,            &
                                                 icomp,          &
                                                 iocc,           &
                                                 ncapos,         &
                                                 trou,           &
                                                 part,           &
                                                 ne,             &
                                                 nd,             &
                                                 ispin,          &
                                                 iorb,           &
                                                 ncf,            &
                                                 metat,          &
                                                 state_ptr,      &
                                                 name_ptr,       &
                                                 nwords,         &
                                                 norb,           &
                                                 nele,           &
                                                 do_dmrg_pt      &
                                                )

  implicit none

  real*8,    intent(in)  :: c(ncf,metat)
  real*8,    intent(out) :: ro4(nwords,nact,nact,nact,nact,metat)
  real*8,    intent(out) :: daaa(metat,nact,nact,nact,nact,nact,nact)
  real*8,    intent(out) :: taa(metat,nact,nact,nact,nact)
  real*8,    intent(out) :: dal(metat,nact,nact)
  real*8,    intent(in)  :: thr
  integer*2, intent(in)  :: ne(*), ispin(*), iorb(*), trou(*), part(*)
  integer,   intent(in)  :: nd(*), icomp(*), iocc(norb,ncapos)
  integer,   intent(in)  :: metat
  integer,   intent(in)  :: state_ptr(metat)
  integer,   intent(in)  :: ncore
  integer,   intent(in)  :: nact
  integer,   intent(in)  :: nwords
  integer,   intent(in)  :: ncf
  integer,   intent(in)  :: nele
  integer,   intent(in)  :: ncapos
  integer,   intent(in)  :: norb
  logical,   intent(in)  :: do_dmrg_pt
  character(len=*), intent(in) :: name_ptr(metat)
  character(len=:), allocatable :: resultfile

  integer                :: iroot, i, j, ip, k, l
  real*8, allocatable    :: tmp(:)

#ifdef DMRG_NEVPT
      if (do_dmrg_pt) then
        if (.not.(rdm_read.or.rdm_distributed)) then
          ! use only the first checkpoint to initialise
          call qcmaquis_interface_init_checkpoint(trim(curr_dir)//'/'//name_ptr(state_ptr(1)))
        end if
      end if
#endif

      select case(nele)
        case(1)
          print*,'just one active electron: are you kidding?'
          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
!             allocate(tmp((2*nact)**2))
!             do iroot = 1, metat
!               tmp = 0
!               call flush(6)
!               call qcmaquis_interface_ctl(                                                           &
!                                           task       = 'imp rdmY',                                   &
!                                           x1         = tmp,                                          &
!                                           ndim       = 2*nact,                                       &
!                                           checkpoint1= name_ptr(state_ptr(iroot)),                   &
!                                           checkpoint2= name_ptr(state_ptr(iroot)),                   &
!                                           msproj     = nspin-1,                                      &
!                                           msprojL    = nspin-1,                                      &
!                                           multiplet  = nspin-1,&! (nspin == 2*S+1) and we need 2*S   &
!                                           multipletL = nspin-1,&! (nspin == 2*S+1) and we need 2*S   &
!                                           rdm1       = .true.,                                       &
!                                           rdm2       = .false.                                       &
!                                          )
!               !> make spinless 1-RDM from alpha-alpha and beta-beta contributions: a-a + b-b
!               do i=1,nact; do j=1,nact
!                 !>                         a-a                               b-b
!                 dal(iroot,i,j) = tmp((2*i-1)+(2*j-2)*(2*nact)) + tmp(1+(2*i-1)+(2*j-1)*(2*nact))
!               end do; end do
!             end do
!             deallocate(tmp)
            do iroot = 1, metat
              ! a call to set_state is required because otherwise the DMRG interface pointer is not initialised
              ! for 3- and 4-RDM measurements that save to a file, this is not required as those function call
              ! set_state already
              call qcmaquis_interface_set_state(int(state_ptr(iroot)-1,4))
              call qcmaquis_interface_get_1rdm(dal(iroot,:,:))
            end do
#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb DEBUG rho1'
      do iroot=1,metat
        print *, 'rho1 for state ',iroot
        do i=1,nact; do j=1,nact
            if(abs(dal(iroot,i,j)) > 1.0d-16)&
            print '(i1,1i2,3x,d19.12)',&
            i,j,dal(iroot,i,j)
        enddo; enddo; call flush(6)
      enddo
#endif
#endif
          else
            call ro1(c,dal,ncore,nact,thr,trou,part,ne,nd,ispin,iorb)
          end if
        case(2)
          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
            ! Unfortunately 2e import still does not work for some reason
            ! It did not work in the previous interface either
            write (6,*) "2-electron case not implemented in DMRG-NEVPT2 yet!"
            stop
!             do iroot = 1, metat
!               call qcmaquis_interface_set_state(int(state_ptr(iroot)-1,4))
!               call qcmaquis_interface_get_2rdm(taa(iroot,:,:,:,:))
!               ! old interface
! !               call qcmaquis_interface_ctl(                                &
! !                                       task  = 'imp rdmX',                 &
! !                                       x2    = taa(iroot,1,1,1,1),         &
! !                                       mdim  = nact**4,                    &
! !                                       state = state_ptr(iroot),           &
! !                                       rdm2  = .true.                      &
! !                                      )
!             end do
#endif
          else
            call ro2(c,taa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb)
          end if
          call bro1(taa,dal,nact,nele,metat)
#ifdef DEBUG_DMRG_NEVPT
          do i=1,nact
             print '(7f12.7)',(dal(1,i,j),j=1,nact)
          enddo
#endif

        case(3)
          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
            do iroot = 1, metat
              if (.not.(rdm_read.or.rdm_distributed)) then
                ! TODO: replace this with the version that returns arrays directly
                ! instead of writing 4-RDM to a file
                call qcmaquis_interface_measure_and_save_3rdm(state_ptr(iroot)-1)
              end if
              write(*,*) "Reading 3-RDM for state "//trim(str(state_ptr(iroot)))//" from "&
                  //"checkpoint file "//trim(name_ptr(state_ptr(iroot)))
              call hdf5_read_rdm(chkpname = trim(curr_dir)//'/'//name_ptr(state_ptr(iroot)), &
                                   rdm3 = daaa(iroot,:,:,:,:,:,:))
! Old interface
!               call qcmaquis_interface_ctl(                                &
!                                       task  = 'imp rdmX',                 &
!                                       x3    = daaa(iroot,1,1,1,1,1,1),    &
!                                       odim  = nact,                       &
!                                       state = state_ptr(iroot),           &
!                                       rdm3  = .true.                      &
!                                      )
            end do
#endif
         else
           call ro3(c,daaa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb)
         endif
         call bro2(daaa,taa,nact,nele,metat)
         call bro1(taa,dal,nact,nele,metat)

        case default ! nele > 3

          if(nele <= 0) stop 'ERROR: negative or zero # active electrons'

          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
            ! If we're skipping the 4-rdm terms, then import 3-RDM from a file, otherwise import the 4-RDM
            if(.not.(no_4rdm_terms)) then
              do iroot = 1, metat
                if (.not.rdm_distributed) then
                  if (.not.rdm_read) then
                    ! evaluate 4rdm

                    ! explicit convertion to integer*4 for the parameter
                    ! TODO: replace this with the version that returns arrays directly
                    ! instead of writing 4-RDM to a file
                    call qcmaquis_interface_measure_and_save_4rdm(state_ptr(iroot)-1)
                  endif
                  write(*,*) "Reading 4-RDM for state "//trim(str(state_ptr(iroot)))//" from "&
                      //"checkpoint file "//trim(name_ptr(state_ptr(iroot)))
                  call hdf5_read_rdm(chkpname = trim(curr_dir)//'/'//trim(name_ptr(state_ptr(iroot))), &
                                    rdm4 = ro4(:,:,:,:,:,iroot))
                else ! rdm_distributed
                  resultfile = trim(molcas_project)// &
                      ".results_state."//trim(str(state_ptr(iroot)-1))//".h5"

                  ! Read paths from rdm_path + '4rdm-scratch.' + iroot
                  call read_distributed_rdm(rdm4=ro4(:,:,:,:,:,iroot), &
                      path=rdm_path//'/4rdm-scratch.'//trim(str(state_ptr(iroot)-1))//'/parts', &
                      result_name=resultfile)
                end if
#ifdef DEBUG_DMRG_NEVPT
      print *, 'imported rho4 for state ',state_ptr(iroot)
      do ip=1,nwords
        do i=1,nact
          do j=1,nact
            do k=1,nact
              do l=1,nact
                  if(abs(ro4(ip,i,j,k,l,iroot)) > 1.0d-16)&
       print '(i3,4i2,3x,d19.12)',ip,i,j,k,l,ro4(ip,i,j,k,l,iroot)
              enddo
            enddo
          enddo
        enddo
      enddo
#endif
              end do
            else
              if (rdm_distributed) write(*,*) "Ignoring distributed 4-RDM flag since no 4-RDM has been requested"

              do iroot = 1, metat
                if (.not.(rdm_read)) then
                  call qcmaquis_interface_measure_and_save_3rdm(state_ptr(iroot)-1)
                end if
                write(*,*) "Reading 3-RDM for state "//trim(str(state_ptr(iroot)))//" from "&
                    //"checkpoint file "//trim(name_ptr(state_ptr(iroot)))
                call hdf5_read_rdm(chkpname = trim(curr_dir)//'/'//name_ptr(state_ptr(iroot)), &
                                    rdm3 = daaa(iroot,:,:,:,:,:,:))
              end do
            end if
#endif
          else
            if(.not.(no_4rdm_terms)) then
              call rofour(c,ro4,ncore,nact,thr,icomp,iocc,ncapos,trou,part,ne,nd,ispin,iorb)
            end if
          end if
          if(.not.(no_4rdm_terms)) then
            call bro3(ro4,daaa,nact,nele,metat)
          end if
          call bro2(daaa,taa,nact,nele,metat)
          call bro1(taa,dal,nact,nele,metat)
      end select

      end subroutine state_specific_densities_driver

!-------------------------------------------------------------------

      subroutine transition_densities_driver(                &
                                             c,              &
                                             daaa,           &
                                             taa,            &
                                             dal,            &
                                             ncore,          &
                                             nact,           &
                                             thr,            &
                                             iocc,           &
                                             ncapos,         &
                                             trou,           &
                                             part,           &
                                             ne,             &
                                             nd,             &
                                             ispin,          &
                                             iorb,           &
                                             ncf,            &
                                             metat,          &
                                             state_ptr,      &
                                             name_ptr,       &
                                             ncoppie,        &
                                             norb,           &
                                             nele,           &
                                             do_dmrg_pt      &
                                            )

  implicit none

  real*8,    intent(in)  :: c(ncf,metat)
  real*8,    intent(out) :: daaa(ncoppie,nact,nact,nact,nact,nact,nact)
  real*8,    intent(out) :: taa(ncoppie,nact,nact,nact,nact)
  real*8,    intent(out) :: dal(ncoppie,nact,nact)
  real*8,    intent(in)  :: thr
  integer*2, intent(in)  :: ne(*), ispin(*), iorb(*), trou(*), part(*)
  integer,   intent(in)  :: nd(*), iocc(norb,ncapos)
  integer,   intent(in)  :: metat
  integer,   intent(in)  :: state_ptr(metat)
  integer,   intent(in)  :: ncore
  integer,   intent(in)  :: nact
  integer,   intent(in)  :: ncoppie
  integer,   intent(in)  :: ncf
  integer,   intent(in)  :: nele
  integer,   intent(in)  :: ncapos
  integer,   intent(in)  :: norb
  logical,   intent(in)  :: do_dmrg_pt
  character(len=*), intent(in) :: name_ptr(metat)
  character(len=:), allocatable :: resultfile

  integer                :: iroot, jroot, i, j, offset_tdm
  integer                :: i1, i2, i3, i4, i5, i6
!   real*8, allocatable    :: tmp(:,:,:,:,:,:)
!   real*8, allocatable    :: tmp1(:)

      offset_tdm = 0
      select case(nele)
        case(1)
          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
            write (*,*) "QD-NEVPT2 for 1 electrons (1-TDM calculation) not implemented yet!"
            stop
!             allocate(tmp1((2*nact)**2))
!             do iroot = 1, metat
!               do jroot = 1, iroot-1
!                 offset_tdm = offset_tdm + 1
!                 !> safety check
!                 if(offset_tdm > ncoppie) exit
!                 tmp1 = 0
!                 call qcmaquis_interface_ctl(                                                           &
!                                             task       = 'imp rdmY',                                   &
!                                             x1         = tmp1,                                         &
!                                             ndim       = 2*nact,                                       &
!                                             checkpoint1= name_ptr(state_ptr(jroot)),                   &
!                                             checkpoint2= name_ptr(state_ptr(iroot)),                   &
!                                             msproj     = nspin-1,                                      &
!                                             msprojL    = nspin-1,                                      &
!                                             multiplet  = nspin-1,&! (nspin == 2*S+1) and we need 2*S   &
!                                             multipletL = nspin-1,&! (nspin == 2*S+1) and we need 2*S   &
!                                             rdm1       = .true.,                                       &
!                                             rdm2       = .false.                                       &
!                                            )
!                 !> make spinless 1-TDM from alpha-alpha and beta-beta contributions: a-a + b-b
!                 do i=1,nact; do j=1,nact
!                   !>                               a-a                                b-b
!                   dal(offset_tdm,i,j) = tmp1((2*i-1)+(2*j-2)*(2*nact)) + tmp1(1+(2*i-1)+(2*j-1)*(2*nact))
!                 end do; end do
!               end do
!             end do
!             deallocate(tmp1)
#endif
          else
            call ro1off(c,dal,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,ncoppie)
          end if
        case(2)
          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
            write (*,*) "QD-NEVPT2 for 2 electrons (2-TDM calculation) not implemented yet!"
            stop
!             do iroot = 1, metat
!               do jroot = 1, iroot-1
!                 offset_tdm = offset_tdm + 1
!                 !> safety check
!                 if(offset_tdm > ncoppie) exit
!                 call qcmaquis_interface_ctl(                                 &
!                                         task   = 'imp rdmY',                 &
!                                         x2     = taa(offset_tdm,1,1,1,1),    &
!                                         mdim   = nact**4,                    &
!                                         state  = state_ptr(jroot),           &
!                                         stateL = state_ptr(iroot),           &
!                                         rdm2   = .true.                      &
!                                        )
!               end do
!             end do
#endif
          else
            call ro2off(c,taa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,ncoppie)
          end if
          call bro1(taa,dal,nact,nele,ncoppie)

        case default ! nele > 2

          if(do_dmrg_pt)then
#ifdef DMRG_NEVPT
            do iroot = 1, metat
              do jroot = 1, iroot-1
                offset_tdm = offset_tdm + 1
                !> safety check
                if(offset_tdm > ncoppie) exit
                  if (.not.rdm_distributed) then
                    if (.not.rdm_read) then
                      call qcmaquis_interface_measure_and_save_trans3rdm( &
                        state_ptr(jroot)-1, state_ptr(iroot)-1)
                    end if
                    resultfile = trim(curr_dir)//'/'//trim(molcas_project)//&
                        ".trans3rdm."//trim(str(state_ptr(jroot)-1))//"_"//trim(str(state_ptr(iroot)-1))//".h5"
                    write(*,*)"Reading 3-TDM for states "//trim(str(state_ptr(jroot)))//" and "//trim(str(state_ptr(iroot)))//&
                      " from file "//trim(resultfile)
                    call hdf5_read_rdm_from_resfile(result_name=trim(resultfile), &
                                                    rdm3 = daaa(offset_tdm,:,:,:,:,:,:),  &
                                                    trans_rdm = .true.)
                  else ! rdm_distributed
                    resultfile = trim(molcas_project)// &
                        ".results_state."//trim(str(state_ptr(jroot)-1))//".h5"

                    ! Read paths from rdm_path + '3rdm-scratch.' + iroot
                    call read_distributed_rdm(rdm3=daaa(offset_tdm,:,:,:,:,:,:), &
                        path=rdm_path//'/3rdm-scratch.'//trim(str(state_ptr(jroot)-1))//"."//trim(str(state_ptr(iroot)-1))// &
                        "/parts", trans_rdm = .true., result_name=resultfile)

                  endif
              end do
            end do

#endif
          else
            call ro3off(c,daaa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,ncoppie)
          end if
          call bro2(daaa,taa,nact,nele,ncoppie)
          call bro1(taa,dal,nact,nele,ncoppie)
      end select

      end subroutine transition_densities_driver

!-------------------------------------------------------------------

      subroutine rofour(c,amat,ncore,nact,thr,icomp,iocc,ncaptot,trou,part,ne,nd,ispin,iorb)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension amat(nwords,nact,nact,nact,nact,metat)
      dimension c(ncf,*)
      integer   icomp(*)
      integer   iocc(norb,ncaptot)
      real*8 , allocatable :: val(:)
      integer, allocatable :: iocca(:),ioccb(:)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      dimension ne(*),nd(*),ispin(*),iorb(*)
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      integer*8 nhits,ijcal
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      integer op,pp,a,b,op1,op2,op3,pp1,pp2,pp3
      integer segno,segno0,segno1,segno2,segno3,segno4,segno5
      nocc=noca-ncore

      allocate(iocca(norb),ioccb(norb),val(metat))

      !print *, 'RHO 4 DEBUG print'
      !print *, '-----------------'
      !print*,'nocc,noca,ncore',nocc,noca,ncore
      !print*,'nact=',nact,' nocc=',noca-ncore,' ncf=',ncf
      !print*,'metat=',metat
      !print*,'ncaptot=',ncaptot
      !print *, '-----------------'
      nfours=0
      ntrips=0
      ndoubs=0
      nsings=0
      nhits=0
      nhitsd=0
      mnext=1
      iconfi=1
      mcapos=0
      znocalc=.false.
      mfin=0
      ijcal=0
      do  mcap=1,ncaptot ! loop over all determinants
         minit=mfin+1
         mfin=minit+icomp(minit)-1
!        if (mod(mcap,100).eq.0)print *,'mcap=',mcap
!        if (mod(mcap,100).eq.0)print *,'nhits,nhitsd',nhits,nhitsd
         call flush(6)
         nec=ne(minit)
         ndm=nd(minit)
         nfin=0
         do  ncap=1,mcap      !do 2
            ninit=nfin+1
            nfin=ninit+icomp(ninit)-1
            nedif=abs(ne(minit)-ne(ninit))
            znocalc=(nedif > 4)
            if(znocalc)nhitsd=nhitsd+1
            if(.not.znocalc)znocalc=(ndiff(iocc(ncore+1,mcap),iocc(ncore+1,ncap),nact) > 8)
            if(znocalc)nhits=nhits+1
            if(znocalc)goto 2
            do m=minit,mfin ! do 3
               nec=ne(m)
               ndm=nd(m)
               call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndm,trou,part,iorb,ispin)
!              if(abs(c(m)).lt.thr) goto 3
               if(mcap.eq.ncap)then
                  nbeg=minit
                  nend=m-1
!                  ci2=cmval*c(m)
!                  val=ci2
                  call case0
               else
                  nbeg=ninit
                  nend=nfin
               endif
               do n=nbeg,nend
!                  if(abs(c(n)).lt.thr) goto 4
                  ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,4,zseg)
                  ijcal=ijcal+1
                  if(ij.eq.0)nbdif=5
                  if(nbdif.gt.4)goto 4
                  if(zseg)then
                     segno=-1
                  else
                     segno=1
                  endif
!                 print *, 'n, m is == ',n,m
                  if(nbdif.eq.1)then
!                    print *, 'diff 1: nv, nw',nv(1), nw(1)
                     nsings=nsings+1
                     n1=iorb(nv(1))
                     ns1=ispin(nv(1))
                     n2=iorb(nw(1))
                     ns2=ispin(nw(1))
                     call case10
                  elseif(nbdif.eq.2)then
                     ndoubs=ndoubs+1
!                    print *, 'diff 2: nv, nw',nv(1:2), nw(1:2)
                     call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorb(nv(1))
                     ns1=ispin(nv(1))
                     n2=iorb(nv(2))
                     ns2=ispin(nv(2))
                     n3=iorb(nw(1))
                     ns3=ispin(nw(1))
                     n4=iorb(nw(2))
                     ns4=ispin(nw(2))
                     call case20
                  elseif(nbdif.eq.3)then
!                    print *, 'diff 3: nv, nw',nv(1:nbdif), nw(1:nbdif)
                     ntrips=ntrips+1
                     call ord3(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorb(nv(1))
                     ns1=ispin(nv(1))
                     n2=iorb(nv(2))
                     ns2=ispin(nv(2))
                     n3=iorb(nv(3))
                     ns3=ispin(nv(3))
                     n4=iorb(nw(1))
                     ns4=ispin(nw(1))
                     n5=iorb(nw(2))
                     ns5=ispin(nw(2))
                     n6=iorb(nw(3))
                     ns6=ispin(nw(3))
                     call case30
                  else          !nbdif.eq.4
!---  renzo debug
                     nfours=nfours+1
!----
!                    print *, 'diff 4: nv, nw',nv(1:nbdif), nw(1:nbdif)
                     call ord4(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorb(nv(1))
                     ns1=ispin(nv(1))
                     n2=iorb(nv(2))
                     ns2=ispin(nv(2))
                     n3=iorb(nv(3))
                     ns3=ispin(nv(3))
                     n4=iorb(nv(4))
                     ns4=ispin(nv(4))
                     n5=iorb(nw(1))
                     ns5=ispin(nw(1))
                     n6=iorb(nw(2))
                     ns6=ispin(nw(2))
                     n7=iorb(nw(3))
                     ns7=ispin(nw(3))
                     n8=iorb(nw(4))
                     ns8=ispin(nw(4))
                  call case40
                  endif
 4             enddo
 3          enddo
 2       enddo
 1    enddo
!---renzo debug
      print*,'number of single    differences=',nsings
      print*,'number of double    differences=',ndoubs
      print*,'number of triple    differences=',ntrips
      print*,'number of quadruple differences=',nfours
      print*,'number of jdifgen calls ',ijcal
!---contributi diagonali
      deallocate(iocca,ioccb,val)
      return
      contains
      subroutine case0
!      print*,'case 0 m=',m
!---  contributo a daaaa
      do istate=1,metat
         valcm=c(m,istate)
         cij=valcm*c(m,istate)
         val(istate)=cij
      enddo
      do i=1,nocca
         ia=iocca(i)
         do j=i+1,nocca
            ic=iocca(j)
            do k=j+1,nocca
               ie=iocca(k)
               do l=k+1,nocca
                  ig=iocca(l)
                  call ordsame(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact,m,m)
               enddo
            enddo
         enddo
      enddo
!---  contributo a dbbbb
      do i=1,noccb
         ia=ioccb(i)
         do j=i+1,noccb
            ic=ioccb(j)
            do k=j+1,noccb
               ie=ioccb(k)
               do l=k+1,noccb
                  ig=ioccb(l)
                  call ordsame(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact,m,m)
               enddo
            enddo
         enddo
      enddo
!---  contributi a daaab,daaba,dabaa,dbaaa
      do i=1,nocca
         ia=iocca(i)
         do j=i+1,nocca
            ic=iocca(j)
            do k=j+1,nocca
               ie=iocca(k)
               do l=1,noccb
                  ig=ioccb(l)
                  call ord31(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact,m,m)
               enddo
            enddo
         enddo
      enddo
!---  contributi a dbbba,dbbab,dbabb,dabbb
      do i=1,noccb
         ia=ioccb(i)
         do j=i+1,noccb
            ic=ioccb(j)
            do k=j+1,noccb
               ie=ioccb(k)
               do l=1,nocca
                  ig=iocca(l)
                  call ord31(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact,m,m)
               enddo
            enddo
         enddo
      enddo
!---  contributo a daabb,dbbaa,dabab,dbaba,dabba,dbaab
      do i=1,nocca
         ia=iocca(i)
         do j=i+1,nocca
            ic=iocca(j)
            do k=1,noccb
               ie=ioccb(k)
               do l=k+1,noccb
                  ig=ioccb(l)
                  call ord22(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact,m,m)
               enddo
            enddo
         enddo
      enddo
      end subroutine case0
      subroutine case10
      n1=n1-ncore
      n2=n2-ncore
      do istate=1,metat
         valcm=c(m,istate)
         cij=valcm*c(n,istate)*segno
         val(istate)=cij
      enddo
      if(ns1.eq.1.and.ns2.eq.1)then
!---a sinistra il piu grande
         if(n1.ge.n2)then
            nl=n1
            nr=n2
         else
            nl=n2
            nr=n1
         endif
!---  contributo a daaaa
         do i=1,nocca
            ia=iocca(i)
            do j=i+1,nocca
               ic=iocca(j)
               do k=j+1,nocca
                  ie=iocca(k)
                  if(ia.ne.n1.and.ic.ne.n1.and.ie.ne.n1.and.ia.ne.n2.and.ic.ne.n2.and.ie.ne.n2)then
                     call ordsame(ia,ic,ie,nl,ia,ic,ie,nr,val,amat,nact,m,n)
                  endif
               enddo
            enddo
         enddo
!---  contributo a dbbba
         do i=1,noccb
            ia=ioccb(i)
            do j=i+1,noccb
               ic=ioccb(j)
               do k=j+1,noccb
                  ie=ioccb(k)
                  call ord31(ia,ic,ie,nl,ia,ic,ie,nr,val,amat,nact,m,n)
               enddo
            enddo
         enddo
!---  contributi a daaab,daaba,dabaa,dbaaa
         do i=1,nocca
            ia=iocca(i)
            do j=i+1,nocca
               ic=iocca(j)
               if(ia.ne.n1.and.ic.ne.n1.and.ia.ne.n2.and.ic.ne.n2)then
                  do k=1,noccb
                     ie=ioccb(k)
!     write (6,*) ia,ic,ie,nl,nr
                     call ord31(ia,ic,nl,ie,ia,ic,nr,ie,val,amat,nact,m,n)
                  enddo
               endif
            enddo
         enddo
!---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
         do i=1,nocca
            ia=iocca(i)
            do j=1,noccb
               ic=ioccb(j)
               do k=j+1,noccb
                  ie=ioccb(k)
                  if(ia.ne.n1.and.ia.ne.n2)then
!     write (6,*) 'calling ord22'
                     call ord22(ia,nl,ic,ie,ia,nr,ic,ie,val,amat,nact,m,n)
                  endif
               enddo
            enddo
         enddo
      elseif(ns1.eq.0.and.ns2.eq.0)then
!---  a sinistra il piu grande
         if(n1.ge.n2)then
            nl=n1
            nr=n2
         else
            nl=n2
            nr=n1
         endif
!---  contributo a dbbbb
         do i=1,noccb
            ia=ioccb(i)
            do j=i+1,noccb
               ic=ioccb(j)
               do k=j+1,noccb
                  ie=ioccb(k)
                  if(ia.ne.n1.and.ic.ne.n1.and.ie.ne.n1.and.ia.ne.n2.and.ic.ne.n2.and.ie.ne.n2)then
!     write (6,*) 'calling ordsame'
                     call ordsame(ia,ic,ie,nl,ia,ic,ie,nr,val,amat,nact,m,n)
                  endif
               enddo
            enddo
         enddo
!---  contributo a daaab
         do i=1,nocca
            ia=iocca(i)
            do j=i+1,nocca
               ic=iocca(j)
               do k=j+1,nocca
                  ie=iocca(k)
!     write (6,*) 'calling ord31'
                  call ord31(ia,ic,ie,nl,ia,ic,ie,nr,val,amat,nact,m,n)
               enddo
            enddo
         enddo
!---  contributi a dbbba,dbbab,dbabb,dabbb
         do i=1,noccb
            ia=ioccb(i)
            do j=i+1,noccb
               ic=ioccb(j)
               if(ia.ne.n1.and.ic.ne.n1.and.ia.ne.n2.and.ic.ne.n2)then
                  do k=1,nocca
                     ie=iocca(k)
!     write (6,*) 'calling ord31'
                     call ord31(ia,ic,nl,ie,ia,ic,nr,ie,val,amat,nact,m,n)
                  enddo
               endif
            enddo
         enddo
!---  contributi a dbbaa,daabb,dbaba,dabab,dbaab,dabba
         do i=1,noccb
            ia=ioccb(i)
            do j=1,nocca
               ic=iocca(j)
               do k=j+1,nocca
                  ie=iocca(k)
                  if(ia.ne.n1.and.ia.ne.n2)then
!     write (6,*) 'calling ord22'
                     call ord22(ic,ie,ia,nl,ic,ie,ia,nr,val,amat,nact,m,n)
                  endif
               enddo
            enddo
         enddo
      endif
!     goto 2
      end subroutine case10
!     due differenze
! 20         continue
!            contains
      subroutine case20
      n1=n1-ncore
      n2=n2-ncore
      n3=n3-ncore
      n4=n4-ncore
      do istate=1,metat
         valcm=c(m,istate)
         cij=valcm*c(n,istate)*segno
         val(istate)=cij
      enddo
!--- n1l,n2l a sinistra
      if((n1+n2).ge.(n3+n4))then
         n1l=n1
         n2l=n2
         n1r=n3
         n2r=n4
      else
         n1l=n3
         n2l=n4
         n1r=n1
         n2r=n2
      endif
      if(ns1.eq.1.and.ns2.eq.1)then
!---  contributo a daaaa
         do i=1,nocca
            ia=iocca(i)
            do j=i+1,nocca
               ic=iocca(j)
               if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4.and.&
                  ic.ne.n1.and.ic.ne.n2.and.ic.ne.n3.and.ic.ne.n4)then
                  call ordsame(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val,amat,nact,m,n)
               endif
            enddo
         enddo
!---  contributo a daaab,daaba,dabaa,dbaaa
         do i=1,nocca
            ia=iocca(i)
            if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4)then
               do j=1,noccb
                  ic=ioccb(j)
                  call ord31(ia,n1l,n2l,ic,ia,n1r,n2r,ic,val,amat,nact,m,n)
               enddo
            endif
         enddo
!---  contributo a dbbaa,daabb,dabab,dbaba,dabba,dbaab
         do i=1,noccb
            ia=ioccb(i)
            do j=i+1,noccb
               ic=ioccb(j)
               call ord22(n1l,n2l,ia,ic,n1r,n2r,ia,ic,val,amat,nact,m,n)
            enddo
         enddo
      elseif(ns1.eq.0.and.ns2.eq.0)then
!---  contributo a dbbbb
         do i=1,noccb
            ia=ioccb(i)
            do j=i+1,noccb
               ic=ioccb(j)
               if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4.and.&
                  ic.ne.n1.and.ic.ne.n2.and.ic.ne.n3.and.ic.ne.n4)then
                  call ordsame(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val,amat,nact,m,n)
               endif
            enddo
         enddo
!---  contributo a dbbba,dbbab,dbabb,dabbb
         do i=1,noccb
            ia=ioccb(i)
            if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4)then
               do j=1,nocca
                  ic=iocca(j)
                  call ord31(ia,n1l,n2l,ic,ia,n1r,n2r,ic,val,amat,nact,m,n)
               enddo
            endif
         enddo
!---  contributo a daabb,dbbaa,dbaba,dabab,dbaab,dabba
         do i=1,nocca
            ia=iocca(i)
            do j=i+1,nocca
               ic=iocca(j)
               call ord22(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val,amat,nact,m,n)
            enddo
         enddo
      elseif(ns1.eq.1.and.ns2.eq.0)then
!---  contributi a daaab,daaba,dabaa,dbaaa
         do i=1,nocca
            ia=iocca(i)
            do j=i+1,nocca
               ic=iocca(j)
               if(ia.ne.n1.and.ia.ne.n3.and.ic.ne.n1.and.ic.ne.n3)then
!---- daaab
                  call ord31(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val,amat,nact,m,n)
               endif
            enddo
         enddo
!---  contributi a dbbba,dbbab,dbabb,dabbb
         do i=1,noccb
            ia=ioccb(i)
            do j=i+1,noccb
               ic=ioccb(j)
               if(ia.ne.n2.and.ia.ne.n4.and.ic.ne.n2.and.ic.ne.n4)then
!---- dbbba
                  call ord31(ia,ic,n2l,n1l,ia,ic,n2r,n1r,val,amat,nact,m,n)
               endif
            enddo
         enddo
!---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
         do i=1,nocca
            ia=iocca(i)
            do j=1,noccb
               ic=ioccb(j)
               if(ia.ne.n1.and.ia.ne.n3.and.ic.ne.n2.and.ic.ne.n4)then
!---- daabb
                  call ord22(ia,n1l,ic,n2l,ia,n1r,ic,n2r,val,amat,nact,m,n)
               endif
            enddo
         enddo
      endif
!     goto 2
      end subroutine case20
!     tre differenze
! 30         continue
!            contains
      subroutine case30
      n1=n1-ncore
      n2=n2-ncore
      n3=n3-ncore
      n4=n4-ncore
      n5=n5-ncore
      n6=n6-ncore
      do istate=1,metat
         valcm=c(m,istate)
         cij=valcm*c(n,istate)*segno
         val(istate)=cij
      enddo
!--- n1l,n2l,n3l a sinistra
      if((n1+n2+n3).ge.(n4+n5+n6))then
         n1l=n1
         n2l=n2
         n3l=n3
         n1r=n4
         n2r=n5
         n3r=n6
      else
         n1l=n4
         n2l=n5
         n3l=n6
         n1r=n1
         n2r=n2
         n3r=n3
      endif
!     casi aaa, bbb, aab, bba
      if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
!---  contributo daaaa
         do i=1,nocca
            ia=iocca(i)
            if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4.and.ia.ne.n5.and.ia.ne.n6)then
               call ordsame(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r,val,amat,nact,m,n)
            endif
         enddo
!---  contributi a daaab,daaba,dabaa,dbaaa
         do i=1,noccb
            ia=ioccb(i)
            call ord31(n1l,n2l,n3l,ia,n1r,n2r,n3r,ia,val,amat,nact,m,n)
         enddo
      elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
!---  contributo dbbbb
         do i=1,noccb
            ia=ioccb(i)
            if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4.and.ia.ne.n5.and.ia.ne.n6)then
               call ordsame(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r,val,amat,nact,m,n)
            endif
         enddo
!---  contributi a dbbba,dbbab,dbabb,dabbb
         do i=1,nocca
            ia=iocca(i)
            call ord31(n1l,n2l,n3l,ia,n1r,n2r,n3r,ia,val,amat,nact,m,n)
         enddo
      elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
!---  contributi a daaab,daaba,dabaa,dbaaa
         do i=1,nocca
            ia=iocca(i)
            if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n4.and.ia.ne.n5)then
!---- daaab
               call ord31(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r,val,amat,nact,m,n)
            endif
         enddo
!---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
         do i=1,noccb
            ia=ioccb(i)
            if(ia.ne.n3.and.ia.ne.n6)then
!---- daabb
               call ord22(n1l,n2l,ia,n3l,n1r,n2r,ia,n3r,val,amat,nact,m,n)
            endif
         enddo
      elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
!---  contributi a dbbba,dbbab,dbabb,dabbb
         do i=1,noccb
            ia=ioccb(i)
            if(ia.ne.n2.and.ia.ne.n3.and.ia.ne.n5.and.ia.ne.n6)then
!---- dbbba
               call ord31(ia,n2l,n3l,n1l,ia,n2r,n3r,n1r,val,amat,nact,m,n)
            endif
         enddo
!---  contributi a dbbaa,daabb,dbaba,dabab,dbaab,dabba
         do i=1,nocca
            ia=iocca(i)
            if(ia.ne.n1.and.ia.ne.n4)then
!---- dbbaa
               call ord22(n1l,ia,n2l,n3l,n1r,ia,n2r,n3r,val,amat,nact,m,n)
            endif
         enddo
      endif
!     goto 2
      end subroutine case30
!     quattro differenze
!     40         continue
!     contains
      subroutine case40
      n1=n1-ncore
      n2=n2-ncore
      n3=n3-ncore
      n4=n4-ncore
      n5=n5-ncore
      n6=n6-ncore
      n7=n7-ncore
      n8=n8-ncore
      do istate=1,metat
         valcm=c(m,istate)
         cij=valcm*c(n,istate)*segno
         val(istate)=cij
      enddo
!---  n1l,n2l,n3l a sinistra
      if((n1+n2+n3+n4).ge.(n5+n6+n7+n8))then
         n1l=n1
         n2l=n2
         n3l=n3
         n4l=n4
         n1r=n5
         n2r=n6
         n3r=n7
         n4r=n8
      else
         n1l=n5
         n2l=n6
         n3l=n7
         n4l=n8
         n1r=n1
         n2r=n2
         n3r=n3
         n4r=n4
      endif
      if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1.and.ns4.eq.1)then
!     contributo a daaaa
         call ordsame(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat,nact,m,n)
      elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0.and.ns4.eq.0)then
!     contributo a dbbbb
         call ordsame(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat,nact,m,n)
      elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1.and.ns4.eq.0)then
!     contributo a daaab
         call ord31(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat,nact,m,n)
      elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0.and.ns4.eq.0)then
!     contributo a dbbba
         call ord31(n2l,n3l,n4l,n1l,n2r,n3r,n4r,n1r,val,amat,nact,m,n)
      elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0.and.ns4.eq.0)then
!     contributo a daabb
         call ord22(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat,nact,m,n)
      endif
      end subroutine case40

      end subroutine rofour
!--------------------------------------------------------------------
      subroutine ro3(c,daaa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daaa(metat,nact,nact,nact,nact,nact,nact)
      dimension c(ncf,*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*)
      integer segno
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      integer, allocatable :: iocca(:), ioccb(:)
      real*8 , allocatable :: cij(:)
      real*8 , allocatable :: dbbb(:,:,:,:,:,:,:), daab(:,:,:,:,:,:,:), dbba(:,:,:,:,:,:,:)

      nocc=noca-ncore
      print*,'nocc,noca',nocc,noca
      print*,'subr. ro3'
      print*,'nact=',nact,' nocc=',nocc,' ncf=',ncf

      allocate(iocca(norb))
      allocate(ioccb(norb))
      allocate(dbbb(metat,nact,nact,nact,nact,nact,nact))
      allocate(daab(metat,nact,nact,nact,nact,nact,nact))
      allocate(dbba(metat,nact,nact,nact,nact,nact,nact))
      allocate(cij(metat))

!---azzeramento
      call zeroe(daaa,metat*nact*6)
      call zeroe(dbbb,metat*nact*6)
      call zeroe(daab,metat*nact*6)
      call zeroe(dbba,metat*nact*6)

      do 1 m=1,ncf ! # determinants
         if (mod(m,100).eq.0)print *,'m=',m
         call flush(6)
!         if(abs(c(m)).lt.thr)goto 1
         nec=ne(m)
         ndm=nd(m)
         call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndm,trou,part,iorb,ispin)
         do 2 n=1,m
!            if(abs(c(n)).lt.thr)goto 2
            do istate=1,metat
               cij(istate)=c(m,istate)*c(n,istate)
            enddo
            if(n.eq.m)goto 3
            ndn=nd(n)
            ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,3,zseg)
            if(ij.eq.0)nbdif=4
            if(nbdif.gt.3)goto 2
            if(zseg)then
               segno=-1
            else
               segno=1
            endif
            if(nbdif.eq.1)then
               n1=iorb(nv(1))
               ns1=ispin(nv(1))
               n2=iorb(nw(1))
               ns2=ispin(nw(1))
               goto 10
            elseif(nbdif.eq.2)then
               call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
               n1=iorb(nv(1))
               ns1=ispin(nv(1))
               n2=iorb(nv(2))
               ns2=ispin(nv(2))
               n3=iorb(nw(1))
               ns3=ispin(nw(1))
               n4=iorb(nw(2))
               ns4=ispin(nw(2))
               goto 20
            else
               call ord3(nv,nw,zseg,segno) !prima gli alfa, poi i beta
               n1=iorb(nv(1))
               ns1=ispin(nv(1))
               n2=iorb(nv(2))
               ns2=ispin(nv(2))
               n3=iorb(nv(3))
               ns3=ispin(nv(3))
               n4=iorb(nw(1))
               ns4=ispin(nw(1))
               n5=iorb(nw(2))
               ns5=ispin(nw(2))
               n6=iorb(nw(3))
               ns6=ispin(nw(3))
               goto 30
            endif
 3          continue
!            ci2=c(m)*c(m)
!     costruzione di daaa
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do k=j+1,nocca
                     ka=iocca(k)
                     do istate=1,metat
                        daaa(istate,ia,ja,ka,ia,ja,ka)=daaa(istate,ia,ja,ka,ia,ja,ka)+cij(istate)
                     enddo
                     call perm3(daaa,ia,ja,ka,ia,ja,ka,nact,metat)
                  enddo
               enddo
            enddo
!     costruzione di dbbb
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do k=j+1,noccb
                     kb=ioccb(k)
                     do istate=1,metat
                        dbbb(istate,ib,jb,kb,ib,jb,kb)=dbbb(istate,ib,jb,kb,ib,jb,kb)+cij(istate)
                     enddo
                     call perm3(dbbb,ib,jb,kb,ib,jb,kb,nact,metat)
                  enddo
               enddo
            enddo
!     costruzione di daab
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do k=1,noccb
                     kb=ioccb(k)
                     do istate=1,metat
                        daab(istate,ia,ja,kb,ia,ja,kb)=daab(istate,ia,ja,kb,ia,ja,kb)+cij(istate)
                     enddo
                     call perm2(daab,ia,ja,kb,ia,ja,kb,nact,metat)
                  enddo
               enddo
            enddo
!     costruzione di dbba
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do k=1,nocca
                     ka=iocca(k)
                     do istate=1,metat
                        dbba(istate,ib,jb,ka,ib,jb,ka)=dbba(istate,ib,jb,ka,ib,jb,ka)+cij(istate)
                     enddo
                     call perm2(dbba,ib,jb,ka,ib,jb,ka,nact,metat)
                  enddo
               enddo
            enddo
            goto 2
!     una differenza
 10         continue
            n1=n1-ncore
            n2=n2-ncore
!            cij=c(m)*c(n)*segno
            if(ns1.eq.1.and.ns2.eq.1)then
!---contributo a daaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ja=iocca(j)
                     if(ia.ne.n1.and.ja.ne.n1.and.ia.ne.n2.and.ja.ne.n2)then
                        do istate=1,metat
                           cijp=cij(istate)*segno
                           daaa(istate,n1,ia,ja,n2,ia,ja)=daaa(istate,n1,ia,ja,n2,ia,ja)+cijp
                           daaa(istate,n2,ia,ja,n1,ia,ja)=daaa(istate,n2,ia,ja,n1,ia,ja)+cijp
                        enddo
                        call perm3(daaa,n1,ia,ja,n2,ia,ja,nact,metat)
                        call perm3(daaa,n2,ia,ja,n1,ia,ja,nact,metat)
                     endif
                  enddo
               enddo
!----contributo a dbba
               do i=1,noccb
                  ib=ioccb(i)
                  do j=i+1,noccb
                     jb=ioccb(j)
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        dbba(istate,ib,jb,n1,ib,jb,n2)=dbba(istate,ib,jb,n1,ib,jb,n2)+cijp
                     enddo
                     call perm2(dbba,ib,jb,n1,ib,jb,n2,nact,metat)
                  enddo
               enddo
!----contributo a daab
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2)then
                     do j=1,noccb
                        jb=ioccb(j)
                        do istate=1,metat
                           cijp=cij(istate)*segno
                           daab(istate,n1,ia,jb,n2,ia,jb)=daab(istate,n1,ia,jb,n2,ia,jb)+cijp
                        enddo
                        call perm2(daab,n1,ia,jb,n2,ia,jb,nact,metat)
                     enddo
                  endif
               enddo
            elseif(ns1.eq.0.and.ns2.eq.0)then
!----contributo a dbbb
               do i=1,noccb
                  ib=ioccb(i)
                  do j=i+1,noccb
                     jb=ioccb(j)
                     if(ib.ne.n1.and.jb.ne.n1.and.ib.ne.n2.and.jb.ne.n2)then
                        do istate=1,metat
                           cijp=cij(istate)*segno
                           dbbb(istate,n1,ib,jb,n2,ib,jb)=dbbb(istate,n1,ib,jb,n2,ib,jb)+cijp
                           dbbb(istate,n2,ib,jb,n1,ib,jb)=dbbb(istate,n2,ib,jb,n1,ib,jb)+cijp
                        enddo
                        call perm3(dbbb,n1,ib,jb,n2,ib,jb,nact,metat)
                        call perm3(dbbb,n2,ib,jb,n1,ib,jb,nact,metat)
                     endif
                  enddo
               enddo
!---contrinuto a daab
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ja=iocca(j)
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        daab(istate,ia,ja,n1,ia,ja,n2)=daab(istate,ia,ja,n1,ia,ja,n2)+cijp
                     enddo
                     call perm2(daab,ia,ja,n1,ia,ja,n2,nact,metat)
                  enddo
               enddo
!----contributo a dbba
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n1.and.ib.ne.n2)then
                     do j=1,nocca
                        ja=iocca(j)
                        do istate=1,metat
                           cijp=cij(istate)*segno
                           dbba(istate,n1,ib,ja,n2,ib,ja)=dbba(istate,n1,ib,ja,n2,ib,ja)+cijp
                        enddo
                        call perm2(dbba,n1,ib,ja,n2,ib,ja,nact,metat)
                     enddo
                  endif
               enddo
            endif
            goto 2
!     due differenze
 20         continue
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
!            cij=c(m)*c(n)*segno
            if(ns1.eq.1.and.ns2.eq.1)then
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4)then
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        daaa(istate,n1,n2,ia,n3,n4,ia)=daaa(istate,n1,n2,ia,n3,n4,ia)+cijp
                        daaa(istate,n3,n4,ia,n1,n2,ia)=daaa(istate,n3,n4,ia,n1,n2,ia)+cijp
                     enddo
                     call perm3(daaa,n1,n2,ia,n3,n4,ia,nact,metat)
                     call perm3(daaa,n3,n4,ia,n1,n2,ia,nact,metat)
                  endif
               enddo
               do i=1,noccb
                  ib=ioccb(i)
                  do istate=1,metat
                     cijp=cij(istate)*segno
                     daab(istate,n1,n2,ib,n3,n4,ib)=daab(istate,n1,n2,ib,n3,n4,ib)+cijp
                  enddo
                  call perm2(daab,n1,n2,ib,n3,n4,ib,nact,metat)
               enddo

            elseif(ns1.eq.0.and.ns2.eq.0)then
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n1.and.ib.ne.n2.and.ib.ne.n3.and.ib.ne.n4)then
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        dbbb(istate,n1,n2,ib,n3,n4,ib)=dbbb(istate,n1,n2,ib,n3,n4,ib)+cijp
                        dbbb(istate,n3,n4,ib,n1,n2,ib)=dbbb(istate,n3,n4,ib,n1,n2,ib)+cijp
                     enddo
                     call perm3(dbbb,n1,n2,ib,n3,n4,ib,nact,metat)
                     call perm3(dbbb,n3,n4,ib,n1,n2,ib,nact,metat)
                  endif
               enddo
               do i=1,nocca
                  ia=iocca(i)
                  do istate=1,metat
                     cijp=cij(istate)*segno
                     dbba(istate,n1,n2,ia,n3,n4,ia)=dbba(istate,n1,n2,ia,n3,n4,ia)+cijp
                  enddo
                  call perm2(dbba,n1,n2,ia,n3,n4,ia,nact,metat)
               enddo
            elseif(ns1.eq.1.and.ns2.eq.0)then
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n3)then
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        daab(istate,n1,ia,n2,n3,ia,n4)=daab(istate,n1,ia,n2,n3,ia,n4)+cijp
                     enddo
                     call perm2(daab,n1,ia,n2,n3,ia,n4,nact,metat)
                  endif
               enddo
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n2.and.ib.ne.n4)then
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        dbba(istate,n2,ib,n1,n4,ib,n3)=dbba(istate,n2,ib,n1,n4,ib,n3)+cijp
                     enddo
                     call perm2(dbba,n2,ib,n1,n4,ib,n3,nact,metat)
                  endif
               enddo
            endif
            goto 2
!     tre differenze
 30         continue
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            n5=n5-ncore
            n6=n6-ncore
!            cij=c(m)*c(n)*segno
!     casi aaa, bbb, aab, bba
            if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
               do istate=1,metat
                  cijp=cij(istate)*segno
                  daaa(istate,n1,n2,n3,n4,n5,n6)=daaa(istate,n1,n2,n3,n4,n5,n6)+cijp
                  daaa(istate,n4,n5,n6,n1,n2,n3)=daaa(istate,n4,n5,n6,n1,n2,n3)+cijp
               enddo
               call perm3(daaa,n1,n2,n3,n4,n5,n6,nact,metat)
               call perm3(daaa,n4,n5,n6,n1,n2,n3,nact,metat)
            elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
               do istate=1,metat
                  cijp=cij(istate)*segno
                  dbbb(istate,n1,n2,n3,n4,n5,n6)=dbbb(istate,n1,n2,n3,n4,n5,n6)+cijp
                  dbbb(istate,n4,n5,n6,n1,n2,n3)=dbbb(istate,n4,n5,n6,n1,n2,n3)+cijp
               enddo
               call perm3(dbbb,n1,n2,n3,n4,n5,n6,nact,metat)
               call perm3(dbbb,n4,n5,n6,n1,n2,n3,nact,metat)
            elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
               do istate=1,metat
                  cijp=cij(istate)*segno
                  daab(istate,n1,n2,n3,n4,n5,n6)=daab(istate,n1,n2,n3,n4,n5,n6)+cijp
               enddo
               call perm2(daab,n1,n2,n3,n4,n5,n6,nact,metat)
            elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
               do istate=1,metat
                  cijp=cij(istate)*segno
                  dbba(istate,n2,n3,n1,n5,n6,n4)=dbba(istate,n2,n3,n1,n5,n6,n4)+cijp
               enddo
               call perm2(dbba,n2,n3,n1,n5,n6,n4,nact,metat)
            endif
 2       continue
 1    continue
!--costruzione matrice ro3 spinless
      do istate=1,metat
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do ip=1,nact
                  do jp=1,nact
                     do kp=1,nact
!                       do istate=1,metat
                           daaa(istate,i,j,k,ip,jp,kp)=daaa(istate,i,j,k &
                                ,ip,jp,kp)+dbbb(istate,i,j,k,ip,jp,kp)   &
                                +daab(istate,i,j,k,ip,jp,kp)+dbba(istate &
                                ,i,j,k,ip,jp,kp)+daab(istate,i,k,j,ip,kp &
                                ,jp)+daab(istate,j,k,i,jp,kp,ip)         &
                                +dbba(istate,i,k,j,ip,kp,jp)+dbba(istate &
                                ,k,j,i,kp,jp,ip)
!                       enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
#ifdef DEBUG_DMRG_NEVPT
       print *, 'blubb DEBUG real constructed rho3'
       do istate=1,metat
       print *, 'rho3 for state ',istate
       do i=1,nact
         do j=1,nact
           do k=1,nact
             do ip=1,nact
               do jp=1,nact
                 do kp=1,nact
                   if(abs(daaa(istate,i,j,k,ip,jp,kp)) > 1.0d-16)&
       print '(i1,5i2,3x,d19.12)',&
       i-1,j-1,k-1,ip-1,jp-1,kp-1,daaa(istate,i,j,k,ip,jp,kp)
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo

       enddo ! istate
#endif

      deallocate(iocca,ioccb,cij,dbbb,daab,dbba)

      end subroutine ro3
!--------------------------------------------------------------------
      subroutine ro3off(c,daaa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer, intent(in) :: ncoppie,nact
      dimension daaa(ncoppie,nact,nact,nact,nact,nact,nact)
      dimension c(ncf,*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*)
      integer segno
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      integer, allocatable :: iocca(:), ioccb(:)
      real*8 , allocatable :: cij(:), cji(:)
      real*8 , allocatable :: dbbb(:,:,:,:,:,:,:), daab(:,:,:,:,:,:,:), dbba(:,:,:,:,:,:,:)
      nocc=noca-ncore
      print*,'nocc,noca',nocc,noca
      print*,'subr. ro3'
      print*,'nact=',nact,' nocc=',nocc,' ncf=',ncf
      allocate(iocca(norb),ioccb(norb),cij(ncoppie),cji(ncoppie))
      allocate(dbbb(ncoppie,nact,nact,nact,nact,nact,nact))
      allocate(daab(ncoppie,nact,nact,nact,nact,nact,nact))
      allocate(dbba(ncoppie,nact,nact,nact,nact,nact,nact))
!---azzeramento
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  do m=1,nact
                     do n=1,nact
                        do ijcou=1,ncoppie
                           daaa(ijcou,i,j,k,l,m,n)=0.d0
                           dbbb(ijcou,i,j,k,l,m,n)=0.d0
                           daab(ijcou,i,j,k,l,m,n)=0.d0
                           dbba(ijcou,i,j,k,l,m,n)=0.d0
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      do 1 m=1,ncf
         if (mod(m,100).eq.0)print *,'m=',m
         call flush(6)
!         if(abs(c(m)).lt.thr)goto 1
         nec=ne(m)
         ndm=nd(m)
         call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndm,trou,part,iorb,ispin)
         do 2 n=1,m
!            if(abs(c(n)).lt.thr)goto 2
            ijcou=0
            do istate=1,metat
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  ijcou=ijcou+1
                  cij(ijcou)=c(m,jstate)*c(n,istate)
                  cji(ijcou)=c(n,jstate)*c(m,istate)
!                  cij(ijcou)=c(m,istate)*c(n,jstate)  !!mind this
!                  cji(ijcou)=c(n,istate)*c(m,jstate)  !!mind this
               enddo
            enddo
            if(n.eq.m)goto 3
            ndn=nd(n)
            ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,3,zseg)
            if(ij.eq.0)nbdif=4
            if(nbdif.gt.3)goto 2
            if(zseg)then
               segno=-1
            else
               segno=1
            endif
            if(nbdif.eq.1)then
               n1=iorb(nv(1))
               ns1=ispin(nv(1))
               n2=iorb(nw(1))
               ns2=ispin(nw(1))
               goto 10
            elseif(nbdif.eq.2)then
               call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
               n1=iorb(nv(1))
               ns1=ispin(nv(1))
               n2=iorb(nv(2))
               ns2=ispin(nv(2))
               n3=iorb(nw(1))
               ns3=ispin(nw(1))
               n4=iorb(nw(2))
               ns4=ispin(nw(2))
               goto 20
            else
               call ord3(nv,nw,zseg,segno) !prima gli alfa, poi i beta
               n1=iorb(nv(1))
               ns1=ispin(nv(1))
               n2=iorb(nv(2))
               ns2=ispin(nv(2))
               n3=iorb(nv(3))
               ns3=ispin(nv(3))
               n4=iorb(nw(1))
               ns4=ispin(nw(1))
               n5=iorb(nw(2))
               ns5=ispin(nw(2))
               n6=iorb(nw(3))
               ns6=ispin(nw(3))
               goto 30
            endif
 3          continue
!            ci2=c(m)*c(m)
!     costruzione di daaa
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do k=j+1,nocca
                     ka=iocca(k)
                     do ijcou=1,ncoppie
                        daaa(ijcou,ia,ja,ka,ia,ja,ka)=daaa(ijcou,ia,ja,ka,ia,ja,ka)+cij(ijcou)
                     enddo
                     call perm3off(daaa,ia,ja,ka,ia,ja,ka,nact,ncoppie)
                  enddo
               enddo
            enddo
!     costruzione di dbbb
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do k=j+1,noccb
                     kb=ioccb(k)
                     do ijcou=1,ncoppie
                        dbbb(ijcou,ib,jb,kb,ib,jb,kb)=dbbb(ijcou,ib,jb,kb,ib,jb,kb)+cij(ijcou)
                     enddo
                     call perm3off(dbbb,ib,jb,kb,ib,jb,kb,nact,ncoppie)
                  enddo
               enddo
            enddo
!     costruzione di daab
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do k=1,noccb
                     kb=ioccb(k)
                     do ijcou=1,ncoppie
                        daab(ijcou,ia,ja,kb,ia,ja,kb)=daab(ijcou,ia,ja,kb,ia,ja,kb)+cij(ijcou)
                     enddo
                     call perm2off(daab,ia,ja,kb,ia,ja,kb,nact,ncoppie)
                  enddo
               enddo
            enddo
!     costruzione di dbba
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do k=1,nocca
                     ka=iocca(k)
                     do ijcou=1,ncoppie
                        dbba(ijcou,ib,jb,ka,ib,jb,ka)=dbba(ijcou,ib,jb,ka,ib,jb,ka)+cij(ijcou)
                     enddo
                     call perm2off(dbba,ib,jb,ka,ib,jb,ka,nact,ncoppie)
                  enddo
               enddo
            enddo
            goto 2
!     una differenza
 10         continue
            n1=n1-ncore
            n2=n2-ncore
!            cij=c(m)*c(n)*segno
            if(ns1.eq.1.and.ns2.eq.1)then
!---contributo a daaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ja=iocca(j)
                     if(ia.ne.n1.and.ja.ne.n1.and.ia.ne.n2.and.ja.ne.n2)then
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           daaa(ijcou,n1,ia,ja,n2,ia,ja)=daaa(ijcou,n1,ia,ja,n2,ia,ja)+cijp
                           daaa(ijcou,n2,ia,ja,n1,ia,ja)=daaa(ijcou,n2,ia,ja,n1,ia,ja)+cjip
                        enddo
                        call perm3off(daaa,n1,ia,ja,n2,ia,ja,nact,ncoppie)
                        call perm3off(daaa,n2,ia,ja,n1,ia,ja,nact,ncoppie)
                     endif
                  enddo
               enddo
!----contributo a dbba
               do i=1,noccb
                  ib=ioccb(i)
                  do j=i+1,noccb
                     jb=ioccb(j)
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           dbba(ijcou,ib,jb,n1,ib,jb,n2)=dbba(ijcou,ib,jb,n1,ib,jb,n2)+cijp
                           dbba(ijcou,ib,jb,n2,ib,jb,n1)=dbba(ijcou,ib,jb,n2,ib,jb,n1)+cjip
                        enddo
                        call perm2off(dbba,ib,jb,n1,ib,jb,n2,nact,ncoppie)
                        call perm2off(dbba,ib,jb,n2,ib,jb,n1,nact,ncoppie)
                     enddo
                  enddo
!----contributo a daab
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2)then
                     do j=1,noccb
                        jb=ioccb(j)
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           daab(ijcou,n1,ia,jb,n2,ia,jb)=daab(ijcou,n1,ia,jb,n2,ia,jb)+cijp
                           daab(ijcou,n2,ia,jb,n1,ia,jb)=daab(ijcou,n2,ia,jb,n1,ia,jb)+cjip
                        enddo
                        call perm2off(daab,n1,ia,jb,n2,ia,jb,nact,ncoppie)
                        call perm2off(daab,n2,ia,jb,n1,ia,jb,nact,ncoppie)
                     enddo
                  endif
               enddo
            elseif(ns1.eq.0.and.ns2.eq.0)then
!----contributo a dbbb
               do i=1,noccb
                  ib=ioccb(i)
                  do j=i+1,noccb
                     jb=ioccb(j)
                     if(ib.ne.n1.and.jb.ne.n1.and.ib.ne.n2.and.jb.ne.n2)then
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           dbbb(ijcou,n1,ib,jb,n2,ib,jb)=dbbb(ijcou,n1,ib,jb,n2,ib,jb)+cijp
                           dbbb(ijcou,n2,ib,jb,n1,ib,jb)=dbbb(ijcou,n2,ib,jb,n1,ib,jb)+cjip
                        enddo
                        call perm3off(dbbb,n1,ib,jb,n2,ib,jb,nact,ncoppie)
                        call perm3off(dbbb,n2,ib,jb,n1,ib,jb,nact,ncoppie)
                     endif
                  enddo
               enddo
!---contrinuto a daab
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ja=iocca(j)
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           daab(ijcou,ia,ja,n1,ia,ja,n2)=daab(ijcou,ia,ja,n1,ia,ja,n2)+cijp
                           daab(ijcou,ia,ja,n2,ia,ja,n1)=daab(ijcou,ia,ja,n2,ia,ja,n1)+cjip
                        enddo
                        call perm2off(daab,ia,ja,n1,ia,ja,n2,nact,ncoppie)
                        call perm2off(daab,ia,ja,n2,ia,ja,n1,nact,ncoppie)
                  enddo
               enddo
!----contributo a dbba
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n1.and.ib.ne.n2)then
                     do j=1,nocca
                        ja=iocca(j)
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           dbba(ijcou,n1,ib,ja,n2,ib,ja)=dbba(ijcou,n1,ib,ja,n2,ib,ja)+cijp
                           dbba(ijcou,n2,ib,ja,n1,ib,ja)=dbba(ijcou,n2,ib,ja,n1,ib,ja)+cjip
                        enddo
                        call perm2off(dbba,n1,ib,ja,n2,ib,ja,nact,ncoppie)
                        call perm2off(dbba,n2,ib,ja,n1,ib,ja,nact,ncoppie)
                     enddo
                  endif
               enddo
            endif
            goto 2
!     due differenze
 20         continue
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
!            cij=c(m)*c(n)*segno
            if(ns1.eq.1.and.ns2.eq.1)then
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4)then
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           daaa(ijcou,n1,n2,ia,n3,n4,ia)=daaa(ijcou,n1,n2,ia,n3,n4,ia)+cijp
                           daaa(ijcou,n3,n4,ia,n1,n2,ia)=daaa(ijcou,n3,n4,ia,n1,n2,ia)+cjip
                        enddo
                        call perm3off(daaa,n1,n2,ia,n3,n4,ia,nact,ncoppie)
                        call perm3off(daaa,n3,n4,ia,n1,n2,ia,nact,ncoppie)
                  endif
               enddo
               do i=1,noccb
                  ib=ioccb(i)
                  do ijcou=1,ncoppie
                     cijp=cij(ijcou)*segno
                     cjip=cji(ijcou)*segno
                     daab(ijcou,n1,n2,ib,n3,n4,ib)=daab(ijcou,n1,n2,ib,n3,n4,ib)+cijp
                     daab(ijcou,n3,n4,ib,n1,n2,ib)=daab(ijcou,n3,n4,ib,n1,n2,ib)+cjip
                  enddo
                  call perm2off(daab,n1,n2,ib,n3,n4,ib,nact,ncoppie)
                  call perm2off(daab,n3,n4,ib,n1,n2,ib,nact,ncoppie)
               enddo

            elseif(ns1.eq.0.and.ns2.eq.0)then
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n1.and.ib.ne.n2.and.ib.ne.n3.and.ib.ne.n4)then
                     do ijcou=1,ncoppie
                        cijp=cij(ijcou)*segno
                        cjip=cji(ijcou)*segno
                        dbbb(ijcou,n1,n2,ib,n3,n4,ib)=dbbb(ijcou,n1,n2,ib,n3,n4,ib)+cijp
                        dbbb(ijcou,n3,n4,ib,n1,n2,ib)=dbbb(ijcou,n3,n4,ib,n1,n2,ib)+cjip
                     enddo
                     call perm3off(dbbb,n1,n2,ib,n3,n4,ib,nact,ncoppie)
                     call perm3off(dbbb,n3,n4,ib,n1,n2,ib,nact,ncoppie)
                  endif
               enddo
               do i=1,nocca
                  ia=iocca(i)
                  do ijcou=1,ncoppie
                     cijp=cij(ijcou)*segno
                     cjip=cji(ijcou)*segno
                     dbba(ijcou,n1,n2,ia,n3,n4,ia)=dbba(ijcou,n1,n2,ia,n3,n4,ia)+cijp
                     dbba(ijcou,n3,n4,ia,n1,n2,ia)=dbba(ijcou,n3,n4,ia,n1,n2,ia)+cjip
                  enddo
                  call perm2off(dbba,n1,n2,ia,n3,n4,ia,nact,ncoppie)
                  call perm2off(dbba,n3,n4,ia,n1,n2,ia,nact,ncoppie)
               enddo
            elseif(ns1.eq.1.and.ns2.eq.0)then
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n3)then
                     do ijcou=1,ncoppie
                        cijp=cij(ijcou)*segno
                        cjip=cji(ijcou)*segno
                        daab(ijcou,n1,ia,n2,n3,ia,n4)=daab(ijcou,n1,ia,n2,n3,ia,n4)+cijp
                        daab(ijcou,n3,ia,n4,n1,ia,n2)=daab(ijcou,n3,ia,n4,n1,ia,n2)+cjip
                     enddo
                  endif
                  call perm2off(daab,n1,ia,n2,n3,ia,n4,nact,ncoppie)
                  call perm2off(daab,n3,ia,n4,n1,ia,n2,nact,ncoppie)
               enddo
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n2.and.ib.ne.n4)then
                     do ijcou=1,ncoppie
                        cijp=cij(ijcou)*segno
                        cjip=cji(ijcou)*segno
                        dbba(ijcou,n2,ib,n1,n4,ib,n3)=dbba(ijcou,n2,ib,n1,n4,ib,n3)+cijp
                        dbba(ijcou,n4,ib,n3,n2,ib,n1)=dbba(ijcou,n4,ib,n3,n2,ib,n1)+cjip
                     enddo
                  endif
                  call perm2off(dbba,n2,ib,n1,n4,ib,n3,nact,ncoppie)
                  call perm2off(dbba,n4,ib,n3,n2,ib,n1,nact,ncoppie)
               enddo
            endif
            goto 2
!     tre differenze
 30         continue
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            n5=n5-ncore
            n6=n6-ncore
!            cij=c(m)*c(n)*segno
!     casi aaa, bbb, aab, bba
            if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
               do ijcou=1,ncoppie
                  cijp=cij(ijcou)*segno
                  cjip=cji(ijcou)*segno
                  daaa(ijcou,n1,n2,n3,n4,n5,n6)=daaa(ijcou,n1,n2,n3,n4,n5,n6)+cijp
                  daaa(ijcou,n4,n5,n6,n1,n2,n3)=daaa(ijcou,n4,n5,n6,n1,n2,n3)+cjip
               enddo
               call perm3off(daaa,n1,n2,n3,n4,n5,n6,nact,ncoppie)
               call perm3off(daaa,n4,n5,n6,n1,n2,n3,nact,ncoppie)
            elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
               do ijcou=1,ncoppie
                  cijp=cij(ijcou)*segno
                  cjip=cji(ijcou)*segno
                  dbbb(ijcou,n1,n2,n3,n4,n5,n6)=dbbb(ijcou,n1,n2,n3,n4,n5,n6)+cijp
                  dbbb(ijcou,n4,n5,n6,n1,n2,n3)=dbbb(ijcou,n4,n5,n6,n1,n2,n3)+cjip
               enddo
               call perm3off(dbbb,n1,n2,n3,n4,n5,n6,nact,ncoppie)
               call perm3off(dbbb,n4,n5,n6,n1,n2,n3,nact,ncoppie)
            elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
               do ijcou=1,ncoppie
                  cijp=cij(ijcou)*segno
                  cjip=cji(ijcou)*segno
                  daab(ijcou,n1,n2,n3,n4,n5,n6)=daab(ijcou,n1,n2,n3,n4,n5,n6)+cijp
                  daab(ijcou,n4,n5,n6,n1,n2,n3)=daab(ijcou,n4,n5,n6,n1,n2,n3)+cjip
               enddo
               call perm2off(daab,n1,n2,n3,n4,n5,n6,nact,ncoppie)
               call perm2off(daab,n4,n5,n6,n1,n2,n3,nact,ncoppie)
            elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
               do ijcou=1,ncoppie
                  cijp=cij(ijcou)*segno
                  cjip=cji(ijcou)*segno
                  dbba(ijcou,n2,n3,n1,n5,n6,n4)=dbba(ijcou,n2,n3,n1,n5,n6,n4)+cijp
                  dbba(ijcou,n5,n6,n4,n2,n3,n1)=dbba(ijcou,n5,n6,n4,n2,n3,n1)+cjip
               enddo
               call perm2off(dbba,n2,n3,n1,n5,n6,n4,nact,ncoppie)
               call perm2off(dbba,n5,n6,n4,n2,n3,n1,nact,ncoppie)
            endif
 2       continue
 1    continue
!--costruzione matrice ro3 spinless
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do ip=1,nact
                  do jp=1,nact
                     do kp=1,nact
                        do ijcou=1,ncoppie
                           daaa(ijcou,i,j,k,ip,jp,kp)=daaa(ijcou,i,j,k    &
                                ,ip,jp,kp)+dbbb(ijcou,i,j,k,ip,jp,kp)     &
                                +daab(ijcou,i,j,k,ip,jp,kp)+dbba(ijcou,i  &
                                ,j,k,ip,jp,kp)+daab(ijcou,i,k,j,ip,kp,jp  &
                                )+daab(ijcou,j,k,i,jp,kp,ip)+dbba(ijcou   &
                                ,i,k,j,ip,kp,jp)+dbba(ijcou,k,j,i,kp,jp   &
                                ,ip)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(iocca,ioccb,cij,cji,dbbb,daab,dbba)
      end subroutine ro3off
!-------------------------------------------------------
      subroutine ro2(c,d2alal,ncore,nact,thr,trou,part,ne,nd,ispin,iorb)
      implicit real*8 (a-h,o-y),logical*1 (z)
      real*8, intent(inout) :: d2alal(metat,nact,nact,nact,nact)
      dimension c(ncf,*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*)
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      logical z1,z2,z3,z4,z5
      integer segno
      integer, allocatable :: iocca(:),ioccb(:)
      real*8 , allocatable :: d2bebe(:,:,:,:,:),d2albe(:,:,:,:,:),cij(:)
      character (len=36) twoRDMfile
! common e dimensionamenti mettere dopo
      allocate(d2bebe(metat,nact,nact,nact,nact),d2albe(metat,nact,nact,nact,nact))
      allocate(cij(metat),iocca(norb),ioccb(norb))
      nocc=noca-ncore
      print*,'nocc,noca',nocc,noca
      call zeroe(d2alal,metat*nact**4)
      call zeroe(d2bebe,metat*nact**4)
      call zeroe(d2albe,metat*nact**4)
      print*,'subr. ro2'
      print*,'nact=',nact,' nocc=',nocc,' ncf=',ncf
      do 1 m=1,ncf
         if (mod(m,100).eq.0)print *,'m=',m
         call flush(6)
!         if(abs(c(m)).lt.thr)goto 1
         ndi=nd(m)
         nec=ne(m)
         call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndi,trou,part,iorb,ispin)
         do 2 n=1,m
!            if(abs(c(n)).lt.thr)goto 2
            do istate=1,metat
               cij(istate)=c(m,istate)*c(n,istate)
            enddo
            if(n.eq.m)goto 3
            ndj=nd(n)
            ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,2,zseg)
            if(ij.eq.0)nbdif=4
            if(nbdif.gt.2)goto 2
               if(zseg)then
                  segno=-1
               else
                  segno=1
               endif
               if(nbdif.eq.1)then
                  n1=iorb(nv(1))
                  ns1=ispin(nv(1))
                  n2=iorb(nw(1))
                  ns2=ispin(nw(1))
                  goto 10
               else
                  call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                  n1=iorb(nv(1))
                  ns1=ispin(nv(1))
                  n2=iorb(nv(2))
                  ns2=ispin(nv(2))
                  n3=iorb(nw(1))
                  ns3=ispin(nw(1))
                  n4=iorb(nw(2))
                  ns4=ispin(nw(2))
!                  print '(a,4i3,l2)','due diff: ',n1,n3,n2,n4,zseg
                  goto 20
               endif
!            goto(10,20),nbdif
 3          continue
!     nessuna differenza
!            ci2=c(m)*c(m)
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do istate=1,metat
                     d2alal(istate,ia,ja,ia,ja)=d2alal(istate,ia,ja,ia,ja)+cij(istate)
                  enddo
                  call permbis(d2alal,ia,ja,ia,ja,nact,metat)
               enddo
            enddo
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do istate=1,metat
                     d2bebe(istate,ib,jb,ib,jb)=d2bebe(istate,ib,jb,ib,jb)+cij(istate)
                  enddo
                  call permbis(d2bebe,ib,jb,ib,jb,nact,metat)
               enddo
            enddo
            do i=1,nocca
               ia=iocca(i)
               do j=1,noccb
                  jb=ioccb(j)
                  do istate=1,metat
                     d2albe(istate,ia,jb,ia,jb)=d2albe(istate,ia,jb,ia,jb)+cij(istate)
                  enddo
               enddo
            enddo
               goto 2
 10            continue
!     una differenza
               n1=n1-ncore
               n2=n2-ncore
!               cij=c(m)*c(n)*segno
               if(ns1.eq.1.and.ns2.eq.1)then
                  do i=1,nocca
                     ia=iocca(i)
                     if(ia.ne.n1.and.ia.ne.n2)then
                        do istate=1,metat
                           cijp=cij(istate)*segno
                           d2alal(istate,n1,ia,n2,ia)=d2alal(istate,n1,ia,n2,ia)+cijp
                        enddo
                        call permbis(d2alal,n1,ia,n2,ia,nact,metat)
                     endif
                  enddo
                  do i=1,noccb
                     ib=ioccb(i)
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        d2albe(istate,n1,ib,n2,ib)=d2albe(istate,n1,ib,n2,ib)+cijp
                        d2albe(istate,n2,ib,n1,ib)=d2albe(istate,n2,ib,n1,ib)+cijp
                     enddo
                  enddo
               else
                  do i=1,noccb
                     ib=ioccb(i)
                     if(ib.ne.n1.and.ib.ne.n2)then
                        do istate=1,metat
                           cijp=cij(istate)*segno
                           d2bebe(istate,n1,ib,n2,ib)=d2bebe(istate,n1,ib,n2,ib)+cijp
                        enddo
                        call permbis(d2bebe,n1,ib,n2,ib,nact,metat)
                     endif
                  enddo
                  do i=1,nocca
                     ia=iocca(i)
                     do istate=1,metat
                        cijp=cij(istate)*segno
                        d2albe(istate,ia,n1,ia,n2)=d2albe(istate,ia,n1,ia,n2)+cijp
                        d2albe(istate,ia,n2,ia,n1)=d2albe(istate,ia,n2,ia,n1)+cijp
                     enddo
                  enddo
               endif
               goto 2
!     due differenze
 20            continue
               n1=n1-ncore
               n2=n2-ncore
               n3=n3-ncore
               n4=n4-ncore
!               cij=c(m)*c(n)*segno
!               print '(a,i1,a,i1,a,i1,a,i1)','Due diff.: ns1=',ns1,
!     $              ' ns2=',ns2,' ns3=',ns3,' ns4=',ns4
               if(ns1.eq.1.and.ns2.eq.1)then
                  do istate=1,metat
                     cijp=cij(istate)*segno
                     d2alal(istate,n1,n2,n3,n4)=d2alal(istate,n1,n2,n3,n4)+cijp
                  enddo
                  call permbis(d2alal,n1,n2,n3,n4,nact,metat)
               elseif(ns1.eq.0.and.ns2.eq.0)then
                  do istate=1,metat
                     cijp=cij(istate)*segno
                     d2bebe(istate,n1,n2,n3,n4)=d2bebe(istate,n1,n2,n3,n4)+cijp
                  enddo
                  call permbis(d2bebe,n1,n2,n3,n4,nact,metat)
               elseif(ns1.eq.1.and.ns2.eq.0)then
                  do istate=1,metat
                     cijp=cij(istate)*segno
                     d2albe(istate,n1,n2,n3,n4)=d2albe(istate,n1,n2,n3,n4)+cijp
                     d2albe(istate,n3,n4,n1,n2)=d2albe(istate,n3,n4,n1,n2)+cijp
                  enddo
               endif
 2          continue
 1       continue
!--costruzione ro2 spinless
         do i=1,nact
            do j=1,nact
               do ip=1,nact
                  do jp=1,nact
                     do istate=1,metat
                        d2alal(istate,i,j,ip,jp)=d2alal(istate,i,j,ip,jp)+d2albe(istate,i,j,ip,jp)&
                                                +d2albe(istate,j,i,jp,ip)+d2bebe(istate,i,j,ip,jp)
                     enddo
                  enddo
               enddo
            enddo
         enddo
#ifdef DEBUG_DMRG_NEVPT
         print *, 'blubb DEBUG rho2'
         do istate=1,metat
         print *, 'rho2 for state ',istate
         do i=1,nact
            do j=1,nact
               do ip=1,nact
                  do jp=1,nact
                        print '(i1,3i2,3x,d19.12)',i,j,ip,jp,d2alal(istate,i,j,ip,jp)
                   enddo
               enddo
            enddo
         enddo
         end do ! loop over states
#endif
         deallocate(d2bebe,d2albe,cij,iocca,ioccb)
         end subroutine ro2
!-------------------------------------------------------
      subroutine ro2off(c,d2alal,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer, intent(in) :: ncoppie
      dimension d2alal(ncoppie,nact,nact,nact,nact)
      dimension c(ncf,*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*)
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      logical z1,z2,z3,z4,z5
      integer segno
      integer, allocatable :: iocca(:),ioccb(:)
      real*8 , allocatable :: cij(:), cji(:)
      real*8 , allocatable :: d2bebe(:,:,:,:,:),d2albe(:,:,:,:,:)

      allocate(d2bebe(ncoppie,nact,nact,nact,nact),d2albe(ncoppie,nact,nact,nact,nact))
      allocate(cij(ncoppie),cji(ncoppie),iocca(norb),ioccb(norb))
      n=nact
      nocc=noca-ncore
      print*,'nocc,noca',nocc,noca
      do i=1,n
         do j=1,n
            do k=1,n
               do l=1,n
                  do ijcou=1,ncoppie
                     d2alal(ijcou,i,j,k,l)=0.d0
                     d2bebe(ijcou,i,j,k,l)=0.d0
                     d2albe(ijcou,i,j,k,l)=0.d0
                  enddo
               enddo
            enddo
         enddo
      enddo
      print*,'subr. ro2'
      print*,'n=',n,' nocc=',nocc,' ncf=',ncf
      do 1 m=1,ncf
         if (mod(m,100).eq.0)print *,'m=',m
         call flush(6)
!         if(abs(c(m)).lt.thr)goto 1
         ndi=nd(m)
         nec=ne(m)
         call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndi,trou,part,iorb,ispin)
         do 2 n=1,m
!            if(abs(c(n)).lt.thr)goto 2
            ijcou=0
            do istate=1,metat
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  ijcou=ijcou+1
                  cij(ijcou)=c(m,jstate)*c(n,istate)
                  cji(ijcou)=c(n,jstate)*c(m,istate)
               enddo
            enddo
            if(n.eq.m)goto 3
            ndj=nd(n)
            ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,2,zseg)
            if(ij.eq.0)nbdif=4
            if(nbdif.gt.2)goto 2
               if(zseg)then
                  segno=-1
               else
                  segno=1
               endif
               if(nbdif.eq.1)then
                  n1=iorb(nv(1))
                  ns1=ispin(nv(1))
                  n2=iorb(nw(1))
                  ns2=ispin(nw(1))
                  goto 10
               else
                  call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                  n1=iorb(nv(1))
                  ns1=ispin(nv(1))
                  n2=iorb(nv(2))
                  ns2=ispin(nv(2))
                  n3=iorb(nw(1))
                  ns3=ispin(nw(1))
                  n4=iorb(nw(2))
                  ns4=ispin(nw(2))
!                  print '(a,4i3,l2)','due diff: ',n1,n3,n2,n4,zseg
                  goto 20
               endif
!            goto(10,20),nbdif
 3          continue
!     nessuna differenza
!            ci2=c(m)*c(m)
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do ijcou=1,ncoppie
                     d2alal(ijcou,ia,ja,ia,ja)=d2alal(ijcou,ia,ja,ia,ja)+cij(ijcou)
                  enddo
                  call permbisoff(d2alal,ia,ja,ia,ja,nact,ncoppie)
               enddo
            enddo
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do ijcou=1,ncoppie
                     d2bebe(ijcou,ib,jb,ib,jb)=d2bebe(ijcou,ib,jb,ib,jb)+cij(ijcou)
                  enddo
                  call permbisoff(d2bebe,ib,jb,ib,jb,nact,ncoppie)
               enddo
            enddo
            do i=1,nocca
               ia=iocca(i)
               do j=1,noccb
                  jb=ioccb(j)
                  do ijcou=1,ncoppie
                     d2albe(ijcou,ia,jb,ia,jb)=d2albe(ijcou,ia,jb,ia,jb)+cij(ijcou)
                  enddo
               enddo
            enddo
               goto 2
 10            continue
!     una differenza
               n1=n1-ncore
               n2=n2-ncore
!               cij=c(m)*c(n)*segno
               if(ns1.eq.1.and.ns2.eq.1)then
                  do i=1,nocca
                     ia=iocca(i)
                     if(ia.ne.n1.and.ia.ne.n2)then
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           d2alal(ijcou,n1,ia,n2,ia)=d2alal(ijcou,n1,ia,n2,ia)+cijp
                           d2alal(ijcou,n2,ia,n1,ia)=d2alal(ijcou,n2,ia,n1,ia)+cjip
                        enddo
                        call permbisoff(d2alal,n1,ia,n2,ia,nact,ncoppie)
                        call permbisoff(d2alal,n2,ia,n1,ia,nact,ncoppie)
                     endif
                  enddo
                  do i=1,noccb
                     ib=ioccb(i)
                     do ijcou=1,ncoppie
                        cijp=cij(ijcou)*segno
                        cjip=cji(ijcou)*segno
                        d2albe(ijcou,n1,ib,n2,ib)=d2albe(ijcou,n1,ib,n2,ib)+cijp
                        d2albe(ijcou,n2,ib,n1,ib)=d2albe(ijcou,n2,ib,n1,ib)+cjip
                     enddo
                  enddo
               else
                  do i=1,noccb
                     ib=ioccb(i)
                     if(ib.ne.n1.and.ib.ne.n2)then
                        do ijcou=1,ncoppie
                           cijp=cij(ijcou)*segno
                           cjip=cji(ijcou)*segno
                           d2bebe(ijcou,n1,ib,n2,ib)=d2bebe(ijcou,n1,ib,n2,ib)+cijp
                           d2bebe(ijcou,n2,ib,n1,ib)=d2bebe(ijcou,n2,ib,n1,ib)+cjip
                        enddo
                        call permbisoff(d2bebe,n1,ib,n2,ib,nact,ncoppie)
                        call permbisoff(d2bebe,n2,ib,n1,ib,nact,ncoppie)
                     endif
                  enddo
                  do i=1,nocca
                     ia=iocca(i)
                     do ijcou=1,ncoppie
                        cijp=cij(ijcou)*segno
                        cjip=cji(ijcou)*segno
                        d2albe(ijcou,ia,n1,ia,n2)=d2albe(ijcou,ia,n1,ia,n2)+cijp
                        d2albe(ijcou,ia,n2,ia,n1)=d2albe(ijcou,ia,n2,ia,n1)+cjip
                     enddo
                  enddo
               endif
               goto 2
!     due differenze
 20            continue
               n1=n1-ncore
               n2=n2-ncore
               n3=n3-ncore
               n4=n4-ncore
!               cij=c(m)*c(n)*segno
!               print '(a,i1,a,i1,a,i1,a,i1)','Due diff.: ns1=',ns1,
!     $              ' ns2=',ns2,' ns3=',ns3,' ns4=',ns4
               if(ns1.eq.1.and.ns2.eq.1)then
                  do ijcou=1,ncoppie
                     cijp=cij(ijcou)*segno
                     cjip=cji(ijcou)*segno
                     d2alal(ijcou,n1,n2,n3,n4)=d2alal(ijcou,n1,n2,n3,n4)+cijp
                     d2alal(ijcou,n3,n4,n1,n2)=d2alal(ijcou,n3,n4,n1,n2)+cjip
                  enddo
                  call permbisoff(d2alal,n1,n2,n3,n4,nact,ncoppie)
                  call permbisoff(d2alal,n3,n4,n1,n2,nact,ncoppie)
               elseif(ns1.eq.0.and.ns2.eq.0)then
                  do ijcou=1,ncoppie
                     cijp=cij(ijcou)*segno
                     cjip=cji(ijcou)*segno
                     d2bebe(ijcou,n1,n2,n3,n4)=d2bebe(ijcou,n1,n2,n3,n4)+cijp
                     d2bebe(ijcou,n3,n4,n1,n2)=d2bebe(ijcou,n3,n4,n1,n2)+cjip
                  enddo
                  call permbisoff(d2bebe,n1,n2,n3,n4,nact,ncoppie)
                  call permbisoff(d2bebe,n3,n4,n1,n2,nact,ncoppie)
               elseif(ns1.eq.1.and.ns2.eq.0)then
                  do ijcou=1,ncoppie
                     cijp=cij(ijcou)*segno
                     cjip=cji(ijcou)*segno
                     d2albe(ijcou,n1,n2,n3,n4)=d2albe(ijcou,n1,n2,n3,n4)+cijp
                     d2albe(ijcou,n3,n4,n1,n2)=d2albe(ijcou,n3,n4,n1,n2)+cjip
                  enddo
               endif
 2          continue
 1       continue
!--costruzione ro2 spinless
         do i=1,nact
            do j=1,nact
               do ip=1,nact
                  do jp=1,nact
                     do ijcou=1,ncoppie
                        d2alal(ijcou,i,j,ip,jp)=d2alal(ijcou,i,j,ip,jp)+d2albe(ijcou,i,j,ip,jp)+d2albe(ijcou,j,i,jp,ip)&
                                               +d2bebe(ijcou,i,j,ip,jp)
                     enddo
                  enddo
               enddo
            enddo
         enddo
#ifdef DEBUG_DMRG_NEVPT
         print *, 'blubb DEBUG trho2'
         do ijcou=1,ncoppie
         print *, 'rho2 for transition ',ijcou
         do i=1,nact
            do j=1,nact
               do ip=1,nact
                  do jp=1,nact
                        print '(i1,3i2,1x,e13.6)',i,j,ip,jp,(-1.0d0)*d2alal(ijcou,i,j,ip,jp)
                   enddo
               enddo
            enddo
         enddo
         end do ! loop over <L|R> state couples
#endif
      deallocate(d2bebe,d2albe,cij,cji,iocca,ioccb)
      end subroutine ro2off
!-----------------------------------------------------
      subroutine ro1(c,d1,ncore,nact,thr,trou,part,ne,nd,ispin,iorb)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d1(metat,nact,nact)
      integer segno
      dimension c(ncf,*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*)
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      real*8 , allocatable :: cij(:)
      allocate(cij(metat))
      n=nact
      nocc=noca-ncore
      nn=(n*n+n)/2
      print*,'subr. ro1'
      print*,'n=',n,' nocc=',nocc,' ncf=',ncf
      print*,'# states = ',metat
      nac=nocc

!     diagonal terms
      do j=1,nac
         do istate=1,metat
            d1(istate,j,j)=2.d0
         enddo
      enddo

      do 1 i=1,ncf
!         if(abs(c(i)).lt.thr)goto 1
         do 2 j=1,i
!            if(abs(c(j)).lt.thr)goto 2
            do istate=1,metat
               cij(istate)=c(i,istate)*c(j,istate)
            enddo
            if(i.eq.j)goto 3
            ij=jdifgen(ne,nd,trou,part,i,j,nv,nw,1,zseg)
            if(ij.eq.0)goto 2
            if(zseg)then
               segno=-1
            else
               segno=1
            endif
            n1=iorb(nv(1))
            ns1=ispin(nv(1))
            n2=iorb(nw(1))
            ns2=ispin(nw(1))
            n2=n2-ncore
            n1=n1-ncore
            do istate=1,metat
               d1(istate,n1,n2)=d1(istate,n1,n2)+cij(istate)*segno
               d1(istate,n2,n1)=d1(istate,n1,n2)
            enddo
            goto 2
 3          nec=ne(i)
            if(nec.eq.0)goto 2
!            cmi2=c(i)*c(i)
            ndi=nd(i)
            do 30 l=1,nec
               lt=trou(l+ndi)
               lp=part(l+ndi)
               lt=iorb(lt)-ncore
               lp=iorb(lp)-ncore
               do istate=1,metat
                  d1(istate,lt,lt)=d1(istate,lt,lt)-cij(istate)
                  d1(istate,lp,lp)=d1(istate,lp,lp)+cij(istate)
               enddo
 30         continue
 2       continue
 1    continue
      deallocate(cij)
      end subroutine ro1
!-----------------------------------------------------
      subroutine ro1off(c,d1,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer, intent(in) :: ncoppie
      dimension d1(ncoppie,nact,nact)
      integer segno
      dimension c(ncf,*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
              dimension ne(*),nd(*),ispin(*),iorb(*)
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      real*8 , allocatable :: cij(:),cji(:)

      allocate(cij(ncoppie),cji(ncoppie))
      n=nact
      nocc=noca-ncore
      nn=(n*n+n)/2
      print*,'subr. ro1off'
      print*,'n=',n,' nocc=',nocc,' ncf=',ncf
      nac=nocc
      do j=1,nac
         do ijcou=1,ncoppie
            d1(ijcou,j,j)=0.d0
         enddo
      enddo
      do 1 i=1,ncf
!         if(abs(c(i)).lt.thr)goto 1
         do 2 j=1,i
!            if(abs(c(j)).lt.thr)goto 2
            ijcou=0
            do istate=1,metat
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  ijcou=ijcou+1
                  cij(ijcou)=c(i,jstate)*c(j,istate)
                  cji(ijcou)=c(j,jstate)*c(i,istate)
               enddo
            enddo
            if(i.eq.j)goto 3
            ij=jdifgen(ne,nd,trou,part,i,j,nv,nw,1,zseg)
            if(ij.eq.0)goto 2
            if(zseg)then
               segno=-1
            else
               segno=1
            endif
            n1=iorb(nv(1))
            ns1=ispin(nv(1))
            n2=iorb(nw(1))
            ns2=ispin(nw(1))
            n2=n2-ncore
            n1=n1-ncore
            do ijcou=1,ncoppie
               d1(ijcou,n1,n2)=d1(ijcou,n1,n2)+cij(ijcou)*segno
               d1(ijcou,n2,n1)=d1(ijcou,n2,n1)+cji(ijcou)*segno
            enddo
            goto 2
 3          nec=ne(i)
            if(nec.eq.0)goto 2
!            cmi2=c(i)*c(i)
            ndi=nd(i)
            do 30 l=1,nec
               lt=trou(l+ndi)
               lp=part(l+ndi)
               lt=iorb(lt)-ncore
               lp=iorb(lp)-ncore
               do ijcou=1,ncoppie
                  d1(ijcou,lt,lt)=d1(ijcou,lt,lt)-cij(ijcou)
                  d1(ijcou,lp,lp)=d1(ijcou,lp,lp)+cij(ijcou)
               enddo
 30         continue
 2       continue
 1    continue
      deallocate(cij,cji)
      end subroutine ro1off
!-------------------------------------------------
      subroutine ord2(nv,nw,zseg,segno)
!--prima gli alfa poi i beta
      logical*1 zseg
      integer segno
      dimension nv(*),nw(*)
      if(nv(2).gt.nv(1))then
         n=nv(1)
         nv(1)=nv(2)
         nv(2)=n
         zseg=.not.zseg
      endif
      if(nw(2).gt.nw(1))then
         n=nw(1)
         nw(1)=nw(2)
         nw(2)=n
         zseg=.not.zseg
      endif
      if(zseg)then
         segno=-1
      else
         segno=1
      endif
      end subroutine ord2
!----------------------------------------------------------------
      subroutine ord3(nv,nw,zseg,segno)
!--prima gli alfa poi i beta
      logical*1 zseg
      integer segno
      dimension nv(*),nw(*)
      if(nv(2).gt.nv(1))then
         n=nv(1)
         nv(1)=nv(2)
         nv(2)=n
         zseg=.not.zseg
      endif
      if(nv(3).gt.nv(1))then
         n=nv(1)
         nv(1)=nv(3)
         nv(3)=n
         zseg=.not.zseg
      endif
      if(nv(3).gt.nv(2))then
         n=nv(2)
         nv(2)=nv(3)
         nv(3)=n
         zseg=.not.zseg
      endif
      if(nw(2).gt.nw(1))then
         n=nw(1)
         nw(1)=nw(2)
         nw(2)=n
         zseg=.not.zseg
      endif
      if(nw(3).gt.nw(1))then
         n=nw(1)
         nw(1)=nw(3)
         nw(3)=n
         zseg=.not.zseg
      endif
      if(nw(3).gt.nw(2))then
         n=nw(2)
         nw(2)=nw(3)
         nw(3)=n
         zseg=.not.zseg
      endif
      if(zseg)then
         segno=-1
      else
         segno=1
      endif
      end subroutine ord3
!----------------------------------------------------------------
      subroutine ord4(nv,nw,zseg,segno)
!--prima gli alfa poi i beta
      logical*1 zseg
      integer segno
      dimension nv(*),nw(*)
      do i=1,4
         imax=nv(i)
         do j=i+1,4
            if(nv(j).gt.imax)then
               n=imax
               imax=nv(j)
               nv(j)=n
               zseg=.not.zseg
            endif
            nv(i)=imax
         enddo
      enddo
      do i=1,4
         imax=nw(i)
         do j=i+1,4
            if(nw(j).gt.imax)then
               n=imax
               imax=nw(j)
               nw(j)=n
               zseg=.not.zseg
            endif
            nw(i)=imax
         enddo
      enddo
      if(zseg)then
         segno=-1
      else
         segno=1
      endif
      end subroutine ord4
!----------------------------------------------------------------
      subroutine giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndm,trou,part,iorb,ispin)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer   iocca, ioccb, jts
      integer*2 trou,part,iorb,ispin
      dimension iocca(*),ioccb(*),trou(*),part(*),iorb(*),ispin(*)
      do i=1,nact
         if(i.le.nocc)then
            iocca(i)=i
            ioccb(i)=i
         else
            iocca(i)=0
            ioccb(i)=0
         endif
      enddo
      do i=1,nec
         it=trou(ndm+i)
         ip=part(ndm+i)
         jts=ispin(it)
         ips=ispin(ip)
         it=iorb(it)-ncore
         ip=iorb(ip)-ncore
         if(jts.eq.1)then
            iocca(it)=0
         else
            ioccb(it)=0
         endif
         if(ips.eq.1)then
            iocca(ip)=ip
         else
            ioccb(ip)=ip
         endif
      enddo
      nocca=0
      noccb=0
      do i=1,nact
         if(iocca(i).ne.0)then
            nocca=nocca+1
            iocca(nocca)=iocca(i)
         endif
         if(ioccb(i).ne.0)then
            noccb=noccb+1
            ioccb(noccb)=ioccb(i)
         endif
      enddo
      end subroutine giveocc2
!--------------------------------------------------------------------
      subroutine perm3(d,i,j,k,l,m,n,nact,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(metat,nact,nact,nact,nact,nact,nact)
!--block 1
      do istate=1,metat
      a=d(istate,i,j,k,l,m,n)
      d(istate,j,i,k,l,m,n)=-a
      d(istate,i,k,j,l,m,n)=-a
      d(istate,k,j,i,l,m,n)=-a
      d(istate,k,i,j,l,m,n)=a
      d(istate,j,k,i,l,m,n)=a
!--block 2
      d(istate,i,j,k,m,l,n)=-a
      d(istate,j,i,k,m,l,n)=a
      d(istate,i,k,j,m,l,n)=a
      d(istate,k,j,i,m,l,n)=a
      d(istate,k,i,j,m,l,n)=-a
      d(istate,j,k,i,m,l,n)=-a
!--block 3
      d(istate,i,j,k,l,n,m)=-a
      d(istate,j,i,k,l,n,m)=a
      d(istate,i,k,j,l,n,m)=a
      d(istate,k,j,i,l,n,m)=a
      d(istate,k,i,j,l,n,m)=-a
      d(istate,j,k,i,l,n,m)=-a
!--block 4
      d(istate,i,j,k,n,m,l)=-a
      d(istate,j,i,k,n,m,l)=a
      d(istate,i,k,j,n,m,l)=a
      d(istate,k,j,i,n,m,l)=a
      d(istate,k,i,j,n,m,l)=-a
      d(istate,j,k,i,n,m,l)=-a
!--block 5
      d(istate,i,j,k,n,l,m)=a
      d(istate,j,i,k,n,l,m)=-a
      d(istate,i,k,j,n,l,m)=-a
      d(istate,k,j,i,n,l,m)=-a
      d(istate,k,i,j,n,l,m)=a
      d(istate,j,k,i,n,l,m)=a
!--block 6
      d(istate,i,j,k,m,n,l)=a
      d(istate,j,i,k,m,n,l)=-a
      d(istate,i,k,j,m,n,l)=-a
      d(istate,k,j,i,m,n,l)=-a
      d(istate,k,i,j,m,n,l)=a
      d(istate,j,k,i,m,n,l)=a
      enddo
      end subroutine perm3
!--------------------------------------------------------------------
      subroutine perm2(d,i,j,k,l,m,n,nact,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(metat,nact,nact,nact,nact,nact,nact)
      do istate=1,metat
!--block 1
      a=d(istate,i,j,k,l,m,n)
      d(istate,j,i,k,l,m,n)=-a
      d(istate,l,m,n,i,j,k)=a
      d(istate,l,m,n,j,i,k)=-a
!--block 2
      d(istate,i,j,k,m,l,n)=-a
      d(istate,j,i,k,m,l,n)=a
      d(istate,m,l,n,i,j,k)=-a
      d(istate,m,l,n,j,i,k)=a
      enddo
      end subroutine perm2
!-----------------------------------------------------------
      subroutine perm3off(d,i,j,k,l,m,n,nact,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(ncoppie,nact,nact,nact,nact,nact,nact)
      do ijcou=1,ncoppie
!--block 1
      a=d(ijcou,i,j,k,l,m,n)
      d(ijcou,j,i,k,l,m,n)=-a
      d(ijcou,i,k,j,l,m,n)=-a
      d(ijcou,k,j,i,l,m,n)=-a
      d(ijcou,k,i,j,l,m,n)=a
      d(ijcou,j,k,i,l,m,n)=a
!--block 2
      d(ijcou,i,j,k,m,l,n)=-a
      d(ijcou,j,i,k,m,l,n)=a
      d(ijcou,i,k,j,m,l,n)=a
      d(ijcou,k,j,i,m,l,n)=a
      d(ijcou,k,i,j,m,l,n)=-a
      d(ijcou,j,k,i,m,l,n)=-a
!--block 3
      d(ijcou,i,j,k,l,n,m)=-a
      d(ijcou,j,i,k,l,n,m)=a
      d(ijcou,i,k,j,l,n,m)=a
      d(ijcou,k,j,i,l,n,m)=a
      d(ijcou,k,i,j,l,n,m)=-a
      d(ijcou,j,k,i,l,n,m)=-a
!--block 4
      d(ijcou,i,j,k,n,m,l)=-a
      d(ijcou,j,i,k,n,m,l)=a
      d(ijcou,i,k,j,n,m,l)=a
      d(ijcou,k,j,i,n,m,l)=a
      d(ijcou,k,i,j,n,m,l)=-a
      d(ijcou,j,k,i,n,m,l)=-a
!--block 5
      d(ijcou,i,j,k,n,l,m)=a
      d(ijcou,j,i,k,n,l,m)=-a
      d(ijcou,i,k,j,n,l,m)=-a
      d(ijcou,k,j,i,n,l,m)=-a
      d(ijcou,k,i,j,n,l,m)=a
      d(ijcou,j,k,i,n,l,m)=a
!--block 6
      d(ijcou,i,j,k,m,n,l)=a
      d(ijcou,j,i,k,m,n,l)=-a
      d(ijcou,i,k,j,m,n,l)=-a
      d(ijcou,k,j,i,m,n,l)=-a
      d(ijcou,k,i,j,m,n,l)=a
      d(ijcou,j,k,i,m,n,l)=a
      enddo
      end subroutine perm3off
!-----------------------------------------------------------
      subroutine perm2off(d,i,j,k,l,m,n,nact,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(ncoppie,nact,nact,nact,nact,nact,nact)
      do ijcou=1,ncoppie
!--block 1
      a=d(ijcou,i,j,k,l,m,n)
      d(ijcou,j,i,k,l,m,n)=-a
!      d(l,m,n,i,j,k)=a
!      d(l,m,n,j,i,k)=-a
!--block 2
      d(ijcou,i,j,k,m,l,n)=-a
      d(ijcou,j,i,k,m,l,n)=a
!      d(m,l,n,i,j,k)=-a
!      d(m,l,n,j,i,k)=a
      enddo
      end subroutine perm2off
!--------------------------------------------------------------------
      subroutine permbis(d,i,j,k,l,nact,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(metat,nact,nact,nact,nact)
      do istate=1,metat
!--block 1
      a=d(istate,i,j,k,l)
      d(istate,j,i,k,l)=-a
      d(istate,k,l,i,j)=a
      d(istate,k,l,j,i)=-a
!--block 2
      d(istate,i,j,l,k)=-a
      d(istate,j,i,l,k)=a
      d(istate,l,k,i,j)=-a
      d(istate,l,k,j,i)=a
      enddo
      end subroutine permbis
!--------------------------------------------------------------------
      subroutine permbisoff(d,i,j,k,l,nact,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(ncoppie,nact,nact,nact,nact)
      do ijcou=1,ncoppie
!--block 1
      a=d(ijcou,i,j,k,l)
      d(ijcou,j,i,k,l)=-a
!      d(ijcou,k,l,i,j)=a
!      d(ijcou,k,l,j,i)=-a
!--block 2
      d(ijcou,i,j,l,k)=-a
      d(ijcou,j,i,l,k)=a
!      d(ijcou,l,k,i,j)=-a
!      d(ijcou,l,k,j,i)=a
      enddo
      end subroutine permbisoff
!--------------------------------------------------------------------
      subroutine permut123(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
!---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      end subroutine permut123
!----------------------------------------------------------
      subroutine permut124(i,j,l,k,ip,jp,lp,kp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
!---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      end subroutine permut124
!--------------------------------------------------------------------
      subroutine permut134(i,l,j,k,ip,lp,jp,kp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
!---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      end subroutine permut134
!--------------------------------------------------------------------
      subroutine permut234(l,i,j,k,lp,ip,jp,kp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
!---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      end subroutine permut234
!--------------------------------------------------------------------
      subroutine permut1234(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
!--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=i
         kp=l
         lp=k
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.3)then
         ip=j
         jp=i
         kp=k
         lp=l
         segno=-1
      elseif(ivolte.eq.4)then
         ip=i
         jp=j
         kp=l
         lp=k
         segno=-1
      endif
      end subroutine permut1234
!--------------------------------------------------------------------
      subroutine permut1324(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
!--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
!         ip=k
!         jp=i
!         kp=l
!         lp=j
         ip=k
         jp=l
         kp=i
         lp=j
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.3)then
!         ip=k
!         jp=i
!         kp=j
!         lp=l
         ip=k
         jp=j
         kp=i
         lp=l
         segno=-1
      elseif(ivolte.eq.4)then
!         ip=i
!         jp=k
!         kp=l
!         lp=j
         ip=i
         jp=l
         kp=k
         lp=j
         segno=-1
      endif
      end subroutine permut1324
!--------------------------------------------------------------------
      subroutine permut1423(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
!--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
!         ip=l
!         jp=i
!         kp=k
!         lp=j
         ip=l
         jp=k
         kp=j
         lp=i
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.3)then
!         ip=l
!         jp=i
!         kp=j
!         lp=k
         ip=l
         jp=j
         kp=k
         lp=i
         segno=-1
      elseif(ivolte.eq.4)then
!         ip=i
!         jp=l
!         kp=k
!         lp=j
         ip=i
         jp=k
         kp=j
         lp=l
         segno=-1
      endif
      end subroutine permut1423
!----------------------------------------------------------------------
      subroutine permut4(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      if(ivolte.eq.1)then
!--permutazioni pari
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
         ip=k
         jp=i
         kp=j
         lp=l
         segno=1
      elseif(ivolte.eq.3)then
         ip=l
         jp=i
         kp=k
         lp=j
         segno=1
      elseif(ivolte.eq.4)then
         ip=j
         jp=k
         kp=i
         lp=l
         segno=1
      elseif(ivolte.eq.5)then
         ip=j
         jp=l
         kp=k
         lp=i
         segno=1
      elseif(ivolte.eq.6)then
         ip=j
         jp=i
         kp=l
         lp=k
         segno=1
      elseif(ivolte.eq.7)then
         ip=l
         jp=j
         kp=i
         lp=k
         segno=1
      elseif(ivolte.eq.8)then
         ip=k
         jp=l
         kp=i
         lp=j
         segno=1
      elseif(ivolte.eq.9)then
         ip=k
         jp=j
         kp=l
         lp=i
         segno=1
      elseif(ivolte.eq.10)then
         ip=l
         jp=k
         kp=j
         lp=i
         segno=1
      elseif(ivolte.eq.11)then
         ip=i
         jp=l
         kp=j
         lp=k
         segno=1
      elseif(ivolte.eq.12)then
         ip=i
         jp=k
         kp=l
         lp=j
         segno=1
!---permutazioni dispari
      elseif(ivolte.eq.13)then
         ip=j
         jp=i
         kp=k
         lp=l
         segno=-1
      elseif(ivolte.eq.14)then
         ip=k
         jp=j
         kp=i
         lp=l
         segno=-1
      elseif(ivolte.eq.15)then
         ip=l
         jp=j
         kp=k
         lp=i
         segno=-1
      elseif(ivolte.eq.16)then
         ip=i
         jp=k
         kp=j
         lp=l
         segno=-1
      elseif(ivolte.eq.17)then
         ip=i
         jp=l
         kp=k
         lp=j
         segno=-1
      elseif(ivolte.eq.18)then
         ip=i
         jp=j
         kp=l
         lp=k
         segno=-1
      elseif(ivolte.eq.19)then
         ip=l
         jp=i
         kp=j
         lp=k
         segno=-1
      elseif(ivolte.eq.20)then
         ip=k
         jp=l
         kp=j
         lp=i
         segno=-1
      elseif(ivolte.eq.21)then
         ip=k
         jp=i
         kp=l
         lp=j
         segno=-1
      elseif(ivolte.eq.22)then
         ip=l
         jp=k
         kp=i
         lp=j
         segno=-1
      elseif(ivolte.eq.23)then
         ip=j
         jp=l
         kp=i
         lp=k
         segno=-1
      elseif(ivolte.eq.24)then
         ip=j
         jp=k
         kp=l
         lp=i
         segno=-1
      endif
      end subroutine permut4
!---------------------------------------------------------------------
      integer FUNCTION Jdifgen(ne,nd,trou,part,II,JJ,nv1,nv2,nma,zsig)
      implicit real*8(a-h,o-y),logical*1(z)
!
!     SSP DE RECHERCHE DES ELEMENTS D'INTERACTION
!     ENTRE DEUX DETERMINANTS, NUMEROTES II ET JJ
!
      integer, parameter :: nmax=4
      dimension ne(*),nd(*),trou(*),part(*)
      INTEGER*2 NE,TROU,PART
      DIMENSION ITA(99),ITB(99),NV1(nmax),NV2(nmax),NATU1(nmax),NATU2(nmax)
!      EQUIVALENCE(NV11,NV1(1)),(NV12,NV1(2)),(NV21,NV2(1)),(NV22,NV2(2))
      common/diff/n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      KR=II
      KL=JJ
!
!     LE DETERMINANT - 1 - EST CHOISI COMME LE
!     DETERMINANT AYANT LE PLUS D'O.M. EXCITEES
!
      zshift=.false.
      IF (NE(KR)-NE(KL)) 5,7,7
    5 KKKK=KR
      KR=KL
      KL=KKKK
      zshift=.true.
    7 NE1=NE(KR)
!     TEST SUR LE FONDAMENTAL
      IF (NE1.EQ.0) GOTO 44
      NE2=NE(KL)
!     NBDIF : EST LE NOMBRE DE SPIN-ORBITALES
!     DIFFERENTES
      NBDIF=NE1-NE2
      IF (NBDIF.GT.nma) GOTO 24
!     CONSTRUCTION POUR LE DETERMINANT -1- ET -2-
!     DE IT ET NS, QUI CONTIENNENT LES NUMEROS
!     DES ORBITALES ET LEURS SPINS
!
      ND1=ND(KR)
      ND2=ND(KL)
      DO 15 J=1,NE1
         ITA(J)=TROU(J+ND1)
         ITA(J+NE1)=PART(J+ND1)
 15   CONTINUE
!
!     ZSIG :  EST LA PARITE DU  NOMBRE DE CROISEMENTS
!     DANS LE DIAGRAMME D INTERACTION
      ZSIG=.FALSE.
!     SI NE2=0 LE DETERMINANT -2- EST LE FONDAMENTAL
      IF (NE2.EQ.0) GOTO 125
 121  CONTINUE
      DO 25 J=1,NE2
         ITB(J)=TROU(J+ND2)
 25   ITB(J+NE2)=PART(J+ND2)
    4 DO 8 K=1,2
         K1=(K-1)*NE1
         K2=(K-1)*NE2
         DO 10 I=1,NE2
            NIT1=ITB(K2+I)
            J1=1+K1
            J2=NE1+K1
            DO 12 J=J1,J2
               NJT2=ITA(J)
               IF (NIT1.NE.NJT2) GOTO 12
               IF (I.EQ.(J-K1)) GOTO 10
               ITA(J)=ITA(I+K1)
               ITA(I+K1)=NJT2
               ZSIG=.NOT.ZSIG
               GOTO 10
 12         CONTINUE
            NBDIF=NBDIF+1
            IF (NBDIF.GT.nma) GOTO 24
            NV1(NBDIF)=NIT1
            NATU1(NBDIF)=K
            NATU2(NBDIF)=I
 10         CONTINUE
 8       CONTINUE
 125     IF (NBDIF.LE.0) GOTO 44
!
!     NAR : EXCITATION SUPLEMENTAIRE DU DETERMINANT -1-
!     PAR RAPPORT AU DETERMINANT -2-
 26      NAR=NE1-NE2
         IF (NAR.LE.0) GOTO 28
 30      DO 32 I=1,NAR
            NET=NE1+1-I
            NV1(I)=ITA(NET)
            NV2(I)=ITA(NET+NE1)
 32      CONTINUE
 28      NBAR=NBDIF-NAR
         IF (NBAR.LE.0) GOTO 34
 36      NAR1=NAR+1
         DO 38 I=NAR1,NBDIF
            K=NATU1(I)
            NI=NATU2(I)
            NU=NV1(I)
            NV2(I)=ITA(NI+NE1*(K-1))
            IF (K.EQ.2) GOTO 38
 40         NV1(I)=NV2(I)
            NV2(I)=NU
            ZSIG=.NOT.ZSIG
 38      CONTINUE
 34      CONTINUE
!     CALCUL DE L ELEMENT DE MATRICE
 2       continue
         if(zshift)then
            do i=1,nbdif
               ndum=nv1(i)
               nv1(i)=nv2(i)
               nv2(i)=ndum
            enddo
         endif
         jdifgen=1
         RETURN
 24      Jdifgen=0
         RETURN
 44      WRITE (6,1006) II,JJ,NE1,NE2
 1006    FORMAT (//5X,'IDENTITE DANS HNTD'/' II,JJ,NE1,NE2 =',4I4)
         WRITE (6,1007) (ITA(K),K=1,NE1)
         WRITE (6,1008) (ITA(K+NE1),K=1,NE1)
         WRITE (6,1007) (ITB(K),K=1,NE2)
         WRITE (6,1008) (ITB(K+NE2),K=1,NE2)
 1007    FORMAT (' TROU',8I4)
 1008    FORMAT (' PART',8I4)
         STOP
         END FUNCTION Jdifgen
!------------------------------------------------------------------
      subroutine accum(amat,ip,jp,kp,lp,mp,np,op,pp,val,segno,metat,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension amat(nwords,nact,nact,nact,nact,metat),val(*)
      integer op,pp,segno
      ind=lindice(ip,jp,kp,lp)
      do istate=1,metat
         amat(ind,mp,np,op,pp,istate)=amat(ind,mp,np,op,pp,istate)+val(istate)*segno
      enddo
!      write (6,'(''accum'',8I3,I6,2f15.8)')
!     * ip,jp,kp,lp,mp,np,op,pp,ind,val(1),amat(ind,mp,np,op,pp)
      end subroutine accum
!------------------------------------------------------------------

      subroutine get_aj_atwo(atwo,nact,ncore,norb,aj)

      use ijkl_utils
      implicit none

      real*8, intent(out) :: aj(*)
      real*8, intent(out) :: atwo(nact,nact,nact,nact)
      integer, intent(in) :: ncore, nact, norb
      integer             :: a,b,ac,bc, k, l, jab
      integer             :: indice,i,j

      !> statement function
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      call dcopy(norb*(norb+1)/2,nevpt_ijkl%oneint,1,aj,1)

!----building Dyall's h eff in aj----------
      do a=1,nact
         ac=a+ncore
         do b=a,nact
            bc=b+ncore
            jab=indice(ac,bc)
            do j=1,ncore
               aj(jab)=aj(jab)+2.d0*ai(ac,bc,j,j)- &
                                    ai(ac,j,j,bc)
            enddo
         enddo
      enddo
!---building matrix atwo
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  atwo(i,j,k,l)=ai(i+ncore,j+ncore,k+ncore,l+ncore)
               enddo
            enddo
         enddo
      enddo

#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb DEBUG atwo'
      do i=1,nact
        do j=1,nact
          do k=1,nact
            do l=1,nact
              if(abs(atwo(i,j,k,l)) > 1.0d-16)&
       print '(i1,3i2,3x,d19.12)',&
       i,j,k,l,atwo(i,j,k,l)
            enddo
          enddo
        enddo
      enddo

      call flush(6)

      print *, 'blubb DEBUG aj'
      do i=1,nact
        ac=i+ncore
        do j=i,nact
            bc=j+ncore
            jab=indice(ac,bc)
              if(abs(aj(jab)) > 1.0d-16)&
       print '(i1,2i2,3x,d19.12)',&
       ac,bc,jab,aj(jab)
        enddo
      enddo

      call flush(6)
#endif


      end subroutine get_aj_atwo
!-------------------------------------------------------------------

      subroutine bamat(amat,atmat,atwo,ro4,daaa,taa,dal,nact,ncore,f,aj,norb,metat,nsym)
      implicit real*8 (a-h,o-y),logical*1 (z)
!---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(metat,nact,nact,nact,nact,nact,nact),taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      dimension ro4(nwords,nact,nact,nact,nact,metat)
      dimension amat(metat,nact,nact,nact,nact,nact,nact)
      dimension atmat(metat,nact,nact,nact,nact,nact,nact)
      dimension atwo(nact,nact,nact,nact)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)

      integer nsym
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)


!--building heff'
!--in AJ is heff
      do a=1,nact
         ac=a+ncore
         do b=a,nact
            bc=b+ncore
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nact
               cc=c+ncore
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum*0.5d0
            enddo
         enddo
      enddo
!-------------------
      call flush(6)
      if (nsym.eq.1) then
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,c,ap,bp,cp,d,istate,ie,dint1,dint2,dint3,dint4,dum1,dum2, &
!$OMP&  dum3,dum4,dum1t,dum2t,dum3t,dum4t,ief) collapse(7)
        do a=1,nact
          do b=1,nact
              do c=1,nact
                do ap=1,nact
                    do bp=1,nact
                      do cp=1,nact
                          do d=1,nact
                            do istate=1,metat
!$OMP ATOMIC
                                amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap   &
                                    ,bp,cp,a,b,c)+f(indice(d+ncore,a+ncore))*eee2(cp,ap,bp,b     &
                                    ,d,c,daaa,taa,dal,nact,istate,metat      &
                                    )-f(indice(c+ncore,d+ncore))*eee2(cp,ap,bp,b,a,d             &
                                    ,daaa,taa,dal,nact,istate,metat          &
                                    )-f(indice(b+ncore,d+ncore))*eee2(cp,ap,bp,d,a,c             &
                                    ,daaa,taa,dal,nact,istate,metat          &
                                    )
!$OMP ATOMIC
                                atmat(istate,ap,bp,cp,a,b,c)=atmat(istate    &
                                    ,ap,bp,cp,a,b,c)+f(indice(d+ncore,a+ncore))*eee2t(cp,ap      &
                                    ,bp,b,d,c,daaa,taa,dal,nact,istate       &
                                    ,metat)-f(indice(c+ncore,d+ncore))*eee2t(cp,ap,bp            &
                                    ,b,a,d,daaa,taa,dal,nact,istate,metat    &
                                    )+f(indice(b+ncore,d+ncore))*eee2t(cp,ap,bp,d,a,c            &
                                    ,daaa,taa,dal,nact,istate,metat          &
                                    )
                            enddo
                            do ie=1,nact
                                dint1=atwo(c,ie,d,a)
                                dint2=atwo(d,ie,ie,a)
                                dint3=atwo(c,d,d,ie)
                                dint4=atwo(b,d,d,ie)
                                do istate=1,metat
                                  dum1=eee2(cp,ap,bp,b,d,ie,daaa,taa,dal,nact,istate,metat)
                                  dum2=eee2(cp,ap,bp,b,d,c,daaa,taa,dal,nact,istate,metat)
                                  dum3=eee2(cp,ap,bp,b,a,ie,daaa,taa,dal,nact,istate,metat)
                                  dum4=eee2(cp,ap,bp,ie,a,c,daaa,taa,dal,nact,istate,metat)
                                  dum1t=-dum1
                                  dum2t=-dum2
                                  dum3t=-dum3
                                  if(bp.eq.b)then
                                      dum1t=dum1t+2.d0*ee2(cp,ap,d,ie,taa,dal,nact,istate,metat)
                                      dum2t=dum2t+2.d0*ee2(cp,ap,d,c,taa,dal,nact,istate,metat)
                                      dum3t=dum3t+2.d0*ee2(cp,ap,a,ie,taa,dal,nact,istate,metat)
                                  endif
!$OMP ATOMIC
                                  atmat(istate,ap,b,cp,a,bp,c)        &
                                      =atmat(istate,ap,b,cp,a,bp,c)   &
                                      +dum1t*dint1-0.5d0*(dum2t*dint2 &
                                      +dum3t*dint3)
                                  dum4t=eee2t(cp,ap,bp,ie,a,c,daaa,taa,dal,nact,istate,metat)
!$OMP ATOMIC
                                  amat(istate,ap,bp,cp,a,b,c)=amat(istate  &
                                      ,ap,bp,cp,a,b,c)+dum1*dint1-0.5d0    &
                                      *(dum2*dint2+dum3*dint3-dum4*dint4)
!$OMP ATOMIC
                                  atmat(istate,ap,bp,cp,a,b,c)=atmat(istate,ap,bp,cp,a,b,c)+0.5d0*dum4t*dint4
                                enddo
                                do ief=1,nact
                                  do istate=1,metat
                                    dint1=atwo(d,ief,ie,a)
                                    dum1=eeee(cp,ap,bp,b,d,ief,ie,c,daaa,taa,dal,ro4,nact,istate,metat)
                                    dum1t=-dum1
                                    if(bp.eq.b)dum1t=dum1t+2.d0*eee2(cp,ap,d,ief,ie,c,daaa,taa,dal,nact,istate,metat)
!$OMP ATOMIC
                                    amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap,bp,cp,a,b,c)+dum1*dint1
!$OMP ATOMIC
                                    atmat(istate,ap,b,cp,a,bp,c)=atmat(istate,ap,b,cp,a,bp,c)+dum1t*dint1
                                    dum2=eeee(cp,ap,bp,b,d,ief,a,ie,daaa,taa,dal,ro4,nact,istate,metat)
                                    dum2t=-dum2
                                    if(bp.eq.b)dum2t=dum2t+2.d0*eee2(cp,ap,d,ief,a,ie,daaa,taa,dal,nact,istate,metat)
                                    dint2=atwo(d,ief,c,ie)
!$OMP ATOMIC
                                    amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap,bp,cp,a,b,c)-dum2*dint2
!$OMP ATOMIC
                                    atmat(istate,ap,b,cp,a,bp,c)=atmat(istate,ap,b,cp,a,bp,c)-dum2t*dint2
                                    dum3=eeee(cp,ap,bp,ie,d,ief,a,c,daaa,taa,dal,ro4,nact,istate,metat)
                                    dum3t=eeeet(cp,ap,bp,ie,d,ief,a,c,daaa,taa,dal,ro4,nact,istate,metat)
                                    dint3=atwo(d,ief,b,ie)
!$OMP ATOMIC
                                    amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap,bp,cp,a,b,c)-dum3*dint3
!$OMP ATOMIC
                                    atmat(istate,ap,bp,cp,a,b,c)=atmat(istate,ap,bp,cp,a,b,c)+dum3t*dint3
                                  enddo
                                enddo
                            enddo
                          enddo
                      enddo
                    enddo
                enddo
              enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      else
        do a=1,nact
          ac=a+ncore
          isa=itsym(ac)
          do b=1,nact
              bc=b+ncore
              isb=itsym(bc)
              isab=its(isa,isb)
              do c=1,nact
                cc=c+ncore
                isc=itsym(cc)
                isabc=its(isab,isc)
                isbc=its(isb,isc)
                isac=its(isa,isc)
                do ap=1,nact
                    apc=ap+ncore
                    isap=itsym(apc)
                    do bp=1,nact
                      bpc=bp+ncore
                      isbp=itsym(bpc)
                      isabp=its(isap,isbp)
                      do cp=1,nact
                          cpc=cp+ncore
                          iscp=itsym(cpc)
                          isabcp=its(isabp,iscp)
                          isatot=its(isabc,isabcp)
                          if (isatot.ne.1) then
                            cycle
                          endif
                          do d=1,nact
                            dc=d+ncore
                            isd=itsym(dc)
                            ida=indice(dc,ac)
                            icd=indice(cc,dc)
                            ibd=indice(bc,dc)
                            do istate=1,metat
                                amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap   &
                                    ,bp,cp,a,b,c)+f(ida)*eee2(cp,ap,bp,b     &
                                    ,d,c,daaa,taa,dal,nact,istate,metat      &
                                    )-f(icd)*eee2(cp,ap,bp,b,a,d             &
                                    ,daaa,taa,dal,nact,istate,metat          &
                                    )-f(ibd)*eee2(cp,ap,bp,d,a,c             &
                                    ,daaa,taa,dal,nact,istate,metat          &
                                    )
                                atmat(istate,ap,bp,cp,a,b,c)=atmat(istate    &
                                    ,ap,bp,cp,a,b,c)+f(ida)*eee2t(cp,ap      &
                                    ,bp,b,d,c,daaa,taa,dal,nact,istate       &
                                    ,metat)-f(icd)*eee2t(cp,ap,bp            &
                                    ,b,a,d,daaa,taa,dal,nact,istate,metat    &
                                    )+f(ibd)*eee2t(cp,ap,bp,d,a,c            &
                                    ,daaa,taa,dal,nact,istate,metat          &
                                    )
                            enddo
                            do ie=1,nact
                                iec=ie+ncore
                                ise=itsym(iec)
                                dint1=atwo(c,ie,d,a)
                                dint2=atwo(d,ie,ie,a)
                                dint3=atwo(c,d,d,ie)
                                dint4=atwo(b,d,d,ie)
                                do istate=1,metat
                                  dum1=eee2(cp,ap,bp,b,d,ie,daaa,taa,dal,nact,istate,metat)
                                  dum2=eee2(cp,ap,bp,b,d,c,daaa,taa,dal,nact,istate,metat)
                                  dum3=eee2(cp,ap,bp,b,a,ie,daaa,taa,dal,nact,istate,metat)
                                  dum4=eee2(cp,ap,bp,ie,a,c,daaa,taa,dal,nact,istate,metat)
                                  dum1t=-dum1
                                  dum2t=-dum2
                                  dum3t=-dum3
                                  if(bp.eq.b)then
                                      dum1t=dum1t+2.d0*ee2(cp,ap,d,ie,taa,dal,nact,istate,metat)
                                      dum2t=dum2t+2.d0*ee2(cp,ap,d,c,taa,dal,nact,istate,metat)
                                      dum3t=dum3t+2.d0*ee2(cp,ap,a,ie,taa,dal,nact,istate,metat)
                                  endif
                                  atmat(istate,ap,b,cp,a,bp,c)        &
                                      =atmat(istate,ap,b,cp,a,bp,c)   &
                                      +dum1t*dint1-0.5d0*(dum2t*dint2 &
                                      +dum3t*dint3)
                                  dum4t=eee2t(cp,ap,bp,ie,a,c,daaa,taa,dal,nact,istate,metat)
                                  amat(istate,ap,bp,cp,a,b,c)=amat(istate  &
                                      ,ap,bp,cp,a,b,c)+dum1*dint1-0.5d0    &
                                      *(dum2*dint2+dum3*dint3-dum4*dint4)
                                  atmat(istate,ap,bp,cp,a,b,c)=atmat(istate,ap,bp,cp,a,b,c)+0.5d0*dum4t*dint4
                                enddo
                                do ief=1,nact
                                  iefc=ief+ncore
                                  isf=itsym(iefc)
                                  isfd=its(isf,isd)
                                  isfde=its(isfd,ise)
                                  istot=its(isabcp,isfde)
                                  do istate=1,metat
                                      if (its(isfde,isa).eq.1.and.its(istot,isbc).eq.1) then
                                        dint1=atwo(d,ief,ie,a)
                                        dum1=eeee(cp,ap,bp,b,d,ief,ie,c,daaa,taa,dal,ro4,nact,istate,metat)
                                        dum1t=-dum1
                                        if(bp.eq.b)dum1t=dum1t+2.d0*eee2(cp,ap,d,ief,ie,c,daaa,taa,dal,nact,istate,metat)
                                        amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap,bp,cp,a,b,c)+dum1*dint1
                                        atmat(istate,ap,b,cp,a,bp,c)=atmat(istate,ap,b,cp,a,bp,c)+dum1t*dint1
                                      endif
                                      if (its(isfde,isc).eq.1.and.its(istot,isab).eq.1) then
                                        dum2=eeee(cp,ap,bp,b,d,ief,a,ie,daaa,taa,dal,ro4,nact,istate,metat)
                                        dum2t=-dum2
                                        if(bp.eq.b)dum2t=dum2t+2.d0*eee2(cp,ap,d,ief,a,ie,daaa,taa,dal,nact,istate,metat)
                                        dint2=atwo(d,ief,c,ie)
                                        amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap,bp,cp,a,b,c)-dum2*dint2
                                        atmat(istate,ap,b,cp,a,bp,c)=atmat(istate,ap,b,cp,a,bp,c)-dum2t*dint2
                                      endif
                                      if (its(isfde,isb).eq.1.and.its(istot,isac).eq.1) then
                                        dum3=eeee(cp,ap,bp,ie,d,ief,a,c,daaa,taa,dal,ro4,nact,istate,metat)
                                        dum3t=eeeet(cp,ap,bp,ie,d,ief,a,c,daaa,taa,dal,ro4,nact,istate,metat)
                                        dint3=atwo(d,ief,b,ie)
                                        amat(istate,ap,bp,cp,a,b,c)=amat(istate,ap,bp,cp,a,b,c)-dum3*dint3
                                        atmat(istate,ap,bp,cp,a,b,c)=atmat(istate,ap,bp,cp,a,b,c)+dum3t*dint3
                                      endif
                                  enddo
                                enddo
                            enddo
                          enddo
                      enddo
                    enddo
                enddo
              enddo
          enddo
        enddo
      end if
      end subroutine bamat
!-------------------------------------------------------------------
      subroutine bbmat(atwo,bmat,btmat,daaa,taa,dal,nact,ncore,f,aj,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
!---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(metat,nact,nact,nact,nact,nact,nact),taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      dimension bmat(metat,nact,nact,nact,nact)
      dimension btmat(metat,nact,nact,nact,nact)
      dimension atwo(nact,nact,nact,nact)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      integer :: indice,i,j
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
!---building matrix heff''
      do a=1,nact
         ac=a+ncore
         do b=a,nact
            bc=b+ncore
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nact
               cc=c+ncore
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum
            enddo
         enddo
      enddo
!-------------------------------------
      do ap=1,nact
         apc=ap+ncore
         do bp=1,nact
            bpc=bp+ncore
            do cp=1,nact
               cpc=cp+ncore
               do a=1,nact
                  ac=a+ncore
                  do istate=1,metat
                     bmat(istate,ap,bp,cp,a)=0.d0
                     btmat(istate,ap,bp,cp,a)=0.d0
                  enddo
                  do c=1,nact
                     cc=c+ncore
                     iac=indice(cc,ac)
                     do istate=1,metat
                        bmat(istate,ap,bp,cp,a)=bmat(istate,ap,bp,cp,a) &
                             -f(iac)*ee2(cp,ap,bp,c,taa,dal,nact,istate &
                             ,metat)
                        btmat(istate,ap,bp,cp,a)=btmat(istate,ap,bp,cp,a&
                             )+aj(iac)*ee2tr(cp,ap,bp,c,taa,dal,nact    &
                             ,istate,metat)
                     enddo
                     do ie=1,nact
                        iec=ie+ncore
                        do ief=1,nact
                           iefc=ief+ncore
                           dint=atwo(a,ie,c,ief)
                           do istate=1,metat
                              dum=eee2(cp,ap,bp,ie,c,ief,daaa,taa,dal,nact,istate,metat)
                              dumt=eee2t(cp,ap,bp,ie,c,ief,daaa,taa,dal,nact,istate,metat)
                              bmat(istate,ap,bp,cp,a)=bmat(istate,ap,bp,cp,a)-dum*dint
                              btmat(istate,ap,bp,cp,a)=btmat(istate,ap,bp,cp,a)+dumt*dint
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end subroutine bbmat
!-------------------------------------------------------------------
      subroutine bcmat(atwo,cmat,ctmat,daaa,taa,dal,nact,ncore,f,aj,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
!---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(metat,nact,nact,nact,nact,nact,nact),taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      dimension cmat(metat,nact,nact,nact,nact)
      dimension ctmat(metat,nact,nact,nact,nact)
      dimension atwo(nact,nact,nact,nact)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      integer :: indice,i,j
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      call flush(6)
!--building heff'
      do a=1,nact
         do b=a,nact
            jab=indice(a+ncore,b+ncore)
            f(jab)=aj(jab)
            do c=1,nact
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum*0.5d0
            enddo
         enddo
      enddo
!-------------------
      cmat=0.0d0
      ctmat=0.0d0
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ap,a,b,c,d,ida,icd,ibd,istate,ie,dint1,dint2,dint3,dint4,     &
!$OMP&  dum1,dum2,dum3,dum4,dum1t,dum2t,dum3t,dum4t,ief) collapse(5)
      do ap=1,nact
         do a=1,nact
            do b=1,nact
               do c=1,nact
                  do d=1,nact
                     ida=indice(d+ncore,a+ncore)
                     icd=indice(c+ncore,d+ncore)
                     ibd=indice(b+ncore,d+ncore)
                     do istate=1,metat
!$OMP ATOMIC
                        cmat(istate,ap,a,b,c)=cmat(istate,ap,a,b,c)   &
                             +f(ida)*ee2(ap,b,d,c,taa,dal,nact,istate &
                             ,metat)-f(icd)*ee2(ap,b,a,d,taa,dal      &
                             ,nact,istate,metat)-f(ibd)*ee2(ap,d      &
                             ,a,c,taa,dal,nact,istate,metat)
!$OMP ATOMIC
                        ctmat(istate,ap,a,b,c)=ctmat(istate,ap,a,b,c) &
                             +f(ida)*ee2t(ap,b,d,c,taa,dal,nact,istate&
                             ,metat)-f(icd)*ee2t(ap,b,a,d,taa         &
                             ,dal,nact,istate,metat)+f(ibd)           &
                             *ee2t(ap,d,a,c,taa,dal,nact,istate,metat)
                     enddo
                     do ie=1,nact
                        dint1=atwo(c,ie,d,a)
                        dint2=atwo(d,ie,ie,a)
                        dint3=atwo(c,d,d,ie)
                        dint4=atwo(b,d,d,ie)
                        do istate=1,metat
                           dum1=ee2(ap,b,d,ie,taa,dal,nact,istate,metat)
                           dum2=ee2(ap,b,d,c,taa,dal,nact,istate,metat)
                           dum3=ee2(ap,b,a,ie,taa,dal,nact,istate,metat)
                           dum4=ee2(ap,ie,a,c,taa,dal,nact,istate,metat)
                           dum1t=ee2t(ap,b,d,ie,taa,dal,nact,istate,metat)
                           dum2t=ee2t(ap,b,d,c,taa,dal,nact,istate,metat)
                           dum3t=ee2t(ap,b,a,ie,taa,dal,nact,istate,metat)
                           dum4t=ee2t(ap,ie,a,c,taa,dal,nact,istate,metat)
!$OMP ATOMIC
                           cmat(istate,ap,a,b,c)=cmat(istate,ap,a,b,c)+dum1*dint1-0.5d0&
                                                *(dum2*dint2+dum3*dint3-dum4*dint4)
!$OMP ATOMIC
                           ctmat(istate,ap,a,b,c)=ctmat(istate,ap,a,b,c)+dum1t*dint1-0.5d0&
                                                *(dum2t*dint2+dum3t*dint3-dum4t*dint4)
                        enddo
                        do ief=1,nact
                           dint1=atwo(d,ief,ie,a)
                           dint2=atwo(d,ief,c,ie)
                           dint3=atwo(d,ief,b,ie)
                           do istate=1,metat
                              dum1=eee2(ap,b,d,ief,ie,c,daaa,taa,dal,nact,istate,metat)
                              dum1t=eee2tl(ap,b,d,ief,ie,c,daaa,taa,dal,nact,istate,metat)
                              dum2=eee2(ap,b,d,ief,a,ie,daaa,taa,dal,nact,istate,metat)
                              dum2t=eee2tl(ap,b,d,ief,a,ie,daaa,taa,dal,nact,istate,metat)
                              dum3=eee2(ap,ie,d,ief,a,c,daaa,taa,dal,nact,istate,metat)
                              dum3t=eee2tl(ap,ie,d,ief,a,c,daaa,taa,dal,nact,istate,metat)
!$OMP ATOMIC
                              cmat(istate,ap,a,b,c)=cmat(istate,ap,a,b,c)+dum1*dint1-dum2*dint2-dum3*dint3
!$OMP ATOMIC
                              ctmat(istate,ap,a,b,c)=ctmat(istate,ap,a,b,c)+dum1t*dint1-dum2t*dint2+dum3t*dint3
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO
      end subroutine bcmat
!-------------------------------------------------------------------
      subroutine bdmat(atwo,dmat,dtmat,taa,dal,nact,ncore,f,aj,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
!---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      dimension dmat(metat,nact,nact)
      dimension dtmat(metat,nact,nact)
      dimension atwo(nact,nact,nact,nact)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      integer :: indice,i,j
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
!---building matrix heff''
      do a=1,nact
         ac=a+ncore
         do b=a,nact
            bc=b+ncore
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nact
               cc=c+ncore
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum
            enddo
         enddo
      enddo
!---------------------------
      do ap=1,nact
         apc=ap+ncore
         do a=1,nact
            ac=a+ncore
            do istate=1,metat
               dmat(istate,ap,a)=0.d0
               dtmat(istate,ap,a)=0.d0
            enddo
            do c=1,nact
               cc=c+ncore
               iac=indice(cc,ac)
               do istate=1,metat
                  dmat(istate,ap,a)=dmat(istate,ap,a)-f(iac)*dal(istate,ap,c)
               dum=-dal(istate,ap,c)
               if(ap.eq.c)dum=dum+2.d0
               dtmat(istate,ap,a)=dtmat(istate,ap,a)+aj(iac)*dum
            enddo
               do ie=1,nact
                  iec=ie+ncore
                  do ief=1,nact
                     iefc=ief+ncore
                     dint=atwo(a,ie,c,ief)
                     do istate=1,metat
                        dum=ee2(ap,ie,c,ief,taa,dal,nact,istate,metat)
                        dumt=ee2t(ap,ie,c,ief,taa,dal,nact,istate,metat)
                        dmat(istate,ap,a)=dmat(istate,ap,a)-dum*dint
                        dtmat(istate,ap,a)=dtmat(istate,ap,a)+dumt*dint
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end subroutine bdmat
!-------------------------------------------------------------------
      subroutine zeroe(a,n)
      real*8, intent(inout) :: a(n)
      a = 0.d0
      end subroutine zeroe
!------------------------------------------------------
      subroutine ord31(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nact,m,n)
      integer segno,segno2,segno3,op,pp,op2,pp2,m,n,mp2,np2,mp,np
      integer ia,ic,ie,ig,ia2,ic2,ie2,ig2,nact,ivolte,jvolte
      integer i1,i2,i3,i4,i5,i6,i7,i8,norm1,norm2
      real*8 amat(nwords,nact,nact,nact,nact,metat),val(*)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      logical*1 ZQDPT,ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG

      norm1=ind_norm%normind(ia,ic,ie,ig)
      norm2=ind_norm%normind(ia2,ic2,ie2,ig2)
      if (norm1.ge.norm2) then
      i1=ia
      i2=ic
      i3=ie
      i4=ig
      i5=ia2
      i6=ic2
      i7=ie2
      i8=ig2
      else
      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig
      endif

      do ivolte=1,6  !aaab
      call permut123(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
!      write (6,'(''ord31 1'',8I3)') mp,np,op,pp
        do jvolte=1,6
        call permut123(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo
      do ivolte=1,6  !aaba
      call permut124(i1,i2,i4,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
!      write (6,'(''ord31 2'',8I3)') mp,np,op,pp
        do jvolte=1,6
        call permut124(i5,i6,i8,i7,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !abaa
      call permut134(i1,i4,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
!      write (6,'(''ord31 3'',8I3)') mp,np,op,pp
        do jvolte=1,6
        call permut134(i5,i8,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !baaa
      call permut234(i4,i1,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6
        call permut234(i8,i5,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      if (norm1.ne.norm2) return
      if (m.eq.n) return

      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig

      do ivolte=1,6  !aaab
      call permut123(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6
        call permut123(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !aaba
      call permut124(i1,i2,i4,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6
        call permut124(i5,i6,i8,i7,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !abaa
      call permut134(i1,i4,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6
        call permut134(i5,i8,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !baaa
      call permut234(i4,i1,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6
        call permut234(i8,i5,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      end subroutine ord31
!-----------------------------------------------------------
      subroutine ord22(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nact,m,n)
      real*8 amat(nwords,nact,nact,nact,nact,metat),val(*)
      integer segno,segno2,segno3,op,pp,op2,pp2,m,n

      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      logical*1 ZQDPT,ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
!     write (6,'(''ord22 in'',8I3)') ia,ic,ie,ig,ia2,ic2,ie2,ig2
      norm1=ind_norm%normind(ia,ic,ie,ig)
      norm2=ind_norm%normind(ia2,ic2,ie2,ig2)
!     write(6,*) 'norme',norm1,norm2
      if (norm1.ge.norm2) then
      i1=ia
      i2=ic
      i3=ie
      i4=ig
      i5=ia2
      i6=ic2
      i7=ie2
      i8=ig2
      else
      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig
      endif
!      write (6,'(''ord22 in'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8

      do ivolte=1,4  !aabb
      call permut1234(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 1'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abab
      call permut1324(i1,i3,i2,i4,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 2'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i5,i7,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abba
      call permut1423(i1,i3,i4,i2,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 3'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i5,i7,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baab
      call permut1423(i3,i1,i2,i4,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 4'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i7,i5,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baba
      call permut1324(i3,i1,i4,i2,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 5'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i7,i5,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !bbaa
      call permut1234(i3,i4,i1,i2,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 6'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i7,i8,i5,i6,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      if (norm1.ne.norm2) return
      if (m.eq.n) return

      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig

!      write (6,'(''ord22 in'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8

      do ivolte=1,4  !aabb
      call permut1234(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 1'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abab
      call permut1324(i1,i3,i2,i4,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 2'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i5,i7,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abba
      call permut1423(i1,i3,i4,i2,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 3'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i5,i7,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baab
      call permut1423(i3,i1,i2,i4,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 4'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i7,i5,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baba
      call permut1324(i3,i1,i4,i2,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 5'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i7,i5,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !bbaa
      call permut1234(i3,i4,i1,i2,mp,np,op,pp,segno,ivolte)
!      write (6,'(''ord22 6'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i7,i8,i5,i6,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      end subroutine ord22
!------------------------------------------------------
      subroutine ordsame(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nact,m,n)
      COMMON /CIP/ IGELS(2),NCFG,NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,&
      ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      logical*1 ZQDPT,ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      real*8 amat(nwords,nact,nact,nact,nact,metat),val(*)
      integer segno,segno2,segno3,op,pp,op2,pp2,n,m
      norm1=ind_norm%normind(ia,ic,ie,ig)
      norm2=ind_norm%normind(ia2,ic2,ie2,ig2)
      if (norm1.ge.norm2) then
      i1=ia
      i2=ic
      i3=ie
      i4=ig
      i5=ia2
      i6=ic2
      i7=ie2
      i8=ig2
      else
      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig
      endif

      do ivolte=1,24
      call permut4(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,24
        call permut4(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo

      if (norm1.ne.norm2) return
      if (m.eq.n) return

      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig

      do ivolte=1,24
      call permut4(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,24
        call permut4(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
!        vals=val*segno*segno2
        segno3=segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,val,segno3,metat,nact)
        enddo
      endif
      enddo


      end subroutine ordsame
!-----------------------------------------------------

      subroutine bro3(ro4,daaa,nact,nele,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daaa(metat,nact,nact,nact,nact,nact,nact)
      dimension ro4(nwords,nact,nact,nact,nact,metat)
!!---controllo su contrazione della ro4
#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb ultra DEBUG ind'
#endif
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do ii=1,nact
                  do jj=1,nact
                     do kk=1,nact
                        do ll=1,nact
                           call ord8(i,j,k,ll,ii,jj,kk,ll,i1,i2,i3,i4,i5,i6,i7,i8)
                           ind=lindice(i1,i2,i3,i4)
                           do istate=1,metat
                              daaa(istate,i,j,k,ii,jj,kk)=daaa(istate,i,j,k,ii,jj,kk)+ro4(ind,i5,i6,i7,i8,istate)
#ifdef DEBUG_DMRG_NEVPT
                             if(k == kk .and. istate == 1 .and. abs(daaa(istate,i,j,k,ii,jj,kk)) > 1.0d-16)then
                               print '(i1,5i2,3x,d19.12,a,i5,2x,d19.12)',&
                               i,j,k,ii,jj,kk,daaa(istate,i,j,k,ii,jj,kk), ' from ind r4 ',ind,ro4(ind,i5,i6,i7,i8,istate)
                             end if
#endif
                           enddo
                        enddo
                        do istate=1,metat
                           daaa(istate,i,j,k,ii,jj,kk)=daaa(istate,i,j,k,ii,jj,kk)/(nele-3)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb DEBUG rho4'
      do istate=1,metat
      print *, 'rho4 for state ',istate
      do ip=1,nwords
        do i=1,nact
          do j=1,nact
            do k=1,nact
              do l=1,nact
                  if(abs(ro4(ip,i,j,k,l,istate)) > 1.0d-16)&
       print '(i3,4i2,3x,d19.12)',&
       ip,i,j,k,l,ro4(ip,i,j,k,l,istate)
              enddo
            enddo
          enddo
        enddo
      enddo

      enddo ! istate

      print *, 'blubb DEBUG rho3'
      do istate=1,metat
      print *, 'rho3 for state ',istate
      do i=1,nact
        do j=1,nact
          do k=1,nact
            do ip=1,nact
              do jp=1,nact
                do kp=1,nact
                  if(abs(daaa(istate,i,j,k,ip,jp,kp)) > 1.0d-16)&
       print '(i1,5i2,3x,d19.12)',&
       i,j,k,ip,jp,kp,daaa(istate,i,j,k,ip,jp,kp)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      enddo ! istate
      call flush(6)
#endif
      end subroutine bro3
!-----------------------------------------------------------------
      subroutine bro2(daaa,taa,nact,nele,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daaa(metat,nact,nact,nact,nact,nact,nact)
      dimension taa(metat,nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            do ii=1,nact
               do jj=1,nact
                  do ll=1,nact
                     do istate=1,metat
                        taa(istate,i,j,ii,jj)=taa(istate,i,j,ii,jj)+daaa(istate,i,j,ll,ii,jj,ll)
                     enddo
                  enddo
                  do istate=1,metat
                     taa(istate,i,j,ii,jj)=taa(istate,i,j,ii,jj)/(nele-2)
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb DEBUG rho2'
      do istate=1,metat
      print *, 'rho2 for state ',istate
      do i=1,nact
        do j=1,nact
          do ip=1,nact
            do jp=1,nact
              if(abs(taa(istate,i,j,ip,jp)) > 1.0d-16)&
       print '(i1,3i2,3x,d19.12)',&
       i,j,ip,jp,taa(istate,i,j,ip,jp)
            enddo
          enddo
        enddo
      enddo

      enddo ! istate
      call flush(6)
#endif
      end subroutine bro2
!-----------------------------------------------------------------
      subroutine bro1(taa,dal,nact,nele,metat)
#ifdef _OPENMP_
      use omp_lib
#endif
      implicit none
      real*8, intent(in)  :: taa(metat,nact,nact,nact,nact)
      real*8, intent(out) :: dal(metat,nact,nact)
      integer, intent(in) :: nact, nele, metat
!     scratch
      integer             :: i, j, ll, istate, nthreads, id

!$omp parallel do shared ( taa, dal, nact, nele, metat )
      do i=1,nact
#ifdef _OPENMP_
!     nthreads = omp_get_num_threads()
!     id       = omp_get_thread_num()
!     print *, 'number of threads = ',nthreads
!     print *, 'my ID is          = ',id
#endif
         do j=1,nact
            do ll=1,nact
               do istate=1,metat
                  dal(istate,i,j)=dal(istate,i,j)+taa(istate,i,ll,j,ll)
               enddo
            enddo
            do istate=1,metat
               dal(istate,i,j)=dal(istate,i,j)/(nele-1)
            enddo
         enddo
      enddo
!$omp end parallel do

#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb DEBUG rho1'
      do istate=1,metat
      print *, 'rho1 for state ',istate
      do i=1,nact
        do j=1,nact
          if(abs(dal(istate,i,j)) > 1.0d-16)&
          print '(i1,1i2,3x,d19.12)',&
          i,j,dal(istate,i,j)
        enddo
      enddo

      enddo ! istate
      call flush(6)
#endif
      end subroutine bro1
!------------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine fdiff(bmat,cmat,btmat,ctmat,nact,metat)
      implicit none

      real*8 , intent(in) :: bmat(metat,nact,nact,nact,nact),cmat(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: btmat(metat,nact,nact,nact,nact),ctmat(metat,nact,nact,nact,nact)
      integer, intent(in) :: nact, metat

      real*8 , parameter  :: threshold = 1.0d-7
      integer             :: istate, i, j, ii, jj
      integer             :: index_ldiff(4)
      real*8              :: bigdiff, val, val2, pdiff

      index_ldiff(1:4) = 0
      do istate=1,metat
         bigdiff=0.d0
         do i =1,nact
            do j=1,nact
               do ii=1,nact
                  do jj=1,nact
                     val=bmat(istate,i,j,ii,jj)
                     val2=cmat(istate,jj,i,j,ii)
                     pdiff=abs(val-val2)
                     if(pdiff.gt.bigdiff)then
                       bigdiff=pdiff
                       index_ldiff(1) = i; index_ldiff(2) = j; index_ldiff(3) = ii; index_ldiff(4) = jj
                     end if
                  enddo
               enddo
            enddo
         enddo
         if(abs(bigdiff) > threshold)then
           write(6,'( a,f14.8,a,i3)') ' maximum difference between bmat and cmat is ',bigdiff, ' for state ',istate
           write(6,'(a,4i4/)') ' matrix indices of bmat(i,j,ii,jj)/cmat(jj,i,j,ii): ',index_ldiff(1:4)
         end if
      enddo

      index_ldiff(1:4) = 0
      do istate=1,metat
         bigdiff=0.d0
         do i =1,nact
            do j=1,nact
               do ii=1,nact
                  do jj=1,nact
                     val=btmat(istate,i,j,ii,jj)
                     val2=ctmat(istate,jj,i,j,ii)
                     pdiff=abs(val-val2)
                     if(pdiff.gt.bigdiff)then
                       bigdiff=pdiff
                       index_ldiff(1) = i; index_ldiff(2) = j; index_ldiff(3) = ii; index_ldiff(4) = jj
                     end if
                  enddo
               enddo
            enddo
         enddo
         if(abs(bigdiff) > threshold)then
           write(6,'( a,f14.8,a,i3)') ' maximum difference between btmat and ctmat is ',bigdiff, ' for state ',istate
           write(6,'(a,4i4/)') ' matrix indices of bmat(i,j,ii,jj)/cmat(jj,i,j,ii): ',index_ldiff(1:4)
         end if
      enddo
      end subroutine fdiff
!------------------------------------------------------------------------------

      subroutine giveocc(iocc,m,nd,ne,trou,part,nocc,norb)
      integer*2 ne(*),trou(*),part(*)
      integer   iocc(*)
      integer nd(*),m,nocc,norb
      do i=1,nocc
         iocc(i)=2
      enddo
      do i=nocc+1,norb
         iocc(i)=0
      enddo
      nstart=nd(m)
!     print *, 'giveocc: m            == ',m
!     print *, 'giveocc: nd(m), ne(m) == ',nd(m), ne(m)
!     print *, 'norb, nocc             == ',norb, nocc
      do i=1,ne(m)
         ib=trou(nstart+i)
         ip=part(nstart+i)
         if(ib.gt.norb)ib=ib-norb
         if(ip.gt.norb)ip=ip-norb
!        print *, 'ib (t), ip (p) of i== ',ib,ip
         iocc(ib)=iocc(ib)-1
         if (ip.le.norb) then
         iocc(ip)=iocc(ip)+1
         endif
      enddo
!     print *, 'giveocc result iocc of det m   == ',iocc(1:norb)
      end subroutine giveocc
!------------------------------------------------------------------
      integer function ndiff(iocc,jocc,norb)
      integer   iocc(*),jocc(*)
      ndiff=0
      do i=1,norb
        ndiff= ndiff + abs(iocc(i)- jocc(i))
      enddo
      end function ndiff
!--------------------------------------------------------
      subroutine esclass(nd,ne,trou,part,icomp,ncf,nocc,norb,nspin)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      integer*2 ne,trou,part
      integer nd,nocc,nspin, icomp
      dimension nd(*),ne(*),trou(*),part(*),icomp(*)
      allocatable iocc(:)
      integer   iocc
      integer :: norb, ncf
      allocate(iocc(norb))
      isz2=nspin-1
      ifirst=1
      ncapos=0
      do idet=1,ncf
         if(idet.eq.ifirst)then
            call giveocc(iocc,idet,nd,ne,trou,part,nocc,norb)
            nsing=0
            do i=1,norb
               if(iocc(i).eq.1)nsing=nsing+1
            enddo
            nha=(nsing-isz2)/2
            nde=noverk(nsing,nha)
            ilast=ifirst+nde-1
            icomp(idet)=nde
            ncapos=ncapos+1
         endif
         if(idet.gt.ifirst.and.idet.le.ilast)then
            icomp(idet)=icomp(ifirst)
            goto 999
         endif
 999     if(idet.eq.ilast)ifirst=ilast+1
      enddo
      print*,'esclass found ',ncapos,' capostipiti'
      deallocate(iocc)
      end subroutine esclass
!------------------------------------------------------------------------
      integer function noverk(n,k)
      if(k.eq.0)then
         noverk=1
      elseif(k.eq.1)then
         noverk=n
      else
         num=n
         iden=1
         do i=2,k
            num=num*(n-i+1)
            iden=iden*i
         enddo
         noverk=num/iden
      endif
      end function noverk
!-------------------------------------------------------------------
      subroutine printtp(ne,nd,trou,part,ncf,norb,c)
      integer*2 ne(*),trou(*),part(*)
      integer nd(*),ncf,norb
      real*8 c(*)
      do i=1,ncf
         ndi=nd(i)
         print '(f10.5,a,$)',c(i),' '
         do j=1,ne(i)
            if(trou(ndi+j).gt.norb)then
               print '(i2,a,$)',trou(ndi+j)-norb,'+'
            else
               print '(i2,a,$)',trou(ndi+j),'-'
            endif
            print '(a,$)',' '
         enddo
         print '(a,$)','==> '
         do j=1,ne(i)
            if(part(ndi+j).gt.norb)then
               print '(i2,a,$)',part(ndi+j)-norb,'+'
            else
               print '(i2,a,$)',part(ndi+j),'-'
            endif
            print '(a,$)',' '
         enddo
         print '(a)',' '
      enddo
      end subroutine printtp

end module koopro4QD
