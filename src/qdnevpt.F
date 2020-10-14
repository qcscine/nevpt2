      module qdnevpt_core

      use nevpt_header
      use hdf5_utils
      use ijkl_utils
      use info_orbital_space
      use info_symmetry
      use nevpt2_cfg ! input parameters
      use info_state_energy ! energies and effective Hamiltonian

      public qdnevpt

      contains

      subroutine qdnevpt(rel_ham,t13,from_molcas)
C-----------------------------------------------------------------------
C
C     This program computes the second order correction to the energy
C     using the partition of the Hamiltonian defined in the
C     ``n-electron valence state perturbation theory'' with the
C     strongly contracted (SC) and partially contracted (PC)
C     versions and using the formalism of the Quasi-Degenerate
C     Perturbation Theory (QD-PT).
C
C-----------------------------------------------------------------------
C
C     INPUT DATA
C
C     NAMELIST FILES
C     ZVERBOSE T intermediate results are printed
C     igelo  = orbital indices of the core (doubly occupied) MO that will be
C              kept frozen in the perturbation expansion.
C
C-----------------------------------------------------------------------


      implicit none

      logical, intent(in) :: rel_ham
      logical, intent(in) :: from_molcas

      !> local variable to control the print level
      integer             :: iprint

      logical*1, pointer :: zgel
      dimension zgel(:)
      logical*1 :: ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD
      logical*1 :: ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     *     NORB,NORB2,NOCB,NOCA,
     *     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     *     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     *     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      real*8, allocatable :: emp(:)
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      logical*1 :: zhandy,zthree, zheff
      COMMON /HND/ ZHANDY,ZTHREE

      common/debug/zdebug
      COMMON /CPU/ TCPU,PCPU
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM(10)
      common/detinf/ntotdet,ntotconf,ntotcap
      common/thresho/thresh
      logical*1          :: zverbose
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*8, allocatable :: fockc(:,:,:),fockv(:,:,:)
      INTEGER*2 NEP,TROUP,PARTP
      integer*1, allocatable :: iocc(:,:)
      real*8, allocatable :: wr(:),wi(:),vr(:,:),vl(:,:),work2(:)
      integer :: idum
      integer :: i, j, ir, irv, ist, is, it, irt, isv, info
      integer :: inita, initb, isym, isz, itrovo, indi, ik, idum2
      integer :: ks, l, kmol, knev, jt, locm
      integer :: k, locmem, memsiz, metat, mmetat
      integer :: nd, ncoppie, ncmax, ncf, ncfg
      integer :: nnrot, nijkl, ngel,m,loffs
      integer :: ntot, nsym, nrot, nref, nr4, norb, nocb, noca, noc
      integer :: jstate, nvmax, nvirt, num, ntrsy, ntotdet, ntotconf
      integer :: ntotcap, norb2, ichange, indice
      real*8  :: stepdb, pcpu, fmpb, tcpu, tau, escf
      real*8  :: brdmp, wtamp, thresh, test, tdegen,een, escf2
      real*8  :: c, convdb, davg, delta, sum_a, sum_c, sum_v
      real*8  :: e2enp, psienep, e2mpp, psimpp
      logical*1 :: zquadru, zdebug
      real*4, intent(out) :: t13
      real*4              :: t1
      integer             :: offset_c, offset_a, offset_v
      real*8, allocatable :: rhotmp(:,:)
      real*8              :: total_sum, total_loc
      real*8              :: nevpt_time_start,nevpt_time_end

C
C     Generalized Koopmans matrices
C
C     gkp1a = generalized Koopmans matrix for alpha electronic
C     affinities
C     gkm1a = generalized Koopmans matrix for alpha ionizzation
C     potentials
C     gkp2aa= generalized Koopmans matrix for alpha double EA
C     ......
C     gk0pa = A matrix used for the 0' class (see report)
C     gk0pd = D matrix used for the 0' class (see report)
C     gk0pf = F matrix used for the 0' class (see report)
C
      real*8 , allocatable :: gkp1a(:,:,:),gkm1a(:,:,:)
      real*8 , allocatable :: gkp2aa(:,:,:,:,:)
      real*8 , allocatable :: gkm2aa(:,:,:,:,:)
      real*8 , allocatable :: gk0pa(:,:,:,:,:),gk0pd(:,:,:,:,:),
     &                        gk0pf(:,:,:)
      real*8 , allocatable :: ro3(:,:,:,:,:,:,:),amat(:,:,:,:,:,:,:),
     &                        bmat(:,:,:,:,:),cmat(:,:,:,:,:),
     &                        dmat(:,:,:)
      real*8 , allocatable :: atmat(:,:,:,:,:,:,:),btmat(:,:,:,:,:),
     *                        ctmat(:,:,:,:,:),dtmat(:,:,:)
      real*8 , allocatable :: rho(:,:,:),w(:),work(:)

Cele QD_new
      real*8 , allocatable :: epsnew(:,:),cvirt(:,:,:),ccore(:,:,:)
      integer, allocatable :: nord_nev2mol(:),nord_mol2nev(:)
      integer, allocatable :: nord_mn_c(:,:),nord_mn_v(:,:)
      real*8 , allocatable :: uc(:,:,:,:),uv(:,:,:,:)
      integer, allocatable :: nmo(:),nc(:),na(:),nv(:),nbe(:)
Cele QD_new
      integer, allocatable :: nord_mn_a(:,:)

      real*8,  allocatable :: vr_tmp(:,:)

C
C     Density matrices SPIN LESS
C
C     da  = one-particle density matrix
C     daa = two-particle density matrix
C     daat= two-hole density matrix
C
      real*8, allocatable :: da(:,:,:)
      real*8, allocatable :: daoff(:,:,:)
      real*8, allocatable :: daa(:,:,:,:,:)
      real*8, allocatable :: daaoff(:,:,:,:,:)
      real*8, allocatable :: ro3off(:,:,:,:,:,:,:)
      real*8, allocatable :: daat(:,:,:,:,:)

C Arrays for the reduced Cholesky vectors, needed for speeding up the calc'n of the Fock matrices with Cholesky
C and also the half-transformed Cholesky vectors for newint_xxxx transforms
      real*8, allocatable :: cho_reduced_J_core(:)
      real*8, allocatable :: cho_reduced_J_act(:,:)
C Arrays for half-transformed newint Cholesky vectors to be initialised as (nchovec,nlvec,nstates)
c TODO: FOR NOW, the transformed vectors for every state are stored at once, which might be quite memory hungry if one wishes to calculate a lot of states
c One should check if the transformation can be performed one state at once, but this can be done at a later stage, once the algorithm works
      real*8, allocatable :: cho_reduced_vec(:,:,:)

C BLAS ddot()
      real*8, external    :: ddot

      integer nos
      DIMENSION NOS(20)
      integer :: ndhb,ndha
      logical*1 :: zdh
      common /deephole/ ndhb,ndha,zdh

      integer :: istate, ii, ib, ia, ias, ibs
      real*8  :: factor

      integer :: muld2h

      !> statement function
      indice(i,j)= max(i,j)*(max(i,j)-1)/2+min(i,j)
      muld2h(i,j)= IEOR(i-1,j-1)+1

      call cpu_time(t1)

      !> minimal print level (TODO: make it adjustable by MOLCAS)
      iprint = print_level

      ntotdet=0
      ntotconf=0
      ntotcap=0
      zverbose=.false.
      zquadru = .false.

      !> make sure the array is allocated for a minimal input
      if(.not.from_molcas)then
        if(.not.allocated(MultGroup%State))
     &  allocate(MultGroup%State(nr_states))
        do istate = 1, nr_states
          MultGroup%State(istate) = istate
        end do
      end if

      !> single-state calculation is always state-specific
      if(nr_states == 1) skip_effective_ham = .true.

      metat = nr_states

      !> initialize HDF5 support and open the files nevpt.h5 and ijkl.h5
      call hdf5_init()
      call hdf5_open(filename, file_id(1))
      call hdf5_open(ijklname, file_id(2))

      if(.not. do_dmrg_pt)then
        call ofiles()
      end if

      !> get information from MOTRA interface: nsym, norb and more...
      norb = 0; nsym = 0
      call initialize_inforb(nsym,norb,from_molcas)

      !> set its and itsym
      call initialize_infsym(nsym,norb)
      ! Leon 29/6/2017: Do not read ref. energies if we're called from molcas
      ! The energies are read from JobIph and provided before the call to
      ! qdnevpt()
      if (.not. from_molcas) then
        if(.not. do_dmrg_pt)then

          read (4) NORB,NOCB,NOCA,METAT,NCF,II,I
          rewind 4
          allocate(ne(ncf))
          allocate(nd(ncf))
          allocate(trou(ii))
          allocate(part(ii))
          allocate(emp(metat))

! moved to info_state_energy module!
!         allocate(e(metat))

          read (4) NORB,NOCB,NOCA,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),
     *       (TROU(J),PART(J),J=1,II),(C,J=1,I),(E(J),J=1,METAT),
     *       (EMP(J),J=1,METAT),(FMPB,J=1,norb*2),
     *       NNROT,(brdmp,J=1,NNROT),MMETAT,(ZHEFF,J=1,MMETAT),
     *       zhandy,zthree,FACTOR,escf

          deallocate(emp)
          deallocate(part)
          deallocate(trou)
          deallocate(nd)
          deallocate(ne)

      !print *, 'initial factor is ...',factor
      !print *, 'escf (core energy) is ... ',escf

Cele 17-01-07
C Attenzione: e(i) sono storate moltiplicate per factor (scritto nel file04)
C ma nel seguito si ricalcola factor con bigint. Questo puo` generare problemi
C (per esempio se lo stesso file04 e` usato con canonizzazioni diverse degli OM
C doppi occupati).
Cele (vedi anche correzione dopo bigint)
          do J=1,metat
          E(J)=E(J)/factor
          enddo
         !print *, ' state energies are ... (original)  ',e(1:metat)

         ! Leon 29/6/2017: This bit is confusing:
         ! apparently e has been read and in the next line it's being read again
         ! from somewhere else? O_o
          call get_state_energies(e,metat,MultGroup%State)
        else
! moved to info_state_energy module!
!         allocate(e(metat)); e = 0
          call get_state_energies(e,metat,MultGroup%State)
        end if
      end if

      thresh=thr_lindep     !threshold for linear dependencies (input parameter)

      !> allocations for pertubation part
      allocate(e2enp(8,metat),psienep(8,metat))
      allocate(e2mpp(8,metat),psimpp(8,metat))
      e2enp = 0; psienep = 0; e2mpp = 0; psimpp = 0

      if(.not. from_molcas)then
        write(6,'(/a)')'--------------------------------------------'
        if(.not. do_dmrg_pt)then
          write(6,'(a)') 'CASSCF   zeroth-order wavefunction'
        else
          write(6,'(a)') 'DMRG-SCF zeroth-order wavefunction'
        endif
        write(6,'(a/)')'--------------------------------------------'
      end if

      if (no_4rdm_terms) then
        write (6,*) 'WARNING: .NO4RDM has been specified. V(+1p) and '
        write (6,*) 'V(-1p) terms will not be calculated.'
      end if

      ncore = inforb%nisht
      nact  = inforb%nasht
      nvirt = norb-ncore-nact

      allocate(nord_mn_c(ncore,nsym),nord_mn_v(nvirt,nsym))
      allocate(nord_mn_a(nact,nsym))
      allocate(nord_nev2mol(norb),nord_mol2nev(norb))
      allocate(nv(nsym),nc(nsym),na(nsym),nbe(nsym),nmo(nsym))

      call set_orbarray_pointers(
     &                           nsym,
     &                           norb,
     &                           ncore,
     &                           nact,
     &                           nvirt,
     &                           itsym,
     &                           nord_mn_c,
     &                           nord_mn_a,
     &                           nord_mn_v,
     &                           nord_nev2mol,
     &                           nord_mol2nev,
     &                           nc,
     &                           na,
     &                           nv,
     &                           nbe,
     &                           nmo,
     &                           ncmax,
     &                           nvmax
     &                          )

      !> integral reading
      call readint(norb,nsym,nord_mol2nev,escf2)

#ifdef IJKL_DEBUG
      do i = 1, NORB*(NORB+1)/2
        print *, 'my nevpt_ijkl%oneint is ...',i,nevpt_ijkl%oneint(i)
      end do
      if (.not.Do_Cholesky) then
        do i = 1, nevpt_ijkl%ni4
        print *, 'my nevpt_ijkl%twoint(i) is ...',i,nevpt_ijkl%twoint(i)
        end do
      end if
#endif

      !> done with basic initialization
      e2mpp  = 0.d0
      psimpp = 0.d0

      call finalize_inforb()


!-----------------------------------------------------------------------
!
!     Building of the core and virtual Fock matrices for each state and symmetry
!     and computing of the canonical core orbitals for each state
!     writen in the basis of the MO on which the bielectronic integrals
!     have been computed.
!
!-----------------------------------------------------------------------
      allocate(da(metat,nact,nact))
      da = 0
      !> get data from file nevpt.h5
      datadim(1)   = metat; datadim(2:3) = nact; datadim_bound = 3

      call hdf5_get_data_dp(file_id(1), "1-RDM"  , datadim, da)

      allocate (uc(ncmax,ncmax,nsym,metat))
      allocate (uv(nvmax,nvmax,nsym,metat))

      allocate (epsnew(norb,metat))
      epsnew = 0

      if (Do_Cholesky) then
        ! initialise arrays needed for calculations with Cholesky decomposition
        if (ncore > 0) then
          allocate(cho_reduced_J_core(nchovec))
          cho_reduced_J_core = 0.0d0
        end if
        allocate(cho_reduced_J_act(nchovec,metat))
        cho_reduced_J_act = 0.0d0
      end if

      do isym=1,nsym
      if (nc(isym).eq.0) cycle
      allocate(w(nc(isym)))
      allocate(work(3*nc(isym)))

      allocate (fockc(nc(isym),nc(isym),metat))

C     Core-Fock operator

! 2016-03-13 Leon: split the Fock operator to state-specific and state-independent part for less computations
      if(iprint > 10) call cpu_time(nevpt_time_start)
!     print *, 'ncore is ...',ncore
! Here comes the state-independent part, saved to the Fock matrix of the ground state

! fill one-electron integrals
      do i=1,nc(isym)
      it=nord_mn_c(i,isym)
       do j=1,i
       jt=nord_mn_c(j,isym)
        indi=indice(it,jt)
        fockc(i,j,1)=nevpt_ijkl%oneint(indi)
       end do
      end do

      !> fill two-electron integrals
      if (Do_Cholesky) then
        ! Simplify the calculation of the Coulomb part of the Fock operator:
        ! F_ij = (i|h|j) = sum_i sum_j sum_k 2(ij|kk) - (ik|kj)
        !                                        ^-- Coulomb part
        ! Since (ij|kk) = sum_J L_ij^J L_kk^J (with L as Cholesky vectors)
        ! sum_k (ij|kk) can be calculated from
        ! sum_J L_ij^J (sum_k L_kk^J), i.e. the sum of L_kk Cholesky vectors
        ! can be stored and used to multiply with L_ij to produce the sum of Coulomb integrals at once
        !
        ! This doesn't work with exchange integrals, unfortunately, so these are added to the Fock operator as before
        ! Update 25-3-2016: TODO: this should also work partially with exchange integrals, implement it if it's worth the computational cost saving (probably it's not)
        do k=1,ncore
          ! element kk in the cholesky array should be at k(k+1)/2
          cho_reduced_J_core=cho_reduced_J_core+
     &        nevpt_ijkl%cholesky_array(:,k*(k+1)/2)
        end do

        !cho_reduced_J_core now contains sum over all diagonal core-core Cholesky vectors

          do i=1,nc(1)
            do j=1,i
              ! Coulomb contribution from the dot product of the sum of Cholesky vectors L_kk
              fockc(i,j,1)=fockc(i,j,1)+2.d0*
     &         ddot(nchovec,nevpt_ijkl%cholesky_array(1,i*(i-1)/2+j),
     &         1,cho_reduced_J_core(1),1)
              ! Exchange contribution
               do k=1,ncore
                fockc(i,j,1)=fockc(i,j,1)-ai(i,k,k,j)
               end do
            end do
          end do
      else
      ! old code
        do i=1,nc(isym)
        it=nord_mn_c(i,isym)
          do j=1,i
          jt=nord_mn_c(j,isym)
            do k=1,ncore
!       print *, 'ai for it,jt, k,k is ...',ai(it,jt,k,k),ai(it,k,k,jt)
              fockc(i,j,1)=fockc(i,j,1)+
     &         2.d0*ai(it,jt,k,k)-ai(it,k,k,jt)
            enddo
          enddo
        enddo
      end if


! And here comes the state-specific part
      do istate=1,metat
      ! copy the state-independent part to the next state if we're not at the last state
      if (istate.lt.metat) then
       fockc(:,:,istate+1) = fockc(:,:,istate)
      end if

      ! Do the same trick for the active Coulomb part as above
      ! except here the integral sum is a sum of products of integrals with density matrix elements
      ! store the integral sum weighted with d.m. elements in cho_reduced_J
      if (Do_Cholesky) then

        do ia=1,nact
          ! ias and ibs have a different meaning here than in the old code, I'm just abusing the variables
          ias = ia+ncore
          do ib=1,nact
            ibs = ib+ncore
            cho_reduced_J_act(:,istate)=
     &        cho_reduced_J_act(:,istate)+
     &        nevpt_ijkl%cholesky_array(:,indice(ias,ibs))
     &        *da(istate,ib,ia)
            ! TODO: is there a way to multiply the density matrix element after summing over the vector elements?
            ! TODO 2: this seems like a good candidate to exploit BLAS DAXPY and/or more parallelism (though the loop in the state-specific part seems good too)
          end do
        end do
        ! now fill the Fock operator
        do i=1,nc(1)
          do j=1,i
          ! Coulomb part
            fockc(i,j,istate)=fockc(i,j,istate)+
     &        ddot(nchovec,nevpt_ijkl%cholesky_array(1,i*(i-1)/2+j),
     &         1,cho_reduced_J_act(1,istate),1)
          ! Exchange part
            do ia=1,nact
              do ib=1,nact
                davg=da(istate,ib,ia)
                fockc(i,j,istate)=fockc(i,j,istate)-
     &             0.5d0*ai(i,ib+ncore,ia+ncore,j)*davg
              end do
            end do
          end do
        end do

      else
!       ! old non-cholesky code

      do i=1,nc(isym)
      it=nord_mn_c(i,isym)
       do j=1,i
        jt=nord_mn_c(j,isym)

        do ia=1,nact
        ias=itsym(ia+ncore)
         do ib=1,nact
         ibs=itsym(ib+ncore)
         if(ias.ne.ibs) cycle
          davg=da(istate,ib,ia)
          fockc(i,j,istate)=fockc(i,j,istate)+
     &      (ai(it,jt,ia+ncore,ib+ncore)-
     &      0.5d0*ai(it,ib+ncore,ia+ncore,jt))*davg
         enddo
        enddo
       enddo
      enddo
      end if

#ifdef _BLUBB_
      print *, 'core fock for symmetry ',isym, 'with dimensions ',
     &          nc(isym)
      do i = 1, nc(isym)
        do j = 1, nc(isym)
          print '(a,i3,a,i3,a,e25.16)','fockc(',i,',',j,') = ',
     &           fockc(i,j,istate)
        end do
      end do
#endif

      info = 0
      call dsyev('V','L',nc(isym),fockc(1,1,istate),nc(isym),w,work,
     &    3*nc(isym),info)
      if(info.ne.0)then
        write (6,*) 'Something wrong in diagonalization of fockc'
        write (6,*) 'symmetry',isym
        write (6,*) 'info value is ',info
        stop -111
      endif


      do i=1,nc(isym)
      epsnew(nord_mn_c(i,isym),istate)=w(i)
      do j=1,nc(isym)
      uc(i,j,isym,istate)=fockc(i,j,istate)
      enddo
      enddo

      enddo !do istate
      deallocate(w)
      deallocate(work)
      deallocate(fockc)
      enddo !do isym

      if(iprint > 10)then
        call cpu_time(nevpt_time_end)
        write (*,'(a,f10.2,a/)'),'Core fock operator constructed in ',
     &  nevpt_time_end-nevpt_time_start,' seconds'
      end if


C     Virtual-Fock operator
! 2016-03-13 Leon: also the virtual Fock operator is split into state-specific and state-independent part
! and the two-electron integrals are treated with Cholesky decomposition

      if(iprint > 10) call cpu_time(nevpt_time_start)

      do isym=1,nsym
      if (nv(isym).eq.0) cycle
      allocate(w(nv(isym)))
      allocate(work(3*nv(isym)))

      allocate (fockv(nv(isym),nv(isym),metat))

! fill one-electron integrals
      do ir=1,nv(isym)
       irt=nord_mn_v(ir,isym)
       do is=1,ir
        ist=nord_mn_v(is,isym)
          fockv(ir,is,1)=nevpt_ijkl%oneint(indice(irt,ist))
        end do
      end do

      if (Do_Cholesky) then
      ! The Coulomb part of the Fock operator is constructed from Cholesky vectors exactly as in the core Fock operator
      ! The diagonal core-core Cholesky vector sum (cho_reduced_J_core) and active-active density-matrix weighted Cholesky vector sum
      ! (cho_reduced_J_act) are reused

      ! 8-4-2016 Leon:  added a ncore>0 safeguard because otherwise the cho_reduced_J_core is requested here which is not allocated

        if (ncore > 0) then
          do ir=1,nv(1)
            irt=nord_mn_v(ir,1)
            do is=1,ir
              ist=nord_mn_v(is,1)
                fockv(ir,is,1)=fockv(ir,is,1)+2.d0*ddot(nchovec,
     &           nevpt_ijkl%cholesky_array(1,irt*(irt-1)/2+ist),
     &           1,cho_reduced_J_core(1),1)
              do j=1,ncore
                fockv(ir,is,1)=fockv(ir,is,1)-ai(irt,j,j,ist)
              enddo
            enddo
          enddo
        end if
      else
      ! old code
        do ir=1,nv(isym)
        irt=nord_mn_v(ir,isym)
        do is=1,ir
          ist=nord_mn_v(is,isym)
          do j=1,ncore
            fockv(ir,is,1)=fockv(ir,is,1)+
     &        2.d0*ai(irt,ist,j,j)-ai(irt,j,j,ist)
          enddo
        enddo
        enddo
      end if

      ! state-specific part of the virtual Fock operator
      do istate=1,metat
      if (istate.lt.metat) then
       fockv(:,:,istate+1) = fockv(:,:,istate)
      end if

      if (Do_Cholesky) then
        ! 8-4-2016 Leon:
        ! if ncore == 0 then cho_reduced_J_act must be set up here
        ! unfortunately if I try to move its setup before the core Fock operator setup, I get a segfault
        if (ncore == 0) then
          do ia=1,nact
            ! ias and ibs have a different meaning here than in the old code, I'm just abusing the variables
            ias = ia+ncore
            do ib=1,nact
              ibs = ib+ncore
              cho_reduced_J_act(:,istate)=
     &          cho_reduced_J_act(:,istate)+
     &          nevpt_ijkl%cholesky_array(:,indice(ias,ibs))
     &         *da(istate,ib,ia)
              ! TODO: is there a way to multiply the density matrix element after summing over the vector elements?
              ! TODO 2: this seems like a good candidate to exploit BLAS DAXPY and/or more parallelism (though the loop in the state-specific part seems good too)
            end do
          end do
        end if

        do ir=1,nv(1)
        irt=nord_mn_v(ir,1)
          do is=1,ir
            ist=nord_mn_v(is,1)
            fockv(ir,is,istate)=fockv(ir,is,istate)+
     &       ddot(nchovec,
     &       nevpt_ijkl%cholesky_array(1,irt*(irt-1)/2+ist),
     &       1,cho_reduced_J_act(1,istate),1)
            do ia=1,nact
              do ib=1,nact
                fockv(ir,is,istate)=fockv(ir,is,istate)-
     &           0.5d0*ai(irt,ib+ncore,ia+ncore,ist)*
     &           da(istate,ib,ia)
              enddo
            enddo
          enddo
        enddo
      ! non-cholesky code
      else
        do ir=1,nv(isym)
        irt=nord_mn_v(ir,isym)
          do is=1,ir
            ist=nord_mn_v(is,isym)
            do ia=1,nact
              ias=itsym(ia+ncore)
              do ib=1,nact
                ibs=itsym(ib+ncore)
                if (ias.ne.ibs) cycle
                davg=da(istate,ib,ia)
                fockv(ir,is,istate)=fockv(ir,is,istate)+
     &           (ai(irt,ist,ia+ncore,ib+ncore)-
     &           0.5d0*ai(irt,ib+ncore,ia+ncore,ist))*davg
              enddo
            enddo
          enddo
        enddo
      end if

      info = 0
      call dsyev('V','L',nv(isym),fockv(1,1,istate),nv(isym),w,work,
     &    3*nv(isym),info)
      if(info.ne.0)then
        write (6,*) 'Something wrong in diagonalization of fockv'
        write (6,*) 'symmetry',isym
        write (6,*) 'info value is ',info
        stop -222
      endif

      do i=1,nv(isym)
      epsnew(nord_mn_v(i,isym),istate)=w(i)
      do j=1,nv(isym)
      uv(i,j,isym,istate)=fockv(i,j,istate)
      enddo
      enddo
      call flush(6)
      enddo !do istate

      if(iprint > 10)then
        call cpu_time(nevpt_time_end)
        write (*,'(a,f10.2,a/)'),'Virt. fock operator constructed in ',
     &  nevpt_time_end-nevpt_time_start,' seconds'
      end if

      deallocate(w)
      deallocate(work)
      deallocate(fockv)
      enddo !do isym

!     !> print the orbital energies.

      write (6,*)
      write (6,*) 'Orbital energies (in CAS terminology)'//
     &            ' used in the NEVPT calculation'
      write (6,*) '-------------------------------------'//
     &            '------------------------------'
      write (6,*)
      write (6,'(''            '',100(A14,I2))')
     &          ('       state ',MultGroup%State(i),i=1,metat)
      if(ncore > 0)then
        write (6,'(/a)') ' inactive'
        write (6,'( a)') ' --------'

        do i=1,ncore
         write(6,'(i10,3x,100f15.8)') i,(epsnew(i,j),j=1,metat)
        enddo
      else
        write (6,'(/a,12x,a)') ' inactive','none'
        write (6,'( a)') ' ---------'
      endif
      write (6,'(/a)') ' secondary'
      write (6,'( a)') ' ---------'
      do i=ncore+nact+1,norb
       write(6,'(i10,3x,100f15.8)') i,(epsnew(i,j),j=1,metat)
      enddo

!     !> set frozen orbitals
      allocate(zgel(norb+norb+1)); zgel = .false.
      ngel=0
      do i=1,nr_frozen_orb
       if (igelo(i).ne.0) then
          ngel=ngel+1
          zgel(igelo(i))=.true.
          zgel(igelo(i)+norb)=.true.
       endif
      enddo
      if (ngel.gt.0) then
         write (6,271) (igelo(i),i=1,ngel)
      endif
 271  format(/'     Frozen orbitals:'/100(5x,18i4/))

      write (6,'(/a)')
     &' -----------------------------------------------'
      write (6,'(a)')
     &' -     Starting the perturbative calculation   -'
      write (6,'(a/)')
     &' -----------------------------------------------'

      call flush(6)

      allocate (cvirt(nvirt,nvirt,metat))

      cvirt = 0

      do istate=1,metat
        do irv=1,nvirt
          ks   = itsym(irv+ncore+nact)
          kmol = nord_nev2mol(irv+ncore+nact)-nbe(ks)-nc(ks)-na(ks)
          do ik=1,nv(ks)
            knev = nord_mol2nev(ik+nbe(ks)+nc(ks)+na(ks))-ncore-nact
            cvirt(knev,irv,istate)=uv(ik,kmol,ks,istate)
          enddo
        enddo
      enddo

      allocate (ccore(ncore,ncore,metat))

      ccore = 0

      do istate=1,metat
        do i=1,ncore
          is   = itsym(i)
          kmol = nord_nev2mol(i)-nbe(is)
          do ik=1,nc(is)
            knev=nord_mol2nev(ik+nbe(is))
            ccore(knev,i,istate)=uc(ik,kmol,is,istate)
          enddo
        enddo
      enddo
c
c
c     Termine V^0 (doppie inattive)
c

      if(Do_Cholesky) then
! Prepare the arrays for storing the half-transformed Cholesky vectors (for the new newint_xxxx transformations)
        allocate(cho_reduced_vec(nchovec,nlvec,metat))
        cho_reduced_vec = 0.0d0

! (half-)transform the Cholesky vectors
        call cholesky_transform(cho_reduced_vec,ccore,cvirt,
     &         ncore,nact,nvirt,metat)
      end if

      call cpu_time(nevpt_time_start)
      ! 27-03-2016 Leon -- extended v0 with Cholesky array as an argument
      ! passing the Cholesky array as an argument even if it hasn't been init'd (ie if we don't do Cholesky)
      ! is an example of a terrible program design
      ! but otherwise the whole file has to be rewritten as a module
      call v0(ncore,nact,nvirt,epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   cvirt,cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a )'),' V(0)   contributions calculated in ',
     & nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)
      if(iprint > 0)then
        write (6,*)
        write (6,*)' Contribution of V(0) to Heff (SC and PC)'
        write (6,*)
        do istate=1,metat
          write (6,'(I2,2x,10F15.10)') MultGroup%State(istate),
     &                              (e2mp(istate,jstate),jstate=1,metat)
        enddo
        write (6,*)
        write (6,*) ' Contribution of V(0) to the norm (SC and PC)'
        write (6,*)
        do istate=1,metat
        write (6,'(''  State '',I2,''  norm = '',F15.10)')
     &           MultGroup%State(istate), psimpp(1,istate)
        enddo
        call flush(6)
      end if

!------------------------------------------------------------------
      !> read in 1-TDM
      ncoppie=metat*(metat-1)/2
      allocate(daoff(ncoppie,nact,nact))

      daoff = 0

      if(ncoppie > 0)then
        !> get data from file
        datadim(1)   = ncoppie; datadim(2:3) = nact; datadim_bound = 3

        call hdf5_get_data_dp(file_id(1), "1-TRDM" , datadim, daoff)
      end if

!------------------------------------------------------------------
!
!     term V^(+1)

      allocate(gkp1a(metat,nact,nact))
      gkp1a = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:3) = nact; datadim_bound = 3

      call hdf5_get_data_dp(file_id(1), "1-EAKO" , datadim, gkp1a)

      call cpu_time(nevpt_time_start)
      call v1k(ncore,nact,nvirt,da,daoff,gkp1a,epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   cvirt,cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a/)'),' V(+1)  contributions calculated in ',
     * nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)

      deallocate (gkp1a)
!------------------------------------------------------------------
!
!     term V^(-1)

      allocate(gkm1a(metat,nact,nact))
      gkm1a = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:3) = nact; datadim_bound = 3
      call hdf5_get_data_dp(file_id(1), "1-IPKO" , datadim, gkm1a)

      call cpu_time(nevpt_time_start)
      call vm1k(ncore,nact,nvirt,da,daoff,gkm1a,epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   ccore,cho_reduced_vec)
      call cpu_time(nevpt_time_end)

      write (6,'(a,f10.2,a/)'),' V(-1)  contributions calculated in ',
     * nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)
      deallocate (gkm1a)
!------------------------------------------------------------------
!
!     term V^(+2)

      allocate(daat(metat,nact,nact,nact,nact))
      allocate(gkp2aa(metat,nact,nact,nact,nact))
      allocate(daaoff(ncoppie,nact,nact,nact,nact))

      daat = 0; gkp2aa = 0; daaoff = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5

      call hdf5_get_data_dp(file_id(1), "2-HRDM" , datadim, daat)
      call hdf5_get_data_dp(file_id(1), "2-EAKO" , datadim, gkp2aa)

      if(ncoppie > 0)then
        datadim(1)   = ncoppie
        call hdf5_get_data_dp(file_id(1), "2-TRDM" , datadim, daaoff)
      end if

      call cpu_time(nevpt_time_start)
      call v2mod(ncore,nact,gkp2aa,daat,daaoff,daoff,epsnew,uv,uc,
     &   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,
     &   nord_mol2nev,cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a/)'),' V(+2)  contributions calculated in ',
     & nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)
      deallocate (gkp2aa)
      deallocate (daat)
!------------------------------------------------------------------
!
!     term V^(-2)

      allocate(daa(metat,nact,nact,nact,nact))
      allocate(gkm2aa(metat,nact,nact,nact,nact))

      daa = 0; gkm2aa = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5

      call hdf5_get_data_dp(file_id(1), "2-RDM"  , datadim, daa)
      call hdf5_get_data_dp(file_id(1), "2-IPKO" , datadim, gkm2aa)

      call cpu_time(nevpt_time_start)
      call vm2mod(ncore,nact,daa,daaoff,gkm2aa,epsnew,uv,uc,
     &   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     &   cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a/)'),' V(-2)  contributions calculated in ',
     & nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)
      deallocate (gkm2aa)

!------------------------------------------------------------------
!
!     term V^(0p)

      allocate(gk0pa(metat,nact,nact,nact,nact))
      allocate(gk0pd(metat,nact,nact,nact,nact))
      allocate(gk0pf(metat,nact,nact))

      gk0pa = 0; gk0pd = 0; gk0pf = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5

      call hdf5_get_data_dp(file_id(1), "gk0pa" , datadim, gk0pa)
      call hdf5_get_data_dp(file_id(1), "gk0pd" , datadim, gk0pd)
      datadim_bound = 3
      call hdf5_get_data_dp(file_id(1), "gk0pf" , datadim, gk0pf)

      call cpu_time(nevpt_time_start)
      call v0pmod(ncore,nact,nvirt,da,daa,daoff,daaoff,gk0pa,gk0pd,
     &   gk0pf,epsnew,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,
     &   nord_mol2nev,cvirt,ccore,compute_rho1st,cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a/)'),' V(0p)  contributions calculated in ',
     & nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)
      deallocate(cvirt, ccore, gk0pa, gk0pd, gk0pf)


!------------------------------------------------------------------
!
!     term V^(-1p)

      if (.not.no_4rdm_terms) then
      allocate(ro3(metat,nact,nact,nact,nact,nact,nact))
      allocate(ro3off(ncoppie,nact,nact,nact,nact,nact,nact))
      allocate(amat(metat,nact,nact,nact,nact,nact,nact))
      allocate(bmat(metat,nact,nact,nact,nact))
      allocate(cmat(metat,nact,nact,nact,nact))
      allocate(dmat(metat,nact,nact))

      ro3 = 0; amat = 0; bmat = 0; cmat = 0; dmat = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:7) = nact; datadim_bound = 7

      call hdf5_get_data_dp(file_id(1), "3-RDM" , datadim, ro3)
      call hdf5_get_data_dp(file_id(1), "amatIP", datadim, amat)
      datadim_bound = 5
      call hdf5_get_data_dp(file_id(1), "bmatIP", datadim, bmat)
      call hdf5_get_data_dp(file_id(1), "cmatIP", datadim, cmat)
      datadim_bound = 3
      call hdf5_get_data_dp(file_id(1), "dmatIP", datadim, dmat)

      if(ncoppie > 0)then
        datadim(1)   = ncoppie; datadim_bound = 7
        call hdf5_get_data_dp(file_id(1), "3-TRDM" , datadim, ro3off)
      end if

      call cpu_time(nevpt_time_start)
      call vmupE(ncore,nact,ro3,ro3off,amat,bmat,cmat,dmat,daa,da,
     &    daaoff,daoff,epsnew,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     &   nord_nev2mol,nord_mol2nev,compute_rho1st,cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a/)'),' V(-1p) contributions calculated in ',
     & nevpt_time_end-nevpt_time_start,' seconds'
      call flush(6)
!------------------------------------------------------------------
!
!     term V^(+1p)

      amat = 0; bmat = 0; cmat = 0; dmat = 0

      !> get data from file
      datadim(1)   = metat; datadim(2:7) = nact; datadim_bound = 7

      call hdf5_get_data_dp(file_id(1), "amatEA", datadim, amat)
      datadim_bound = 5
      call hdf5_get_data_dp(file_id(1), "bmatEA", datadim, bmat)
      call hdf5_get_data_dp(file_id(1), "cmatEA", datadim, cmat)
      datadim_bound = 3
      call hdf5_get_data_dp(file_id(1), "dmatEA", datadim, dmat)

      call cpu_time(nevpt_time_start)
      call vpupE(ncore,nact,ro3,ro3off,amat,bmat,cmat,dmat,daa,da,daaoff
     &     ,daoff,epsnew,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     &   nord_nev2mol,nord_mol2nev,compute_rho1st,cho_reduced_vec)
      call cpu_time(nevpt_time_end)
      write (6,'(a,f10.2,a/)'),' V(+1p) contributions calculated in ',
     & nevpt_time_end-nevpt_time_start,' seconds'
      deallocate(amat, bmat, cmat, dmat, ro3, ro3off, daa, daaoff)
      call flush(6)
      end if
!------------------------------------------------------------------
!
!     Write the final results.

      if(skip_effective_ham)then
        write (6,'(//a//)') ' state-specific NEVPT2: FINAL RESULTS'
      else
        write (6,'(//a//)') ' QD-NEVPT2: FINAL RESULTS'
      end if
      write (6,'(a)')' Second order contribution of each class to the'//
     &            ' energies'
      if(.not.skip_effective_ham)
     &write (6,'(a)') ' (to the diagonal elements of Heff)'
      do istate=1,metat
      write(6,'(/a,i4)') ' state number:', MultGroup%State(istate)
      write(6,'(a)')    ' -----------------'
      write(6,'(a)')    ' class        norm(sc)           energy(sc)'//
     &'           norm(pc)           energy(pc)'
       write (6,123) '(0)  ',psimpp(1,istate),e2mpp(1,istate)
     $      ,psienep(1,istate),e2enp(1,istate)
       write (6,123) '(+1) ',psimpp(2,istate),e2mpp(2,istate)
     $      ,psienep(2,istate),e2enp(2,istate)
       write (6,123) '(-1) ',psimpp(3,istate),e2mpp(3,istate)
     $      ,psienep(3,istate),e2enp(3,istate)
       write (6,123) '(+2) ',psimpp(4,istate),e2mpp(4,istate)
     $      ,psienep(4,istate),e2enp(4,istate)
       write (6,123) '(-2) ',psimpp(5,istate),e2mpp(5,istate)
     $      ,psienep(5,istate),e2enp(5,istate)
       write (6,123) '(+1)''',psimpp(6,istate),e2mpp(6,istate)
     $      ,psienep(6,istate),e2enp(6,istate)
       write (6,123) '(-1)''',psimpp(7,istate),e2mpp(7,istate)
     $      ,psienep(7,istate),e2enp(7,istate)
       write (6,123) '(0)'' ',psimpp(8,istate),e2mpp(8,istate)
     $      ,psienep(8,istate),e2enp(8,istate)
       print '(1x,85a1)',('-',i=1,85)
       write (6,123) 'Total',psimp(istate),e2mp(istate,istate)
     $      ,psien(istate),e2en(istate,istate)
      enddo
c     123 format (2x,A,4f20.15)
 123  format (1x,A,4f20.15)
      allocate(wr(metat),wi(metat),vr(metat,metat),vl(metat,metat))
      wr = 0; wi = 0; vr =0; vl = 0
      allocate(work2(10*metat))

      do i=1,metat
       e2mp(i,i)=e2mp(i,i)+e(i)
      enddo

      if(.not.skip_effective_ham)then

       print '(/a)','*******************************************'
       print '(a )','* Second order Effective Hamiltonian (SC) *'
       print '(a/)','*******************************************'
       call matout(metat,metat,e2mp,metat,1.d0)

       print '(a )','**************************************************'
       print '(a )','* Zero + second order Effective Hamiltonian (SC) *'
       print '(a/)','**************************************************'
       call matout(metat,metat,e2mp,metat,1.d0)
       write (6,*)
       call dgeev('V','V',metat,e2mp,metat,wr,wi,vl,metat,vr,metat,
     &            work2,size(work2),info)
       if(info.ne.0)print*,'result of diagonalization was',info
c---   riordinamento
       do i=1,metat
        do j=i+1,metat
         if(wr(j).lt.wr(i))then
            wtamp=wr(i)
            wr(i)=wr(j)
            wr(j)=wtamp
            wtamp=wi(i)
            wi(i)=wi(j)
            wi(j)=wtamp
            do k=1,metat
             wtamp=vr(k,i)
             vr(k,i)=vr(k,j)
             vr(k,j)=wtamp
             wtamp=vl(k,i)
             vl(k,i)=vl(k,j)
             vl(k,j)=wtamp
            enddo
         endif
        enddo
       enddo
c---------------------
       write (6,*)
       print*,'Eigenvalues of Heff (non-Hermitian):'
       write (6,*)
       print*,'n.       Real          Imag     '
       do istate=1,metat
        print '(i2,1x,f15.8,1x,f12.8)',MultGroup%State(istate),
     &                                 wr(istate),wi(istate)
       enddo
       write (6,*)
       print *,'Right eigenvectors'
       call matout(metat,metat,vr,metat,1.d0)
       write (6,*)
       print *,'Left eigenvectors'
       call matout(metat,metat,vl,metat,1.d0)
       write (6,*)

       !> put the eigenvalues on the diagonal of a square matrix
       call dcopy(metat,wr,1,heff_sc,metat+1)

       !> on output: eigenvectors are on heff_sc, Heff is in e2mp
       call hermitian_heff('SC',heff_sc,vr,wr,metat,e2mp,
     &                     MultGroup%State)

       !> save eigenvalues (overwrite norm - no longer needed)
       call dcopy(metat,wr,1,psimp,1)


c------calculation quadratic distances--------------------------
        if(zquadru)then
          call quadrut(metat,vr,norb,ncore,nact,da,daoff,nsym,'SC')
        endif
!------calculation quadratic distances--------------------------
      else
        !> save eigenvalues (overwrite norm - no longer needed)
        call dcopy(metat,e2mp,metat+1,psimp,1)
      end if

      if(no_pc) goto 999

      do i=1,metat
       e2en(i,i)=e2en(i,i)+e(i)
      enddo

      if(.not.skip_effective_ham)then
       print '(/a)','*******************************************'
       print '(a )','* Second order effective Hamiltonian (PC) *'
       print '(a/)','*******************************************'
       call matout(metat,metat,e2en,metat,1.d0)


       print '(/a)','**************************************************'
       print '(a )','* Zero + second order effective Hamiltonian (PC) *'
       print '(a/)','**************************************************'
       call matout(metat,metat,e2en,metat,1.d0)

       call dgeev('V','V',metat,e2en,metat,wr,wi,vl,metat,vr,metat,
     &            work2,size(work2),info)
       if(info.ne.0)print*,'result of diagonalization was',info
c---  riordinamento
       do i=1,metat
         do j=i+1,metat
          if(wr(j).lt.wr(i))then
             wtamp=wr(i)
             wr(i)=wr(j)
             wr(j)=wtamp
             wtamp=wi(i)
             wi(i)=wi(j)
             wi(j)=wtamp
             do k=1,metat
              wtamp=vr(k,i)
              vr(k,i)=vr(k,j)
              vr(k,j)=wtamp
              wtamp=vl(k,i)
              vl(k,i)=vl(k,j)
              vl(k,j)=wtamp
             enddo
          endif
         enddo
        enddo
c---------------------
        write (6,*)
        print*,'Eigenvalues of Heff (non-Hermitian):'
        write (6,*)
        print*,'n.       Real          Imag     '
        do istate=1,metat
         print '(i2,1x,f15.8,1x,f12.8)',MultGroup%State(istate),
     &                                  wr(istate),wi(istate)
        enddo
        write (6,*)
        print *,'Right eigenvectors'
        call matout(metat,metat,vr,metat,1.d0)
        write (6,*)
        print *,'Left eigenvectors'
        call matout(metat,metat,vl,metat,1.d0)
        write (6,*)

        !> put the eigenvalues on the diagonal of a square matrix
        call dcopy(metat,wr,1,heff_pc,metat+1)

        !> on output: eigenvectors are on heff_pc, Heff is in e2en
        call hermitian_heff('PC',heff_pc,vr,wr,metat,e2en,
     &                      MultGroup%State)
        !> save eigenvalues (overwrite norm - no longer needed)
        call dcopy(metat,wr,1,psien,1)

!------calculation quadratic distances--------------------------
        if(zquadru)then
           call quadrut(metat,vr,norb,ncore,nact,da,daoff,nsym,'PC')
        endif
!------calculation quadratic distances--------------------------
      else
        !> save eigenvalues (overwrite norm - no longer needed)
        call dcopy(metat,e2en,metat+1,psien,1)
      end if

999   continue
      WRITE (6,'(/)')
      if(skip_effective_ham)then
        WRITE (6,1061)
      else
        WRITE (6,1060)
      end if
      WRITE (6,2000)
      WRITE (6,'(a/)') ' state   diagonalization     '//
     &'strongly-contracted     partially-contracted'
      DO M=1,METAT
        if(.not.no_pc)then
          WRITE (6,1100) ' @', MultGroup%State(m),E(m),psimp(m),psien(m)
        else
          WRITE (6,1110) ' @', MultGroup%State(m),E(m),psimp(m),
     &                         ' not requested'
        end if
      end do
      WRITE (6,2000)
      WRITE (6,'(/)')
 1060 FORMAT ('                 Summary of the QD-NEVPT2 energies')
 1061 FORMAT ('                   Summary of the NEVPT2 energies')
 1100 FORMAT (a2,I4,F18.7,6x,F18.7,7x,F18.7)
 1110 FORMAT (a2,I4,2F19.7,9x,a)
 2000 FORMAT (1X,72('-'))

      if(allocated(wr)) deallocate(wr,wi,vr,vl,work2)
      deallocate(daoff)
      deallocate(e2enp,psienep,e2mpp,psimpp)
      deallocate(zgel)
      if(allocated(igelo)) deallocate(igelo)
      if(allocated(cho_reduced_vec)) deallocate(cho_reduced_vec)


      if(.not.compute_rho1st)then
        deallocate(nord_mn_c,nord_mn_a,nord_mn_v)
        deallocate(nord_nev2mol,nord_mol2nev)
        deallocate(nv,nc,na,nmo,nbe)
        deallocate(da)
      else
        allocate(rho(norb,norb,metat))
        rho = 0.0d0

        !> active-active
        do istate=1,metat
          do i=1,nact
            do j=1,nact
              rho(ncore+i,ncore+j,istate)=da(istate,i,j)
            enddo
          enddo
        enddo
        !> inactive-inactive
        do istate=1,metat
          do i=1,ncore
           rho(i,i,istate)=2.0d0
          enddo
        enddo

        !> inactive-virtual
        rewind 27
        do istate=1,metat
          read(27)((rho(i,j,istate),i=1,ncore),j=ncore+nact+1,norb)
          do i=1,ncore
            do j=ncore+nact+1,norb
              rho(j,i,istate)=rho(i,j,istate)
            enddo
          enddo
        enddo

        !> active-virtual
        do istate=1,metat
          read(27)((rho(i,j,istate),i=ncore+1,ncore+nact),
     &                              j=ncore+nact+1,norb)
          do i=ncore+1,ncore+nact
            do j=ncore+nact+1,norb
             rho(j,i,istate)=rho(i,j,istate)
            enddo
          enddo
        enddo

        !> inactive-active
        do istate=1,metat
          read(27)((rho(i,j,istate),i=1,ncore),j=ncore+1,ncore+nact)
          do i=1,ncore
            do j=ncore+1,ncore+nact
              rho(j,i,istate)=rho(i,j,istate)
            enddo
          enddo
        enddo

        do istate=1,metat
          print*,'Density matrix 0   for state ', istate
          call matout(norb,norb,rho(1,1,istate),norb,1.d0)
        end do

        close(27,status='keep')

        do istate=1,metat

          if(compute_nooccn)then
            total_sum = 0.0d0
            print '(/a,i3)',
     &      ' Natural orbital occupation numbers for state ',
     &        MultGroup%State(istate)
            print '(a    )',
     &      ' !NOTE! Approximate numbers --> no i-i/v-v terms!'
            print '(a/   )',
     &      ' ------------------------------------------------'
            do isym=1,nsym
              total_loc = 0.0d0
              if(nmo(isym) == 0)then
               print '(a,i3,a/)', ' irrep ',isym,': no orbitals'
              else
               print '(a,i3,a/)', ' irrep ',isym,': orbital occupations'
              end if
              allocate(w(nmo(isym)))
              allocate(work(3*nmo(isym)))

              allocate(rhotmp(nmo(isym),nmo(isym)))
              rhotmp = 0; w = 0; work = 0

              offset_c = 0
              offset_a = nc(isym)
              offset_v = nc(isym) + na(isym)

              !> symmetry: inactive - inactive block
              do ir=1,nc(isym)
                irt=nord_mn_c(ir,isym)
                do is=1,ir
                  ist=nord_mn_c(is,isym)
                  rhotmp(offset_c+ir,offset_c+is)=rho(irt,ist,istate)
                  rhotmp(offset_c+is,offset_c+ir)=rho(ist,irt,istate)
                end do
              end do

              !> symmetry: active - active block
              do ir=1,na(isym)
                irt=nord_mn_a(ir,isym)
                do is=1,ir
                  ist=nord_mn_a(is,isym)
                  rhotmp(offset_a+ir,offset_a+is)=rho(irt,ist,istate)
                  rhotmp(offset_a+is,offset_a+ir)=rho(ist,irt,istate)
                end do
              end do

              !> symmetry: virtual - virtual block
              do ir=1,nv(isym)
                irt=nord_mn_v(ir,isym)
                do is=1,ir
                  ist=nord_mn_v(is,isym)
                  rhotmp(offset_v+ir,offset_v+is)=rho(irt,ist,istate)
                  rhotmp(offset_v+is,offset_v+ir)=rho(ist,irt,istate)
                end do
              end do

#ifdef DEBUG_RHOFIRST
              print*,'Density matrix 0   for state ',
     &               MultGroup%State(istate), ' subblock',isym
              call matout(nmo(isym),nmo(isym),rhotmp,nmo(isym),1.d0)
#endif

              !> symmetry: inactive - active block
              do ir=1,nc(isym)
                irt=nord_mn_c(ir,isym)
                do is=1,na(isym)
                  ist=nord_mn_a(is,isym)
                  rhotmp(offset_c+ir,offset_a+is)=rho(irt,ist,istate)
                  rhotmp(offset_a+is,offset_c+ir)=rho(ist,irt,istate)
                end do
              end do

              !> symmetry: active - virtual block
              do ir=1,na(isym)
                irt=nord_mn_a(ir,isym)
                do is=1,nv(isym)
                  ist=nord_mn_v(is,isym)
                  rhotmp(offset_a+ir,offset_v+is)=rho(irt,ist,istate)
                  rhotmp(offset_v+is,offset_a+ir)=rho(ist,irt,istate)
                end do
              end do

              !> symmetry: inactive - virtual block
              do ir=1,nc(isym)
                irt=nord_mn_c(ir,isym)
                do is=1,nv(isym)
                  ist=nord_mn_v(is,isym)
                  rhotmp(offset_c+ir,offset_v+is)=rho(irt,ist,istate)
                  rhotmp(offset_v+is,offset_c+ir)=rho(ist,irt,istate)
                end do
              end do

#ifdef DEBUG_RHOFIRST
              print*,'Density matrix 0+1 for state ',
     &               MultGroup%State(istate), ' subblock',isym
              call matout(nmo(isym),nmo(isym),rhotmp,nmo(isym),1.d0)
#endif

              if(nmo(isym) <= 0) cycle

              call dsyev('V','L',nmo(isym),rhotmp,nmo(isym),w,work,
     *              3*nmo(isym),info)
                if(info.ne.0)then
                  write (6,*) ' failed to diagonalize rho-1st'
                  write (6,*) 'symmetry',isym
                  stop -75
                endif

                call flush(6)
                print '(6f12.8)',((w(i)),i=1,nmo(isym))

                call order_eigvalvec(rhotmp,w,nmo(isym),
     &                                        nmo(isym),nmo(isym))

                sum_c = 0.0d0; sum_a = 0.0d0; sum_v = 0.0d0
                do i = 1, nc(isym)
                  sum_c = sum_c + (w(offset_c+i))
                end do
                do i = 1, na(isym)
                  sum_a = sum_a + (w(offset_a+i))
                end do
                do i = 1, nv(isym)
                  sum_v = sum_v + (w(offset_v+i))
                end do

                do i = 1, nmo(isym)
                  total_loc = total_loc + w(i)
                  total_sum = total_sum + w(i)
                end do

                print '(/a   )',
     &          ' ------------------------------------------------'
                print '( a,f15.6/)',
     &               ' Total sum of occupation numbers =',total_loc
                print '( a,f15.6 )',
     &            ' = core    occupation numbers    :',sum_c
                print '( a,f15.6 )',
     &            ' + active  occupation numbers    :',sum_a
                print '( a,f15.6/)',
     &            ' + virtual occupation numbers    :',sum_v

              deallocate(w)
              deallocate(work)
              deallocate(rhotmp)
            enddo !do isym
            print '( a,f15.6/)',
     &            ' Final sum of occupation numbers =',total_sum

            !call order_eigvalvec(rho(1,1,istate),w,norb,norb,norb)
            !print*,'core occupation numbers:'
            !print '(6f12.8)',(w(i),i=norb-ncore+1,norb)
            !print*,'active occupation numbers:'
            !print '(6f12.8)',(w(i),i=norb-ncore-nact+1,norb-ncore)
            !print*,'virtual occupation numbers:'
            !print '(6f12.8)',(w(i),i=1,norb-ncore-nact)
          end if
        end do

        deallocate(da,rho)
        deallocate(nord_mn_c,nord_mn_a,nord_mn_v)
        deallocate(nord_nev2mol,nord_mol2nev)
        deallocate(nv,nc,na,nmo,nbe)
      endif
      if(allocated(cho_reduced_J_core)) deallocate(cho_reduced_J_core)
      if(allocated(cho_reduced_J_act)) deallocate(cho_reduced_J_act)

      call finalize_ijkl()
      call finalize_infsym()

      !> close the files nevpt.h5+ijkl.h5 and turn-off HDF5 support
      call hdf5_close(file_id(1)); call hdf5_close(file_id(2))
      call hdf5_exit()

      call cpu_time(t13)
      t13 = t13 - t1

      end subroutine qdnevpt
c----------------------------------------------------------------

      subroutine order_eigvalvec(eigenvectors,eigenvalues,ldim,
     &                           nvalues,nevec)

      implicit none

      real*8 , intent(inout) :: eigenvectors(*)
      real*8 , intent(inout) :: eigenvalues(*)
      integer, intent(in)    :: ldim
      integer, intent(in)    :: nvalues
      integer, intent(in)    :: nevec

      integer                :: offset, i, j, imin
      real*8                 :: emin

      offset = 1
      do i = 1, nvalues-1
         emin = eigenvalues(i)
         imin = i
         do j = i+1, nvalues
           if(abs(eigenvalues(j)) >= abs(emin))then
             emin = eigenvalues(j)
             imin = j
           endif
         end do
         if(imin /= i)then
           eigenvalues(imin)= eigenvalues(i)
           eigenvalues(i)   = emin
           if(nevec > 0)then
             call dswap(nevec,eigenvectors(offset),1,
     &                        eigenvectors((imin-1)*ldim+1),1)
           endif
         endif
         offset = offset + ldim
      end do

      end subroutine order_eigvalvec
c----------------------------------------------------------------

      subroutine quadrut(metat,v,norb,ncore,nact,ro,rooff,nsym
     $     ,mode)
      use info_symmetry
      implicit real*8 (a-h,o-z)
      dimension v(metat,*),ro(metat,nact,nact),rooff(metat*(metat-1)/2
     $     ,nact,nact)
      character*2 mode
      integer a,b
      dimension qn(6),q0(6),q1(6),q2(6)
      character*8 word(6)
      data word /'XXSECMOM','XYSECMOM','XZSECMOM','YYSECMOM','YZSECMOM',
     $     'ZZSECMOM'/
      allocatable q(:,:),cmo(:,:),f(:),qq(:,:)
      allocate(q(norb,norb))
c      allocate(cmo(norb,norb))
      allocate(f(norb*(norb+1)/2))
      allocate(qq(metat,metat))
c--file17 is the file produced by prop02 (file DIPQUAD)
      rewind 17
      read(17)
      read(17)
      read(17)
      do iword=1,6
       read(17)(f(i),i=1,norb*(norb+1)/2)
       ij=0
       do i=1,norb
        do j=1,i
         ij=ij+1
         q(i,j)=f(ij)
         q(j,i)=f(ij)
        enddo
       enddo
       call trasforma(norb,q,nsym)
c---- calculating all traces among metat states--
c----in MOLCAS the quadrupoles are calculated wrt the center of mass
c----notice that only electronic components are calculated !!!!!!
c----DALTON instead calculates wrt the current origin
       call traces(qq,q,ro,rooff,ncore,nact,norb,metat)
       write (6,*)
       write (6,*)
       print'(3a)','Value of ',word(iword)
     $      ,' in the primitive zero order space: '
       write (6,*)
       call matout(metat,metat,qq,metat,1.d0)
c---- transformation to the corrected space
       write (6,*)
        print'(3a)','Value of ',word(iword)
     $       ,' in the corrected zero order space: '
       do mu=1,metat
        tamp=0.d0
        do istate=1,metat
         do jstate=1,metat
          tamp=tamp+qq(istate,jstate)*v(istate,mu)*v(jstate,mu)
         enddo
        enddo
        print '(a,i2,a,f16.8)','state #',mu,' ',tamp
       enddo
      enddo
      return
      end subroutine quadrut
c------------------------------------------------
      subroutine traces(qq,q,ro,rooff,ncore,nact,norb,metat)
      implicit real*8 (a-h,o-z)
      integer a,b
      dimension qq(metat,*),q(norb,*),ro(metat,nact,nact),
     $     rooff(metat*(metat-1)/2,nact,nact)
      call zeroe(qq,metat**2)
      ijcou=0
      do istate=1,metat
       do i=1,ncore
        qq(istate,istate)=qq(istate,istate)+2.d0*q(i,i)
       enddo
       do a=1,nact
        ia=a+ncore
        do b=1,nact
         ib=b+ncore
         qq(istate,istate)=qq(istate,istate)+q(ia,ib)*ro(istate,b,a)
        enddo
       enddo
       do jstate=1,istate
        if(istate.eq.jstate)cycle
        ijcou=ijcou+1
        do a=1,nact
         ia=a+ncore
         do b=1,nact
          ib=b+ncore
          qq(istate,jstate)=qq(istate,jstate)+q(ia,ib)*rooff(ijcou,b,a)
          qq(jstate,istate)=qq(jstate,istate)+q(ia,ib)*rooff(ijcou,a,b)
         enddo
        enddo
       enddo
      enddo
      return
      end subroutine traces
c-------------------------------------------------------------
      subroutine readdip(xdip,n,key)
      implicit real*8 (a-h,o-z)
      character*8 mlabel(4),key,ast,label
      dimension xdip(*)
      ast='********'
      nn=n*(n+1)/2
      rewind 17
 10   read(17,end=100)mlabel
      if(.not.(mlabel(1).eq.ast.and.mlabel(4).eq.key))goto10
      read(17)(xdip(i),i=1,nn)
      print '(2a)','dipole  matrix ',key
 100  continue
      return
      end subroutine readdip
c---------------------------------------
      SUBROUTINE OFILES()
      use nevpt2_cfg ! input parameters
      implicit none
      open(4,file=file04,form='UNFORMATTED',status='UNKNOWN')
      rewind 4
      END subroutine ofiles
C***********************************************************************
      SUBROUTINE MATOUT(N,N2,A,NN,factor)
      real*8 A,factor
      DIMENSION A(NN,*)
C
   20 FORMAT(7X,5(8X,I2,7X))
   30 FORMAT(1X,I2,4X,5(F17.10))
   35 FORMAT(8X,11(5X,I3,3X))
   37 FORMAT(8X,11(4X,A4,3X))
   36 FORMAT(/)
C
      DO 40 M=1,N2,5
      K=M+4
      IF(K.GE.N2) K=N2
   10 WRITE(6,20) (J,J=M,K)
      DO 40 I=1,N
      WRITE(6,30) I,(A(I,J)/factor,J=M,K)
 40   CONTINUE
      call flush(6)
      RETURN
      END subroutine matout
C***********************************************************************
      subroutine v0(ncore,nact,nvirt,epsnew,uv,uc,
     &   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     &   cvirt,cho_reduced_vec)
      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec

      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     &     NORB,NORB2,NOCB,NOCA,
     &     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     &     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     &     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap

Cele QD_new
      real *8 epsnew(norb,metat),cvirt(nvirt,nvirt,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable dumai(:,:,:),traint(:,:,:)
      real*8 cho_reduced_vec(nchovec,nlvec,metat)

      if (.not.Do_Cholesky) then
        allocate (dumai(nvirt,nvirt,metat))
      end if
        allocate (traint(nvirt,nvirt,metat))
Cele QD_new

      do 800 ii=1,ncore
         if (zgel(ii)) goto 800
         isy=itsym(ii)
         do 700 ij=ii,ncore
            if (zgel(ij)) goto 700
            jsy=itsym(ij)
            ijsy=its(isy,jsy)
#ifdef _DEBUG_
            if (zverbose) then
            write (6,*) 'Excitations ',ii,ij,' ->'
            call flush(6)
            endif
#endif


        traint = 0
         if (Do_Cholesky) then
           call newint_irsj_cholesky(cho_reduced_vec,traint,
     &    metat,ij,ii,ncore,nact,nvirt)
         else
           dumai  = 0
           call newint_irsjnew(metat,norb,ncore,nact,nvirt,
     *      factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *      nord_nev2mol,nord_mol2nev,ij,ii,cvirt,dumai,traint)
         end if
            do 600 ir=ncore+nact+1,norb
               call flush(6)
               irv=ir-ncore-nact
               irsy=itsym(ir)
               irijsy=its(ijsy,irsy)

               do 500 is=ir,norb
                  isv=is-ncore-nact
                  issy=itsym(is)
                  irsijsy=its(irijsy,issy)
                  if(irsijsy.ne.1)goto 500

                  zrs=ir.eq.is
                  zij=ii.eq.ij

                  do istate=1,metat
                  ab2 = 0.0d0
! !                   if (Do_Cholesky) then
!                     ab1 = cholesky_traint(cho_reduced_vec,
!      &                   ii,ir,ij,is,istate,metat)
!                     if(.not.(zrs.or.zij)) then
!                       ab2 = cholesky_traint(cho_reduced_vec,
!      &                   ii,is,ij,ir,istate,metat)
!                     end if
!                   else
                  ab1 = traint(irv,isv,istate)
                  if(.not.(zrs.or.zij)) then
                    ab2  = traint(isv,irv,istate)
                  end if
!                   end if
                  if(zrs.and.zij)       dint = ab1
                  if(zrs.neqv.zij)      dint = sqrt(2.d0)*ab1
                  if(.not.(zrs.or.zij)) dint = 2.d0*sqrt(ab1**2+ab2**2
     $               -ab1*ab2)
                  co=-dint/(epsnew(ir,istate)+epsnew(is,istate)-
     *                      epsnew(ii,istate)-epsnew(ij,istate))
                  e2mp(istate,istate)=e2mp(istate,istate)+dint*co
                  e2en(istate,istate)=e2mp(istate,istate)
                  psimp(istate)=psimp(istate)+co*co
                  psien(istate)=psimp(istate)
                  e2mpp(1,istate)=e2mpp(1,istate)+dint*co
                  psimpp(1,istate)=psimpp(1,istate)+co*co
                  e2enp(1,istate)=e2mpp(1,istate)
                  psienep(1,istate)=psimpp(1,istate)
                  enddo
 500           enddo
 600        enddo

 700     enddo

 800  enddo

      if(allocated(dumai)) deallocate (dumai)
      if(allocated(traint)) deallocate (traint)

      return
      end subroutine v0
c----------------------------------------------------------------
      subroutine v1k(ncore,nact,nvirt,da,daoff,gkp1a,epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   cvirt,cho_reduced_vec)
      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/thresho/thresh
      real*4 tarray(2)
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension da(metat,nact,nact),gkp1a(metat,nact,nact),
     $     daoff(metat*(metat-1)/2,nact,nact)
      allocatable s(:,:),coef(:,:,:),w(:,:),work(:),caux(:,:),
     $     dnormk(:,:),epsimk(:),den(:),e2(:,:),psi2(:)

Cele QD_new
      real *8 epsnew(norb,metat),cvirt(nvirt,nvirt,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable cint(:),dumai(:,:,:),traint(:,:,:),traintk(:,:,:)

      real*8 cho_reduced_vec(nchovec,nlvec,metat)

      allocate(cint(metat))
      if (.not.Do_Cholesky) then
        allocate (dumai(nact,nvirt,metat))
        dumai=0.0d0
      end if

      allocate (traint(nact,nvirt,metat))
      allocate (traintk(nact,nvirt,metat))
      traint=0.0d0
      traintk=0.0d0
Cele QD_new

c
      allocate(dnormk(metat,metat))
      allocate(epsimk(metat))
      allocate(den(metat))
      allocate(psi2(metat))
      allocate(e2(metat,metat))
      e2 = 0.0d0
      call cpu_time(tini)
      do  ii=1,ncore
         if (zgel(ii)) goto 800
         isy=itsym(ii)

         do  ij=ii,ncore
            if (zgel(ij)) goto 700
            jsy=itsym(ij)
            ijsy=its(isy,jsy)

            if (zverbose) then
            write (6,*) 'Excitations ',ii,ij,' ->'
            call flush(6)
            endif

      if (Do_Cholesky) then
       call newint_airj_cholesky(cho_reduced_vec,traint,
     &    metat,ii,ij,ncore,nact,nvirt)
      if (ij.ne.ii)
     & call newint_airj_cholesky(cho_reduced_vec,traintk,
     &    metat,ij,ii,ncore,nact,nvirt)
      else
      call newint_airjnew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ii,ij,traint,dumai,cvirt)
      if (ij.ne.ii)
     *call newint_airjnew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ij,ii,traintk,dumai,cvirt)
      end if
            do  is=ncore+nact+1,norb
               issy=itsym(is)
               isijsy=its(ijsy,issy)
               isv=is-ncore-nact
c      write (6,*) 'Eccitazione con gli indici',ii,ij,' -> ',is

               do istate=1,metat
                  dnormk(istate,istate)=0.d0
                  epsimk(istate)=0.d0
                  do jstate=1,istate-1
                     dnormk(istate,jstate)=0.d0
                     dnormk(jstate,istate)=0.d0
                  enddo
               enddo

               do ia=1,nact
                  iaa=ia+ncore
                  iasy=itsym(iaa)
                  if(its(isijsy,iasy).ne.1)cycle

                  do iap=1,nact
                     iaap=iap+ncore
                     iapsy=itsym(iaap)
                     if(its(isijsy,iapsy).ne.1)cycle

                     do istate=1,metat
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,ia,ii,isv,ij,istate,contj)
        contj=traint(ia,isv,istate)
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,iap,ii,isv,ij,istate,contjp)
        contjp=traint(iap,isv,istate)
                     if (ij.ne.ii) then
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,ia,ij,isv,ii,istate,contk)
        contk=traintk(ia,isv,istate)
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,iap,ij,isv,ii,istate,contkp)
        contkp=traintk(iap,isv,istate)
                     endif
                     if (ij.eq.ii) then
                        cint(istate)=0.5d0*contj*contjp
                     else
                        cint(istate)=contj*contjp+contk*contkp-
     *                       0.5d0*contj*contkp-0.5d0*contk*contjp
                     endif
                        if (ia.eq.iap) then
                           ro=2.d0-da(istate,ia,iap)
                        else
                           ro=-da(istate,ia,iap)
                        endif
                        dnormk(istate,istate)=dnormk(istate,istate)+2.d0
     $                       *ro*cint(istate)
                        epsimk(istate)=epsimk(istate)+gkp1a(istate,ia
     $                       ,iap)*cint(istate)
                        do jstate=1,istate-1
                           ijcou=(istate-1)*(istate-2)/2+jstate
                           ro=-daoff(ijcou,iap,ia)
                           ro2=-daoff(ijcou,ia,iap)
                           dnormk(jstate,istate)=dnormk(jstate,istate)+2
     $                          .d0*ro*cint(istate)
                           dnormk(istate,jstate)=dnormk(istate,jstate)+2
     $                          .d0*ro2*cint(jstate)
                        enddo
                     enddo
                  enddo
               enddo
               do istate=1,metat
                  den(istate)=0.d0
                  if (abs(dnormk(istate,istate)).gt.0.d0)then
                     epsimk(istate)=epsimk(istate)/dnormk(istate,istate)
                     den(istate)=epsnew(ii,istate)+epsnew(ij,istate)-
     *               epsnew(is,istate)-epsimk(istate)
                     co=dnormk(istate,istate)/den(istate)
                     con=co/den(istate)
                     e2mp(istate,istate)=e2mp(istate,istate)+co
                     psimp(istate)=psimp(istate)+con
                     e2mpp(2,istate)=e2mpp(2,istate)+co
                     psimpp(2,istate)=psimpp(2,istate)+con
                  endif
               enddo
               do istate=1,metat
                  e2(istate,istate)=e2mpp(2,istate)
                  do jstate=1,istate-1
                     if(den(istate).ne.0.d0)e2mp(jstate,istate)
     $                    =e2mp(jstate,istate)+dnormk(jstate,istate)
     $                    /den(istate)
                     if(den(jstate).ne.0.d0)e2mp(istate,jstate)
     $                    =e2mp(istate,jstate)+dnormk(istate,jstate)
     $                    /den(jstate)
                     if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate
     $                    ,istate)+dnormk(jstate,istate)/den(istate)
                     if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate
     $                    ,jstate)+dnormk(istate,jstate)/den(jstate)
                  enddo
               enddo
 600        enddo

 700     enddo

 800  enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*)' Contribution of V(+1) to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(+1) to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(2,istate)
      enddo
      end if
      deallocate(cint)
      if(allocated(traint)) deallocate(traint)
      if(allocated(traintk)) deallocate(traintk)
      if(allocated(dumai)) deallocate(dumai)
      deallocate(dnormk)
      deallocate(epsimk)
      deallocate(den)
      write (6,*)
      call flush(6)

      if(no_pc) goto 999
c----Partially contracted NEV-PT---------------------------------
      allocate(s(1:nact,1:nact))
      allocate(coef(nact,nact,metat))
      allocate(w(nact,metat))
      allocate(work(3*nact))
      do istate=1,metat
         do i=1,nact
            do j=1,nact
               if(i.ne.j)then
                  s(i,j)=-da(istate,i,j)
               else
                  s(i,j)=2.0d0-da(istate,i,j)
               endif
               coef(i,j,istate)=gkp1a(istate,i,j)*0.5d0
              enddo
             enddo
         info = 0
         call dsyev('V','U',nact,s,nact,w(1,istate),work,3*nact,info)
      if(info.ne.0)then
         print *,'something wrong in diagonalization: sub v1k'
         print *,' info = ',info
         stop -999
      endif
      thr=thresh
      ndep=0
      do i=1,nact
       if(w(i,istate).lt.0.d0)then
        ndep=ndep+1
        write (6,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $           ,w(i,istate)
        elseif(w(i,istate).lt.thr)then
            ndep=ndep+1
        endif
      enddo
      if (ndep.gt.0) then
      write (6,*)'Linear dependencies analysis:'
      write (6,'(a,i4,a)')' There are ',ndep,' linear dependencies'
      endif
      nnew=nact-ndep
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact,1:nnew))
      call caux0p(caux,s,w(1,istate),nact,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact,nnew,1.d0,coef(1,1,istate),nact,caux,nact
     $     ,0.d0,s,nact)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact,1.d0,caux,nact,s,nact,0.d0
     $     ,coef(1,1,istate),nact)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      info = 0
      call dsyev('V','U',nnew,coef(1,1,istate),nact,w(1,istate),work,3
     $     *nact,info)
      if(info.ne.0)then
         print *,'something wrong in diagonalization: coef'
         print *,' info = ',info
         stop -999
      endif
c--   back transformation
c--   moltiplico CAUX(nact,nnew) per coef(nnew,nnew)
c--   e metto in S(nact,nnew)
      call dgemm('N','N',nact,nnew,nnew,1.d0,caux,nact,coef(1,1,istate)
     $     ,nact,0.d0,s,nact)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call zeroe(coef(1,1,istate),nact*nact)
      call copia(coef(1,1,istate),s,nact,nnew)
      deallocate(caux)
      ioffw=1
      if(zverbose)then
         print*,'eigenvalues are:'
         do i=1,nact
            print '(a,i2,a,f12.8)','eps(',i,')=',w(i,istate)
         enddo
      endif
      enddo
      deallocate(s)
      call e2p1cont(e2,psi2,coef,w,da,daoff,epsnew,uv,uc,factor,
     *      nact,ncore,nvirt,norb,nnew,metat,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,zgel,
     *      cvirt,cho_reduced_vec)
      if(print_level > 0)then
      write (6,*)
      write (6,*)' Contribution of V(+1) to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(+1) to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),psi2(istate)
      enddo
      end if
      do istate=1,metat
         e2enp(2,istate)=e2(istate,istate)
         psienep(2,istate)=psi2(istate)
         psien(istate)=psien(istate)+psienep(2,istate)
         e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
         do jstate=1,istate-1
            e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
            e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
         enddo
      enddo
      deallocate(coef)
      deallocate(w)
      deallocate(work)
999   deallocate(psi2)
      deallocate(e2)
      return
      end subroutine v1k
c----------------------------------------------------------------
      subroutine fill2 (s,da,nact)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension s(nact,nact),da(nact,nact)
      do i=1,nact
         do j=1,nact
            if(i.ne.j)then
               s(i,j)=-da(i,j)
            else
               s(i,i)=2.d0-da(i,i)
            endif
         enddo
      enddo
      return
      end subroutine fill2
c-----------------------------------------------------------
      subroutine vm1k(ncore,nact,nvirt,da,daoff,gkm1a,epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   ccore,cho_reduced_vec)
      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
!
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      common/thresho/thresh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension da(metat,nact,nact),daoff(metat*(metat-1)/2,nact,nact)
     $     ,gkm1a(metat,nact,nact)
      allocatable s(:,:),coef(:,:,:),w(:,:),work(:),caux(:,:),dnormk(:,:
     $     ),epsimk(:),den(:),e2(:,:),psi2(:)
Cele QD_new
      real *8 epsnew(norb,metat),ccore(ncore,ncore,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable cint(:),contj(:),contk(:),contjp(:),contkp(:),
     *            traint(:,:,:),traintk(:,:,:),dumai(:,:,:)

      ! Cholesky half-transformed vectors
      real*8 cho_reduced_vec(nchovec,nlvec,metat)

      allocate(cint(metat))
      allocate(contj(metat))
      allocate(contk(metat))
      allocate(contjp(metat))
      allocate(contkp(metat))
      if (.not.Do_Cholesky) allocate (dumai(ncore,nact,metat))
      allocate (traint(ncore,nact,metat))
      allocate (traintk(ncore,nact,metat))
Cele QD_new

c
      allocate(dnormk(metat,metat))
      allocate(e2(metat,metat))
      allocate(epsimk(metat))
      allocate(den(metat))
      allocate(psi2(metat))
      e2 = 0.0d0

      do 800 ir=ncore+nact+1,norb
      irv=ir-ncore-nact
      irsy=itsym(ir)

      do 700 is=ir,norb
      isv=is-ncore-nact
      issy=itsym(is)
      irssy=its(irsy,issy)

      if (zverbose) then
      write (6,*) 'Excitations ',' ->',ir,is
      call flush(6)
      endif

      if (Do_Cholesky) then
        call newint_ris_a_cholesky(cho_reduced_vec,traint,
     &    metat,ir,is,ncore,nact,nvirt)
        if (ir.ne.is) then
          call newint_ris_a_cholesky(cho_reduced_vec,traintk,
     &    metat,is,ir,ncore,nact,nvirt)
        end if
      else

      call newint_ris_anew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,irv,isv,traint,dumai,ccore)
      if (ir.ne.is)
     *call newint_ris_anew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,isv,irv,traintk,dumai,ccore)
      end if
      do 600 ii=1,ncore
      if (zgel(ii)) goto 600
      isy=itsym(ii)
      irsisy=its(irssy,isy)



c      write (6,*) 'Eccitazione con gli indici',ii,' -> ',ir,is

      do istate=1,metat
         dnormk(istate,istate)=0.d0
         epsimk(istate)=0.d0
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            dnormk(istate,jstate)=0.d0
            dnormk(jstate,istate)=0.d0
         enddo
      enddo

c      call newint_ris_a(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,irv,ii,isv,dintj)
c      call newint_ris_a(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,isv,ii,irv,dintk)

      do ia=1,nact
         iaa=ia+ncore
         iasy=itsym(iaa)
         if(its(irsisy,iasy).ne.1)cycle

      do istate=1,metat
         contj(istate)=traint(ii,ia,istate)
         if (ir.ne.is)  contk(istate)=traintk(ii,ia,istate)
c         write (6,*) contj(istate),contk(istate)
      enddo

         do iap=1,nact
            iaap=iap+ncore
            iapsy=itsym(iaap)
            if(its(irsisy,iapsy).ne.1)cycle


      do istate=1,metat
       contjp(istate)=traint(ii,iap,istate)
       if (ir.ne.is) contkp(istate)=traintk(ii,iap,istate)
c         write (6,*) contjp(istate),contkp(istate)
          if (ir.eq.is) then
             cint(istate)=.5d0*contj(istate)*contjp(istate)
          else
             cint(istate)=contj(istate)*contjp(istate)+
     *          contk(istate)*contkp(istate)-
     *          0.5d0*contj(istate)*contkp(istate)-
     *          0.5d0*contk(istate)*contjp(istate)
          endif
               dnormk(istate,istate)=dnormk(istate,istate)+2.d0
     $              *da(istate,ia,iap)*cint(istate)
               epsimk(istate)=epsimk(istate)+
     *               gkm1a(istate,ia,iap)*cint(istate)
               do jstate=1,istate-1
                  ijcou=(istate-1)*(istate-2)/2+jstate
                  dnormk(jstate,istate)=dnormk(jstate,istate)+2.d0
     $                 *daoff(ijcou,ia,iap)*cint(istate)
                  dnormk(istate,jstate)=dnormk(istate,jstate)+2.d0
     $                 *daoff(ijcou,iap,ia)*cint(jstate)
               enddo
            enddo
         enddo
      enddo
      do istate=1,metat
         den(istate)=0.d0
         if (abs(dnormk(istate,istate)).ge.1.d-14) then
            epsimk(istate)=-epsimk(istate)/dnormk(istate,istate)
            den(istate)=epsimk(istate)+epsnew(ii,istate)-
     *                  epsnew(ir,istate)-epsnew(is,istate)
            co=dnormk(istate,istate)/den(istate)
            con=co/den(istate)
            e2mp(istate,istate)=e2mp(istate,istate)+co
            psimp(istate)=psimp(istate)+con
            e2mpp(3,istate)=e2mpp(3,istate)+co
            psimpp(3,istate)=psimpp(3,istate)+con
         endif
      enddo
      do istate=1,metat
         e2(istate,istate)=e2mpp(3,istate)
         do jstate=1,istate-1
            if(den(istate).ne.0.d0)e2mp(jstate,istate)=e2mp(jstate
     $           ,istate)+dnormk(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2mp(istate,jstate)=e2mp(istate
     $           ,jstate)+dnormk(istate,jstate)/den(jstate)
            if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate,istate)
     $           +dnormk(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate,jstate)
     $           +dnormk(istate,jstate)/den(jstate)
         enddo
      enddo



 600  enddo

 700  enddo

 800  enddo

      if(print_level > 0)then
      write (6,*)
      write (6,*)' Contribution of V(-1) to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(-1) to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(3,istate)
      enddo
      end if
      if(allocated(traint))deallocate(traint)
      if(allocated(traintk))deallocate(traintk)
      if(allocated(dumai))deallocate(dumai)
      deallocate(cint)
      deallocate(contj)
      deallocate(contk)
      deallocate(contjp)
      deallocate(contkp)
      deallocate(dnormk)
      deallocate(epsimk)
      deallocate(den)
      call flush(6)

      if(no_pc) goto 999
c----Partially contracted NEV-PT-------------
      allocate(s(1:nact,1:nact))
      allocate(coef(nact,nact,metat))
      allocate(w(nact,metat))
      allocate(work(3*nact))
      do istate=1,metat
         do i=1,nact
            do j=1,nact
               s(i,j)=da(istate,i,j)
               coef(i,j,istate)=gkm1a(istate,i,j)*0.5d0
            enddo
         enddo
         info = 0
         call dsyev('V','U',nact,s,nact,w(1,istate),work,3*nact,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: s (v(-1))'
           print *,' info = ',info
           stop -888
         endif
         thr=thresh
         ndep=0
         do i=1,nact
            if(w(i,istate).lt.0.d0)then
               ndep=ndep+1
        if(zverbose) write (6,'(a,i3,a,f12.8)')
     $              'negative eigenvalue: eps(',i,')=',w(i,istate)
            elseif(w(i,istate).lt.thr)then
               ndep=ndep+1
            endif
         enddo
         if (ndep.gt.0.and.zverbose) then
            write (6,*)'Linear dependencies analysis:'
            write (6,'(a,i4,a)')' There are ',ndep
     $           ,' linear dependencies'
         endif
         nnew=nact-ndep
c--   costruzione della matrice di trasformazione
         allocate (caux(1:nact,1:nnew))
         call caux0p(caux,s,w(1,istate),nact,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
         call dsymm('L','U',nact,nnew,1.d0,coef(1,1,istate),nact,caux
     $        ,nact,0.d0,s,nact)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
         call dgemm('T','N',nnew,nnew,nact,1.d0,caux,nact,s,nact,0.d0
     $        ,coef(1,1,istate),nact)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact,w(1,istate),work
     $        ,3*nact,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: coef (v(-1))'
           print *,' info = ',info
           stop -888
         endif
c--   back transformation
c--   moltiplico CAUX(nact,nnew) per coef(nnew,nnew)
c--   e metto in S(nact,nnew)
         call dgemm('N','N',nact,nnew,nnew,1.d0,caux,nact,coef(1,1
     $        ,istate),nact,0.d0,s,nact)
c--   metto in coef i coefficienti trasformati (for clarity sake)
         call zeroe(coef(1,1,istate),nact*nact)
         call copia(coef(1,1,istate),s,nact,nnew)
         deallocate(caux)
         if(zverbose)then
            print*,'eigenvalues are:'
            do i=1,nact
               print '(a,i2,a,f12.8)','eps(',i,')=',w(i,istate)
            enddo
         endif
      enddo
      deallocate(s)
      call e2m1cont(e2,psi2,coef,w,da,daoff,nact,nnew,ncore
     $     ,norb,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,
     *      ccore,cho_reduced_vec)
      if(print_level > 0)then
      write (6,*)
      write (6,*)' Contribution of V(-1) to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(-1) to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psi2(istate)
      enddo
      end if
      call flush(6)
      do istate=1,metat
         e2enp(3,istate)=e2(istate,istate)
         psienep(3,istate)=psi2(istate)
         psien(istate)=psien(istate)+psienep(3,istate)
         e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
            e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
         enddo
      enddo
      deallocate(coef)
      deallocate(w)
      deallocate(work)
999   deallocate(psi2)
      return
      end subroutine vm1k
c-------------------------------------------------------------
      subroutine e2m1cont(e2,psi2,coef,w,da,daoff,nact,nnew,ncore
     $     ,norb,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,
     *      ccore,cho_reduced_vec)

      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec

      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(nact,*),coef(nact,nact,metat),da(metat,nact
     $     ,nact),daoff(metat*(metat-1)/2,nact,nact)
      dimension zgel(*),e2(metat,*),psi2(*)
      allocatable s(:,:,:,:),rsmui(:,:),rsmuip(:,:),deno(:)
      integer a,ap
Cele QD_new
      real *8 epsnew(norb,metat),ccore(ncore,ncore,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable traint(:,:,:),traintk(:,:,:),dumai(:,:,:)

      ! Cholesky half-transformed vectors
      real*8 cho_reduced_vec(nchovec,nlvec,metat)

      nvirt=norb-nact-ncore
      if (.not.Do_Cholesky)allocate (dumai(ncore,nact,metat))
      allocate (traint(ncore,nact,metat))
      allocate (traintk(ncore,nact,metat))
Cele QD_new

      allocate(s(nact,nact,metat,metat))
      allocate(rsmui(metat,metat))
      allocate(rsmuip(metat,metat))
      allocate(deno(metat))
      ncoppie=metat*(metat-1)/2
      s = 0.0d0
c---costruzione matrice S
      do a=1,nact
         do mu=1,nact
            do ap=1,nact
               ijcou=0
               do istate=1,metat
                  s(a,mu,istate,istate)=s(a,mu,istate,istate)+da(istate
     $                 ,ap,a)*coef(ap,mu,istate)
                  do jstate=1,istate
                     if(jstate.eq.istate)cycle
                     ijcou=ijcou+1
                     coe=coef(ap,mu,istate)
                     if(coe.ne.0.d0)s(a,mu,jstate,istate)=s(a,mu,jstate
     $                    ,istate)+daoff(ijcou,ap,a)*coe
                     coe=coef(ap,mu,jstate)
                     if(coe.ne.0.d0)s(a,mu,istate,jstate)=s(a,mu,istate
     $                    ,jstate)+daoff(ijcou,a,ap)*coe
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      do ir=ncore+nact+1,norb
      irv=ir-ncore-nact
      irs=itsym(ir)
         do is=ir,norb
         isv=is-ncore-nact
         iss=itsym(is)
         irss=its(irs,iss)

         if (Do_Cholesky) then
           call newint_ris_a_cholesky(cho_reduced_vec,traint,
     &       metat,ir,is,ncore,nact,nvirt)
           if (ir.ne.is) then
              call newint_ris_a_cholesky(cho_reduced_vec,traintk,
     &        metat,is,ir,ncore,nact,nvirt)
           end if
         else
           call newint_ris_anew(metat,norb,ncore,nact,nvirt,
     *      factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *      nord_nev2mol,nord_mol2nev,irv,isv,traint,dumai,ccore)
            if(ir.ne.is)
     *        call newint_ris_anew(metat,norb,ncore,nact,nvirt,
     *         factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *         nord_nev2mol,nord_mol2nev,isv,irv,traintk,dumai,ccore)
         end if

         do i=1,ncore
           if (zgel(i)) goto 1
           iis=itsym(i)
           iirss=its(irss,iis)
c      call newint_ris_a(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,irv,i,isv,dintj)
c      call newint_ris_a(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,isv,i,irv,dintk)
               do mu=1,nact
                  call zeroe(rsmui,metat**2)
                  call zeroe(rsmuip,metat**2)
                  do a=1,nact
                  ias=itsym(a+ncore)
                  if (ias.ne.iirss) goto 2
                     if(ir.ne.is)then
                        do istate=1,metat
                        cint=traint(i,a,istate)+traintk(i,a,istate)
                        bint=traint(i,a,istate)-traintk(i,a,istate)
                           rsmui(istate,istate)=rsmui(istate,istate)
     $                          +cint*s(a,mu,istate,istate)
                           rsmuip(istate,istate)=rsmuip(istate,istate)
     $                          +bint*s(a,mu,istate,istate)
                           do jstate=1,metat
                              if(jstate.eq.istate)cycle
                              rsmui(jstate,istate)=rsmui(jstate,istate)
     $                             +cint*s(a,mu,jstate,istate)
                              rsmuip(jstate,istate)=rsmuip(jstate,istate
     $                             )+bint*s(a,mu,jstate,istate)
                           enddo
                        enddo
                     else
                        do istate=1,metat
                         cint=traint(i,a,istate)
                           rsmui(istate,istate)=rsmui(istate,istate)
     $                          +cint*s(a,mu,istate,istate)
                           do jstate=1,metat
                              if(jstate.eq.istate)cycle
                              rsmui(jstate,istate)=rsmui(jstate,istate)
     $                             +cint*s(a,mu,jstate,istate)
                           enddo
                        enddo
                     endif
  2               enddo ! do a
                  do istate=1,metat
                     deno(istate)=w(mu,istate)+epsnew(ir,istate)+
     *                        epsnew(is,istate)-epsnew(i,istate)
                  enddo
                  if(ir.ne.is)then
                     do istate=1,metat
                        contri=(0.5d0*rsmui(istate,istate)**2+1.5d0
     $                       *rsmuip(istate,istate)**2)/deno(istate)
                     e2(istate,istate)=e2(istate,istate)-contri
                     psi2(istate)=psi2(istate)+contri/deno(istate)
                     do jstate=1,istate
                        if(jstate.eq.istate)cycle
                        contri=(0.5d0*rsmui(jstate,istate)*rsmui(istate
     $                       ,istate)+1.5d0*rsmuip(jstate,istate)
     $                       *rsmuip(istate,istate))/deno(istate)
                        e2(jstate,istate)=e2(jstate,istate)-contri
                        contri=(0.5d0*rsmui(istate,jstate)*rsmui(jstate
     $                       ,jstate)+1.5d0*rsmuip(istate,jstate)
     $                       *rsmuip(jstate,jstate))/deno(jstate)
                        e2(istate,jstate)=e2(istate,jstate)-contri
                     enddo
                  enddo
                  else
                     do istate=1,metat
                        contri=rsmui(istate,istate)**2/deno(istate)
                        e2(istate,istate)=e2(istate,istate)-contri
                        psi2(istate)=psi2(istate)+contri/deno(istate)
                        do jstate=1,istate
                           if(jstate.eq.istate)cycle
                           contri=rsmui(jstate,istate)*rsmui(istate
     $                          ,istate)/deno(istate)
                           e2(jstate,istate)=e2(jstate,istate)-contri
                           contri=rsmui(istate,jstate)*rsmui(jstate
     $                          ,jstate)/deno(jstate)
                           e2(istate,jstate)=e2(istate,jstate)-contri
                        enddo
                     enddo
                  endif
               enddo
 1          enddo
         enddo
      enddo
      deallocate(s)
      deallocate(rsmui)
      deallocate(rsmuip)
      deallocate(deno)
      deallocate(traint)
      deallocate(traintk)
      if (allocated(dumai))deallocate(dumai)
      return
      end subroutine e2m1cont
c-----------------------------------------------------------
      subroutine e2m2cont(e2,psi2,coef,w,daa,daaoff,nact,ncore,norb
     $     ,nnew,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(nact**2,*),coef(nact**2,nact**2,metat),
     $     daa(metat,nact,nact,nact,nact),
     $     daaoff(metat*(metat-1)/2,nact,nact,nact,nact)
      dimension zgel(*),e2(metat,*),psi2(*)
      allocatable s(:,:,:,:),rsmu(:,:),deno(:)
      integer a,ap,b,bp,ab,abp

      ! Cholesky half-transformed vectors
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable dint1(:,:,:)
      allocate (dint1(nact,nact,metat))
      nvirt=norb-nact-ncore
Cele QD_new

      ncoppie=metat*(metat-1)/2
      allocate(s(nact**2,nact**2,metat,metat))
      allocate(rsmu(metat,metat))
      allocate(deno(metat))
      s = 0.0d0
      nact2=nact**2
c--caso r#s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            do mu=1,nact2
               abp=0
               do ap=1,nact
                  do bp=1,nact
                     abp=abp+1
                     ijcou=0
                     do istate=1,metat
                        coe=coef(abp,mu,istate)
                        if(coe.ne.0.d0)s(ab,mu,istate,istate)=s(ab,mu
     $                       ,istate,istate)+daa(istate,ap,bp,a,b)
     $                       *coe
                        do jstate=1,istate
                           if(jstate.eq.istate)cycle
                           ijcou=ijcou+1
                           coe=coef(abp,mu,istate)
                           if(coe.ne.0.d0)s(ab,mu,jstate,istate)=s(ab,mu
     $                         ,jstate,istate)+daaoff(ijcou,bp,ap,b,a)
     $                          *coe
                           coe=coef(abp,mu,jstate)
                           if(coe.ne.0.d0)s(ab,mu,istate,jstate)=s(ab,mu
     $                          ,istate,jstate)+daaoff(ijcou,b,a,bp,ap)
     $                          *coe
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      do ir=ncore+nact+1,norb
      irv=ir-ncore-nact
      irs=itsym(ir)
         do is=ir+1,norb
         isv=is-ncore-nact
         iss=itsym(is)
         irss=its(irs,iss)
          if (Do_Cholesky) then
            call newint_rasb_cholesky(cho_reduced_vec,dint1,
     &    metat,ir,is,ncore,nact,nvirt)
          else
      call newint_rs_ab(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,irv,isv,dint1)
          end if
            do mu=1,nact2
               call zeroe(rsmu,metat**2)
               ab=0
               do a=1,nact
               ias=itsym(a+ncore)
               irsas=its(irss,ias)
                  do b=1,nact
                     ab=ab+1
                     ibs=itsym(b+ncore)
                     irsabs=its(irsas,ibs)
                     if (irsabs.ne.1) cycle
                     do istate=1,metat
                      rsmu(istate,istate)=rsmu(istate,istate)+
     $                 dint1(b,a,istate)*s(ab,mu,istate,istate)
                        do jstate=1,metat
                           if(jstate.eq.istate)cycle
                           rsmu(jstate,istate)=rsmu(jstate,istate)+
     $                        dint1(b,a,istate)*s(ab,mu,jstate,istate)
                        enddo
                     enddo
                  enddo
               enddo
               do istate=1,metat
                  deno(istate)=w(mu,istate)+
     *            epsnew(ir,istate)+epsnew(is,istate)
               enddo
               do istate=1,metat
                  contri=rsmu(istate,istate)**2/deno(istate)
                  e2(istate,istate)=e2(istate,istate)-contri
                  psi2(istate)=psi2(istate)+contri/deno(istate)
                  do jstate=1,istate
                     if(jstate.eq.istate)cycle
                     contri=rsmu(jstate,istate)*rsmu(istate,istate)
     $                    /deno(istate)
                     e2(jstate,istate)=e2(jstate,istate)-contri
                     contri=rsmu(istate,jstate)*rsmu(jstate,jstate)
     $                    /deno(jstate)
                     e2(istate,jstate)=e2(istate,jstate)-contri
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(s)
      deallocate(dint1)
      return
      end subroutine e2m2cont
c-----------------------------------------------------------
      subroutine e2m2contp(e2,psi2,coef,w,daa,daaoff,nact,ncore,norb
     $     ,nnew,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(nact*(nact+1)/2,*),coef(nact*(nact+1)/2,nact
     $     *(nact+1)/2,metat),daa(metat,nact,nact,nact,nact),
     $     daaoff(metat*(metat-1)/2,nact,nact,nact,nact)
      dimension zgel(*),e2(metat,*),psi2(*)
      allocatable s(:,:,:,:),rsmu(:,:),deno(:)
      integer a,ap,b,bp,ab,abp
      ! Cholesky half-transformed vectors
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable dint1(:,:,:)
      allocate (dint1(nact,nact,metat))
      nvirt=norb-nact-ncore
Cele QD_new

      ncoppie=metat*(metat-1)/2
      nn1=nact*(nact+1)/2
      allocate(s(nn1,nn1,metat,metat))
      allocate(rsmu(metat,metat))
      allocate(deno(metat))
      s = 0.0d0
c--caso r=s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=a,nact
            ab=ab+1
            do mu=1,nn1
               abp=0
               do ap=1,nact
                  do bp=ap,nact
                     abp=abp+1
                     ijcou=0
                     do istate=1,metat
                        coe=coef(abp,mu,istate)
                        if(coe.ne.0.d0)s(ab,mu,istate,istate)=s(ab,mu
     $                       ,istate,istate)+0.5d0*(daa(istate,ap,bp,a,b
     $                       )+daa(istate,ap,bp,b,a))*coe
                        do jstate=1,istate
                           if(jstate.eq.istate)cycle
                           ijcou=ijcou+1
                           coe=coef(abp,mu,istate)
                           if(coe.ne.0.d0)s(ab,mu,jstate,istate)=s(ab,mu
     $                          ,jstate,istate)+0.5d0*(daaoff(ijcou,bp
     $                          ,ap,b,a)+daaoff(ijcou,ap,bp,b,a))
     $                          *coe
                           coe=coef(abp,mu,jstate)
                           if(coe.ne.0.d0)s(ab,mu,istate,jstate)=s(ab,mu
     $                          ,istate,jstate)+0.5d0*(daaoff(ijcou,b
     $                          ,a,ap,bp)+daaoff(ijcou,b,a,bp,ap))*coe
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      do ir=ncore+nact+1,norb
      irv=ir-ncore-nact
      if (Do_Cholesky) then
        call newint_rasb_cholesky(cho_reduced_vec,dint1,
     &    metat,ir,ir,ncore,nact,nvirt)
      else
        call newint_rs_ab(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,irv,irv,dint1)
      end if
       do mu=1,nn1
        call zeroe(rsmu,metat**2)
        ab=0
        do a=1,nact
         ias=itsym(a+ncore)
         do b=a,nact
          ab=ab+1
          ibs=itsym(b+ncore)
          if (ias.ne.ibs) cycle
          do istate=1,metat
          dint2=dint1(a,b,istate)
          if(a.ne.b)dint2=dint2*2.d0
           rsmu(istate,istate)=rsmu(istate,istate)+dint2
     $          *s(ab,mu,istate,istate)
           do jstate=1,metat
            if(jstate.eq.istate)cycle
            rsmu(jstate,istate)=rsmu(jstate,istate)+
     $           dint2*s(ab,mu,jstate,istate)
           enddo
          enddo
         enddo
        enddo
        do istate=1,metat
         deno(istate)=w(mu,istate)+epsnew(ir,istate)+epsnew(ir,istate)
        enddo
        do istate=1,metat
         contri=rsmu(istate,istate)**2/deno(istate)
         e2(istate,istate)=e2(istate,istate)-contri
         psi2(istate)=psi2(istate)+contri/deno(istate)
         do jstate=1,istate
          if(jstate.eq.istate)cycle
          contri=rsmu(jstate,istate)*rsmu(istate,istate)
     $         /deno(istate)
          e2(jstate,istate)=e2(jstate,istate)-contri
          contri=rsmu(istate,jstate)*rsmu(jstate,jstate)
     $         /deno(jstate)
          e2(istate,jstate)=e2(istate,jstate)-contri
         enddo
        enddo
       enddo
      enddo
      deallocate (dint1)
      return
      end subroutine e2m2contp
c-----------------------------------------------------------
      subroutine e2p2cont(e2,psi2,coef,w,daa,daaoff,nact,ncore,
     $  norb,nnew,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *  ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      use info_symmetry
      ! Leon -- for Cholesky
      use nevpt2_cfg, only : Do_Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec

      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(nact**2,*),coef(nact**2,nact**2,metat),
     $     daa(metat,nact,nact,nact,nact),zgel(*),daaoff(metat*(metat-1)
     $     /2,nact,nact,nact,nact),e2(metat,*),psi2(*)
      integer a,ap,b,bp,ab,abp
      allocatable s(:,:,:,:),aijmu(:,:),deno(:),traint(:,:,:)
      ! Leon
      ! traint as in v2mod()
      ! Leon -- for Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      nvirt=norb-nact-ncore
Cele QD_new

      ncoppie=metat*(metat-1)/2
      allocate(s(nact**2,nact**2,metat,metat))
      allocate(aijmu(metat,metat))
      allocate(deno(metat))
      s = 0.0d0

      if(Do_Cholesky) then
        allocate(traint(nact,nact,metat))
        traint = 0.0d0
      end if
c--caso r#s
c---costruzione matrice S
      nact2=nact**2
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            do mu=1,nact2
               abp=0
               do ap=1,nact
                  do bp=1,nact
                     abp=abp+1
                     ijcou=0
                     do istate=1,metat
                        coe=coef(abp,mu,istate)
                        if(coe.ne.0.d0)s(ab,mu,istate,istate)=s(ab,mu
     $                       ,istate,istate)+daa(istate,ap,bp,a,b)
     $                       *coe
                        do jstate=1,istate
                           if(jstate.eq.istate)cycle
                           ijcou=ijcou+1
                           coe=coef(abp,mu,istate)
                           if(coe.ne.0.d0)s(ab,mu,jstate,istate)=s(ab,mu
     $                          ,jstate,istate)+daaoff(ijcou,b,a,bp,ap)
     $                          *coe
                           coe=coef(abp,mu,jstate)
                           if(coe.ne.0.d0)s(ab,mu,istate,jstate)=s(ab,mu
     $                          ,istate,jstate)+daaoff(ijcou,bp,ap,b,a)
     $                          *coe
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      do i=1,ncore
        if (zgel(i)) goto 1
        iis=itsym(i)
         do j=i+1,ncore
           if (zgel(j)) goto 2
            ijs=itsym(j)
            iijs=its(iis,ijs)

            if(Do_Cholesky) then
              call newint_iajb_cholesky(cho_reduced_vec,traint,
     &    metat,i,j,ncore,nact,nvirt)
            end if
            do mu=1,nact2
               call zeroe(aijmu,metat**2)
               ab=0
               do a=1,nact
               ias=itsym(a+ncore)
               iijas=its(iijs,ias)
                  do b=1,nact
                     ab=ab+1
                     ibs=itsym(b+ncore)
                     if (ibs.ne.iijas) cycle
                     do istate=1,metat
                      if (Do_Cholesky) then
                        dint = traint(b,a,istate)
                      else
      call newint_iajb(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,b,j,a,istate,dint)
                      end if
                        aijmu(istate,istate)=aijmu(istate,istate)+
     $                       dint*s(ab,mu,istate,istate)
                        do jstate=1,metat
                           if(jstate.eq.istate)cycle
                        aijmu(jstate,istate)=aijmu(jstate,istate)+
     $                       dint*s(ab,mu,jstate,istate)
                     enddo
                  enddo
               enddo
            enddo
            do istate=1,metat
             deno(istate)=w(mu,istate)-epsnew(i,istate)-epsnew(j,istate)
            enddo
               do istate=1,metat
                  contri=aijmu(istate,istate)**2/deno(istate)
                  e2(istate,istate)=e2(istate,istate)-contri
                  psi2(istate)=psi2(istate)+contri/deno(istate)
                  do jstate=1,istate
                     if(jstate.eq.istate)cycle
                     contri=aijmu(jstate,istate)*aijmu(istate,istate)
     $                    /deno(istate)
                     e2(jstate,istate)=e2(jstate,istate)-contri
                     contri=aijmu(istate,jstate)*aijmu(jstate,jstate)
     $                    /deno(jstate)
                     e2(istate,jstate)=e2(istate,jstate)-contri
                  enddo
               enddo
            enddo
 2       enddo
 1    enddo

      if(allocated(traint))deallocate(traint)
      return
      end subroutine e2p2cont
c-----------------------------------------------------------
      subroutine e2p2contp(e2,psi2,coef,w,daa,daaoff,nact,ncore,norb
     $     ,nnew,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      use info_symmetry
      ! Leon
      ! Cholesky stuff as in e2p2cont()
      use nevpt2_cfg, only : Do_Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(nact*(nact+1)/2,*),coef(nact*(nact+1)/2,nact
     $     *(nact+1)/2,metat),daa(metat,nact,nact
     $     ,nact,nact),zgel(*),daaoff(metat*(metat-1)/2,nact,nact,nact
     $     ,nact),e2(metat,*),psi2(*)
      integer a,ap,b,bp,ab,abp
      allocatable s(:,:,:,:),aijmu(:,:),deno(:),traint(:,:,:)
      ! Leon -- for Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      nvirt=norb-nact-ncore
Cele QD_new

      ncoppie=metat*(metat-1)/2
      nn1=nact*(nact+1)/2
      allocate(s(nn1,nn1,metat,metat))
      allocate(aijmu(metat,metat))
      allocate(deno(metat))
      s = 0.0d0
      if (Do_Cholesky) then
        allocate(traint(nact,nact,metat))
        traint = 0.0d0
      end if

c--caso r=s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=a,nact
            ab=ab+1
            do mu=1,nn1
               abp=0
               do ap=1,nact
                  do bp=ap,nact
                     abp=abp+1
                     ijcou=0
                     do istate=1,metat
                        coe=coef(abp,mu,istate)
                        if(coe.ne.0.d0)s(ab,mu,istate,istate)=s(ab,mu
     $                       ,istate,istate)+0.5d0*(daa(istate,ap,bp,a,b
     $                       )+daa(istate,ap,bp,b,a))*coe
                        do jstate=1,istate
                           if(jstate.eq.istate)cycle
                           ijcou=ijcou+1
                           coe=coef(abp,mu,istate)
                           if(coe.ne.0.d0)s(ab,mu,jstate,istate)=s(ab,mu
     $                          ,jstate,istate)+0.5d0*(daaoff(ijcou,b
     $                          ,a,bp,ap)+daaoff(ijcou,b,a,ap,bp))*coe
                           coe=coef(abp,mu,jstate)
                           if(coe.ne.0.d0)s(ab,mu,istate,jstate)=s(ab,mu
     $                          ,istate,jstate)+0.5d0*(daaoff(ijcou,bp
     $                          ,ap,b,a)+daaoff(ijcou,ap,bp,b,a))
     $                          *coe
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      do i=1,ncore
        if (zgel(i)) goto 1
         if (Do_Cholesky) then
           call newint_iajb_cholesky(cho_reduced_vec,traint,
     &    metat,i,i,ncore,nact,nvirt)
         end if
         do mu=1,nn1
            call zeroe(aijmu,metat**2)
            ab=0
            do a=1,nact
               ab=ab+1
               ias=itsym(a+ncore)
               do istate=1,metat
                if (Do_Cholesky) then
                  dint = traint(a,a,istate)
                else
      call newint_iajb(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,a,i,a,istate,dint)
                end if
                  aijmu(istate,istate)=aijmu(istate,istate)+
     $                 dint*s(ab,mu,istate,istate)
                  do jstate=1,metat
                     if(jstate.eq.istate)cycle
                     aijmu(jstate,istate)=aijmu(jstate,istate)+
     $                    dint*s(ab,mu,jstate,istate)
                  enddo
               enddo
               do b=a+1,nact
                  ab=ab+1
                  ibs=itsym(b+ncore)
                  if (ias.ne.ibs) cycle
                  do istate=1,metat
                    if(Do_Cholesky) then
                      dint = traint(b,a,istate)
                    else
      call newint_iajb(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,b,i,a,istate,dint)
                    end if
                  dint=dint*2.d0
                     aijmu(istate,istate)=aijmu(istate,istate)+
     $                    dint*s(ab,mu,istate,istate)
                     do jstate=1,metat
                        if(jstate.eq.istate)cycle
                        aijmu(jstate,istate)=aijmu(jstate,istate)+
     $                       dint*s(ab,mu,jstate,istate)
                     enddo
                  enddo
               enddo
            enddo
            do istate=1,metat
               deno(istate)=w(mu,istate)-2.d0*epsnew(i,istate)
            enddo
            do istate=1,metat
               contri=aijmu(istate,istate)**2/deno(istate)
               e2(istate,istate)=e2(istate,istate)-contri
               psi2(istate)=psi2(istate)+contri/deno(istate)
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  contri=aijmu(jstate,istate)*aijmu(istate,istate)
     $                 /deno(istate)
                  e2(jstate,istate)=e2(jstate,istate)-contri
                  contri=aijmu(istate,jstate)*aijmu(jstate,jstate)
     $                 /deno(jstate)
                  e2(istate,jstate)=e2(istate,jstate)-contri
               enddo
            enddo
         enddo
 1    enddo
      if(allocated(traint))deallocate(traint)
      return
      end subroutine e2p2contp
c----------------------------------------------------------------
! Leon: introduced also cho_reduced_vec here for Cholesky
! might be subject to a larger restructuring in order not to call the integral transformation again
      subroutine e2p1cont(e2,psi2,coef,w,da,daoff,epsnew,uv,uc,
     *   factor,nact,ncore,nvirt,norb,nnew,metat,nmo,nc,na,nv,nbe,
     *   ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,zgel,
     *   cvirt,cho_reduced_vec)
      use info_symmetry
! Leon -- for Cholesky
      use nevpt2_cfg, only : Do_Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
!
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(nact,*),coef(nact,nact,metat),da(metat,nact
     $     ,nact),zgel(*),daoff(metat*(metat-1)/2,nact,nact),e2(metat,*)
     $     ,psi2(*)
      allocatable s(:,:,:,:),rjimu(:,:),rjimup(:,:),deno(:)
      integer a,ap
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat),
     *       cvirt(nvirt,nvirt,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable traint(:,:,:),traintk(:,:,:),dumai(:,:,:)

! Leon -- for Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
!
      nvirt=norb-nact-ncore
      allocate (traint(nact,nvirt,metat))
      allocate (traintk(nact,nvirt,metat))
      if (.not.Do_Cholesky)allocate (dumai(nact,nvirt,metat))
Cele QD_new

      allocate(s(nact,nact,metat,metat))
      allocate(rjimu(metat,metat))
      allocate(rjimup(metat,metat))
      allocate(deno(metat))
      ncoppie=metat*(metat-1)/2
      s = 0.0d0
c---  costruzione matrice S
      do a=1,nact
       do mu=1,nact
        do ap=1,nact
         ijcou=0
         do istate=1,metat
          if(ap.ne.a)then
             ro=-da(istate,ap,a)
          else
             ro=2.d0-da(istate,ap,a)
          endif
          coe=coef(ap,mu,istate)
          if(coe.ne.0.d0)s(a,mu,istate,istate)=s(a,mu,istate
     $         ,istate)+ro*coe
          do jstate=1,istate-1
           ijcou=ijcou+1
           coe=coef(ap,mu,istate)
           if(coe.ne.0.d0)s(a,mu,jstate,istate)=s(a,mu,jstate
     $          ,istate)-daoff(ijcou,a,ap)*coe
           coe=coef(ap,mu,jstate)
           if(coe.ne.0.d0)s(a,mu,istate,jstate)=s(a,mu,istate
     $          ,jstate)-daoff(ijcou,ap,a)*coe
          enddo
         enddo
        enddo
       enddo
      enddo
c--------------------------------
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      denomin=1.D40
      do i=1,ncore
       if (zgel(i)) cycle
       iis=itsym(i)
       do j=i,ncore
        if (zgel(j)) cycle
        ijs=itsym(j)
        iijs=its(iis,ijs)
        ! Leon
        ! TODO: Restructure the code so that newint_xxxx() do NOT have to be called the 2nd time here
        ! but traints are reused from v_xx()

          if (Do_Cholesky) then
            call newint_airj_cholesky(cho_reduced_vec,traint,
     &        metat,i,j,ncore,nact,nvirt)
            if (j.ne.i)
     &        call newint_airj_cholesky(cho_reduced_vec,traintk,
     &          metat,j,i,ncore,nact,nvirt)
          else
              call newint_airjnew(metat,norb,ncore,nact,nvirt,
     *          factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *          nord_nev2mol,nord_mol2nev,i,j,traint,dumai,cvirt)
              if (j.ne.i)
     *          call newint_airjnew(metat,norb,ncore,nact,nvirt,
     *            factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *            nord_nev2mol,nord_mol2nev,j,i,traintk,dumai,cvirt)
          end if

         do ir=ncore+nact+1,norb
         irv=ir-ncore-nact
         irs=itsym(ir)
         irijs=its(irs,iijs)


         do mu=1,nact
          call zeroe(rjimu,metat**2)
          call zeroe(rjimup,metat**2)
          do a=1,nact
          ias=itsym(a+ncore)
          if (ias.ne.irijs) cycle
           if(i.ne.j)then
              do istate=1,metat
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,a,i,irv,j,istate,traint)
          contj=traint(a,irv,istate)
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,a,j,irv,i,istate,traint1)
          contk=traintk(a,irv,istate)
              cint=contj+contk
              bint=contj-contk
               rjimu(istate,istate)=rjimu(istate,istate)
     $              +cint*s(a,mu,istate,istate)
               rjimup(istate,istate)=rjimup(istate,istate)
     $              +bint*s(a,mu,istate,istate)
               do jstate=1,metat
               if (istate.eq.jstate) cycle
                rjimu(jstate,istate)=rjimu(jstate,istate)
     $               +cint*s(a,mu,jstate,istate)
                rjimup(jstate,istate)=rjimup(jstate,istate
     $               )+bint*s(a,mu,jstate,istate)
               enddo
              enddo
           else
              do istate=1,metat
c      call newint_airj(metat,norb,ncore,nact,nvirt,
c     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
c     *   nord_nev2mol,nord_mol2nev,a,i,irv,i,istate,traint)
              cint=traint(a,irv,istate)
               rjimu(istate,istate)=rjimu(istate,istate)
     $              +cint*s(a,mu,istate,istate)
               do jstate=1,metat
               if (istate.eq.jstate) cycle
                rjimu(jstate,istate)=rjimu(jstate,istate)
     $               +cint*s(a,mu,jstate,istate)
               enddo
              enddo
           endif
          enddo !do a
          do istate=1,metat
           deno(istate)=w(mu,istate)+epsnew(ir,istate)-
     *                  epsnew(i,istate)-epsnew(j,istate)
           if(abs(deno(istate)).lt.denomin)denomin=abs(deno(istate))
          enddo
          if(i.ne.j)then
             do istate=1,metat
              contri=(0.5d0*rjimu(istate,istate)**2+1.5d0
     $             *rjimup(istate,istate)**2)/deno(istate)
              e2(istate,istate)=e2(istate,istate)-contri
              psi2(istate)=psi2(istate)+contri/deno(istate)
              do jstate=1,istate-1
               contri=(0.5d0*rjimu(jstate,istate)*rjimu(istate
     $              ,istate)+1.5d0*rjimup(jstate,istate)
     $              *rjimup(istate,istate))/deno(istate)
               e2(jstate,istate)=e2(jstate,istate)-contri
               contri=(0.5d0*rjimu(istate,jstate)*rjimu(jstate
     $              ,jstate)+1.5d0*rjimup(istate,jstate)
     $              *rjimup(jstate,jstate))/deno(jstate)
               e2(istate,jstate)=e2(istate,jstate)-contri
              enddo
             enddo
          else
             do istate=1,metat
              contri=rjimu(istate,istate)**2/deno(istate)
              e2(istate,istate)=e2(istate,istate)-contri
              psi2(istate)=psi2(istate)+contri/deno(istate)
              do jstate=1,istate-1
               contri=rjimu(jstate,istate)*rjimu(istate
     $              ,istate)/deno(istate)
               e2(jstate,istate)=e2(jstate,istate)-contri
               contri=rjimu(istate,jstate)*rjimu(jstate
     $              ,jstate)/deno(jstate)
               e2(istate,jstate)=e2(istate,jstate)-contri
              enddo
             enddo
          endif
         enddo
 3      enddo
 2     enddo
      enddo
      deallocate(s)
      deallocate(rjimu)
      deallocate(rjimup)
      deallocate(deno)
      deallocate(traint)
      deallocate(traintk)
      if(allocated(dumai))deallocate(dumai)
      return
      end subroutine e2p1cont
c------------------------------------------------------
      subroutine v2mod(ncore,nact,gkp2aa,daat,daaoff,daoff,
     *   epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   cho_reduced_vec)
      use info_symmetry
      ! Leon -- for Cholesky
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common/thresho/thresh
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension gkp2aa(metat,nact,nact,nact,nact)
      dimension daat(metat,nact,nact,nact,nact)
      dimension daaoff(metat*(metat-1)/2,nact,nact,nact,nact)
      dimension daoff(metat*(metat-1)/2,nact,nact)
      allocatable coef(:,:,:),s(:,:),w(:,:),work(:),caux(:,:)
      allocatable daatoff(:,:,:,:,:)
      allocatable dn(:,:),dk(:),den(:),e2(:,:),psi2(:)
      ! Leon -- for Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable cint(:),dint1(:),dint2(:),traint(:,:,:)
      ! Leon
      ! added traint array to call newint_iajb() in batches as in other contributions
      ! for now, it will be only implemented for Cholesky
      nvirt=norb-nact-ncore
      allocate(cint(metat))
      allocate(dint1(metat))
      allocate(dint2(metat))
Cele QD_new

      allocate(dn(metat,metat))
      allocate(dk(metat))
      allocate(den(metat))
      allocate(e2(metat,metat))
      allocate(psi2(metat))
      e2 = 0.0d0

      if(Do_Cholesky) then
        allocate(traint(nact,nact,metat))
        traint=0.0d0
      end if

c
c     t ne u
c

      ncoppie=metat*(metat-1)/2
      do 1800 ir=1,ncore
      if (zgel(ir)) goto 1800
      irsy=itsym(ir)

      do 1700 is=ir,ncore
      if (zgel(is)) goto 1700
      issy=itsym(is)
      irssy=its(irsy,issy)

      if (zverbose) then
      write (6,*) 'Excitations ',ir,is,' ->'
      call flush(6)
      endif

c
c     Si inizia il ciclo sugli attivi
c

      do istate=1,metat
         dn(istate,istate)=0.d0
         dk(istate)=0.d0
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            dn(istate,jstate)=0.d0
            dn(jstate,istate)=0.d0
         enddo
      enddo

      ijcou=0
      if (Do_Cholesky) then
        call newint_iajb_cholesky(cho_reduced_vec,traint,
     &    metat,ir,is,ncore,nact,nvirt)
      end if
      do ic=ncore+1,ncore+nact
         ica=ic-ncore
         do id=ncore+1,ncore+nact
            ida=id-ncore
            icsy=itsym(ic)
            idsy=itsym(id)
            icdsy=its(icsy,idsy)
            isimm=its(icdsy,irssy)
            if (isimm.ne.1) goto 1121
            if (Do_Cholesky) then
              dint2(:)=traint(ica,ida,:)
            else
              do istate=1,metat
      call newint_iajb(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ir,ica,is,ida,istate,dint2(istate))
              enddo
            end if
            do icp=ncore+1,ncore+nact
               icpa=icp-ncore
               do idp=ncore+1,ncore+nact
                  idpa=idp-ncore
                  icpsy=itsym(icp)
                  idpsy=itsym(idp)
                  icpdpsy=its(icpsy,idpsy)
                  isimm=its(icpdpsy,irssy)
                  if (isimm.ne.1) goto 1120
                  do istate=1,metat
                     dum=daat(istate,icpa,idpa,ica,ida)
                     dumk=gkp2aa(istate,icpa,idpa,ica,ida)
                     ! indentation nightmare!!!!!
      if (Do_Cholesky) then
                        dint1(istate)=traint(icpa,idpa,istate)
      else
      call newint_iajb(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ir,icpa,is,idpa,istate,
     *   dint1(istate))
      end if
                     cint(istate)=dint1(istate)*dint2(istate)
                     if (ir.eq.is) then
                        cint(istate)=cint(istate)*0.25d0
                        dum=dum+daat(istate,icpa,idpa,ida,ica)
                        dumk=dumk+gkp2aa(istate,idpa,icpa,ica,ida)
                     endif
                    dn(istate,istate)=dn(istate,istate)+cint(istate)*dum
                     dk(istate)=dk(istate)+cint(istate)*dumk
                     do jstate=1,istate-1
                     ijcou=(istate-1)*(istate-2)/2+jstate
                        dum=ro2toff(icpa,idpa,ica,ida,daaoff,daoff,nact
     $                       ,ijcou,ncoppie)
                        if(ir.eq.is)dum=dum+ro2toff(idpa,icpa,ica
     $                       ,ida,daaoff,daoff,nact,ijcou,ncoppie)
                        dn(jstate,istate)=dn(jstate,istate)+
     *                                    dum*cint(istate)
                        dum=ro2toff(ica,ida,icpa,idpa,daaoff,daoff,nact
     $                       ,ijcou,ncoppie)
                        if(ir.eq.is)dum=dum+ro2toff(ica,ida,idpa
     $                       ,icpa,daaoff,daoff,nact,ijcou,ncoppie)
                        dn(istate,jstate)=dn(istate,jstate)+
     $                                    dum*cint(jstate)
                     enddo
                  enddo
 1120          enddo
            enddo
 1121    enddo
      enddo

      do istate=1,metat
         den(istate)=0.d0
         if(abs(dn(istate,istate)).ge.1.d-14)then
            epsimk=dk(istate)/dn(istate,istate)
            den(istate)=epsnew(ir,istate)+epsnew(is,istate)-epsimk
            co=dn(istate,istate)/den(istate)
            con=co/den(istate)
            e2mp(istate,istate)=e2mp(istate,istate)+co
            psimp(istate)=psimp(istate)+con
            e2mpp(4,istate)=e2mpp(4,istate)+co
            psimpp(4,istate)=psimpp(4,istate)+con
         endif
      enddo
      do istate=1,metat
         e2(istate,istate)=e2mpp(4,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            if(den(istate).ne.0.d0)e2mp(jstate,istate)=e2mp(jstate
     $           ,istate)+dn(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2mp(istate,jstate)=e2mp(istate
     $           ,jstate)+dn(istate,jstate)/den(jstate)
            if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate
     $           ,istate)+dn(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate
     $           ,jstate)+dn(istate,jstate)/den(jstate)
         enddo
      enddo

 1700 enddo
 1800 enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*)' Contribution of V(+2) to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(+2) to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(4,istate)
      enddo
      end if
      deallocate(cint)
      deallocate(dint1)
      deallocate(dint2)
      deallocate(dn)
      deallocate(den)
      if(allocated(traint)) deallocate(traint)

      if(no_pc) goto 999
c---Partially contracted NEV-PT-----------------------------------
      allocate(s(nact**2,nact**2))
      allocate(coef(nact**2,nact**2,metat))
      allocate(w(nact**2,metat))
      allocate(work(3*nact**2))
      allocate(daatoff(ncoppie,nact,nact,nact,nact))
      do icpa=1,nact
         do idpa=1,nact
            do ica=1,nact
               do ida=1,nact
                  do ijcou=1,ncoppie
                     dum=ro2toff(icpa,idpa,ica,ida,daaoff,daoff,nact
     $                    ,ijcou,ncoppie)
                     daatoff(ijcou,icpa,idpa,ica,ida)=dum
                  enddo
               enddo
            enddo
         enddo
      enddo
c---------------
      do istate=1,metat
         call fillcm2(s,daat,nact,istate,metat)
         call fillcm2(coef(1,1,istate),gkp2aa,nact,istate,metat)
      nnew=nact**2
      nact3=nact**2
      info = 0
      call dsyev('V','U',nact3,s,nact3,w(1,istate),work,3*nact3,info)
      if(info.ne.0)then
         print *,'something wrong in diagonalization: s (v(+2))'
         print *,' info = ',info
         stop -777
      endif

      if(zverbose)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i,istate)
         enddo
      print*,'Analisi delle dipendenze lineari:'
      endif
      thr=thresh
      ndep=0
      do i=1,nact3
         if(w(i,istate).lt.0.d0)then
            ndep=ndep+1
      if(zverbose)print '(a,i3,a,f12.8)','negative eigenvalue: eps(',
     *            i,')=',w(i,istate)
         elseif(w(i,istate).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if(zverbose)print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w(1,istate),nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef(1,1,istate),nact3,caux
     $     ,nact3,0.d0,s,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef(1,1,istate),nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      info = 0
      call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work,3
     $     *nact3,info)
      if(info.ne.0)then
        print *,'something wrong in diagonalization: coef (v(+2))'
        print *,' info = ',info
        stop -776
      endif
      if(zverbose)then
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i,istate)
         enddo
         call flush(6)
      endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef(1,1,istate
     $     ),nact3,0.d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call zeroe(coef(1,1,istate),nact3*nact3)
      call copia(coef(1,1,istate),s,nact3,nnew)
      deallocate(caux)
      enddo
      deallocate(s)
      call e2p2cont(e2,psi2,coef,w,daat,daatoff,nact,ncore,norb
     $     ,nnew,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      if(zverbose)then
         print*,'Partially contracted contribution: class (+2): i.ne.j'
         print '(a,(5f20.14))','E2=',(e2(istate,istate),istate=1
     $        ,metat)
         print '(a,(5f20.14))','psi2=',(psi2(istate),istate=1,metat)
      endif
      nact2=nact*(nact+1)/2
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(nact2,nact2,metat))
      allocate(w(nact2,metat))
      allocate(work(1:3*nact2))
      do istate=1,metat
         call fillcm2p(s,daat,nact,istate,metat)
         call fillcm2p(coef(1,1,istate),gkp2aa,nact,istate,metat)
         nact3=nact2
         nnew=nact2
         info = 0
         call dsyev('V','U',nact3,s,nact3,w(1,istate),work,3*nact3,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: s (v(+2))'
           print *,' info = ',info
           stop -666
         endif
         if(zverbose)then
            print*,'matrice metrica diagonalizzata con info=',info
            print*,'autovalori matrice metrica:'
            do i=1,nact3
               print '(i4,f14.8)',i,w(i,istate)
            enddo
         print*,'Analisi delle dipendenze lineari:'
         endif
         thr=thresh
         ndep=0
         do i=1,nact3
            if(w(i,istate).lt.0.d0)then
               ndep=ndep+1
         if(zverbose)print '(a,i3,a,f12.8)','negative eigenvalue: eps(',
     *               i,')=',w(i,istate)
            elseif(w(i,istate).lt.thr)then
               ndep=ndep+1
            endif
         enddo
        if(zverbose)print'(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
         nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
         allocate (caux(1:nact3,1:nnew))
         call caux0p(caux,s,w(1,istate),nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
         call dsymm('L','U',nact3,nnew,1.d0,coef(1,1,istate),nact3,caux
     $        ,nact3,0.d0,s,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
         call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $        ,coef(1,1,istate),nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $        ,3*nact3,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: coef (v(+2))'
           print *,' info = ',info
           stop -665
         endif
         if(zverbose)then
            print*,'diagonalization got info=',info
            print*,'eigenvalues:'
            do i=1,nnew
               print '(a,i4,a,f15.8)','eps(',i,')=',w(i,istate)
            enddo
            call flush(6)
         endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
         call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef(1,1
     $        ,istate),nact3,0.d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
         call zeroe(coef(1,1,istate),nact2*nact2)
         call copia(coef(1,1,istate),s,nact3,nnew)
         deallocate(caux)
      enddo
      call e2p2contp(e2,psi2,coef,w,daat,daatoff,nact,ncore,norb
     $     ,nnew,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      if(zverbose)then
         print*,'Partially contracted contribution: class (-1): i.eq.j'
         print '(a,(5f20.14))','E2=',(e2(istate,istate),istate=1
     $        ,metat)
         print '(a,(5f20.14))','psi2=',(psi2(istate),istate=1,metat)
      endif
      do istate=1,metat
         e2enp(4,istate)=e2enp(4,istate)+e2(istate,istate)
         psienep(4,istate)=psienep(4,istate)+psi2(istate)
         psien(istate)=psien(istate)+psienep(4,istate)
         e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
            e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
         enddo
      enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*)' Contribution of V(+2) to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(+2) to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psienep(4,istate)
      enddo
      end if
      call flush(6)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      deallocate(daatoff)
999   deallocate(psi2)
      return
      end subroutine v2mod
c----------------------------------------------------------------
      subroutine vm2mod(ncore,nact,daa,daaoff,gkm2aa,epsnew,uv,uc,
     *   nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,
     *   cho_reduced_vec)
      use info_symmetry
      ! Leon -- for Cholesky
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common/thresho/thresh
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension daa(metat,nact,nact,nact,nact),daaoff(metat*(metat-1)/2
     $     ,nact,nact,nact,nact)
      dimension gkm2aa(metat,nact,nact,nact,nact)
      allocatable coef(:,:,:),s(:,:),w(:,:),work(:),caux(:,:)
      allocatable dn(:,:),dk(:),den(:),e2(:,:),psi2(:)

      ! Leon -- Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable cint(:),dint1(:,:,:)
      nvirt=norb-nact-ncore
      allocate(cint(metat))
      allocate(dint1(nact,nact,metat))
Cele QD_new

      allocate(dn(metat,metat))
      allocate(dk(metat))
      allocate(den(metat))
      allocate(e2(metat,metat))
      allocate(psi2(metat))
      e2 = 0.0d0
c
      do 1800 ir=ncore+nact+1,norb
      irv=ir-ncore-nact
      irsy=itsym(ir)

      do 1700 is=ir,norb
      isv=is-ncore-nact
      issy=itsym(is)
      irssy=its(irsy,issy)

      if (zverbose) then
      write (6,*) 'Excitations ',' ->',ir,is
      call flush(6)
      endif

c      write (6,*) 'Eccitazione con gli indici',' -> ',ir,is

c
c     Si inizia il ciclo sugli attivi
c
      do istate=1,metat
         dn(istate,istate)=0.d0
         dk(istate)=0.d0
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            dn(istate,jstate)=0.d0
            dn(jstate,istate)=0.d0
         enddo
      enddo

      if(Do_Cholesky) then
       ! TODO: Maybe batching the integrals over virtuals instead of actives would bring better parallelisation
       call newint_rasb_cholesky(cho_reduced_vec,dint1,
     &    metat,is,ir,ncore,nact,nvirt)
      else
        call newint_rs_ab(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,isv,irv,dint1)
      end if

      do ic=ncore+1,ncore+nact
         ica=ic-ncore
         do id=ncore+1,ncore+nact
            ida=id-ncore
            icsy=itsym(ic)
            idsy=itsym(id)
            icdsy=its(icsy,idsy)
            isimm=its(icdsy,irssy)
            if (isimm.ne.1) goto 1121
            do icp=ncore+1,ncore+nact
               icpa=icp-ncore
               do idp=ncore+1,ncore+nact
                  idpa=idp-ncore
                  icpsy=itsym(icp)
                  idpsy=itsym(idp)
                  icpdpsy=its(icpsy,idpsy)
                  isimm=its(icpdpsy,irssy)
                  if (isimm.ne.1) goto 1120
                  ijcou=0
                  do istate=1,metat
                     dum=daa(istate,icpa,idpa,ica,ida)
                     dumk=gkm2aa(istate,icpa,idpa,ica,ida)
                     cint(istate)=dint1(ica,ida,istate)*
     *                            dint1(icpa,idpa,istate)
                     if (ir.eq.is) then
                        cint(istate)=cint(istate)*0.25d0
                        dum=dum+daa(istate,icpa,idpa,ida,ica)
                        dumk=dumk+gkm2aa(istate,idpa,icpa,ica,ida)
                     endif
                    dn(istate,istate)=dn(istate,istate)+dum*cint(istate)
                     dk(istate)=dk(istate)+dumk*cint(istate)
                     do jstate=1,istate-1
                        ijcou=ijcou+1
                        dum=daaoff(ijcou,ica,ida,icpa,idpa)
                        if(ir.eq.is)dum=dum+daaoff(ijcou,ica,ida,idpa
     $                       ,icpa)
                        dn(jstate,istate)=dn(jstate,istate)+
     *                                    dum*cint(istate)
                        dum=daaoff(ijcou,icpa,idpa,ica,ida)
                        if(ir.eq.is)dum=dum+daaoff(ijcou,idpa,icpa,ica
     $                       ,ida)
                        dn(istate,jstate)=dn(istate,jstate)+
     *                                    dum*cint(jstate)
                     enddo
                  enddo
 1120          enddo
            enddo
 1121    enddo
      enddo

      do istate=1,metat
         den(istate)=0.d0
         if(abs(dn(istate,istate)).ge.1.d-14)then
            epsimk=-dk(istate)/dn(istate,istate)
     $
            den(istate)=epsimk-epsnew(ir,istate)-epsnew(is,istate)
            co=dn(istate,istate)/den(istate)
            con=co/den(istate)
            e2mp(istate,istate)=e2mp(istate,istate)+co
            psimp(istate)=psimp(istate)+con
            e2mpp(5,istate)=e2mpp(5,istate)+co
            psimpp(5,istate)=psimpp(5,istate)+con
         endif
      enddo
      do istate=1,metat
         e2(istate,istate)=e2mpp(5,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            if(den(istate).ne.0.d0)e2mp(jstate,istate)=e2mp(jstate
     $           ,istate)+dn(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2mp(istate,jstate)=e2mp(istate
     $           ,jstate)+dn(istate,jstate)/den(jstate)
            if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate
     $           ,istate)+dn(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate
     $           ,jstate)+dn(istate,jstate)/den(jstate)
         enddo
      enddo

 1700 enddo
 1800 enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(-2) to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(-2) to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(5,istate)
      enddo
      end if
      deallocate(cint)
      deallocate(dint1)
      deallocate(dn)
      deallocate(den)

      if(no_pc) goto 999
c---Partially contracted NEV-PT-----------------------------------
      allocate(s(nact**2,nact**2))
      allocate(coef(nact**2,nact**2,metat))
      allocate(w(nact**2,metat))
      allocate(work(3*nact**2))
      do istate=1,metat
         call fillcm2(s,daa,nact,istate,metat)
         call fillcm2(coef(1,1,istate),gkm2aa,nact,istate,metat)
         nnew=nact**2
         nact3=nact**2
         info = 0
         call dsyev('V','U',nact3,s,nact3,w(1,istate),work,3*nact3,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: s (v(-2))'
           print *,' info = ',info
           stop -555
         endif
         if(zverbose)then
            print*,'matrice metrica diagonalizzata con info=',info
            print*,'autovalori matrice metrica:'
            do i=1,nact3
               print '(i4,f14.8)',i,w(i,istate)
            enddo
         print*,'Analisi delle dipendenze lineari:'
         endif
         thr=thresh
         ndep=0
         do i=1,nact3
            if(w(i,istate).lt.0.d0)then
               ndep=ndep+1
         if(zverbose)print '(a,i3,a,f12.8)','negative eigenvalue: eps(',
     *          i,')=',w(i,istate)
            elseif(w(i,istate).lt.thr)then
               ndep=ndep+1
            endif
         enddo
       if(zverbose)print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
         nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
         allocate (caux(1:nact3,1:nnew))
         call caux0p(caux,s,w(1,istate),nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
         call dsymm('L','U',nact3,nnew,1.d0,coef(1,1,istate),nact3,caux
     $        ,nact3,0.d0,s,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
         call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $        ,coef(1,1,istate),nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $        ,3*nact3,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: coef (v(-2))'
           print *,' info = ',info
           stop -554
         endif
         if(zverbose)then
            print*,'diagonalization got info=',info
            print*,'eigenvalues:'
            do i=1,nnew
               print '(a,i4,a,f15.8)','eps(',i,')=',w(i,istate)
            enddo
            call flush(6)
         endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
         call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef(1,1
     $        ,istate),nact3,0.d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
         call zeroe(coef(1,1,istate),nact3*nact3)
         call copia(coef(1,1,istate),s,nact3,nnew)
c     stop 'vm2mod'
         deallocate(caux)
        enddo
        deallocate(s)
      call e2m2cont(e2,psi2,coef,w,daa,daaoff,nact,ncore,norb,nnew
     $     ,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,
     *      cho_reduced_vec)
      if(zverbose)then
         print*,'Partially contracted contribution: class (-2): r.ne.s'
         print '(a,(5f20.14))','E2=',(e2(istate,istate),istate=1
     $        ,metat)
         print '(a,(5f20.14))','psi2=',(psi2(istate),istate=1,metat)
       call flush(6)
      endif
      nact2=nact*(nact+1)/2
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(nact2,nact2,metat))
      allocate(w(nact2,metat))
      allocate(work(1:3*nact2))
      do istate=1,metat
         call fillcm2p(s,daa,nact,istate,metat)
         call fillcm2p(coef(1,1,istate),gkm2aa,nact,istate,metat)
         nact3=nact2
         nnew=nact2
         info = 0
         call dsyev('V','U',nact3,s,nact3,w(1,istate),work,3*nact3,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: s (v(-2))'
           print *,' info = ',info
           stop -553
         endif
         if(zverbose)then
            print*,'matrice metrica diagonalizzata con info=',info
            print*,'autovalori matrice metrica:'
            do i=1,nact3
               print '(i4,f14.8)',i,w(i,istate)
            enddo
         print*,'Analisi delle dipendenze lineari:'
         endif
         thr=thresh
         ndep=0
         do i=1,nact3
            if(w(i,istate).lt.0.d0)then
               ndep=ndep+1
         if(zverbose)print '(a,i3,a,f12.8)','negative eigenvalue: eps(',
     *         i,')=',w(i,istate)
            elseif(w(i,istate).lt.thr)then
               ndep=ndep+1
            endif
         enddo
       if(zverbose)print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
         nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
         allocate (caux(1:nact3,1:nnew))
         call caux0p(caux,s,w(1,istate),nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
         call dsymm('L','U',nact3,nnew,1.d0,coef(1,1,istate),nact3,caux
     $        ,nact3,0.d0,s,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
         call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $        ,coef(1,1,istate),nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $        ,3*nact3,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: coef (v(-2))'
           print *,' info = ',info
           stop -552
         endif
         if(zverbose)then
            print*,'diagonalization got info=',info
            print*,'eigenvalues:'
            do i=1,nnew
               print '(a,i4,a,f15.8)','eps(',i,')=',w(i,istate)
            enddo
            call flush(6)
         endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
         call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef(1,1
     $        ,istate),nact3,0.d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
         call zeroe(coef(1,1,istate),nact3*nact3)
         call copia(coef(1,1,istate),s,nact3,nnew)
         deallocate(caux)
      enddo
      deallocate(s)
      call e2m2contp(e2,psi2,coef,w,daa,daaoff,nact,ncore,norb,nnew
     $     ,zgel,metat,epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,cho_reduced_vec)
      if(zverbose)then
         print*,'Partially contracted contribution: class (-2): r.eq.s'
         print '(a,(5f20.14))','E2=',(e2(istate,istate),istate=1
     $        ,metat)
         print '(a,(5f20.14))','psi2=',(psi2(istate),istate=1,metat)
      endif
      do istate=1,metat
         e2enp(5,istate)=e2enp(5,istate)+e2(istate,istate)
         psienep(5,istate)=psienep(5,istate)+psi2(istate)
         psien(istate)=psien(istate)+psienep(5,istate)
         e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
            e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
         enddo
      enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(-2) to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(-2) to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psienep(5,istate)
      enddo
      end if
      call flush(6)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
999   deallocate(psi2)
      return
      end subroutine vm2mod
c-------------------------------------------------------------
      subroutine provasymm(daa,nact,discrepm,discrep2m)
      implicit real*8 (a-h,o-z)
      real*8 daa(nact,nact,nact,nact)
c--controllo hermitianita
      discrepm=0.d0
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  discrep=abs(daa(i,j,k,l)-daa(k,l,i,j))
                  if(discrep.gt.discrepm)discrepm=discrep
               enddo
            enddo
         enddo
      enddo
      discrep2m=0.d0
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  discrep=abs(daa(i,j,k,l)-daa(j,i,l,k))
                  if(discrep.gt.discrep2m)discrep2m=discrep
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine provasymm
c-------------------------------------------------------------
      subroutine provasymm2(daa,nact,discrepm)
      implicit real*8 (a-h,o-z)
      real*8 daa(nact,nact,nact,nact)
c--controllo hermitianita
      discrepm=0.d0
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  discrep=abs(daa(j,i,k,l)-daa(l,k,i,j))
                  if(discrep.gt.discrepm)discrepm=discrep
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine provasymm2
c-------------------------------------------------------------
      subroutine provasymm3(daa,nact,discrepm)
      implicit real*8 (a-h,o-z)
      real*8 daa(nact,nact,nact,nact,nact,nact)
c--controllo hermitianita
      discrepm=0.d0
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  do m=1,nact
                     do n=1,nact
                        discrep=abs(daa(i,j,k,l,m,n)-daa(l,m,n,i,j,k))
                        if(discrep.gt.discrepm)discrepm=discrep
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine provasymm3
c-------------------------------------------------------------
      subroutine provasymm4(daa,nact,discrepm)
      implicit real*8 (a-h,o-z)
      real*8 daa(nact,nact)
c--controllo hermitianita
      discrepm=0.d0
      do i=1,nact
         do j=1,nact
            discrep=abs(daa(i,j)-daa(j,i))
            if(discrep.gt.discrepm)discrepm=discrep
         enddo
      enddo
      return
      end subroutine provasymm4
c----------------------------------------------------------------
      subroutine fillcm2(s,daa,nact,istate,metat)
      real*8 s(nact**2,nact**2),daa(metat,nact,nact,nact,nact)
c---caso r#s
      ij=0
      do i=1,nact
         do j=1,nact
            ij=ij+1
            kl=0
            do k=1,nact
               do l=1,nact
                  kl=kl+1
                  s(ij,kl)=daa(istate,i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine fillcm2
c----------------------------------------------------------------
      subroutine fillcm2p(s,daa,nact,istate,metat)
      real*8 s(nact*(nact+1)/2,nact*(nact+1)/2),daa(metat,nact,nact,nact
     $     ,nact)
c---caso r=s
      ij=0
      do i=1,nact
         do j=i,nact
            ij=ij+1
            kl=0
            do k=1,nact
               do l=k,nact
                  kl=kl+1
                  s(ij,kl)=daa(istate,i,j,k,l)+daa(istate,i,j,l,k)
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine fillcm2p
c-------------------------------------------------------------
      subroutine fillc0p(s,coef,ga,gd,daa,da,nact,istate,metat)
      real*8 s(2*nact**2,2*nact**2),daa(metat,nact,nact,nact,nact)
      real*8 coef(2*nact**2,2*nact**2),da(metat,nact,nact)
      real*8 ga(metat,nact,nact,nact,nact),gd(metat,nact,nact,nact,nact)
!       real*8 ee2
      integer a,b,ap,bp,ab,abp
      nact2=nact**2
c---parte nord-ovest
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=0
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=ga(istate,b,a,ap,bp)
                  s(ab,abp)=2.d0*ee2(bp,ap,a,b,daa,da,nact,istate,metat)
               enddo
            enddo
         enddo
      enddo
c---parte nord-est
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=nact2
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=-0.5d0*ga(istate,b,a,ap,bp)
                  s(ab,abp)=-ee2(bp,ap,a,b,daa,da,nact,istate,metat)
               enddo
            enddo
         enddo
      enddo
c---parte sud-ovest
      ab=nact2
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=0
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=-0.5d0*ga(istate,b,a,ap,bp)
                  s(ab,abp)=-ee2(bp,ap,a,b,daa,da,nact,istate,metat)
               enddo
            enddo
         enddo
      enddo
c---parte sud-est
      ab=nact2
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=nact2
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=gd(istate,b,a,ap,bp)
                  s(ab,abp)=-daa(istate,bp,a,b,ap)
                  if(a.eq.ap)s(ab,abp)=s(ab,abp)+2.d0*da(istate,bp,b)
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine fillc0p
c----------------------------------------------------------------
      subroutine fillcm1p(s,coef,amat,ro3,daa,da,nact,nact3,istate,metat
     $     )
      real*8 s(nact3,nact3),daa(metat,nact,nact,nact,nact)
      real*8 coef(nact3,nact3),da(metat,nact,nact)
      real*8 ro3(metat,nact,nact,nact,nact,nact,nact)
      real*8 amat(metat,nact,nact,nact,nact,nact,nact)
!       real*8 eee2
      integer a,b,c,ap,bp,cp,abc,abcp
      abc=0
      do a=1,nact
         do b=1,nact
            do c=1,nact
               abc=abc+1
               abcp=0
               do ap=1,nact
                  do bp=1,nact
                     do cp=1,nact
                        abcp=abcp+1
                        coef(abc,abcp)=amat(istate,ap,bp,cp,a,b,c)
                        s(abc,abcp)=eee2(cp,ap,bp,b,a,c,ro3,daa,da,nact
     $                       ,istate,metat)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine fillcm1p
c----------------------------------------------------------------
      subroutine fillcp1p(s,coef,amat,ro3,daa,da,nact,nact3,istate,metat
     $     )
      real*8 s(nact3,nact3),daa(metat,nact,nact,nact,nact)
      real*8 coef(nact3,nact3),da(metat,nact,nact)
      real*8 ro3(metat,nact,nact,nact,nact,nact,nact)
      real*8 amat(metat,nact,nact,nact,nact,nact,nact)
!       real*8 eee2t
      integer a,b,c,ap,bp,cp,abc,abcp
      abc=0
      do a=1,nact
         do b=1,nact
            do c=1,nact
               abc=abc+1
               abcp=0
               do ap=1,nact
                  do bp=1,nact
                     do cp=1,nact
                        abcp=abcp+1
                        coef(abc,abcp)=amat(istate,ap,bp,cp,a,b,c)
                        s(abc,abcp)=eee2t(cp,ap,bp,b,a,c,ro3,daa,da,nact
     $                       ,istate,metat)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine fillcp1p
c---------------------------------------------------------------
      subroutine v0pmod(ncore,nact,nvirt,da,daa,daoff,daaoff,gk0pa,
     *    gk0pd,gk0pf,epsnew,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *    nord_nev2mol,nord_mol2nev,cvirt,ccore,compute_rho1st,
     *    cho_reduced_vec)

      use info_symmetry
      ! Leon -- for Cholesky
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec

      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      logical, intent(in) :: compute_rho1st

      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     *     NORB,NORB2,NOCB,NOCA,
     *     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     *     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     *     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common/thresho/thresh
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension gk0pa(metat,nact,nact,nact,nact),gk0pd(metat,nact,nact
     $     ,nact,nact),gk0pf(metat,nact,nact)
      dimension daa(metat,nact,nact,nact,nact),daaoff(metat*(metat-1)/2
     $     ,nact,nact,nact,nact)
      dimension da(metat,nact,nact),daoff(metat*(metat-1)/2,nact,nact)
      allocatable s(:,:),coef(:,:,:),w(:,:),work(:),caux(:,:),
     $     dn(:,:),dkoo(:),dummy(:,:),dummyk(:),den(:
     $     ),e2(:,:),psi2(:)
Cele QD_new
      allocatable heffd(:,:)
      real *8 epsnew(norb,metat),cvirt(nvirt,nvirt,metat),
     *        ccore(ncore,ncore,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable djm(:,:,:),dkm(:,:,:)

      ! Leon -- Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      !> statement function
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)


      real*8, external    :: ddot
! inline (statement) functions

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &   indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
      nvirt=norb-nact-ncore


      allocate(heffd(norb*(norb+1)/2,metat))
      allocate (djm(nact,nact,metat))
      allocate (dkm(nact,nact,metat))
Cele QD_new

      allocate(dn(metat,metat))
      allocate(dummy(metat,metat))
      allocate(dkoo(metat))
      allocate(dummyk(metat))
      allocate(den(metat))
      allocate(e2(metat,metat))
      allocate(psi2(metat))
      e2 = 0.0d0
      ncoppie=metat*(metat-1)/2
c---- building Dyall's h eff----------

      do ik=1,ncore
         do ir=ncore+nact+1,norb
         irv=ir-ncore-nact
         jkr=indice(ik,ir)
            do istate=1,metat
      call newint_ir(metat,norb,ncore,nact,nevpt_ijkl%oneint,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ik,irv,istate,traint)
            heffd(jkr,istate)=traint
            do j=1,ncore
c      write (6,*) 'vop',ik,ir,istate,j
c      call flush(6)
      if (Do_Cholesky) then
        traint = cholesky_traint(j,j,ik,ir,istate)
        traint1 = cholesky_traint(j,ik,j,ir,istate)
      else
      call newint_ijkr(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,j,j,ik,irv,istate,traint)
      call newint_ijkr(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,j,ik,j,irv,istate,traint1)
      end if
            oint=2.d0*traint-traint1
             heffd(jkr,istate)=heffd(jkr,istate)+oint
            enddo
            enddo
         enddo
      enddo

c

      do 1800 ii=1,ncore
         if (zgel(ii)) goto 1800
         iis=itsym(ii)

         do 1700 ir=ncore+nact+1,norb
            irs=itsym(ir)
            iirs=its(iis,irs)
            irv=ir-nact-ncore

         if (zverbose) then
         write (6,*) 'Excitations ',ii,' ->',ir
         call flush(6)
         endif

c      write (6,*) 'Eccitazione con gli indici',ii,' -> ',ir
c      call flush(6)

c
c     Si inizia il ciclo sugli attivi
c
            do istate=1,metat
               dn(istate,istate)=0.d0
               dkoo(istate)=0.d0
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  dn(istate,jstate)=0.d0
                  dn(jstate,istate)=0.d0
               enddo
            enddo

c      write (6,*) 'Entro J'
c      call flush(6)
       if (Do_Cholesky) then
         call newint_irab_cholesky(cho_reduced_vec,djm,
     &    metat,ii,ir,ncore,nact,nvirt)
         call newint_iarb_cholesky(cho_reduced_vec,dkm,
     &    metat,ii,ir,ncore,nact,nvirt)
       else
        call newint_ir_ab_J(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ii,irv,djm,cvirt,ccore)
      call newint_ir_ab_K(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ii,irv,dkm,cvirt,ccore)
       end if
c      write (6,*) 'Entro K'
c      call flush(6)

c
c

            do ia=ncore+1,ncore+nact
               iaa=ia-ncore
               ias=itsym(ia)
               do ib=ncore+1,ncore+nact
                  iba=ib-ncore
                  ibs=itsym(ib)
                  iabsy=its(ias,ibs)
                  if (its(iabsy,iirs).ne.1) goto 1121
                  do iap=ncore+1,ncore+nact
                     iapa=iap-ncore
                     iaps=itsym(iap)
                     do ibp=ncore+1,ncore+nact
                        ibpa=ibp-ncore
                        ibps=itsym(ibp)
                        iabps=its(iaps,ibps)
                        if (its(iabps,iirs).ne.1) goto 1120

                        do istate=1,metat
                        dj=djm(iaa,iba,istate)
                        djp=djm(iapa,ibpa,istate)
                        dk=dkm(iaa,iba,istate)
                        dkp=dkm(iapa,ibpa,istate)

                        dint1=2.d0*djp*dj-djp*dk-dkp*dj
                        dint2=dkp*dk
                        dint11=djp*dj-0.5d0*djp*dk-0.5d0*dj*dkp
                        dint22=dkp*dk
                           dum1=ee2(ibpa,iapa,iaa,iba,daa,da,nact,istate
     $                          ,metat)
                           dum2=ee2tr(ibpa,iba,iapa,iaa,daa,da,nact
     $                          ,istate,metat)
                           if(iaa.eq.iba) dum2=dum2+da(istate,ibpa,iapa)

                           dn(istate,istate)=dn(istate,istate)+dint1
     $                          *dum1+dint2*dum2

                           dumk1=gk0pa(istate,ibpa,iapa,iaa,iba)
                           dumk2=gk0pd(istate,ibpa,iapa,iaa,iba)
                           dkoo(istate)=dkoo(istate)+dint11*dumk1+dint22
     $                          *dumk2
                           do jstate=1,metat
                              if(jstate.eq.istate)cycle

                  if (jstate.lt.istate) then
                  ijcou=(istate-1)*(istate-2)/2+jstate
                  else
                  ijcou=(jstate-1)*(jstate-2)/2+istate
                  endif
                              dum1=ee2(ibpa,iapa,iaa,iba,daaoff,daoff
     $                             ,nact,ijcou,ncoppie)
                              dum2=ee2tr(ibpa,iba,iapa,iaa,daaoff,daoff
     $                             ,nact,ijcou,ncoppie)
                              if(iaa.eq.iba) dum2=dum2+
     $                             daoff(ijcou ,iapa,ibpa)
                              dn(jstate,istate)=dn(jstate,istate)+dint1
     $                             *dum1+dint2*dum2
                           enddo
                        enddo

 1120                enddo
                  enddo
 1121          enddo
            enddo

            do istate=1,metat
               dummy(istate,istate)=0.d0
               dummyk(istate)=0.d0
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  dummy(istate,jstate)=0.d0
                  dummy(jstate,istate)=0.d0
               enddo
            enddo

            do ia=ncore+1,ncore+nact
               iaa=ia-ncore
               do ib=ncore+1,ncore+nact
                  iba=ib-ncore
                  do istate=1,metat
                  dint1=2.d0*djm(iaa,iba,istate)-dkm(iaa,iba,istate)
                  dint2=dint1*0.5d0
                     dum=da(istate,iba,iaa)
                     dummy(istate,istate)=dummy(istate,istate)+dum*dint1
                     do jstate=1,metat
                        if(jstate.eq.istate)cycle

                  if (jstate.lt.istate) then
                  ijcou=(istate-1)*(istate-2)/2+jstate
                  else
                  ijcou=(jstate-1)*(jstate-2)/2+istate
                  endif
                        dum=daoff(ijcou,iba,iaa)+daoff(ijcou,iaa,iba)
                        dummy(jstate,istate)=dummy(jstate,istate)+dum
     $                       *dint1
                     enddo
                     dumk=gk0pf(istate,iaa,iba)
                     dummyk(istate)=dummyk(istate)+dumk*dint2
                  enddo
               enddo
            enddo

            iir=indice(ir,ii)
            ijcou=0
            do istate=1,metat
              dn(istate,istate)=dn(istate,istate)+2.d0*heffd(iir,istate)
     $              *dummy(istate,istate)+2.d0*heffd(iir,istate)**2
              dkoo(istate)=dkoo(istate)+heffd(iir,istate)*dummyk(istate)
               do jstate=1,istate-1
                  ijcou=ijcou+1
                  dn(jstate,istate)=dn(jstate,istate)+heffd(iir,istate)
     $                 *dummy(jstate,istate)
                  dn(istate,jstate)=dn(istate,jstate)+heffd(iir,jstate)
     $                 *dummy(istate,jstate)
               enddo
            enddo


            do istate=1,metat
               den(istate)=0.d0
               if(abs(dn(istate,istate)).gt.1.d-14)then
                  epsimk=dkoo(istate)/dn(istate,istate)
     $
                  den(istate)=epsnew(ii,istate)-epsnew(ir,istate)-epsimk
                  co=dn(istate,istate)/den(istate)
                  con=co/den(istate)
                  e2mp(istate,istate)=e2mp(istate,istate)+co
                  psimp(istate)=psimp(istate)+con
                  e2mpp(8,istate)=e2mpp(8,istate)+co
                  psimpp(8,istate)=psimpp(8,istate)+con
               endif
            enddo
            do istate=1,metat
               e2(istate,istate)=e2mpp(8,istate)
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  if(den(istate).ne.0.d0)e2mp(jstate,istate)=e2mp(jstate
     $                 ,istate)+dn(jstate,istate)/den(istate)
                  if(den(jstate).ne.0.d0)e2mp(istate,jstate)=e2mp(istate
     $                 ,jstate)+dn(istate,jstate)/den(jstate)
                  if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate
     $                 ,istate)+dn(jstate,istate)/den(istate)
                  if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate
     $                 ,jstate)+dn(istate,jstate)/den(jstate)
               enddo
            enddo
 1700    enddo
 1800 enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(0)'' to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(0)'' to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(8,istate)
      enddo
      end if
      deallocate(djm)
      deallocate(dkm)
      deallocate(dn)
      deallocate(den)

      if(no_pc) goto 999
c---  Partially contracted NEV-PT
      nact2=2*nact**2
      allocate(s(nact2,nact2))
      allocate(coef(nact2,nact2,metat))
      allocate(w(nact2,metat))
      allocate(work(3*nact2))
      call zeroe(w,nact2*metat)
      call zeroe(coef,nact2*nact2*metat)
      do istate=1,metat
         call fillc0p(s,coef(1,1,istate),gk0pa,gk0pd,daa,da,nact,istate
     $        ,metat)
         info = 0
         call dsyev('V','U',nact2,s,nact2,w(1,istate),work,3*nact2,info)
         if(info.ne.0)then
           print *,'something wrong in diagonalization: s (v(0p))'
           print *,' info = ',info
           stop -444
         endif
         if(zverbose)then
            print*,'matrice metrica diagonalizzata con info=',info
            print*,'autovalori matrice metrica:'
            do i=1,nact2
               print '(i2,f14.8)',i,w(i,istate)
            enddo
         print*,'Analisi delle dipendenze lineari:'
         endif
         ndep=0
         thr=thresh
         do i=1,nact2
            if(w(i,istate).lt.0.d0)then
               ndep=ndep+1
         if(zverbose)print '(a,i3,a,f12.8)','negative eigenvalue: eps(',
     *         i,')=',w(i,istate)
            elseif(w(i,istate).lt.thr)then
               ndep=ndep+1
            endif
         enddo
        if(zverbose)print'(a,i3,a)','Ci sono',ndep,' dipendenze lineari'
         nnew=nact2-ndep
c--   costruzione della matrice di trasformazione
         allocate (caux(1:nact2,1:nnew))
         call caux0p(caux,s,w(1,istate),nact2,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
         call dsymm('L','U',nact2,nnew,1.d0,coef(1,1,istate),nact2,caux
     $        ,nact2,0.d0,s,nact2)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
         call dgemm('T','N',nnew,nnew,nact2,1.d0,caux,nact2,s,nact2,0.d0
     $        ,coef(1,1,istate),nact2)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact2,w(1,istate),work
     $        ,3*nact2,info)
         if(info.ne.0)then
            print*,'something wrong in diagonalization: sub v0pmod'
            stop -443
         endif
         if(zverbose)then
            print*,'diagonalization got info=',info
            print*,'eigenvalues:'
            do i=1,nnew
               print '(a,i2,a,f15.8)','eps(',i,')=',w(i,istate)
            enddo
         endif
         call flush(6)
c--   back transformation
c--   moltiplico CAUX(nact2,nnew) per coef(nnew,nnew)
c--   e metto in S(nact2,nnew)
         call dgemm('N','N',nact2,nnew,nnew,1.d0,caux,nact2
     $        ,coef(1,1,istate),nact2,0.d0,s,nact2)
c     print*,'nuovi coefficienti'
c     call matout(nact2,nnew,s,nact2,1.d0)
c--   metto in coef i coefficienti trasformati (for clarity sake)
         call zeroe(coef(1,1,istate),nact2*nact2)
         call copia(coef(1,1,istate),s,nact2,nnew)
c--   inizio il calcolo perturbativo
         deallocate(caux)
      enddo
      deallocate(s)
      call e20pcont(e2,psi2,heffd,coef,w,daa,da,daaoff,daoff
     $     ,nact,ncore,nvirt,norb,nnew,zgel,metat,
     *      epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,
     *      cvirt,ccore,compute_rho1st,cho_reduced_vec)
      do istate=1,metat
         e2enp(8,istate)=e2(istate,istate)
         psienep(8,istate)=psi2(istate)
         psien(istate)=psien(istate)+psienep(8,istate)
         e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
            e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
         enddo
      enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(0)'' to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(0)'' to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psienep(8,istate)
      enddo
      end if
      deallocate(coef)
      deallocate(w)
      deallocate(work)
999   deallocate(psi2)
      deallocate(heffd)
      return
      end subroutine v0pmod
c----------------------------------------------------------------
      subroutine caux0p(caux,s,w,nact2,nnew,thr)
      implicit real*8 (a-h,o-z)
      dimension caux(nact2,nnew),s(nact2,nact2),w(*)
      jnew=0
      do j=1,nact2
         if(w(j).ge.thr)then
            jnew=jnew+1
            val=sqrt(w(j))
            do i=1,nact2
               caux(i,jnew)=s(i,j)/val
            enddo
         endif
      enddo
      return
      end subroutine caux0p
c----------------------------------------------------------------
      subroutine copia(coef,s,nact2,nnew)
      implicit real*8 (a-h,o-z)
      dimension coef(nact2,nact2),s(nact2,nact2)
      do mu=1,nnew
         do i=1,nact2
            coef(i,mu)=s(i,mu)
         enddo
      enddo
      return
      end subroutine copia
c----------------------------------------------------------------
      subroutine e20pcont(e2,psi2,heffd,coef,w,daa,da,daaoff,daoff
     $     ,nact,ncore,nvirt,norb,nnew,zgel,metat,
     *      epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,
     *      cvirt,ccore,compute_rho1st,cho_reduced_vec)

      use info_symmetry
      ! Leon -- for Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      use nevpt2_cfg, only : Do_Cholesky


      implicit real*8 (a-h,o-z)
      logical, intent(in) :: compute_rho1st
      logical*1 zgel
      allocatable rhoir(:,:,:),vaux(:,:),rimu(:,:),rimup(:,:),deno(:)
      allocatable s(:,:,:,:),sp(:,:,:,:),spp(:,:,:,:)
      dimension coef(2*nact**2,2*nact**2,metat),
     $     w(2*nact**2,*),e2(metat,*),psi2(*),
     $     daa(metat,nact,nact,nact,nact),da(metat,nact,nact),
     *     heffd(norb*(norb+1)/2,metat),
     $     zgel(*),daaoff(metat*(metat-1)/2,nact,nact,nact,nact),
     $     daoff(metat*(metat-1)/2,nact,nact)
      real*8 cho_reduced_vec(nchovec,nlvec,metat) ! where to store the vectors
Cele QD_new
      real *8 epsnew(norb,metat),cvirt(nvirt,nvirt,metat),
     *        ccore(ncore,ncore,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable hir(:),raib(:),rabi(:),djm(:,:,:),dkm(:,:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nvirt=norb-nact-ncore
      allocate (hir(metat))
      allocate (raib(metat))
      allocate (rabi(metat))
      allocate (djm(nact,nact,metat))
      allocate (dkm(nact,nact,metat))

Cele QD_new

      nn      =   nact**2
      nn2     = nn*2
      ncoppie = metat*(metat-1)/2
      nvirt   = norb-ncore-nact
      allocate(s(nact**2,nn2,metat,metat))
      allocate(sp(nact**2,nn2,metat,metat))
      allocate(spp(nact**2,nn2,metat,metat))
      allocate(rimu(metat,metat))
      allocate(rimup(metat,metat))
      allocate(deno(metat))
      if(compute_rho1st)then
         open(27,file='RHOFILE',form='UNFORMATTED',status='UNKNOWN')
         rewind 27
         allocate(rhoir(ncore,nvirt,metat))
         allocate(vaux(nn2,metat))
         call zeroe(rhoir,ncore*nvirt*metat)
         call zeroe(vaux,nn2*metat)
      endif
      call zeroe(s,nn*nn2*metat*metat)
      call zeroe(sp,nn*nn2*metat*metat)
      call zeroe(spp,nn*nn2*metat*metat)
c--costruiamo le matrici ausiliarie S(ab,mu),S'(ab,mu) e S''(ab,mu)
c--e il vettore vaux(mu) se si vuole la rho perturbata al primo ordine
      iab=0
      do ia=1,nact
         do ib=1,nact
            iab=iab+1
            iabp=0
            do iap=1,nact
               do ibp=1,nact
                  iabp=iabp+1
                  ijcou=0
                  do istate=1,metat
                     val=ee2(ibp,iap,ia,ib,daa,da,nact,istate,metat)
                     val2=-daa(istate,iap,ib,ia,ibp)
                     if(iap.eq.ia)val2=val2+2.d0*da(istate,ib,ibp)
                     call daxpy(nn2,val,coef(iabp,1,istate),nn2,s(iab,1
     $                    ,istate,istate),nn)
                     call daxpy(nn2,val,coef(iabp+nn,1,istate),nn2
     $                    ,sp(iab,1,istate,istate),nn)
                     call daxpy(nn2,val2,coef(iabp+nn,1,istate),nn2
     $                    ,spp(iab,1,istate,istate),nn)
                     if(compute_rho1st .and. iabp == 1)then
                       do j = 1, nn2
                       vaux(j,istate) =
     &                 vaux(j,istate) + da(istate,ia,ib)
     &                                * (2.0d0*coef(iab,j,istate)
     &                                        -coef(iab+nn,j,istate))
                       end do
                     end if
                     do jstate=1,istate
                        if(jstate.eq.istate)cycle
                        ijcou=ijcou+1
                        val=ee2(ibp,iap,ia,ib,daaoff,daoff,nact,ijcou
     $                       ,ncoppie)
                        val2=-daaoff(ijcou,iap,ib,ia,ibp)
                        if(iap.eq.ia)val2=val2+2.d0*daoff(ijcou,ib
     $                       ,ibp)
                        call daxpy(nn2,val,coef(iabp,1,jstate),nn2
     $                       ,s(iab,1,istate,jstate),nn)
                        call daxpy(nn2,val,coef(iabp+nn,1,jstate),nn2
     $                       ,sp(iab,1,istate,jstate),nn)
                        call daxpy(nn2,val2,coef(iabp+nn,1,jstate),nn2
     $                       ,spp(iab,1,istate,jstate),nn)
                        val=ee2(ib,ia,iap,ibp,daaoff,daoff,nact,ijcou
     $                       ,ncoppie)
                        val2=-daaoff(ijcou,ia,ibp,iap,ib)
                        if(iap.eq.ia)val2=val2+2.d0*daoff(ijcou,ibp
     $                       ,ib)
                        call daxpy(nn2,val,coef(iabp,1,istate),nn2
     $                       ,s(iab,1,jstate,istate),nn)
                        call daxpy(nn2,val,coef(iabp+nn,1,istate),nn2
     $                       ,sp(iab,1,jstate,istate),nn)
                        call daxpy(nn2,val2,coef(iabp+nn,1,istate),nn2
     $                       ,spp(iab,1,jstate,istate),nn)
                     enddo
                  enddo
               enddo
         enddo
      enddo
      enddo
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      do i=1,ncore              !loop sui core
      iis=itsym(i)
         if (zgel(i)) goto 1
         do ir=ncore+nact+1,norb !loop sui virtuali
            irv=ir-ncore-nact
            irs=itsym(ir)
            iirs=its(iis,irs)
            iir=indice(ir,i)
       if (Do_Cholesky) then
         call newint_irab_cholesky(cho_reduced_vec,djm,
     &    metat,i,ir,ncore,nact,nvirt)
         call newint_iarb_cholesky(cho_reduced_vec,dkm,
     &    metat,i,ir,ncore,nact,nvirt)
       else
         call newint_ir_ab_J(metat,norb,ncore,nact,nvirt,
     *    factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *    nord_nev2mol,nord_mol2nev,i,irv,djm,cvirt,ccore)
         call newint_ir_ab_K(metat,norb,ncore,nact,nvirt,
     *    factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *    nord_nev2mol,nord_mol2nev,i,irv,dkm,cvirt,ccore)
        end if
            do mu=1,nn2
               call zeroe(rimu,metat**2)
               call zeroe(rimup,metat**2)
               iab=0
               do ia=1,nact
               ias=itsym(ia+ncore)
               iiras=its(iirs,ias)
                  do ib=1,nact
                     iab=iab+1
                     ibs=itsym(ib+ncore)
                     if (iiras.ne.ibs) cycle

                     do istate=1,metat
                     hir(istate)=heffd(iir,istate)
                        raib(istate)=djm(ia,ib,istate)
                        rabi(istate)=dkm(ia,ib,istate)
                        rimu(istate,istate)=rimu(istate,istate)+(2.d0
     $         *raib(istate)-rabi(istate))*s(iab,mu,istate,istate)+2.d0
     $                *da(istate,ib,ia)*hir(istate)*coef(iab,mu,istate)
                        rimup(istate,istate)=rimup(istate,istate)-
     *                raib(istate)*sp(iab,mu,istate,istate)+rabi(istate)
     *           *spp(iab,mu,istate,istate)-hir(istate)*da(istate,ib,ia)
     $                       *coef(iab+nn,mu,istate)
                        do jstate=1,istate-1
                        ijcou=(istate-1)*(istate-2)/2+jstate
                           rimu(jstate,istate)=rimu(jstate,istate)+(2.d0
     $          *raib(istate)-rabi(istate))*s(iab,mu,jstate,istate)+2.d0
     $                     *daoff(ijcou,ib,ia)*hir(istate)*coef(iab,mu
     $                          ,istate)
                           rimup(jstate,istate)=rimup(jstate,istate)
     $               -raib(istate)*sp(iab,mu,jstate,istate)+rabi(istate)
     $                          *spp(iab,mu,jstate,istate)-hir(istate)
     $                          *daoff(ijcou,ib,ia)*coef(iab+nn,mu
     $                          ,istate)
                           rimu(istate,jstate)=rimu(istate,jstate)+(2.d0
     $          *raib(jstate)-rabi(jstate))*s(iab,mu,istate,jstate)+2.d0
     $              *daoff(ijcou,ia,ib)*hir(jstate)*coef(iab,mu,jstate)
                           rimup(istate,jstate)=rimup(istate,jstate)
     $               -raib(jstate)*sp(iab,mu,istate,jstate)+rabi(jstate)
     $                          *spp(iab,mu,istate,jstate)-hir(jstate)
     $                          *daoff(ijcou,ia,ib)*coef(iab+nn,mu
     $                          ,jstate)

                        enddo
                     enddo
                  enddo
               enddo
               do istate=1,metat
                  deno(istate)=w(mu,istate)+
     *             epsnew(ir,istate)-epsnew(i,istate)
!                   if (dabs(deno(istate))<0.01) then
!                     print *,"warning, denominator for psien < 0.01 (",
!      *             deno(istate),")"
!                     print *,"w(",mu,",",istate,")=",w(mu,
!      *             istate),"epsnew(,",ir,",",istate,")=",epsnew(ir,
!      *             istate),"epsnew(,",i,",",istate,")=",epsnew(i,
!      *             istate)
!                     print *,"rimu=",rimu(istate,istate)," rimup=",
!      *             rimup(istate,istate)
!                   end if
               enddo
               do istate=1,metat
                  co=rimu(istate,istate)+rimup(istate,istate)
                  cirmu=co/deno(istate)
                  cont=co*cirmu
                  e2(istate,istate)=e2(istate,istate)-cont

!                   if (dabs(deno(istate))<0.01) then
!                     print *,"psi^2 before:",psi2(istate)
!                   end if
                  psi2(istate)=psi2(istate)+cirmu**2

!                   if (dabs(deno(istate))<0.01) then
!                     print *,"psi^2 after:",psi2(istate)
!                   end if

                  if(compute_rho1st)then
                     ivirt=ir-ncore-nact
                     rhoir(i,ivirt,istate)=rhoir(i,ivirt,istate)
     &                                    -cirmu*vaux(mu,istate)
                  endif

                  do jstate=1,istate-1
                     contri=(rimu(jstate,istate)+rimup(jstate,istate))
     $                    *(rimu(istate,istate)+rimup(istate,istate))
     $                    /deno(istate)
                     e2(jstate,istate)=e2(jstate,istate)-contri
                     contri=(rimu(istate,jstate)+rimup(istate,jstate))
     $                    *(rimu(jstate,jstate)+rimup(jstate,jstate))
     $                    /deno(jstate)
                     e2(istate,jstate)=e2(istate,jstate)-contri
                  enddo
               enddo
            enddo
         enddo
 1    enddo
      if(compute_rho1st)then
        do istate = 1, metat
          write(27)((rhoir(i,ivirt,istate),i=1,ncore),ivirt=1,nvirt)
        end do
        deallocate(rhoir,vaux)
      end if
      deallocate(s)
      deallocate(sp)
      deallocate(spp)
      deallocate(djm)
      deallocate(dkm)
      deallocate(rabi)
      deallocate(raib)
      return
      end subroutine e20pcont
c----------------------------------------------------------------
      subroutine e21pcont(e2,psi2,heffd,coef,w,ro3,daa,da,ro3off,daaoff
     $     ,daoff,nact,ncore,norb,nnew,zgel,metat,epsnew,uv,uc,factor,
     * nmo,nc,na,nv,nbe,ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,
     * compute_rho1st,cho_reduced_vec)
      use info_symmetry
      ! Leon -- for Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      use nevpt2_cfg, only : Do_Cholesky
      implicit real*8 (a-h,o-z)
      logical*1 zgel
      logical, intent(in) :: compute_rho1st
      dimension coef(nact**3,nact**3,metat),
     $     w(nact**3,*),ro3(nact,nact,nact,nact,nact
     $     ,nact),daa(nact,nact,nact,nact),da(nact,nact)
     *     ,heffd(norb*(norb+1)/2,metat),zgel(*),e2(metat,*),psi2(*)
      real*8, intent(in) :: ro3off(metat*(metat-1)/2,nact,nact,nact,
     &                                               nact,nact,nact)
      real*8, intent(in) :: daaoff(metat*(metat-1)/2,nact,nact,
     &                                               nact,nact)
      real*8, intent(in) :: daoff(metat*(metat-1)/2,nact,nact)

      real*4 tarray(2)
      logical*1 zcond
      allocatable rhoid(:,:,:),vaux(:,:,:),s(:,:,:,:),sp(:,:,:,:),
     &            aimu(:,:),deno(:)
      real*8 ::   myval
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable baic(:),hia(:),dint1(:,:,:,:)

      real*8 cho_reduced_vec(nchovec,nlvec,metat) ! Cholesky vectors
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nvirt=norb-nact-ncore
      allocate (baic(metat))
      allocate (hia(metat))
      allocate (dint1(nact,nact,nact,metat))
Cele QD_new

      nact3=nact**3
      nact3s=nact3**2
      metat2=metat**2
      ncoppie=metat*(metat-1)/2
      allocate(s(nact3,nact3,metat,metat))
      allocate(sp(nact3,nact,metat,metat))
      allocate(deno(metat))
      allocate(aimu(metat,metat))
      s  = 0.0d0; sp = 0.0d0
      if(compute_rho1st)then
        allocate(rhoid(ncore,nact,metat))
        allocate(vaux(nact,nact3,metat))
        rhoid = 0.0d0; vaux = 0.0d0
      end if
      iabc=0
      iabcp=0
      ijcou=0
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ia,ib,ic,iap,ibp,icp,iabc,
C$OMP& iabcp,zcond,ijcou,istate,jstate) COLLAPSE(6)
      do ia=1,nact
         do ib=1,nact
            do ic=1,nact
               do iap=1,nact
                  do ibp=1,nact
                     do icp=1,nact
                        ! Leon
                        ! iabc, iabcp and ijcou must be calculated from the loop indices, otherwise the loop cannot be parallelised!
                        iabc=ic+nact*((ib-1)+nact*(ia-1))
                        iabcp=icp+nact*((ibp-1)+nact*(iap-1))
                        zcond=ib.eq.1.and.ic.eq.1
                        ijcou=0
                        do istate=1,metat
                           s(:,iabc,istate,istate)=
     $  s(:,iabc,istate,istate)+eee2t(icp,iap,ibp,ib,ia,ic
     $        ,ro3,daa,da,nact,istate,metat)*coef(:,iabcp,istate)
                           if (zcond) then
                             sp(:,ia,istate,istate)=
     $ sp(:,ia,istate,istate)+ee2tr(icp,iap
     $ ,ibp,ia,daa,da,nact,istate,metat)*coef(:,iabcp,istate)
                           end if
                           do jstate=1,istate-1
                              ijcou=istate*(istate-3)/2+jstate+1
c---val(jstate,istate) wants RO(jstate,istate)!!!
!                               val(istate,jstate)=eee2t(icp,iap,ibp,ib,ia
!      $                             ,ic,ro3off,daaoff,daoff,nact,ijcou
!      $                             ,ncoppie)
                              s(:,iabc,istate,jstate)=
     $ s(:,iabc,istate,jstate)+eee2t(icp,iap,ibp,ib,ia,ic,ro3off,daaoff,
     $   daoff,nact,ijcou,ncoppie)*coef(:,iabcp,jstate)

                              s(:,iabc,jstate,istate)=
     $ s(:,iabc,jstate,istate)+eee2t(ic,ia,ib,ibp,iap,icp,ro3off,daaoff,
     $          daoff,nact,ijcou,ncoppie)*coef(:,iabcp,istate)
                              if (zcond) then
                                sp(:,ia,istate,jstate)=
     $ sp(:,ia,istate,jstate)+ee2tr(icp,iap,ibp,ia,daaoff,daoff,nact,
     $    ijcou,ncoppie)*coef(:,iabcp,jstate)

                                sp(:,ia,jstate,istate)=
     $ sp(:,ia,jstate,istate)+ee2tl(ia,ibp,iap,icp,daaoff,daoff,nact,
     $    ijcou,ncoppie)*coef(:,iabcp,istate)
                              end if
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
      if(compute_rho1st)then
         myval = 0.0d0
         do id=1,nact
            do mu=1,nact3
               iabc=0
               do ia=1,nact
                  do ib=1,nact
                     do ic=1,nact
                        iabc=iabc+1
                        do istate=1,metat
                       myval=ee2tl(id,ib,ia,ic,daa,da,nact,istate,metat)
                       vaux(id,mu,istate) = vaux(id,mu,istate)
     &                                    + coef(mu,iabc,istate)
     &                                    * myval
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c-------------------------------------------------------
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      do i=1,ncore   !loop sugli orbitali di core
        iis=itsym(i)
        if (zgel(i)) goto 1
          if (Do_Cholesky) then
            call newint_iabc_cholesky(cho_reduced_vec,dint1,
     &         metat,i,ncore,nact,nvirt)
          else
            call newint_i_abc(metat,norb,ncore,nact,nvirt,
     *        factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *        nord_nev2mol,nord_mol2nev,i,dint1)
          end if
         do mu=1,nact3
            call zeroe(aimu,metat**2)
            iabc=0
            do ia=1,nact
            ias=itsym(ia+ncore)
               if (iis.eq.ias) then
               ira=indice(i,ia+ncore)
               do istate=1,metat
               hia(istate)=heffd(ira,istate)
                  aimu(istate,istate)=aimu(istate,istate)+hia(istate)*
     *                sp(mu,ia,istate,istate)
                  do jstate=1,istate-1
                    aimu(jstate,istate)=aimu(jstate,istate)+hia(istate)*
     *                sp(mu,ia,jstate,istate)
                     aimu(istate,jstate)=aimu(istate,jstate)+
     *                 hia(jstate)*sp(mu,ia,istate,jstate)
                  enddo
               enddo
               endif !if (iis.eq.ias)
               do ib=1,nact
               ibs=itsym(ib+ncore)
               iabsy=its(ias,ibs)
                  do ic=1,nact
                     ics=itsym(ic+ncore)
                     iabcs=its(iabsy,ics)
                     iabc=iabc+1
                     if (iabcs.ne.iis)cycle
                     do istate=1,metat
                      baic(istate)=dint1(ia,ic,ib,istate)
                        aimu(istate,istate)=aimu(istate,istate)+
     *                   baic(istate)*s(mu,iabc,istate,istate)
                        do jstate=1,istate-1
                           aimu(jstate,istate)=aimu(jstate,istate)+
     *                   baic(istate)*s(mu,iabc,jstate,istate)
                           aimu(istate,jstate)=aimu(istate,jstate)+
     *                   baic(jstate)*s(mu,iabc,istate,jstate)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do istate=1,metat
               deno(istate)=w(mu,istate)-epsnew(i,istate)
            enddo
            do istate=1,metat
               co=aimu(istate,istate)/deno(istate)
               cont=aimu(istate,istate)*co
               e2(istate,istate)=e2(istate,istate)-cont
               psi2(istate)=psi2(istate)+co**2
               if(compute_rho1st)then
                 do id=1,nact
                   rhoid(i,id,istate) = rhoid(i,id,istate)
     &                                - co
     &                                * vaux(id,mu,istate)
                 end do
               end if
               do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  co=aimu(jstate,istate)/deno(istate)
                  cont=aimu(istate,istate)*co
                  e2(jstate,istate)=e2(jstate,istate)-cont
                  co=aimu(istate,jstate)/deno(jstate)
                  cont=aimu(jstate,jstate)*co
                  e2(istate,jstate)=e2(istate,jstate)-cont
               enddo
            enddo
         enddo
 1    enddo
      if(compute_rho1st)then
        do istate = 1, metat
          write(27)((rhoid(i,id,istate),i=1,ncore),id=1,nact)
        end do
        deallocate(rhoid,vaux)
      end if
      deallocate(s)
      deallocate(sp)
      deallocate(deno)
      deallocate(aimu)
      deallocate(dint1)
      return
      end subroutine e21pcont
c------------------------------------------------------------------
      subroutine e2m1pcont(e2,psi2,heffd,coef,w,ro3,daa,da,ro3off,daaoff
     $     ,daoff,nact,ncore,norb,nnew,zgel,metat,
     *      epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *    nord_nev2mol,nord_mol2nev,nsym,compute_rho1st,cho_reduced_vec)
      use info_symmetry
      ! Leon -- for Cholesky
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      use nevpt2_cfg, only : Do_Cholesky
      implicit real*8 (a-h,o-z)
      logical, intent(in) :: compute_rho1st
      logical*1 zgel
      dimension coef(nact**3,nact**3,metat),
     $     w(nact**3,*),ro3(metat,nact,nact
     $     ,nact,nact,nact,nact),daa(metat,nact,nact,nact,nact),da(metat
     $     ,nact,nact),heffd(norb*(norb+1)/2,metat)
     *     ,zgel(*),ro3off(metat*(metat-1)/2,nact
     $     ,nact,nact,nact,nact,nact),daaoff(metat*(metat-1)/2,nact,nact
     $     ,nact,nact),daoff(metat*(metat-1)/2,nact,nact),e2(metat,*)
     $     ,psi2(*)
      real*8, allocatable :: rhodr(:,:,:),vaux(:,:,:),
     &                       s(:,:,:,:),sp(:,:,:,:),rmu(:,:),deno(:)
      real*8              :: myval
      real*4 tarray(2)
      logical*1 zcond

      real*8 cho_reduced_vec(nchovec,nlvec,metat) ! Cholesky vectors
Cele QD_new
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      allocatable rabc(:),hra(:),dint1(:,:,:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nvirt=norb-nact-ncore
      allocate (rabc(metat))
      allocate (hra(metat))
      allocate (dint1(nact,nact,nact,metat))
Cele QD_new

      nact3=nact**3
      nact3s=nact3**2
      nvirt=norb-ncore-nact
      ncoppie=metat*(metat-1)/2
      metat2=metat**2
      allocate(s(nact3,nact3,metat,metat))
      allocate(sp(nact3,nact,metat,metat))
      allocate(deno(metat))
      allocate(rmu(metat,metat))
      s = 0.0d0; sp = 0.0d0

      if(compute_rho1st)then
        allocate(rhodr(nact,nvirt,metat))
        allocate(vaux(nact,nact3,metat))
        rhodr = 0.0d0; vaux = 0.0d0
      end if

      iabc=0
      iabcp=0
      ijcou=0
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ia,ib,ic,iap,ibp,icp,iabc,
C$OMP& iabcp,zcond,ijcou,istate,jstate) COLLAPSE(6)
      do ia=1,nact
       do ib=1,nact
        do ic=1,nact
         do iap=1,nact
          do ibp=1,nact
           do icp=1,nact
            iabc=ic+nact*((ib-1)+nact*(ia-1))
            iabcp=icp+nact*((ibp-1)+nact*(iap-1))
            ijcou=0
            zcond=ib.eq.1.and.ic.eq.1
            do istate=1,metat
!              val(istate,istate)=eee2(icp,iap,ibp,ib,ia,ic
!      $            ,ro3,daa,da,nact,istate,metat)
!              if(zcond)val2(istate,istate)=ee2(icp,iap,ibp
!      $            ,ia,daa,da,nact,istate,metat)
             s(:,iabc,istate,istate)=s(:,iabc,istate,istate)+
     $  eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact,istate,metat)
     $  *coef(:,iabcp,istate)
             if (zcond) then
               sp(:,ia,istate,istate)=sp(:,ia,istate,istate)+
     $  ee2(icp,iap,ibp,ia,daa,da,nact,istate,metat)
     $  *coef(:,iabcp,istate)
             end if
             do jstate=1,istate-1
!               if(jstate.eq.istate)cycle
              ijcou=istate*(istate-3)/2+jstate+1
              s(:,iabc,istate,jstate)=
     $ s(:,iabc,istate,jstate)+eee2(icp,iap,ibp,ib,ia,ic,ro3off,daaoff,
     $   daoff,nact,ijcou,ncoppie)*coef(:,iabcp,jstate)

              s(:,iabc,jstate,istate)=
     $ s(:,iabc,jstate,istate)+eee2(ic,ia,ib,ibp,iap,icp,ro3off,daaoff,
     $          daoff,nact,ijcou,ncoppie)*coef(:,iabcp,istate)
              if (zcond) then
                sp(:,ia,istate,jstate)=
     $ sp(:,ia,istate,jstate)+ee2(icp,iap,ibp,ia,daaoff,daoff,nact,
     $    ijcou,ncoppie)*coef(:,iabcp,jstate)

                sp(:,ia,jstate,istate)=
     $ sp(:,ia,jstate,istate)+ee2(ia,ibp,iap,icp,daaoff,daoff,nact,
     $    ijcou,ncoppie)*coef(:,iabcp,istate)
              end if
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
C$OMP END PARALLEL DO
      if(compute_rho1st)then
        myval = 0.0d0
        do id=1,nact
          do mu=1,nact3
            iabc=0
            do ia=1,nact
              do ib=1,nact
                do ic=1,nact
                  iabc=iabc+1
                  do istate = 1, metat
                    !myval=ee2_dypc(id,ib,ia,ic,daa,da,nact,istate,metat)
                    myval=ee2(id,ib,ia,ic,daa,da,nact,istate,metat)
                    vaux(id,mu,istate) = vaux(id,mu,istate)
     &                                 + coef(mu,iabc,istate)
     &                                 * myval
                  end do
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
c-------------------------------------------------------
      call cpu_time(t4)
      if(print_level > 5)
     &print '(a,f8.2,a)','time to build S',t4-t3,' sec.'
      call zeroe(e2,metat**2)
      call zeroe(psi2,metat)
      do ir=ncore+nact+1,norb   !loop sui virtuali
      irv=ir-ncore-nact
      irs=itsym(ir)

      if (Do_Cholesky) then
        call newint_abcr_cholesky(cho_reduced_vec,dint1,
     &    metat,ir,ncore,nact,nvirt)
      else
      call newint_r_abc(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,irv,dint1)
      end if
       do mu=1,nact3
        call zeroe(rmu,metat**2)
        iabc=0

        do ia=1,nact
         ias=itsym(ia+ncore)
         ira=indice(ir,ia+ncore)

         if (ias.eq.irs) then
         do istate=1,metat
         hra(istate)=heffd(ira,istate)
      rmu(istate,istate)=rmu(istate,istate)+hra(istate)*sp(mu,ia
     $                 ,istate,istate)
      do jstate=1,istate-1
      rmu(jstate,istate)=rmu(jstate,istate)+hra(istate)*sp(mu,ia
     $                    ,jstate,istate)
      rmu(istate,jstate)=rmu(istate,jstate)+hra(jstate)*sp(mu,ia
     $                    ,istate,jstate)
      enddo
      enddo
      endif !if (ias.eq.irs) then

         do ib=1,nact
         ibs=itsym(ib+ncore)
         iabsy=its(ibs,ias)
          do ic=1,nact
          ics=itsym(ic+ncore)
          iabcs=its(iabsy,ics)
           iabc=iabc+1
           if (iabcs.ne.irs) cycle
           do istate=1,metat
       rabc(istate)=dint1(ia,ic,ib,istate)
      rmu(istate,istate)=rmu(istate,istate)+rabc(istate)*s(mu
     $                       ,iabc,istate,istate)
      do jstate=1,istate-1
      rmu(jstate,istate)=rmu(jstate,istate)+rabc(istate)
     $                          *s(mu,iabc,jstate,istate)
      rmu(istate,jstate)=rmu(istate,jstate)+rabc(jstate)
     $                          *s(mu,iabc,istate,jstate)
      enddo
      enddo
          enddo
         enddo
        enddo
        do istate=1,metat
         deno(istate)=w(mu,istate)+epsnew(ir,istate)
        enddo
        do istate=1,metat
         co=rmu(istate,istate)/deno(istate)
         cont=rmu(istate,istate)*co
         e2(istate,istate)=e2(istate,istate)-cont
         psi2(istate)=psi2(istate)+co**2

         if(compute_rho1st)then
           do id=1,nact
             rhodr(id,irv,istate) = rhodr(id,irv,istate)
     &                            - co
     &                            * vaux(id,mu,istate)
           end do
         endif

         do jstate=1,istate
          if(jstate.eq.istate)cycle
          co=rmu(jstate,istate)/deno(istate)
          cont=rmu(istate,istate)*co
          e2(jstate,istate)=e2(jstate,istate)-cont
          co=rmu(istate,jstate)/deno(jstate)
          cont=rmu(jstate,jstate)*co
          e2(istate,jstate)=e2(istate,jstate)-cont
         enddo
        enddo
       enddo
      enddo
      if(compute_rho1st)then
        do istate = 1, metat
          write(27)((rhodr(i,ivirt,istate),i=1,nact),ivirt=1,nvirt)
        end do
        deallocate(rhodr,vaux)
      end if
      deallocate(s)
      deallocate(sp)
      deallocate(deno)
      deallocate(rmu)
      deallocate(hra)
      deallocate(rabc)
      deallocate(dint1)
      return
      end subroutine e2m1pcont
c----------------------------------------------------------------
      subroutine vpupE(ncore,nact,ro3,ro3off,amat,bmat,cmat,dmat,daa,da
     $     ,daaoff,daoff,epsnew,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,compute_rho1st,cho_reduced_vec)

      use info_symmetry
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level

      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      logical, intent(in) :: compute_rho1st

      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     *     NORB,NORB2,NOCB,NOCA,
     *     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     *     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     *     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common/thresho/thresh
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension ro3(metat,nact,nact,nact,nact,nact,nact),
     *     amat(metat,nact,nact,nact,nact,nact,nact),
     *     bmat(metat,nact,nact,nact,nact),
     *     cmat(metat,nact,nact,nact,nact),
     *     daa(metat,nact,nact,nact,nact),
     *     dmat(metat,nact,nact),
     *     da(metat,nact,nact)
      dimension ro3off(metat*(metat-1)/2,nact,nact,nact,nact,nact,nact),
     $     daaoff(metat*(metat-1)/2,nact,nact,nact,nact),
     $     daoff(metat*(metat-1)/2,nact,nact)
      allocatable s(:,:),coef(:,:,:),caux(:,:),w(:,:),work(:),
     $     dumm(:,:),dumk(:),den(:),e2(:,:),psi2(:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
Cele QD_new
      allocatable heffd(:,:),cint(:),cint2(:),dint1(:,:,:,:)
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax

      real*8 cho_reduced_vec(nchovec,nlvec,metat) ! Cholesky vectors
      allocatable traint_arr(:,:,:), traint1_arr(:,:,:)

      nvirt=norb-nact-ncore
      allocate(heffd(norb*(norb+1)/2,metat))
      allocate(dint1(nact,nact,nact,metat))
      allocate(cint(metat))
      allocate(cint2(metat))
Cele QD_new

      allocate(dumm(metat,metat))
      allocate(dumk(metat))
      allocate(den(metat))
      allocate(e2(metat,metat))
      allocate(psi2(metat))

      if(Do_Cholesky) then
        allocate(traint_arr(ncore,nact,metat))
        allocate(traint1_arr(ncore,nact,metat))
      end if
      ncoppie=metat*(metat-1)/2
      e2 = 0.0d0
c---- building Dyall's h eff----------
      do ik=ncore+1,ncore+nact
      ika=ik-ncore
      iks=itsym(ik)
        do ir=1,ncore
        irs=itsym(ir)
        jkr=indice(ik,ir)
          do istate=1,metat
            heffd(jkr,istate)=0.d0
            if (irs.ne.iks) cycle
            call newint_ia(metat,norb,ncore,nact,nevpt_ijkl%oneint,
     *       factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *       nord_nev2mol,nord_mol2nev,ir,ika,istate,traint)
            heffd(jkr,istate)=traint
          end do
        end do
      end do

      do j=1,ncore
        if (Do_Cholesky) then
          call newint_iija_cholesky(cho_reduced_vec,traint_arr,
     &    metat,j,ncore,nact,nvirt)
          call newint_ijia_cholesky(cho_reduced_vec,traint1_arr,
     &    metat,j,ncore,nact,nvirt)
        end if
        do ik=ncore+1,ncore+nact
          ika=ik-ncore
          iks=itsym(ik)
          do ir=1,ncore
            irs=itsym(ir)
            jkr=indice(ik,ir)
            do istate=1,metat
              if (Do_Cholesky) then
                traint=traint_arr(ir,ika,istate)
                traint1=traint1_arr(ir,ika,istate)
              else
                call newint_ijka(metat,norb,ncore,nact,nvirt,
     *       factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *       nord_nev2mol,nord_mol2nev,j,j,ir,ika,istate,traint)
                call newint_ijka(metat,norb,ncore,nact,nvirt,
     *       factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *       nord_nev2mol,nord_mol2nev,ir,j,j,ika,istate,traint1)
              end if
              heffd(jkr,istate)=heffd(jkr,istate)+2.d0*traint-traint1
            enddo
          enddo
        enddo
      enddo

      if (allocated(traint_arr)) deallocate(traint_arr)
      if (allocated(traint1_arr)) deallocate(traint1_arr)

      do 800 ii=1,ncore
       if (zgel(ii)) goto 800
       iis=itsym(ii)

        if (zverbose) then
        write (6,*) 'Excitations ',ii,' ->'
        call flush(6)
        endif
      if (Do_Cholesky) then
        call newint_iabc_cholesky(cho_reduced_vec,dint1,
     &    metat,ii,ncore,nact,nvirt)
      else
        call newint_i_abc(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ii,dint1)
      end if
c      write (6,*) 'Eccitazione con gli indici',ii,' -> '

       do istate=1,metat
        dumm(istate,istate)=0.d0
        dumk(istate)=0.d0
        do jstate=1,istate
         if(jstate.eq.istate)cycle
         dumm(istate,jstate)=0.d0
         dumm(jstate,istate)=0.d0
        enddo
       enddo

       do ia=1,nact
        iac=ia+ncore
        ias=itsym(iac)
        iia=indice(ii,iac)
        do iap=1,nact
         iapc=iap+ncore
         iaps=itsym(iapc)
         iiap=indice(ii,iapc)
         ijcou=0
         do istate=1,metat
          rodum=-da(istate,iap,ia)
          if (iap.eq.ia) rodum=rodum+2.d0
          dumm(istate,istate)=dumm(istate,istate)+
     *      heffd(iiap,istate)*rodum*heffd(iia,istate)
          dumk(istate)=dumk(istate)+heffd(iiap,istate)*
     *         dmat(istate,iap,ia)*heffd(iia,istate)
          do jstate=1,istate-1
           ijcou=ijcou+1
           rodum=-daoff(ijcou,iap,ia)
           dumm(jstate,istate)=dumm(jstate,istate)+
     *         heffd(iiap,istate)*rodum*heffd(iia,istate)
           rodum=-daoff(ijcou,ia,iap)
           dumm(istate,jstate)=dumm(istate,jstate)+
     *         heffd(iiap,jstate)*rodum*heffd(iia,jstate)
          enddo
         enddo
         do ibp=1,nact
          ibpc=ibp+ncore
          ibps=itsym(ibpc)
          iabps=its(iaps,ibps)
          do icp=1,nact
           icpc=icp+ncore
           icps=itsym(icpc)
           iabcps=its(iabps,icps)
           if (iabcps.ne.iis) cycle
           ijcou=0
           do istate=1,metat
            cint(istate)=dint1(iap,icp,ibp,istate)
            dumm(istate,istate)=dumm(istate,istate)+2.d0*cint(istate)*
     *         ee2tr(icp,iap,ibp,ia,daa,da,nact,istate,metat)*
     *         heffd(iia,istate)
            dumk(istate)=dumk(istate)+cint(istate)*
     *           bmat(istate,iap,ibp,icp,ia)*heffd(iia,istate)
            dumk(istate)=dumk(istate)+cint(istate)*
     *           cmat(istate,ia,iap,ibp,icp)*heffd(iia,istate)
            do jstate=1,istate-1
             ijcou=ijcou+1
             dumm2=ee2tr(icp,iap,ibp,ia,daaoff,daoff,nact,ijcou
     $            ,ncoppie)+ee2tl(ia,ibp,iap,icp,daaoff,daoff,nact
     $            ,ijcou,ncoppie)
             dumm(jstate,istate)=dumm(jstate,istate)+cint(istate)
     $            *dumm2*heffd(iia,istate)
             dumm(istate,jstate)=dumm(istate,jstate)+cint(jstate)
     $            *dumm2*heffd(iia,jstate)
            enddo
           enddo
           do ib=1,nact
            ibc=ib+ncore
            ibs=itsym(ibc)
            iabsy=its(ias,ibs)
            do ic=1,nact
             icc=ic+ncore
             ics=itsym(icc)
             iabcs=its(iabsy,ics)
             if (iabcs.ne.iis) cycle
             ijcou=0
             do istate=1,metat
              cint2(istate)=dint1(ia,ic,ib,istate)
              dumm(istate,istate)=dumm(istate,istate)+cint(istate)*
     *          eee2t(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact,istate,metat)
     $             *cint2(istate)
              dumk(istate)=dumk(istate)+cint(istate)*
     *          amat(istate,iap,ibp,icp,ia,ib,ic)*cint2(istate)
              do jstate=1,istate-1
               ijcou=ijcou+1
               dumm(jstate,istate)=dumm(jstate,istate)+cint(istate)
     $              *eee2t(icp,iap,ibp,ib,ia,ic,ro3off,daaoff,daoff
     $              ,nact,ijcou,ncoppie)*cint2(istate)
               dumm(istate,jstate)=dumm(istate,jstate)+cint(jstate)
     $              *eee2t(ic,ia,ib,ibp,iap,icp,ro3off,daaoff,daoff
     $              ,nact,ijcou,ncoppie)*cint2(jstate)
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       do istate=1,metat
        den(istate)=0.d0
        if(dumm(istate,istate).ne.0.d0)then
           epsimk=dumk(istate)/dumm(istate,istate)
           den(istate)=epsnew(ii,istate)-epsimk
           co=dumm(istate,istate)/den(istate)
           con=co/den(istate)
           e2mp(istate,istate)=e2mp(istate,istate)+co
           psimp(istate)=psimp(istate)+con
           e2mpp(6,istate)=e2mpp(6,istate)+co
           psimpp(6,istate)=psimpp(6,istate)+con
        endif
       enddo
       do istate=1,metat
        e2(istate,istate)=e2mpp(6,istate)
        do jstate=1,istate
         if(jstate.eq.istate)cycle
         if(den(istate).ne.0.d0)e2mp(jstate,istate)=e2mp(jstate
     $        ,istate)+dumm(jstate,istate)/den(istate)
         if(den(jstate).ne.0.d0)e2mp(istate,jstate)=e2mp(istate
     $        ,jstate)+dumm(istate,jstate)/den(jstate)
         if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate,istate)
     $        +dumm(jstate,istate)/den(istate)
         if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate,jstate)
     $        +dumm(istate,jstate)/den(jstate)
        enddo
       enddo

 800  enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(+1)'' to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(+1)'' to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(6,istate)
      enddo
      end if
      deallocate(dint1)
      deallocate(cint)
      deallocate(cint2)
      call flush(6)

      if(no_pc) goto 999

      nact3=nact**3
      allocate(s(1:nact3,1:nact3))
      allocate(coef(1:nact3,1:nact3,metat))
      allocate(w(nact3,metat))
      allocate(work(1))
      w = 0; coef = 0
      do istate=1,metat
       call fillcp1p(s,coef(1,1,istate),amat,ro3,daa,da,nact,nact3
     $      ,istate,metat)
       info = 0
       call dsyev('V','U',nact3,s,nact3,w(1,istate),work,-1,info)
       if(info.ne.0)then
           print *,'something wrong in diagonalization: s (v(+1p))'
           print *,' info = ',info
           stop -333
       endif
       lwork=work(1)
       deallocate(work)
       allocate(work(lwork))
       info = 0
       call dsyev('V','U',nact3,s,nact3,w(1,istate),work,lwork,info)
       if(info.ne.0)then
          print*,'something wrong in diagonalization: sub vpupe'
          stop -332
       endif
       if(zverbose)then
          print*,'matrice metrica diagonalizzata con info=',info
          print*,'autovalori matrice metrica:'
          do i=1,nact3
           print '(i4,f14.8)',i,w(i,istate)
          enddo
       print '(a,f8.2,a)','time for diagonalization of S',t-tini
     $      ,' sec.'
       print*,'Analisi delle dipendenze lineari:'
       endif
       thr=thresh
       ndep=0
       do i=1,nact3
        if(w(i,istate).lt.0.d0)then
           ndep=ndep+1
       if(zverbose)  print '(a,i3,a,f12.8)','negative eigenvalue: eps(',
     *      i,')=',w(i,istate)
        elseif(w(i,istate).lt.thr)then
           ndep=ndep+1
        endif
       enddo
       if(zverbose)print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
       nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
       allocate (caux(1:nact3,1:nnew))
       call caux0p(caux,s,w(1,istate),nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
       call dsymm('L','U',nact3,nnew,1.d0,coef(1,1,istate),nact3,caux
     $      ,nact3,0.d0,s,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
       call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $      ,coef(1,1,istate),nact3)
c       t2=etime(tarray)
      if(zverbose) then
      call cpu_time(t2)
       print '(a,f8.2,a)','time for matrix transf.',t2-t,' sec.'
      endif
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
       call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $      ,-1,info)
       lwork=work(1)
       deallocate(work)
       allocate(work(lwork))
       info = 0
       call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $      ,lwork,info)
       if(info.ne.0)then
          print*,'something wrong in diagonalization: sub vpupe'
          stop -331
       endif
       if(zverbose)then
      call cpu_time(t3)
       print '(a,f8.2,a)','time for diagonalization of H ',t3-t2
     $      ,' sec.'
          print*,'diagonalization got info=',info
          print*,'eigenvalues:'
          do i=1,nnew
           print '(a,i4,a,f15.8)','eps(',i,')=',w(i,istate)
          enddo
       endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
       call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef(1,1
     $      ,istate),nact3,0.d0,s,nact3)
c     print*,'nuovi coefficienti'
c     call matout(nact3,nnew,s,nact3,1.d0)
c--   metto in coef i coefficienti trasformati (for clarity sake)
c     call copia(coef,s,nact3,nnew)
c     do mu=1,nnew
       do mu=1,nact3
        do i=1,nact3
         if(mu.le.nnew)then
            coef(mu,i,istate)=s(i,mu)
         else
            coef(mu,i,istate)=0.d0
         endif
        enddo
       enddo
c--   inizio il calcolo perturbativo
c--   costruiamo le matrici ausiliarie S(abc,mu) e S'(a,mu)
       deallocate(caux)
c     allocate(s(1:nnew,1:nact3))
c     allocate(sp(1:nnew,1:nact3))
      enddo
      call e21pcont(e2,psi2,heffd,coef,w,ro3,daa,da,ro3off,daaoff,daoff
     $     ,nact,ncore,norb,nnew,zgel,metat,
     *      epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,compute_rho1st,
     *      cho_reduced_vec)
      do istate=1,metat
       e2enp(6,istate)=e2(istate,istate)
       psienep(6,istate)=psi2(istate)
       psien(istate)=psien(istate)+psienep(6,istate)
       e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
       do jstate=1,istate
        if(jstate.eq.istate)cycle
        e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
        e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
       enddo
      enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(+1)'' to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(+1)'' to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psienep(6,istate)
      enddo
      end if
      call flush(6)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
999   deallocate(heffd)
      deallocate(dumm)
      deallocate(dumk)
      return
      end subroutine vpupE
c----------------------------------------------------------------
      subroutine vmupE(ncore,nact,ro3,ro3off,amat,bmat,cmat,dmat,daa,da
     $     ,daaoff,daoff,epsnew,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,compute_rho1st,cho_reduced_vec)

      use info_symmetry
      use nevpt2_cfg, only : Do_Cholesky, no_pc, MultGroup, print_level
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec

      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      logical, intent(in) :: compute_rho1st

      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      pointer nd,ne,trou,part
      dimension nd(:),ne(:),trou(:),part(:)
      COMMON /DET/ ND,NE,TROU,PART
      pointer :: e2enp, psienep, e2mpp, psimpp
      dimension e2enp(:,:),psienep(:,:),e2mpp(:,:),psimpp(:,:)
      COMMON /PERTURBER/e2enp,psienep,e2mpp,psimpp,zverbose
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common/thresho/thresh
      common /deephole/ ndhb,ndha,zdh
      common/detinf/ntotdet,ntotconf,ntotcap
      dimension ro3(metat,nact,nact,nact,nact,nact,nact),
     *     amat(metat,nact,nact,nact,nact,nact,nact),
     *     bmat(metat,nact,nact,nact,nact),
     *     cmat(metat,nact,nact,nact,nact),
     *     daa(metat,nact,nact,nact,nact),
     *     dmat(metat,nact,nact),
     *     da(metat,nact,nact),
     $     daoff(metat*(metat-1)/2,nact,nact),
     $     daaoff(metat*(metat-1)/2,nact,nact,nact,nact),
     $     ro3off(metat*(metat-1)/2,nact,nact,nact,nact,nact,nact)
      allocatable s(:,:),coef(:,:,:),w(:,:),work(:),caux(:,:),
     $     dumm(:,:),dumk(:),den(:),e2(:,:),psi2(:),provmat(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      ! Leon -- Cholesky
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
Cele QD_new
      allocatable heffd(:,:),cint(:),cint2(:)
      real *8 epsnew(norb,metat)
      integer nsym,nord_nev2mol(norb),nord_mol2nev(norb)
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax

      real*8 dummii,dummij,dummji,dumki
      ! newint_xxxx() is batched here, we use traintX_arr arrays corresponding to traint and traint1 without batches
      allocatable dint1(:,:,:,:),traint_arr(:,:,:),traint1_arr(:,:,:)
      allocatable traint_abbr(:,:,:)
      nvirt=norb-nact-ncore
      allocate(heffd(norb*(norb+1)/2,metat))
      allocate(dint1(nact,nact,nact,metat))
      allocate(cint(metat))
      allocate(cint2(metat))

      if(Do_Cholesky) then
        allocate(traint_arr(nact,nvirt,metat))
        allocate(traint1_arr(nact,nvirt,metat))
        allocate(traint_abbr(nact,nact,metat))
      end if
Cele QD_new

      allocate(dumm(metat,metat))
      allocate(dumk(metat))
      allocate(den(metat))
      allocate(e2(metat,metat))
      allocate(provmat(metat,metat))
      allocate(psi2(metat))
      e2 = 0.0d0
      ncoppie=metat*(metat-1)/2
      metat2=metat**2
c----building Dyall's h eff----------
      ! Leon -- Split the H_Dyall building into three loops to ease batching of integrals
      do ik=ncore+1,ncore+nact
      ia=ik-ncore
      iks=itsym(ik)
         do ir=ncore+nact+1,norb
            irs=itsym(ir)
            irv=ir-ncore-nact
            jkr=indice(ik,ir)
            do istate=1,metat
            if (iks.ne.irs) cycle
      call newint_ar(metat,norb,ncore,nact,nevpt_ijkl%oneint,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ia,irv,istate,traint)
            heffd(jkr,istate)=traint
            end do
         end do
      end do

      do j=1,ncore
        if (Do_Cholesky) then
          call newint_arij_cholesky(cho_reduced_vec,traint_arr,
     &    metat,j,j,ncore,nact,nvirt)
          call newint_airj_cholesky(cho_reduced_vec,traint1_arr,
     &    metat,j,j,ncore,nact,nvirt)
        end if
        do ik=ncore+1,ncore+nact
          ia=ik-ncore
          iks=itsym(ik)
          do ir=ncore+nact+1,norb
            irs=itsym(ir)
            irv=ir-ncore-nact
            jkr=indice(ik,ir)
            do istate=1,metat
              if (iks.ne.irs) cycle
                if (Do_Cholesky) then
                  traint = traint_arr(ia,irv,istate)
                  traint1 = traint1_arr(ia,irv,istate)
                else
                  call newint_arij(metat,norb,ncore,nact,nvirt,
     *             factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *             nord_nev2mol,nord_mol2nev,ia,irv,j,j,istate,traint)
                  call newint_airj(metat,norb,ncore,nact,nvirt,
     *             factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *            nord_nev2mol,nord_mol2nev,ia,j,irv,j,istate,traint1)
                end if
            heffd(jkr,istate)=heffd(jkr,istate)+
     *        2.d0*traint-traint1
            end do
          end do
        end do
      end do

      do ir=ncore+nact+1,norb
        irs=itsym(ir)
        irv=ir-ncore-nact
        if (Do_Cholesky) then
          call newint_abbr_cholesky(cho_reduced_vec,traint_abbr,
     &    metat,ir,ncore,nact,nvirt)
        ! TODO: there should be a way to avoid the newint_abbr_cholesky call here and use the same batch for building the Dyall Hamiltonian here and calculating the SC contributions below
        end if
        do ik=ncore+1,ncore+nact
          ia=ik-ncore
          iks=itsym(ik)
          jkr=indice(ik,ir)
          do istate=1,metat
            if (iks.ne.irs) cycle
            do j=ncore+1,ncore+nact
              ib=j-ncore
              if (Do_Cholesky) then
                traint = traint_abbr(ia,ib,istate)
              else
                call newint_abcr(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ia,ib,ib,irv,istate,traint)
              end if
              heffd(jkr,istate)=heffd(jkr,istate)-traint
            enddo
          enddo
        enddo
      enddo

      if (allocated(traint_arr)) deallocate(traint_arr)
      if (allocated(traint1_arr)) deallocate(traint1_arr)
      if (allocated(traint_abbr)) deallocate(traint_abbr)

      do 800 ir=ncore+nact+1,norb
      irs=itsym(ir)
      irv=ir-ncore-nact

      if (zverbose) then
      write (6,*) 'Excitations ',' ->',ir
      call flush(6)
      endif
      if (Do_Cholesky) then
        call newint_abcr_cholesky(cho_reduced_vec,dint1,
     &    metat,ir,ncore,nact,nvirt)
      else
        call newint_r_abc(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,irv,dint1)
      end if

      dumm=0.0d0
      dumk=0.0d0

      if (nsym.eq.1) then

      ! Parallelised code without symmetry support


       do istate=1,metat
         dumki=0
         dummii=0
         do jstate=1,istate
           dummij=0
           dummji=0
           ijcou=istate*(istate-3)/2+jstate+1
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ib,ic,ibc,icc,ia,iap,iac,
C$OMP& ira,iapc,irap,ibp,ibpc,icp,icpc,dumm2) COLLAPSE(2)
C$OMP& REDUCTION(+:dummii,dummij,dummji,dumki)
           do ia=1,nact
             do iap=1,nact
               iac=ia+ncore
               ira=indice(ir,iac)
               iapc=iap+ncore
               irap=indice(ir,iapc)
               if (istate.eq.jstate) then
                 dummii=dummii+
     *          heffd(irap,istate)*da(istate,iap,ia)*heffd(ira,istate)
                 dumki=dumki+heffd(irap,istate)
     *            *dmat(istate,iap,ia)*heffd(ira,istate)
               else
                 dummji=dummji+
     $            heffd(irap,istate)*daoff(ijcou,ia,iap)*
     $            heffd(ira,istate)
                 dummij=dummij+
     $            heffd(irap,jstate)*daoff(ijcou,iap,ia)*
     $            heffd(ira,jstate)
               endif

               do ibp=1,nact
                 ibpc=ibp+ncore
                 do icp=1,nact
                   icpc=icp+ncore
                   if (istate.eq.jstate) then
                     dummii=dummii+2.d0*
     *                 dint1(iap,icp,ibp,istate)*
     *                 ee2(icp,iap,ibp,ia,daa,da,nact,istate,metat)*
     *                 heffd(ira,istate)
                     dumki=dumki+dint1(iap,icp,ibp,
     *                 istate)*bmat(istate,iap,ibp,icp,ia)*
     *                 heffd(ira,istate)+dint1(iap,icp,ibp,istate)*
     *                 cmat(istate,ia,iap,ibp,icp)*heffd(ira,istate)
                   else
                     dumm2=ee2(icp,iap,ibp,ia,daaoff,daoff,nact,ijcou,
     $                 ncoppie)+ee2(ia,ibp,iap,icp,daaoff,daoff,nact,
     $                 ijcou,ncoppie)
                     dummji=dummji+
     $                dint1(iap,icp,ibp,istate)*dumm2*heffd(ira,istate)
                     dummij=dummij+
     $                dint1(iap,icp,ibp,jstate)*dumm2*heffd(ira,jstate)
                   endif

                   do ib=1,nact
                     do ic=1,nact
                       ibc=ib+ncore
                       icc=ic+ncore
                       if (istate.eq.jstate) then
                         dummii=dummii+
     *                    dint1(iap,icp,ibp,istate)*
     *                    eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,
     *                    da,nact,istate,metat)*
     *                    dint1(ia,ic,ib,istate)
                         dumki=dumki+
     *                    dint1(iap,icp,ibp,istate)*
     *                    amat(istate,iap,ibp,icp,ia,ib,ic)*
     *                    dint1(ia,ic,ib,istate)
                       else
                         dummji=dummji+
     $                    dint1(iap,icp,ibp,istate)*
     $                    eee2(icp,iap,ibp,ib,ia,ic,ro3off,
     $                    daaoff,daoff,nact,ijcou,ncoppie)*
     $                    dint1(ia,ic,ib,istate)
                         dummij=dummij+
     $                    dint1(iap,icp,ibp,jstate)*
     $                    eee2(ic,ia,ib,ibp,iap,icp,ro3off,
     $                    daaoff,daoff,nact,ijcou,ncoppie)*
     $                    dint1(ia,ic,ib,jstate)
                       endif
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
C$OMP END PARALLEL DO
           dumm(istate,jstate)=dummij
           dumm(jstate,istate)=dummji
         enddo
         dumk(istate)=dumki
         dumm(istate,istate)=dummii
       enddo
      else

      ! If symmetry is present, use the old code

        do ia=1,nact
        iac=ia+ncore
        ias=itsym(iac)
        ira=indice(ir,iac)
        do iap=1,nact
          iapc=iap+ncore
          iaps=itsym(iapc)
          irap=indice(ir,iapc)
          ijcou=0
          if ((ias.eq.irs).and.(iaps.eq.irs)) then
          do istate=1,metat
            dumm(istate,istate)=dumm(istate,istate)+heffd(irap,istate)
     *          *da(istate,iap,ia)*heffd(ira,istate)
            dumk(istate)=dumk(istate)+heffd(irap,istate)
     *          *dmat(istate,iap,ia)*heffd(ira,istate)
            do jstate=1,istate
                if(jstate.eq.istate)cycle
                ijcou=ijcou+1
                dumm(jstate,istate)=dumm(jstate,istate)+
     $             heffd(irap,istate)*
     $             daoff(ijcou,ia,iap)*heffd(ira,istate)
                dumm(istate,jstate)=dumm(istate,jstate)+
     $             heffd(irap,jstate)
     $             *daoff(ijcou,iap,ia)*heffd(ira,jstate)
            enddo
          enddo
          endif !if ((ias.eq.irs).and.(iaps.eq.irs))

          do ibp=1,nact
          ibpc=ibp+ncore
          ibps=itsym(ibpc)
          iabps=its(iaps,ibps)
          do icp=1,nact
            icpc=icp+ncore
            icps=itsym(icpc)
            iabcps=its(iabps,icps)
            if (iabcps.ne.irs) cycle
            ijcou=0
            do istate=1,metat
            cint(istate)=dint1(iap,icp,ibp,istate)
              dumm(istate,istate)=dumm(istate,istate)+2.d0*cint(istate)*
     *           ee2(icp,iap,ibp,ia,daa,da,nact,istate,metat)*
     *           heffd(ira,istate)
              dumk(istate)=dumk(istate)+cint(istate)*
     *           bmat(istate,iap,ibp,icp,ia)*heffd(ira,istate)
              dumk(istate)=dumk(istate)+cint(istate)*
     *           cmat(istate,ia,iap,ibp,icp)*heffd(ira,istate)
              do jstate=1,istate
                  if(jstate.eq.istate)cycle
                  ijcou=ijcou+1
                  dumm2=ee2(icp,iap,ibp,ia,daaoff,daoff,nact,ijcou
     $               ,ncoppie)+ee2(ia,ibp,iap,icp,daaoff,daoff,nact
     $               ,ijcou,ncoppie)
                  dumm(jstate,istate)=dumm(jstate,istate)+cint(istate)
     $               *dumm2*heffd(ira,istate)
                  dumm(istate,jstate)=dumm(istate,jstate)+cint(jstate)
     $               *dumm2*heffd(ira,jstate)
              enddo
            enddo

            do ib=1,nact
            ibc=ib+ncore
            ibs=itsym(ibc)
            iabsy=its(ias,ibs)
            do ic=1,nact
              icc=ic+ncore
              ics=itsym(icc)
              iabcs=its(iabsy,ics)
              if (iabcs.ne.irs) cycle
              ijcou=0
              do istate=1,metat
          cint2(istate)=dint1(ia,ic,ib,istate)
                dumm(istate,istate)=dumm(istate,istate)+cint(istate)*
     *          eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact,istate,metat)*
     *          cint2(istate)
                dumk(istate)=dumk(istate)+cint(istate)*
     *          amat(istate,iap,ibp,icp,ia,ib,ic)*cint2(istate)
                do jstate=1,istate-1
                    ijcou=ijcou+1
                    dumm(jstate,istate)=dumm(jstate,istate)+cint(istate)
     $                 *eee2(icp,iap,ibp,ib,ia,ic,ro3off,daaoff,daoff
     $                 ,nact,ijcou,ncoppie)*cint2(istate)
                    dumm(istate,jstate)=dumm(istate,jstate)+cint(jstate)
     $                 *eee2(ic,ia,ib,ibp,iap,icp,ro3off,daaoff,daoff
     $                 ,nact,ijcou,ncoppie)*cint2(jstate)
                enddo
              enddo
          enddo
        enddo
        enddo
        enddo
        enddo
        enddo

      end if

      do istate=1,metat
         den(istate)=0.d0
         if(dumm(istate,istate).ne.0.d0)then
            epsimk=-dumk(istate)/dumm(istate,istate)
     $
            den(istate)=epsimk-epsnew(ir,istate)
            co=dumm(istate,istate)/den(istate)
            con=co/den(istate)
            e2mp(istate,istate)=e2mp(istate,istate)+co
            psimp(istate)=psimp(istate)+con
            e2mpp(7,istate)=e2mpp(7,istate)+co
            psimpp(7,istate)=psimpp(7,istate)+con
         endif
      enddo
      do istate=1,metat
         e2(istate,istate)=e2mpp(7,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            if(den(istate).ne.0.d0)e2mp(jstate,istate)=e2mp(jstate
     $           ,istate)+dumm(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2mp(istate,jstate)=e2mp(istate
     $           ,jstate)+dumm(istate,jstate)/den(jstate)
            if(den(istate).ne.0.d0)e2(jstate,istate)=e2(jstate,istate)
     $           +dumm(jstate,istate)/den(istate)
            if(den(jstate).ne.0.d0)e2(istate,jstate)=e2(istate,jstate)
     $           +dumm(istate,jstate)/den(jstate)
c------------------
         enddo
      enddo

 800  enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(-1)'' to Heff (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(-1)'' to the norm (SC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psimpp(7,istate)
      enddo
      end if
      call flush(6)
      deallocate(provmat)
      deallocate(dint1)
      deallocate(cint)
      deallocate(cint2)
      if(no_pc) goto 999
c--------------------
c--   Partially contracted NEV-PT class (-1')
      nact3=nact**3
      allocate(s(nact3,nact3))
      allocate(coef(1:nact3,1:nact3,metat))
      allocate(w(nact3,metat))
      allocate(work(1))
      w = 0.0d0; coef = 0.0d0
      do istate=1,metat
         call fillcm1p(s,coef(1,1,istate),amat,ro3,daa,da,nact,nact3
     $        ,istate,metat)
         info = 0
         call dsyev('V','U',nact3,s,nact3,w(1,istate),work,-1,info)
         lwork=work(1)
         deallocate(work)
         allocate(work(lwork))
         info = 0
         call dsyev('V','U',nact3,s,nact3,w(1,istate),work,lwork,info)
         if(info.ne.0)then
            print*,'something wrong in diagonalization: sub vmupe'
            stop -330
         endif
         if(zverbose)then
            print*,'matrice metrica diagonalizzata con info=',info
            print*,'autovalori matrice metrica:'
            do i=1,nact3
               print '(i4,f14.8)',i,w(i,istate)
            enddo
         print '(a,f8.2,a)','time for diagonalization of S',t-tini
     $        ,' sec.'
         print*,'Analisi delle dipendenze lineari:'
         endif
         thr=thresh
         ndep=0
         do i=1,nact3
            if(w(i,istate).lt.0.d0)then
               ndep=ndep+1
       if(zverbose)print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,
     *      ')=',w(i,istate)
            elseif(w(i,istate).lt.thr)then
               ndep=ndep+1
            endif
         enddo
       if(zverbose)print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
         nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
         allocate (caux(1:nact3,1:nnew))
         call caux0p(caux,s,w(1,istate),nact3,nnew,thr)
c         if(zverbose)print*,'check: final jnew is',jnew
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
         call dsymm('L','U',nact3,nnew,1.d0,coef(1,1,istate),nact3,caux
     $        ,nact3,0.d0,s,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
         call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $        ,coef(1,1,istate),nact3)
c         t2=etime(tarray)
       if(zverbose) then
      call cpu_time(t2)
         print '(a,f8.2,a)','time for matrix transf.',t2-t,' sec.'
      endif
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $        ,-1,info)
         lwork=work(1)
         deallocate(work)
         allocate(work(lwork))
         info = 0
         call dsyev('V','U',nnew,coef(1,1,istate),nact3,w(1,istate),work
     $        ,lwork,info)
         if(info.ne.0)then
            print*,'something wrong in diagonalization: sub vmupe'
            stop -330
         endif
c         t3=etime(tarray)
         if(zverbose)then
      call cpu_time(t3)
         print '(a,f8.2,a)','time for diagonalization of H ',t3-t2
     $        ,' sec.'
            print*,'diagonalization got info=',info
            print*,'eigenvalues:'
            do i=1,nnew
               print '(a,i4,a,f15.8)','eps(',i,')=',w(i,istate)
            enddo
            call flush(6)
         endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
         call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef(1,1
     $        ,istate),nact3,0.d0,s,nact3)
c     print*,'nuovi coefficienti'
c     call matout(nact3,nnew,s,nact3,1.d0)
c--   metto in coef i coefficienti trasformati (for clarity sake)
c     call copia(coef,s,nact3,nnew)
c         do mu=1,nnew
         do mu=1,nact3
            do i=1,nact3
               if(mu.le.nnew)then
                  coef(mu,i,istate)=s(i,mu)
               else
                  coef(mu,i,istate)=0.d0
               endif
            enddo
         enddo
c--   inizio il calcolo perturbativo
c--   costruiamo le matrici ausiliarie S(abc,mu) e S'(a,mu)
         deallocate(caux)
      enddo
      deallocate(s)
      call e2m1pcont(e2,psi2,heffd,coef,w,ro3,daa,da,ro3off
     $     ,daaoff,daoff,nact,ncore,norb,nnew,zgel,metat,
     *      epsnew,uv,uc,factor,nmo,nc,na,nv,nbe,
     *      ncmax,nvmax,nord_nev2mol,nord_mol2nev,nsym,compute_rho1st,
     *      cho_reduced_vec)
      do istate=1,metat
         e2enp(7,istate)=e2(istate,istate)
         psienep(7,istate)=psi2(istate)
         psien(istate)=psien(istate)+psienep(7,istate)
         e2en(istate,istate)=e2en(istate,istate)+e2(istate,istate)
         do jstate=1,istate
            if(jstate.eq.istate)cycle
            e2en(istate,jstate)=e2en(istate,jstate)+e2(istate,jstate)
            e2en(jstate,istate)=e2en(jstate,istate)+e2(jstate,istate)
         enddo
      enddo
      if(print_level > 0)then
      write (6,*)
      write (6,*) ' Contribution of V(-1)'' to Heff (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(I2,2x,10F15.10)')  MultGroup%State(istate),
     &                    (e2(istate,jstate),
     *                    jstate=1,metat)
      enddo
      write (6,*)
      write (6,*) ' Contribution of V(-1)'' to the norm (PC)'
      write (6,*)
      do istate=1,metat
      write (6,'(''  State '',I2,''  norm ='',F15.10)')
     &             MultGroup%State(istate),
     *             psienep(7,istate)
      enddo
      end if
      call flush(6)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
999   deallocate(heffd)
      deallocate(dumm)
      deallocate(dumk)
      return
      end subroutine vmupE
!-------------------------------------------------------------------

      function eee2(a,b,c,d,e,f,daaa,taa,dal,nact,istate,metat) ! index swap compared to eexx_utils
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f
      dimension daaa(metat,nact,nact,nact,nact,nact,nact),taa(metat,nact
     $     ,nact,nact,nact),dal(metat,nact,nact)
c--the same as eee but with the spinless density matrices as input
      eee2=daaa(istate,b,d,f,a,c,e)
      if(e.eq.d)eee2=eee2+taa(istate,b,f,a,c)
      if(e.eq.b)eee2=eee2+taa(istate,f,d,a,c)
      if(b.eq.c)eee2=eee2+taa(istate,d,f,a,e)
      if(b.eq.c.and.e.eq.d)eee2=eee2+dal(istate,f,a)
      return
      end function eee2
c-------------------------------------------------------------------
      function ee2(a,b,c,d,taa,dal,nact,istate,metat) ! index swap compared to eexx_utils
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      ee2=taa(istate,b,d,a,c)
      if(b.eq.c)ee2=ee2+dal(istate,d,a)
      return
      end function ee2
c-------------------------------------------------------------------
      function ee2_dypc(a,b,c,d,taa,dal,nact,istate,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,istate
      dimension taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      ee2_dypc=taa(istate,a,c,b,d)
      if(b.eq.c)ee2_dypc=ee2_dypc+dal(istate,a,d)
      return
      end function ee2_dypc
c-------------------------------------------------------------------
      function eee2t(a,b,c,d,e,f,daaa,taa,dal,nact,istate,metat) ! same as in eexxx_utils
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f
      dimension daaa(metat,nact,nact,nact,nact,nact,nact),taa(metat,nact
     $     ,nact,nact,nact),dal(metat,nact,nact)
      eee2t=-eee2(a,b,d,c,e,f,daaa,taa,dal,nact,istate,metat)
      if(c.eq.d)eee2t=eee2t+2.d0*ee2(a,b,e,f,taa,dal,nact,istate,metat)
      return
      end function eee2t
c-------------------------------------------------------------------
      function ee2tr(a,b,c,d,taa,dal,nact,istate,metat)  ! index swap compared to eexx_utils
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      ee2tr=-ee2(a,b,d,c,taa,dal,nact,istate,metat)
      if(c.eq.d)ee2tr=ee2tr+2.d0*dal(istate,b,a)
      return
      end function ee2tr
c-------------------------------------------------------------------
      function ee2tl(a,b,c,d,taa,dal,nact,istate,metat)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      ee2tl=-ee2(b,a,c,d,taa,dal,nact,istate,metat)
      if(a.eq.b)ee2tl=ee2tl+2.d0*dal(istate,d,c)
      return
      end function ee2tl
!-------------------------------------------------------------------
      function ro2toff(ap,bp,a,b,d2,d1,nact,ijcou,ncoppie)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,ap,bp,cp
      dimension d2(ncoppie,nact,nact,nact,nact),d1(ncoppie
     $     ,nact,nact)
      ro2toff=d2(ijcou,ap,bp,a,b)
      if(ap.eq.a)ro2toff=ro2toff-2.d0*d1(ijcou,bp,b)
      if(bp.eq.b)ro2toff=ro2toff-2.d0*d1(ijcou,ap,a)
      if(ap.eq.b)ro2toff=ro2toff+d1(ijcou,bp,a)
      if(bp.eq.a)ro2toff=ro2toff+d1(ijcou,ap,b)
      return
      end function ro2toff
!-------------------------------------------------------------------

      subroutine zeroe(a,n)
      real*8 a(*)
      do i=1,n
         a(i)=0.d0
      enddo
      return
      end subroutine zeroe
!-------------------------------------------------------------------

      subroutine trasforma(norb,rho,nsym)
      use info_symmetry
      implicit real*8 (a-h,o-z)
      dimension rho(norb,norb)
      dimension isymtot(nsym),nbeg(nsym),iperm(norb),jperm(norb)
      dimension a(norb,norb)
      isymtot=0
      do i=1,norb
       do isy=1,nsym
        if(itsym(i).eq.isy)then
           isymtot(isy)=isymtot(isy)+1
           exit
        endif
       enddo
      enddo
c      print*,'isymtot',isymtot
      nbeg(1)=0
      do isy=2,nsym
       nbeg(isy)=nbeg(isy-1)+isymtot(isy-1)
      enddo
c      print*,'nbeg',nbeg
      do isy=1,nsym
       j=0
       do i=1,norb
        if(itsym(i).eq.isy)then
           j=j+1
           jj=nbeg(isy)+j
           iperm(i)=jj
           jperm(jj)=i  ! inverse of iperm
        endif
       enddo
      enddo
      do i=1,norb
       k=iperm(i)
       do j=1,norb
        l=iperm(j)
        a(i,j)=rho(k,l)
       enddo
      enddo
      rho=a
      return
      end subroutine trasforma
Cele QD_new
c----------------------------------------------------
      subroutine newint_irsjnew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,j,cvirt,dumai,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none
      integer metat,i,j,ii,ij,ik,il,irv,isv,it,itv,ir
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,kmol,lmol,is,js,ks,ls
      integer inev,jnev,knev,lnev
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 factor,dum1,dum12,dum123
      real*8 dumai(nvirt,nvirt,metat),traint(nvirt,nvirt,metat),
     *       cvirt(nvirt,nvirt,metat)

!     !> Transformation of integrals of the type <ir|sj> [in the other formalism (ir,sj)]

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      js=itsym(j)
      jmol=nord_nev2mol(j)-nbe(js)



      do istate=1,metat
      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      dum1=uc(ii,imol,is,istate)
      do ij=1,nc(js)
      jnev=nord_mol2nev(ij+nbe(js))
      dum12=dum1*uc(ij,jmol,js,istate)
      do irv=1,nvirt
      ir=irv+ncore+nact
      do itv=1,nvirt
      it=itv+ncore+nact
      traint(itv,irv,istate)=traint(itv,irv,istate)+
     *    dum12*ai(inev,it,ir,jnev)
      enddo
      enddo
      enddo
      enddo
      enddo

      do istate=1,metat
      call dgemm('N','N',nvirt,nvirt,nvirt,1.d0,traint(1,1,istate),
     *  nvirt,cvirt(1,1,istate),nvirt,0.d0,dumai(1,1,istate),nvirt)
      call dgemm('T','N',nvirt,nvirt,nvirt,1.d0,cvirt(1,1,istate),nvirt,
     *     dumai(1,1,istate),nvirt,0.d0,traint(1,1,istate),nvirt)
      enddo
      return
      end subroutine newint_irsjnew
c----------------------------------------------------
      subroutine newint_airj(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ia,i,ir,j,istate,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,j,ii,ij,il,jt,kt,lt,ir,ia
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,lmol,is,js,ks,ls
      integer inev,jnev,lnev
      real*8  factor,traint,dum12
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <ar|ij> [in the other formalism (ai,rj)]
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      js=itsym(j)
      jmol=nord_nev2mol(j)-nbe(js)

      kt=ia+ncore
      ks=itsym(kt)

      lt=ir+ncore+nact
      ls=itsym(lt)
      lmol=nord_nev2mol(lt)-nbe(ls)-nc(ls)-na(ls)

      traint=0.d0

      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      do ij=1,nc(js)
      jnev =nord_mol2nev(ij+nbe(js))
      dum12=uc(ii,imol,is,istate)*uc(ij,jmol,js,istate)
      do il=1,nv(ls)
      lnev=nord_mol2nev(il+nbe(ls)+nc(ls)+na(ls))
      traint=traint+dum12*uv(il,lmol,ls,istate)*ai(kt,inev,lnev,jnev)
      enddo
      enddo
      enddo

      return
      end subroutine newint_airj
c----------------------------------------------------
      subroutine newint_airjnew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,j,traint,dumai,cvirt)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,j,ii,ij,il,jt,kt,lt,ir,ia
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,lmol,is,js,ks,ls
      integer inev,jnev,lnev
      real*8  factor,dum12
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 cvirt(nvirt,nvirt,metat),traint(nact,nvirt,metat),
     *       dumai(nact,nvirt,metat)


c     Transformation of integrals of the type <ar|ij> [in the other formalism (ai,rj)]

      do istate=1,metat
      do ia=1,nact
      do ir=1,nvirt
      traint(ia,ir,istate)=0.d0
      dumai(ia,ir,istate)=0.d0
      enddo
      enddo
      enddo


      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      js=itsym(j)
      jmol=nord_nev2mol(j)-nbe(js)

      do ia=1,nact
      do ir=1,nvirt
      lt=ir+ncore+nact
      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      do ij=1,nc(js)
      jnev=nord_mol2nev(ij+nbe(js))
      do istate=1,metat
      dum12=uc(ii,imol,is,istate)*uc(ij,jmol,js,istate)
      dumai(ia,ir,istate)=dumai(ia,ir,istate)+
     *       dum12*ai(ia+ncore,inev,lt,jnev)
      enddo
      enddo
      enddo
      enddo
      enddo

      do istate=1,metat
      call dgemm('N','N',nact,nvirt,nvirt,1.d0,dumai(1,1,istate),
     *  nact,cvirt(1,1,istate),nvirt,0.d0,traint(1,1,istate),nact)
      enddo

      return
      end subroutine newint_airjnew
c----------------------------------------------------
      subroutine newint_ris_a(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ir,i,iu,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,ii,ij,ik,it,kt,lt,ir,ia,iu,ijs,ijks
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,kmol,is,js,ks,ls
      integer inev,jnev,knev
      real*8  factor,traint(nact,metat),dum12
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <rs|ia> [in the other formalism (ri,sa)]
c

      it=ir+ncore+nact
      is=itsym(it)
      imol=nord_nev2mol(it)-nbe(is)-nc(is)-na(is)

      js=itsym(i)
      jmol=nord_nev2mol(i)-nbe(js)
      ijs=its(is,js)

      kt=iu+ncore+nact
      ks=itsym(kt)
      kmol=nord_nev2mol(kt)-nbe(ks)-nc(ks)-na(ks)
      ijks=its(ijs,ks)

      do ia=1,nact
      lt=ia+ncore
      ls=itsym(lt)

      do istate=1,metat

      traint(ia,istate)=0.d0
      if (ls.ne.ijks) cycle

      do ii=1,nv(is)
      inev=nord_mol2nev(ii+nbe(is)+nc(is)+na(is))
      do ij=1,nc(js)
      jnev=nord_mol2nev(ij+nbe(js))
      dum12=uv(ii,imol,is,istate)*uc(ij,jmol,js,istate)
      do ik=1,nv(ks)
      knev=nord_mol2nev(ik+nbe(ks)+nc(ks)+na(ks))
      traint(ia,istate)=traint(ia,istate)+dum12*
     *           uv(ik,kmol,ks,istate)*ai(inev,jnev,knev,lt)
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end subroutine newint_ris_a
c----------------------------------------------------
      subroutine newint_ris_anew(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ir,iu,traint,dumai,ccore)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,ii,ij,ik,it,kt,lt,ir,ia,iu,ijs,ijks,j
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,kmol,is,js,ks,ls
      integer inev,jnev,knev
      real*8  factor,traint(ncore,nact,metat),dum12
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat),
     *       ccore(ncore,ncore,metat),dumai(ncore,nact,metat)

c
c     Transformation of integrals of the type <rs|ia> [in the other formalism (ri,sa)]
c
      do istate=1,metat
      do ia=1,nact
      do i=1,ncore
      traint(i,ia,istate)=0.d0
      dumai(i,ia,istate)=0.d0
      enddo
      enddo
      enddo

      it=ir+ncore+nact
      is=itsym(it)
      imol=nord_nev2mol(it)-nbe(is)-nc(is)-na(is)

      kt=iu+ncore+nact
      ks=itsym(kt)
      kmol=nord_nev2mol(kt)-nbe(ks)-nc(ks)-na(ks)

      do j=1,ncore
      js=itsym(j)
      ijs=its(is,js)
      ijks=its(ijs,ks)

      do ia=1,nact
      lt=ia+ncore
      ls=itsym(lt)

      if (ls.ne.ijks) cycle

      do istate=1,metat

      do ii=1,nv(is)
      inev=nord_mol2nev(ii+nbe(is)+nc(is)+na(is))
      do ik=1,nv(ks)
      dum12=uv(ii,imol,is,istate)*uv(ik,kmol,ks,istate)
      knev=nord_mol2nev(ik+nbe(ks)+nc(ks)+na(ks))
      dumai(j,ia,istate)=dumai(j,ia,istate)+
     *        dum12*ai(inev,j,knev,lt)
      enddo
      enddo
      enddo
      enddo
      enddo

      do istate=1,metat
      if(ncore == 0) cycle
      call dgemm('T','N',ncore,nact,ncore,1.d0,ccore(1,1,istate),
     *  ncore,dumai(1,1,istate),ncore,0.d0,traint(1,1,istate),ncore)
      enddo

#ifdef _CHOL_DEBUG_
       do istate=1,metat
       do i=1,ncore
       do ia=1,nact
       write (6,*) 'dumai,traint',i,ia,istate,dumai(i,ia,istate),
     *              traint(i,ia,istate)
       enddo
       enddo
       enddo
#endif


      return
      end subroutine newint_ris_anew
c----------------------------------------------------
      subroutine newint_iajb(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,ia,j,ib,istate,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,j,ii,ik,jt,lt,ia,ib
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,kmol,is,js,ks,ls
      integer inev,knev
      real*8  factor,traint
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <ij|ab> [in the other formalism (ia,jb)]
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      jt=ia+ncore
      js=itsym(jt)

      ks=itsym(j)
      kmol=nord_nev2mol(j)-nbe(ks)

      lt=ib+ncore
      ls=itsym(lt)

      traint=0.d0

      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      do ik=1,nc(ks)
      knev=nord_mol2nev(ik+nbe(ks))
      traint=traint+uc(ii,imol,is,istate)*uc(ik,kmol,ks,istate)*
     *        ai(inev,jt,knev,lt)
      enddo
      enddo

      return
      end subroutine newint_iajb
c----------------------------------------------------
      subroutine newint_rs_ab(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,r,s,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,ii,ik,it,jt,kt,lt,ia,ib,r,s,iks,jls
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,kmol,is,js,ks,ls,inev,knev
      real*8  factor,traint(nact,nact,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <rs|ab> [in the other formalism (ra,sb)]
c

      it=r+ncore+nact
      is=itsym(it)
      imol=nord_nev2mol(it)-nbe(is)-nc(is)-na(is)

      kt=s+ncore+nact
      ks=itsym(kt)
      kmol=nord_nev2mol(kt)-nbe(ks)-nc(ks)-na(ks)
      iks=its(is,ks)

      do istate=1,metat

      do ia=1,nact
      jt=ia+ncore
      js=itsym(jt)

      do ib=1,nact
      lt=ib+ncore
      ls=itsym(lt)
      jls=its(js,ls)

      traint(ia,ib,istate)=0.d0
      if (iks.ne.jls) cycle

      do ii=1,nv(is)
      inev=nord_mol2nev(ii+nbe(is)+nc(is)+na(is))
      do ik=1,nv(ks)
      knev=nord_mol2nev(ik+nbe(ks)+nc(ks)+na(ks))
      traint(ia,ib,istate)=traint(ia,ib,istate)+uv(ii,imol,is,istate)*
     *        uv(ik,kmol,ks,istate)*ai(inev,jt,knev,lt)
      enddo
      enddo
      enddo
      enddo
      enddo


      return
      end subroutine newint_rs_ab
c----------------------------------------------------
      subroutine newint_ir(metat,norb,ncore,nact,onel,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,irv,istate,traint)
      use info_symmetry
      implicit none
      integer metat,i,j,ii,ij,irv,jt
      integer norb,nact,ncore,indice,istate,indi
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,is,js
      integer inev,jnev
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 onel(norb*(norb+1)/2),factor,traint
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

c
c     Transformation of the monoelectronic integral with indices i,r
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      jt=irv+ncore+nact
      js=itsym(jt)
      jmol=nord_nev2mol(jt)-nbe(js)-nc(js)-na(js)

      traint=0.d0

      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      do ij=1,nv(js)
      jnev=nord_mol2nev(ij+nbe(js)+nc(js)+na(js))

      indi=indice(inev,jnev)
      traint=traint+uc(ii,imol,is,istate)*uv(ij,jmol,js,istate)*
     *        onel(indi)

      enddo
      enddo

      return
      end subroutine newint_ir
c----------------------------------------------------
      subroutine newint_ijkr(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,j,k,irv,istate,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,j,k,ii,ij,ik,il,irv,lt
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,kmol,lmol,is,js,ks,ls
      integer inev,jnev,knev,lnev
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 factor,traint,dum12,dum123

c
c     Transformation of integrals of the type <ik|jr> [in the other formalism (ij,kr)]
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)
      js=itsym(j)
      jmol=nord_nev2mol(j)-nbe(js)
      ks=itsym(k)
      kmol=nord_nev2mol(k)-nbe(ks)
      lt=irv+ncore+nact
      ls=itsym(lt)
      lmol=nord_nev2mol(lt)-nbe(ls)-nc(ls)-na(ls)

      traint=0.d0


      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      do ij=1,nc(js)
      jnev=nord_mol2nev(ij+nbe(js))
      dum12=uc(ii,imol,is,istate)*uc(ij,jmol,js,istate)
      do ik=1,nc(ks)
      knev=nord_mol2nev(ik+nbe(ks))
      dum123=dum12*uc(ik,kmol,ks,istate)
      do il=1,nv(ls)
      lnev=nord_mol2nev(il+nbe(ls)+nc(ls)+na(ls))
      traint=traint+dum123*uv(il,lmol,ls,istate)*
     *        ai(inev,jnev,knev,lnev)
      enddo
      enddo
      enddo
      enddo

      return
      end subroutine newint_ijkr
c----------------------------------------------------
      subroutine newint_ir_ab_J(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,r,traint,cvirt,ccore)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,ii,ij,jt,kt,lt,ia,ib,r,ijs,kls,j
      integer norb,nact,ncore,nvirt,indice,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,is,js,ks,ls
      integer inev,jnev
      real*8  factor,traint(nact,nact,metat),dumai,dumint
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 cvirt(nvirt,nvirt,metat),ccore(ncore,ncore,metat)
      allocatable dumai(:,:),dumint(:,:)
      allocate (dumai(ncore,nvirt))
      allocate (dumint(ncore,metat))

c
c     Transformation of integrals of the type <ia|rb> [in the other formalism (ir,ab)]
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      jt=r+ncore+nact
      js=itsym(jt)
      jmol=nord_nev2mol(jt)-nbe(js)-nc(js)-na(js)
      ijs=its(is,js)


      do ia=1,nact
      kt=ia+ncore
      ks=itsym(kt)

      do ib=1,nact
      lt=ib+ncore
      ls=itsym(lt)
      kls=its(ks,ls)

      do j=1,ncore
      do ij=1,nvirt
      ii=ij+ncore+nact
      dumai(j,ij)=ai(j,ii,kt,lt)
      enddo
      enddo

      do istate=1,metat
      call dgemm('N','N',ncore,1,nvirt,1.d0,dumai,
     *  ncore,cvirt(1,r,istate),nvirt,0.d0,dumint(1,istate),ncore)
      call dgemm('T','N',1,1,ncore,1.d0,ccore(1,i,istate),ncore,
     *     dumint(1,istate),ncore,0.d0,traint(ia,ib,istate),1)
      enddo

      enddo
      enddo


      return
      end subroutine newint_ir_ab_J
c----------------------------------------------------
      subroutine newint_ir_ab_K(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,r,traint,cvirt,ccore)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,ii,ik,jt,kt,lt,ia,ib,r,iks,jls,j,ij
      integer norb,nact,ncore,nvirt,indice,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,kmol,is,js,ks,ls
      integer inev,knev
      real *8 cvirt(nvirt,nvirt,metat),ccore(ncore,ncore,metat)
      real*8  factor,traint(nact,nact,metat),dumai,dumint
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      allocatable dumai(:,:),dumint(:,:)
      allocate (dumai(ncore,nvirt))
      allocate (dumint(ncore,metat))

c
c     Transformation of integrals of the type <ir|ab> [in the other formalism (ia,rb)]
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      kt=r+ncore+nact
      ks=itsym(kt)
      kmol=nord_nev2mol(kt)-nbe(ks)-nc(ks)-na(ks)
      iks=its(is,ks)


      do ia=1,nact
      jt=ia+ncore
      js=itsym(jt)

      do ib=1,nact
      lt=ib+ncore
      ls=itsym(lt)
      jls=its(js,ls)

      do j=1,ncore
      do ij=1,nvirt
      ii=ij+ncore+nact
      dumai(j,ij)=ai(j,jt,ii,lt)
      enddo
      enddo

      do istate=1,metat
      call dgemm('N','N',ncore,1,nvirt,1.d0,dumai,
     *  ncore,cvirt(1,r,istate),nvirt,0.d0,dumint(1,istate),ncore)
      call dgemm('T','N',1,1,ncore,1.d0,ccore(1,i,istate),ncore,
     *     dumint(1,istate),ncore,0.d0,traint(ia,ib,istate),1)
      enddo

      enddo
      enddo


      return
      end subroutine newint_ir_ab_K
c----------------------------------------------------
      subroutine newint_ar(metat,norb,ncore,nact,onel,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ia,irv,istate,traint)
      use info_symmetry
      implicit none
      integer metat,i,j,ij,irv,jt,ia,it
      integer norb,nact,ncore,indice,istate,indi
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer jmol,js
      integer jnev
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 onel(norb*(norb+1)/2),factor,traint
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

c
c     Transformation of the monoelectronic integral with indices a,r
c

      it=ia+ncore

      jt=irv+ncore+nact
      js=itsym(jt)
      jmol=nord_nev2mol(jt)-nbe(js)-nc(js)-na(js)

      traint=0.d0

      do ij=1,nv(js)
      jnev=nord_mol2nev(ij+nbe(js)+nc(js)+na(js))

      indi=indice(it,jnev)
      traint=traint+uv(ij,jmol,js,istate)*
     *        onel(indi)

      enddo

      return
      end subroutine newint_ar
c----------------------------------------------------
      subroutine newint_arij(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ia,ir,i,j,istate,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,j,ij,ik,il,it,jt,ir,ia
      integer norb,nact,ncore,nvirt,indice,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer jmol,kmol,lmol,js,ks,ls
      integer jnev,knev,lnev
      real*8  factor,traint,dum12
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <ai|rj> [in the other formalism (ar,ij)]
c

      it=ia+ncore

      jt=ir+ncore+nact
      js=itsym(jt)
      jmol=nord_nev2mol(jt)-nbe(js)-nc(js)-na(js)

      ks=itsym(i)
      kmol=nord_nev2mol(i)-nbe(ks)

      ls=itsym(j)
      lmol=nord_nev2mol(j)-nbe(ls)

      traint=0.d0

      do ij=1,nv(js)
      jnev=nord_mol2nev(ij+nbe(js)+nc(js)+na(js))
      do ik=1,nc(ks)
      knev=nord_mol2nev(ik+nbe(ks))
      dum12=uv(ij,jmol,js,istate)*uc(ik,kmol,ks,istate)
      do il=1,nc(ls)
      lnev=nord_mol2nev(il+nbe(ls))
      traint=traint+dum12*uc(il,lmol,ls,istate)*ai(it,jnev,knev,lnev)
      enddo
      enddo
      enddo

      return
      end subroutine newint_arij
c----------------------------------------------------
      subroutine newint_abcr(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,ia,ib,ic,r,istate,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,il,it,jt,kt,lt,ia,ib,r,ic
      integer norb,nact,ncore,nvirt,indice,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer lmol,is,js,ks,ls,lis,lijs
      integer inev,jnev,knev,lnev
      real*8  factor,traint
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 dumJ,dumK,dumM

c
c     Transformation of integrals of the type <ac|br> [in the other formalism (ab,cr)]
c
      lt   = r+ncore+nact
      ls   = itsym(lt)
      lmol = nord_nev2mol(lt)-nbe(ls)-nc(ls)-na(ls)

      it   = ia+ncore
      is   = itsym(it)
      lis  = its(ls,is)

      jt   = ib+ncore
      js   = itsym(jt)
      lijs = its(lis,js)

      kt   = ic+ncore
      ks   = itsym(kt)

      traint=0.d0

      if (lijs.ne.ks) return

      do il=1,nv(ls)
      lnev=nord_mol2nev(il+nbe(ls)+nc(ls)+na(ls))
      traint=traint+uv(il,lmol,ls,istate)*ai(it,jt,kt,lnev)
      enddo

!     print *, 'traint is for i, j, k, r --> ',it,jt,kt,lt, ' --> ',
!    &          traint


      return
      end subroutine newint_abcr
c----------------------------------------------------
      subroutine newint_r_abc(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,r,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,il,it,jt,kt,lt,ia,ib,r,ic
      integer norb,nact,ncore,nvirt,indice,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer lmol,is,js,ks,ls,lis,lijs
      integer lnev
      real*8  factor,traint(nact,nact,nact,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <ac|br> [in the other formalism (ab,cr)]
c

      lt=r+ncore+nact
      ls=itsym(lt)
      lmol=nord_nev2mol(lt)-nbe(ls)-nc(ls)-na(ls)

      do ia=1,nact
      it=ia+ncore
      is=itsym(it)
      lis=its(ls,is)

      do ib=1,nact
      jt=ib+ncore
      js=itsym(jt)
      lijs=its(lis,js)

      do ic=1,nact
      kt=ic+ncore
      ks=itsym(kt)

      do istate=1,metat

      traint(ia,ib,ic,istate)=0.d0

      if (lijs.ne.ks) cycle

      do il=1,nv(ls)
      lnev=nord_mol2nev(il+nbe(ls)+nc(ls)+na(ls))
      traint(ia,ib,ic,istate)=traint(ia,ib,ic,istate)+
     *        uv(il,lmol,ls,istate)*ai(it,jt,kt,lnev)
      enddo

      enddo
      enddo
      enddo
      enddo

      return
      end subroutine newint_r_abc
c----------------------------------------------------
      subroutine newint_i_abc(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,ijs,ijks,il,it,jt,kt,ia,ib,ic
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer is,js,ks,ls
      integer lnev,lmol
      real*8  factor,traint(nact,nact,nact,metat)
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)

c
c     Transformation of integrals of the type <ac|bi> [in the other formalism (ab,ci)]
c


      ls=itsym(i)
      lmol=nord_nev2mol(i)-nbe(ls)

      do ia=1,nact
      it=ia+ncore
      is=itsym(it)

      do ib=1,nact
      jt=ib+ncore
      js=itsym(jt)
      ijs=its(is,js)

      do ic=1,nact
      kt=ic+ncore
      ks=itsym(kt)
      ijks=its(ijs,ks)

      do istate=1,metat

      traint(ia,ib,ic,istate)=0.d0
      if (ijks.ne.ls) cycle

      do il=1,nc(ls)
      lnev=nord_mol2nev(il+nbe(ls))
      traint(ia,ib,ic,istate)=traint(ia,ib,ic,istate)+
     *        uc(il,lmol,ls,istate)*ai(it,jt,kt,lnev)
      enddo
      enddo
      enddo
      enddo
      enddo


      return
      end subroutine newint_i_abc
c----------------------------------------------------
      subroutine newint_ia(metat,norb,ncore,nact,onel,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,ia,istate,traint)
      use info_symmetry
      implicit none
      integer metat,i,j,ii,jt,ia,is,js
      integer norb,nact,ncore,indice,istate,indi
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol
      integer inev
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 onel(norb*(norb+1)/2),factor,traint
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

c
c     Transformation of the monoelectronic integral with indices i,a
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)

      jt=ia+ncore
      js=itsym(jt)

      traint=0.d0

      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))

      indi=indice(inev,jt)
      traint=traint+uc(ii,imol,is,istate)*onel(indi)

      enddo

      return
      end subroutine newint_ia
c----------------------------------------------------
      subroutine newint_ijka(metat,norb,ncore,nact,nvirt,
     *   factor,nsym,uv,uc,nmo,nc,na,nv,nbe,ncmax,nvmax,
     *   nord_nev2mol,nord_mol2nev,i,j,k,ia,istate,traint)

      use ijkl_utils, only : ai
      use info_symmetry

      implicit none

      integer metat,i,j,k,ii,ij,ik,lt,ia
      integer norb,nact,ncore,nvirt,istate
      integer nsym,isym,nord_nev2mol(norb),
     *        nord_mol2nev(norb)
      integer imol,jmol,kmol,is,js,ks
      integer inev,jnev,knev
      integer nmo(nsym),nc(nsym),na(nsym),nv(nsym),nbe(nsym),ncmax,nvmax
      real*8 uv(nvmax,nvmax,nsym,metat),uc(ncmax,ncmax,nsym,metat)
      real*8 factor,traint,dum12

c
c     Transformation of integrals of the type <ik|ja> [in the other formalism (ij,ka)]
c

      is=itsym(i)
      imol=nord_nev2mol(i)-nbe(is)
      js=itsym(j)
      jmol=nord_nev2mol(j)-nbe(js)
      ks=itsym(k)
      kmol=nord_nev2mol(k)-nbe(ks)

      lt=ia+ncore

      traint=0.d0


      do ii=1,nc(is)
      inev=nord_mol2nev(ii+nbe(is))
      do ij=1,nc(js)
      jnev=nord_mol2nev(ij+nbe(js))
      dum12=uc(ii,imol,is,istate)*uc(ij,jmol,js,istate)
      do ik=1,nc(ks)
      knev=nord_mol2nev(ik+nbe(ks))
      traint=traint+dum12*uc(ik,kmol,ks,istate)*ai(inev,jnev,knev,lt)
      enddo
      enddo
      enddo

      return
      end subroutine newint_ijka
c------------------------------------
!     prepare the half-transformed Cholesky vectors
      subroutine cholesky_transform(cho_reduced_vec,ccore,cvirt,
     &         ncore,nact,nvirt,metat)
      use ijkl_utils
      implicit none
      integer,intent(in) :: ncore,nact,nvirt,metat
      real*8 cho_reduced_vec(nchovec,nlvec,metat) ! where to store the vectors
      real*8 ccore(ncore,ncore,metat) ! Core and virtual Fock operators
      real*8 cvirt(nvirt,nvirt,metat) !

      real*8,allocatable :: cho_quarter_transform(:,:,:) ! store the result of a quarter-transform here

      integer i,j,r,ii,rr,ir,irr,a,ia,istate
      integer indice

      indice(i,r)=max(i,r)*(max(i,r)-1)/2+min(i,r)

      ! TODO: maybe convert it to a matrix form



      ! this routine calculates the Cholesky half-transformed vectors
      ! the code is derived from an O(n^5) "smart" integral transform, though for a Cholesky vector it's only O(n^3)
      ! Core-core |ij) vectors

      allocate(cho_quarter_transform(nchovec,ncore,ncore))

      do istate=1,metat
        cho_quarter_transform = 0.0d0
        do i=1,ncore
          do j=1,ncore
            do ii=1,ncore
              cho_quarter_transform(:,i,j)=
     &  cho_quarter_transform(:,i,j)+ccore(ii,i,istate)*
     &  nevpt_ijkl%cholesky_array(:,indice(ii,j))
            end do
          end do
        end do

        do i=1,ncore
          do j=i,ncore ! this must be j=i because cho_reduced_vec is symmetric!
            do ii=1,ncore
              cho_reduced_vec(:,indice(i,j),istate) =
     &  cho_reduced_vec(:,indice(i,j),istate)+
     &  ccore(ii,j,istate)*cho_quarter_transform(:,i,ii)
            end do
          end do
        end do
      end do

      if (allocated(cho_quarter_transform)) then
        deallocate(cho_quarter_transform)
      end if

      ! Core-virtual |ir) vectors
      allocate(cho_quarter_transform(nchovec,ncore,nvirt))

      do istate=1,metat
        cho_quarter_transform = 0.0d0
        do i=1,ncore
          do r=1,nvirt
            ir = r + ncore + nact
            do ii=1,ncore
              cho_quarter_transform(:,i,r)=
     &  cho_quarter_transform(:,i,r)+ccore(ii,i,istate)*
     &  nevpt_ijkl%cholesky_array(:,indice(ii,ir))
            end do
          end do
        end do

        do i=1,ncore
          do r=1,nvirt
            ir = r + ncore + nact
            do rr=1,nvirt
              cho_reduced_vec(:,indice(i,ir),istate) =
     &  cho_reduced_vec(:,indice(i,ir),istate)+
     &  cvirt(rr,r,istate)*cho_quarter_transform(:,i,rr)
            end do
          end do
        end do
      end do


      if (allocated(cho_quarter_transform)) then
        deallocate(cho_quarter_transform)
      end if

      !! stupid transformation code (O(n^4)), you can enable it for testing (though it's untested and may be buggy)
!       do istate=1,metat
!         do i=1,ncore
!           do r=1,nvirt
!             ir = r + ncore + nact
!             do ii=1,ncore
!               do rr=1,nvirt
!                 irr = rr + ncore + nact
!                 cho_reduced_vec(:,indice(i,ir),istate) =
!      &            cho_reduced_vec(:,indice(i,ir),istate) +
!      &            ccore(ii,i,metat)*cvirt(rr,r,istate)*
!      &           nevpt_ijkl%cholesky_array(:,indice(ii,irr))
!               end do
!             end do
!           end do
!         end do
!       end do


      ! core-active |ai) Cholesky vectors
      do istate=1,metat
        do i=1,ncore
          do a=1,nact
            ia = a + ncore
            do ii=1,ncore
              cho_reduced_vec(:,indice(i,ia),istate) =
     &  cho_reduced_vec(:,indice(i,ia),istate)+
     &  ccore(ii,i,istate)*nevpt_ijkl%cholesky_array(:,indice(ii,ia))
            end do
          end do
        end do
      end do

      ! active-virtual |ar) Cholesky vectors
      do istate=1,metat
        do r=1,nvirt
          ir = r + ncore + nact
          do a=1,nact
            ia = a + ncore
            do rr=1,nvirt
              irr = rr + ncore + nact
              cho_reduced_vec(:,indice(ir,ia),istate) =
     &  cho_reduced_vec(:,indice(ir,ia),istate)+
     &  cvirt(rr,r,istate)*nevpt_ijkl%cholesky_array(:,indice(irr,ia))
            end do
          end do
        end do
      end do

      ! active-active |ab) vectors need not be transformed.
      ! for now, they're copied as-is to the new vectors. If we need to
      ! save some memory, we may omit this step and copy the
      ! vectors directly from cholesky_array on integral reconstruction.
      do istate=1, metat
        do i=ncore+1,ncore+nact
          do j=i,ncore+nact
            cho_reduced_vec(:,indice(i,j),istate) =
     &   nevpt_ijkl%cholesky_array(:,indice(i,j))
          end do
        end do
      end do

      ! TODO: other integrals...
      end subroutine cholesky_transform





      ! The following routines return a batch of transformed integrals, (i,x,x,j) analogously to newint_xxxx() calculated
      ! with cholesky_traint() from the transformed Cholesky vectors


      ! one could avoid this subroutine and use cholesky_traint() directly in v0(),v1k() etc.
      ! but it should be easier to parallelise if it's in a separate subroutine

      subroutine newint_irsj_cholesky(cho_reduced_vec,traint,
     &    metat,i,j,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,j,ncore,nact,nvirt,metat
      integer r,s,istate,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nvirt,nvirt,metat)

C BLAS ddot()
      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &   indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)


C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(istate,r,s) collapse(3)
      do istate=1,metat
        do r=1,nvirt
          do s=1,nvirt
            ! TODO: is the matrix symmetric? check it! maybe it can add additional computational saving
            traint(r,s,istate) = cholesky_traint(
     &       i,r+ncore+nact,s+ncore+nact,j,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_irsj_cholesky

      subroutine newint_airj_cholesky(cho_reduced_vec,traint,
     &    metat,i,j,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,j,ncore,nact,nvirt,metat
      integer r,a,istate,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nvirt,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(r,a,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do r=1,nvirt
            traint(a,r,istate) = cholesky_traint(
     &       a + ncore,i,r + ncore + nact,j,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_airj_cholesky

      subroutine newint_arij_cholesky(cho_reduced_vec,traint,
     &    metat,i,j,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,j,ncore,nact,nvirt,metat
      integer r,a,istate,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nvirt,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(r,a,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do r=1,nvirt
            traint(a,r,istate) = cholesky_traint(
     &       a + ncore,r + ncore + nact,i,j,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_arij_cholesky

      subroutine newint_ris_a_cholesky(cho_reduced_vec,traint,
     &    metat,r,s,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: r,s,ncore,nact,nvirt,metat
      integer i,a,istate,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(ncore,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,a,istate) collapse(3)
      do istate=1,metat
        do i=1,ncore
          do a=1,nact
            traint(i,a,istate) = cholesky_traint(
     &       r,i,s,a + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_ris_a_cholesky

      subroutine newint_iajb_cholesky(cho_reduced_vec,traint,
     &    metat,i,j,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,j,ncore,nact,nvirt,metat
      integer a,b,istate,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
            traint(a,b,istate) = cholesky_traint(
     &       i,a + ncore,j,b + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_iajb_cholesky

      subroutine newint_rasb_cholesky(cho_reduced_vec,traint,
     &    metat,r,s,ncore,nact,nvirt) ! this subroutine is technically totally analogous to newint_iajb_cholesky, but may be changed in the future for better parallelisation
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: r,s,ncore,nact,nvirt,metat
      integer a,b,istate,i,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
            traint(a,b,istate) = cholesky_traint(
     &       r,a + ncore,s,b + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_rasb_cholesky

      subroutine newint_irab_cholesky(cho_reduced_vec,traint,
     &    metat,i,r,ncore,nact,nvirt) ! corresponds to newint_ir_ab_J but keeps to the common nomenclature
      ! again, traints are calculated in active x active batches, possibly changing the batching (eg over virtuals) might improve parallelism
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,r,ncore,nact,nvirt,metat
      integer a,b,istate,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
            traint(a,b,istate) = cholesky_traint(
     &       i,r,a + ncore,b + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_irab_cholesky
      subroutine newint_iarb_cholesky(cho_reduced_vec,traint,
     &    metat,i,r,ncore,nact,nvirt) ! corresponds to newint_ir_ab_K
      ! but keeps to the common nomenclature
      ! again, traints are calculated in active x active batches,
      ! possibly changing the batching (eg over virtuals) might improve
      ! parallelism
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,r,ncore,nact,nvirt,metat
      integer a,b,istate,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
            traint(a,b,istate) = cholesky_traint(
     &       i,a + ncore,r,b + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_iarb_cholesky

      subroutine newint_abcr_cholesky(cho_reduced_vec,traint,
     &    metat,r,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: r,ncore,nact,nvirt,metat
      integer a,b,c,istate,i,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,c,istate) collapse(4)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
            do c=1,nact
              traint(a,b,c,istate) = cholesky_traint(
     &        a + ncore,b + ncore,c + ncore,r,istate)
            end do
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_abcr_cholesky


      subroutine newint_abbr_cholesky(cho_reduced_vec,traint,
     &    metat,r,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: r,ncore,nact,nvirt,metat
      integer a,b,istate,i,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,istate) collapse(3)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
              traint(a,b,istate) = cholesky_traint(
     &        a + ncore,b + ncore,b + ncore,r,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_abbr_cholesky

      subroutine newint_iija_cholesky(cho_reduced_vec,traint,
     &    metat,i,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,ncore,nact,nvirt,metat
      integer j,a,istate,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(ncore,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,j,istate) collapse(3)
      do istate=1,metat
        do j=1,ncore
          do a=1,nact
              traint(j,a,istate) = cholesky_traint(
     &        i,i,j,a + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_iija_cholesky

      subroutine newint_ijia_cholesky(cho_reduced_vec,traint,
     &    metat,i,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,ncore,nact,nvirt,metat
      integer j,a,istate,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(ncore,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,j,istate) collapse(3)
      do istate=1,metat
        do j=1,ncore
          do a=1,nact
              traint(j,a,istate) = cholesky_traint(
     &        i,j,i,a + ncore,istate)
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_ijia_cholesky

      subroutine newint_iabc_cholesky(cho_reduced_vec,traint,
     &    metat,i,ncore,nact,nvirt)
      use ijkl_utils, only : nevpt_ijkl,nchovec,nlvec
      implicit none
      integer,intent(in) :: i,ncore,nact,nvirt,metat
      integer a,b,c,istate,j,k,l
      real*8 cho_reduced_vec(nchovec,nlvec,metat)
      real*8 traint(nact,nact,nact,metat)

      real*8, external    :: ddot
! inline (statement) functions
      integer indice
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      real*8 cholesky_traint
      cholesky_traint(i,j,k,l,istate)=ddot(nchovec,cho_reduced_vec(1,
     &  indice(i,j),istate),1,cho_reduced_vec(1,indice(k,l),istate),1)
C$OMP PARALLEL DO DEFAULT(shared) PRIVATE(a,b,c,istate) collapse(4)
      do istate=1,metat
        do a=1,nact
          do b=1,nact
            do c=1,nact
              traint(a,b,c,istate) = cholesky_traint(
     &        a + ncore,b + ncore,c + ncore,i,istate)
            end do
          end do
        end do
      end do
C$OMP END PARALLEL DO
      end subroutine newint_iabc_cholesky

!-------------------------------------------------------------------------------
      subroutine hermitian_heff(key,heff_evc,vr,wr,metat,heff,state_ptr)

      !> similarity transformation of Heff to make it Hermitian (transform also the right eigenvectors)
      !> see eq. (19) in JCP, 121, 4043 (2004)
      !>          ~                 ~                              ~     ~
      !> Heff' = Psi S^{-1/2} _^_ (Psi S^{-1/2})+ with S_{ij} = <Psi_i|Psi_j> and
      !>  ~
      !> Psi_i is the i-th right eigenvector and
      !> _^_ is the diagonal matrix of eigenvalues of Heff
!-------------------------------------------------------------------------------
      implicit none

      real*8,  intent(inout) :: vr(metat,metat)
      real*8,  intent(inout) :: heff_evc(metat,metat)
      real*8,  intent(inout) :: heff(metat,metat)
      real*8,  intent(inout) :: wr(metat)
      integer, intent(in)    :: metat
      integer, intent(in)    :: state_ptr(metat)
      character(len=2)       :: key

      real*8, allocatable    :: s(:,:), kvec(:,:), kvec2(:,:)
      real*8, allocatable    :: vr_tmp(:,:), work(:)
      real*8                 :: dfact, tmp
      real*8, external       :: ddot
      integer                :: istate, jstate, k, info, lwork
!-------------------------------------------------------------------------------

        allocate(s(metat,metat),vr_tmp(metat,metat),
     &           work(max(1,3*metat-1)),kvec(metat,metat),
     &           kvec2(metat,metat))
        s = 0; vr_tmp = 0; work = 0; kvec = 0; kvec2 = 0

        !> calculate overlap matrix S
        do istate = 1, metat
          do jstate = 1, metat
            s(istate,jstate) = ddot(metat,vr(1,istate),1,vr(1,jstate),1)
          end do
        end do

       !print *,'overlap matrix of right eigenvectors --> ',key
       !call matout(metat,metat,s,metat,1.d0)

        !> form S^{-1/2} in s
        call rtmatinv(s, metat, wr, kvec, kvec2)

        !> transform eigenvectors C' = CS^{-1/2}
        call dgemm('n','n',metat,metat,metat,1.0d0,
     &             vr    ,metat,!0.0d0,
     &             s     ,metat, 0.0d0,
     &             vr_tmp,metat)

        !> build Heff' as C' _^_ C'+
        call dgemm('n','n',metat,metat,metat,1.0d0,
     &             vr_tmp  ,metat,!0.0d0,
     &             heff_evc,metat, 0.0d0,
     &             s       ,metat)

        call dgemm('n','t',metat,metat,metat,1.0d0,
     &             s       ,metat,!0.0d0,
     &             vr_tmp  ,metat, 0.0d0,
     &             heff_evc,metat)

        print *,' Heff~ (Hermitian) --> ',key
        call matout(metat,metat,heff_evc,metat,1.d0)

        !> save Hermitian Heff
        call dcopy(metat**2,heff_evc,1,heff,1)

        !> diagonalize Hermitian Heff~ matrix
        info = 0; lwork = max(1,3*metat-1)
        call dsyev('v','u',metat,heff_evc,metat,wr,work,lwork,info)
        if (info.ne.0)
     &  write(6,*) " hermitian_heff: error in dsyev call! info = ",info

        print *,'eigenvectors of Heff~ (Hermitian) --> ',key
        call matout(metat,metat,heff_evc,metat,1.d0)
        print *,'eigenvalues of Heff~ matrix (Hermitian) --> ',key
        call matout(metat,1,wr,metat,1.d0)

        !> calculate overlap matrix S
        do istate = 1, metat
          do jstate = 1, metat
            s(istate,jstate) =
     &      ddot(metat,heff_evc(1,istate),1,heff_evc(1,jstate),1)
          end do
        end do

       !print *,'overlap matrix of Hermitian right eigenvectors --> ',
     & !         key
       !call matout(metat,metat,s,metat,1.d0)

        deallocate(s,vr_tmp,work,kvec,kvec2)


      end subroutine hermitian_heff
!-------------------------------------------------------------------------------

      subroutine hermitian_right_evc(key,vr,metat)

      !> transform the right eigenvectors of Heff
      !> see eq. (19) in JCP, 121, 4043 (2004)
      !>                                                ~     ~          ~
      !> Heff' = S^{-1/2} Heff S^{1/2} with S_{ij} = <Psi_i|Psi_j> and Psi_i is the i-th right eigenvector
!-------------------------------------------------------------------------------
      implicit none

      real*8,  intent(inout) :: vr(metat,metat)
      integer, intent(in)    :: metat
      character(len=2)       :: key

      real*8, allocatable    :: s(:,:), eval(:), kvec(:,:), kvec2(:,:)
      real*8, allocatable    :: vr_tmp(:,:), work(:)
      real*8                 :: dfact
      real*8, external       :: ddot
      integer                :: istate, jstate, info, lwork
!-------------------------------------------------------------------------------

        allocate(s(metat,metat),vr_tmp(metat,metat),eval(metat),
     &           work(max(1,3*metat-1)),kvec(metat,metat),
     &           kvec2(metat,metat))
        s = 0; vr_tmp = 0; eval = 0; work = 0; kvec = 0; kvec2 = 0

        !> calculate overlap matrix S
        do istate = 1, metat
          do jstate = 1, metat
            s(istate,jstate) = ddot(metat,vr(1,istate),1,vr(1,jstate),1)
          end do
        end do

        call dcopy(metat**2,vr,1,vr_tmp,1)

!       print *,'overlap matrix of right eigenvectors'
!       call matout(metat,metat,s,metat,1.d0)

!       !> form S^{-1/2}
        call rtmatinv(s, metat, eval, kvec, kvec2)
!
!        print *,'trafo matrix'
!       call matout(metat,metat,kvec,metat,1.d0)

        !> transform eigenvectors
        call dgemm('n','n',metat,metat,metat,1.0d0,
     &             vr_tmp,metat,!0.0d0,
     &             s     ,metat, 0.0d0,
     &             vr,    metat)

        print *,'transformed right eigenvectors (Hermitian) ',key
        call matout(metat,metat,vr,metat,1.d0)

        !> calculate overlap matrix S
!       do istate = 1, metat
!         do jstate = 1, metat
!           s(istate,jstate) = ddot(metat,vr(1,istate),1,vr(1,jstate),1)
!         end do
!       end do

!       print *,'overlap matrix of transformed right eigenvectors'
!       call matout(metat,metat,s,metat,1.d0)

!       write(6,'(/a)') 'norm of right eigenvectors (Hermitian) -> ',key
!       write(6,'( a)')'-----------------------------------------------'
!       do istate=1,metat
!         eval(1:metat) = vr(1:metat,istate)
!        write(6,'(a,i3,a,f16.10)') ' # state ',istate,
!    &  ': ',norm2(eval)
!       end do

        deallocate(s,vr_tmp,eval,work,kvec,kvec2)


      end subroutine hermitian_right_evc
!-------------------------------------------------------------------------------

      SUBROUTINE rtmat(A, p, eig, G, GD)
!     Returns A as its pxp symetric matrix square-root.
      INTEGER p
      real*8 A(p,p), eig(p), G(p,p), GD(p,p)
!
      INTEGER i,j, lwork, info
      real*8, allocatable :: work(:)

      !> diagonalize matrix
      allocate(work(max(1,3*p-1))); work = 0
      info = 0; lwork = max(1,3*p-1)
      call dsyev('v','u',p,a,p,eig,work,lwork,info)
      if (info.ne.0)
     &write(6,*) " rmatinv: error in dsyev call! info = ",info
      call dcopy(p**2,a,1,g,1)
      deallocate(work)

!     print *,'S^{1/2}> eigenvectors of A'
!     call matout(p,p,g,p,1.d0)

!     print *,'S^{1/2}> eigenvalues of A'
!     call matout(p,1,eig,p,1.d0)

!     eig now contains the eigenvalues of A, and the corresponding
!     eigenvectors are now the columns of G.
!     Set A as the diagonal matrix of eigenvalue square roots:
      do 200 i=1,p
       do 100 j=1,p
        A(i,j)=0.0d0
 100   continue
       A(i,i)=dsqrt(eig(i))
 200  continue

      call dgemm('n','n',p,p,p,
     &           1.0d0,g,p,a,p,
     &           0.0d0,gd,p)

      call dgemm('n','t',p,p,p,
     &           1.0d0,gd,p,g,p,
     &           0.0d0,a,p)

      END SUBROUTINE rtmat

      SUBROUTINE rtmatinv(A, p, eig, G, GD)
!     Returns A as its pxp symetric matrix inv square-root.
      INTEGER p
      real*8 A(p,p), eig(p), G(p,p), GD(p,p)
!
      INTEGER i,j, lwork, info
      real*8, allocatable :: work(:)

      !> diagonalize matrix
      allocate(work(max(1,3*p-1))); work = 0
      info = 0; lwork = max(1,3*p-1)
      call dsyev('v','u',p,a,p,eig,work,lwork,info)
      if (info.ne.0)
     &write(6,*) " rmatinv: error in dsyev call! info = ",info
      call dcopy(p**2,a,1,g,1)
      deallocate(work)

!     print *,'S^{-1/2}> eigenvectors of A'
!     call matout(p,p,g,p,1.d0)

!     print *,'S^{-1/2}> eigenvalues of A'
!     call matout(p,1,eig,p,1.d0)

!     eig now contains the eigenvalues of A, and the corresponding
!     eigenvectors are now the columns of G.
!     Set A as the diagonal matrix of eigenvalue square roots:
      do 200 i=1,p
       do 100 j=1,p
        A(i,j)=0.0d0
 100   continue
       A(i,i)=1.0d0/dsqrt(eig(i))
 200  continue

      call dgemm('n','n',p,p,p,
     &           1.0d0,g,p,a,p,
     &           0.0d0,gd,p)

      call dgemm('n','t',p,p,p,
     &           1.0d0,gd,p,g,p,
     &           0.0d0,a,p)
      END SUBROUTINE rtmatinv

      end module qdnevpt_core
