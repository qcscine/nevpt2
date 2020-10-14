module ijkl_utils

!> add IJKL_DEBUG flag to the cmake compilation "-DIJKL_DEBUG" if debug print is what you want
!#define IJKL_DEBUG

#ifdef _OPENMP_
  use omp_lib
#endif
  use hdf5_utils
  use info_symmetry
  use info_orbital_space
  use nevpt2_cfg

  implicit none

  public initialize_ijkl
  public finalize_ijkl
  public determine_n2int
  public reijkl
  public readint
  public ext1int2nev1int
  public ext2int2nev2int
  public ai
#ifdef unused
  public itijklt
#endif

  type integrals
     integer*8            :: ntwoint   =  0
     integer              :: noneint   =  0
     integer*8            :: ni4
     integer, allocatable :: nad(:)
     integer, allocatable :: kt(:)
     integer, allocatable :: ndeb(:)
     integer, allocatable :: num(:)
     integer, allocatable :: indic(:)
     integer, allocatable :: jndic(:)
     integer, allocatable :: lndic(:)
     real*8 , allocatable :: twoint(:)
     real*8 , allocatable :: oneint(:)
     real*8 , allocatable :: cholesky_array(:,:)
  end type integrals
  type (integrals), save, public :: nevpt_ijkl

  ! Cholesky set sizes/dimensions
  ! For symmetry supprt, those have to be dimension(8) arrays

    integer,public :: nlvec ! number of vectors in the Cholesky set
    integer,public :: nchovec ! number of Cholesky sets (should be norb(norb+1)/2 or something
    !maybe these two are better off in nevpt_ijkl?

  private
  ! HDF5 handles
    integer(HID_T) :: dset_id
    integer(HID_T) :: dspace_id
!     integer :: norb_ ! number of orbitals
!     ! it's required here for calculating the integrals from the Cholesky vectors
!     ! it should be extended to an array with norb for each symmetry species if we want to support symmetry
  contains

!------------------------------------------------------------------------------

      subroutine initialize_ijkl(n1int,norb)

        integer, intent(in) :: norb
        integer, intent(in) :: n1int

        integer             :: ndim

        ndim                = norb*(norb+1)/2
        nevpt_ijkl%noneint  = n1int
        nevpt_ijkl%ntwoint  = nevpt_ijkl%ni4 ! hopefully that's still valid with Cholesky

!         norb_ = norb
        if (.not.(Do_Cholesky)) then
          allocate(nevpt_ijkl%twoint(nevpt_ijkl%ntwoint))
          nevpt_ijkl%twoint   = 0
        end if
        allocate(nevpt_ijkl%oneint(nevpt_ijkl%noneint))
        allocate(nevpt_ijkl%num(ndim),nevpt_ijkl%indic(ndim),nevpt_ijkl%jndic(ndim),nevpt_ijkl%lndic(ndim))

        nevpt_ijkl%oneint   = 0
        nevpt_ijkl%num      = 0
        nevpt_ijkl%indic    = 0
        nevpt_ijkl%jndic    = 0
        nevpt_ijkl%lndic    = 0

        ! Cholesky vector reading init is done in readint()

      end subroutine initialize_ijkl
!------------------------------------------------------------------------

      subroutine finalize_ijkl()

        if(Do_Cholesky) call hdf5_close_cholesky(dset_id,dspace_id)

        if(allocated(nevpt_ijkl%oneint)) deallocate(nevpt_ijkl%oneint)
        if(allocated(nevpt_ijkl%twoint)) deallocate(nevpt_ijkl%twoint)
        if(allocated(nevpt_ijkl%indic))  deallocate(nevpt_ijkl%indic)
        if(allocated(nevpt_ijkl%jndic))  deallocate(nevpt_ijkl%jndic)
        if(allocated(nevpt_ijkl%lndic))  deallocate(nevpt_ijkl%lndic)
        if(allocated(nevpt_ijkl%num))    deallocate(nevpt_ijkl%num)
        if(allocated(nevpt_ijkl%ndeb))   deallocate(nevpt_ijkl%ndeb)
        if(allocated(nevpt_ijkl%nad))    deallocate(nevpt_ijkl%nad)
        if(allocated(nevpt_ijkl%kt))     deallocate(nevpt_ijkl%kt)
        if(allocated(nevpt_ijkl%cholesky_array)) deallocate(nevpt_ijkl%cholesky_array)

      end subroutine finalize_ijkl
!------------------------------------------------------------------------

      SUBROUTINE determine_n2int(norb,nsym)
!
!     !brief>: compute the number of two-electron integrals
!
      integer, intent(in)  :: norb, nsym
      integer              :: ndegen, noc, ij, i, j, ijs, nsm1, k, kp1, njp, is, js
      integer, allocatable :: nos(:)

      allocate(nos(nsym)); nos = 0

      !> determine dimensions for nad, kt, ndeb

      allocate(nevpt_ijkl%ndeb(nsym)); nevpt_ijkl%ndeb = 0

      NJP=0
      do IS=1,NSYM
        do JS=1,IS
          NJP=NJP+1
        end do
      end do

      allocate(nevpt_ijkl%nad(njp),nevpt_ijkl%kt(njp)); nevpt_ijkl%nad = 0; nevpt_ijkl%kt = 0

      !> start the real stuff...

      DO 802 I=1,NORB
        IS=ITSYM(I)
        DO 803 J=1,I
          JS=ITSYM(J)
          IJS=ITS(IS,JS)
  803   NOS(IJS)=NOS(IJS)+1
  802 CONTINUE
      !print *, 'nos array is ',NOS(1:nsym)

      DO 804 K=1,NSYM-1
      KP1=K+1
      !print *, 'integrals for sym ',k,nevpt_ijkl%NDEB(k)+(NOS(K)*(NOS(K)+1))/2
  804 nevpt_ijkl%NDEB(KP1) = nevpt_ijkl%NDEB(K)+(NOS(K)*(NOS(K)+1))/2
!     print *, 'total integrals ', nevpt_ijkl%NDEB(NSYM),(NOS(NSYM)*(NOS(NSYM)+1))/2

      nevpt_ijkl%ni4 = nevpt_ijkl%NDEB(NSYM)+(NOS(NSYM)*(NOS(NSYM)+1))/2

      WRITE (6,74) nevpt_ijkl%ni4
   74 FORMAT (/' total number of 4-index 2e-integrals:',I12)

      NJP=0
      DO 820 IS=1,NSYM
      DO 820 JS=1,IS
      NJP=NJP+1
      IF(JS.EQ.IS)GOTO 821
      nevpt_ijkl%NAD(NJP)=0
      nevpt_ijkl%KT(NJP)=0
      GOTO 820
  821 nevpt_ijkl%NAD(NJP)=IS
      nevpt_ijkl%KT(NJP)=1
  820 CONTINUE

      deallocate(nos)

      end subroutine determine_n2int
!-------------------------------------------------------------------------------

      subroutine reijkl(norb,nsym)
!
!     PREPARATION POUR LE STOCKAGE DES INTEGRALES BIELECTONIQUES
!
      integer, intent(in)  :: norb, nsym
      integer, allocatable :: nos(:)

      integer   :: ndegen, noc, ij, i, j, ijs, nsm1, k, kp1, njp, is, js

      allocate(nos(nsym)); nos = 0

      IJ=0
      nevpt_ijkl%num(1) = 0
      DO I = 2, norb
        DO J = 1, I
          IJ=IJ+1
          nevpt_ijkl%num(IJ+1)=nevpt_ijkl%num(IJ)+IJ
        end do
      end do

      DO 802 I = 1, NORB
        IS = ITSYM(I)
        DO 803 J = 1, I
                    JS       = ITSYM(J)
                    IJS      = ITS(IS,JS)
                    NOS(IJS) = NOS(IJS) + 1
                         IJ  = nevpt_ijkl%num(I)+J
        nevpt_ijkl%indic(IJ) = IJS
  803   nevpt_ijkl%jndic(IJ) = NOS(IJS)
  802 CONTINUE


      deallocate(nos)

      end subroutine reijkl
!-------------------------------------------------------------------------------

      subroutine readint(norb,nsym,nord_mol2nev,ecore)

      integer, intent(in)  :: norb
      integer, intent(in)  :: nsym
      integer, intent(in)  :: nord_mol2nev(*)
      real*8,  intent(out) :: ecore

      !> read core energy
      ecore = 0.0d0
      datadim(1)    = 1; datadim_bound = 1
      !call hdf5_get_data(file_id(2),"ecore ",datadim,ecore)

      !> compute all necessary offsets, dimensions, etc
      call determine_n2int(norb,nsym)

      !> initialize integral storage arrays
      call initialize_ijkl(norb*(norb+1)/2,norb)

      !> initialize indices array
      call reijkl(norb,nsym)

      call ext1int2nev1int(                    &
                           nevpt_ijkl%oneint,  &
                           (norb*(norb+1))/2,  &
                           inforb%nosh,        &
                           nsym,               &
                           norb,               &
                           nord_mol2nev,       &
                           6                   &
                          )
      if (.not.(Do_Cholesky)) then
        call ext2int2nev2int(                    &
                             norb,               &
                             inforb%nosh,        &
                             inforb%iosh,        &
                             nsym,               &
                             nevpt_ijkl%ndeb,    &
                             nord_mol2nev,       &
                             nevpt_ijkl%twoint,  &
                             6                   &
                            )
        else
          if (nsym > 1) stop "Cholesky decomposition with symmetry not supported yet."
          ! initialise Cholesky reading routines
          call hdf5_init_rd_cholesky(file_id(2),dset_id,dspace_id,nsym,nlvec,nchovec)

          ! ... and the Cholesky array
          allocate(nevpt_ijkl%cholesky_array(nchovec,nlvec))

          ! read the whole Cholesky set into memory
          call hdf5_read_all_cholesky(dset_id,dspace_id,nevpt_ijkl%cholesky_array)
        end if
      end subroutine readint
!-------------------------------------------------------------------

      subroutine ext2int2nev2int(norb,                &
                                 nr_orb_per_sym,      &
                                 iosh,                &
                                 nsym,                &
                                 ndeb,                &
                                 ireost,              &
                                 twoint,              &
                                 lupri)
!
!     !> brief: read integrals from external file and store on nevpt integral array
!     !> author: S. Knecht
!
      !> @params
      integer, intent(in) :: norb
      integer, intent(in) :: nr_orb_per_sym(*)
      integer, intent(in) :: iosh(*)
      integer, intent(in) :: nsym
      integer, intent(in) :: ndeb(*)
      integer, intent(in) :: ireost(*)
      integer, intent(in) :: lupri
      real*8 , intent(out):: twoint(*)
      !> local varaibales
      integer              :: itamp, itamp2, indn, i, j, k, l, ij, kl, indtot
      integer              :: lbatch
      integer              :: nsp, nsr, nsq, nss, nssm, nspq, nspqr, nspqrs
      integer              :: ii, jj, kk, ll, itmp, ktmp, istart, jstart
      real*8, allocatable  :: batch(:)
      character(len=7)     :: nevpt_host_program = 'molcas '
      integer              :: noccend, noccendi, noccendj
      integer*8            :: mbatch

      indtot = 0

      !> loop over symmetry batches for P Q R S
      do nsp = 1, nsym

       if(nr_orb_per_sym(nsp) == 0) cycle

        do nsq = 1, nsp

          if(nr_orb_per_sym(nsp) == 0) cycle

          nspq = ieor(nsp-1,nsq-1)+1

          do nsr = 1, nsp

            if(nr_orb_per_sym(nsr) == 0) cycle

            nspqr = ieor(nspq-1,nsr-1)+1

            nssm = nsr ; if(nsr == nsp) nssm = nsq

            do nss = 1, nssm

              if(nr_orb_per_sym(nss) == 0) cycle

              nspqrs = ieor(nspqr-1,nss-1)+1

              if(nspqrs /= 1) cycle

#ifdef IJKL_DEBUG
              write(lupri,'(a,4i1,a/)') ' sorting ints for integral class: (',nsp,nsq,nsr,nss,')'
              call flush(lupri)
#endif

              !> non-zero (by symmetry) batch of integrals
              write(datatag,'(a1,i3,a1,i3,a1,i3,a1,i3)') 'p',nsp,'q',nsq,'r',nsr,'s',nss
              tagx = 16
              datadim(1)        = 1; datadim_bound = 1
              call hdf5_get_data(file_id(2),"XXXXXX",datadim,mbatch)
              lbatch            = mbatch

#ifdef IJKL_DEBUG
              write(lupri,*) ' # integrals ==> ',lbatch
              call flush(lupri)
#endif

              allocate(batch(lbatch)); batch = 0

              datadim(1)        = lbatch
              call hdf5_get_data(file_id(2),"XXXXXX",datadim,batch)

#ifdef IJKL_DEBUG
             !write(lupri,*) ' batch is ... ',batch
              call flush(lupri)
#endif

              indn = 0

              do

#ifdef IJKL_DEBUG
                write(lupri,*) 'indn is ... ',indn
                call flush(lupri)
#endif
                if(indn >= lbatch) exit

                !> determine external indices i,j,k,l
                do k  = 1, nr_orb_per_sym(nsr)

                                 ktmp = nr_orb_per_sym(nss)
                  if(nss == nsr) ktmp = k

                  do l = 1, ktmp

                                   istart = 1
                    if(nsr == nsp) istart = k

                    if(nsp == nsr)then
                      noccendi = k
                    else
                      noccendi = nr_orb_per_sym(nsp)
                    end if

                    do i = istart, nr_orb_per_sym(nsp)
!                   do i = 1, noccendi

!                                      itmp = nr_orb_per_sym(nsq)
!                     if(nsq == nsp)   itmp = i
!                                    jstart = 1
!                     if(nsq == nss) jstart = l

!                     prevent double counting of symmetry redundant indices: set upper summation index
                      jstart = 1
!                     if IJKL in same irrep, restrict j to at most i
                      if(nsq == nsp .and. nsr == nss)then !.and.ISYM == KSYM
                        noccendj = i
!                       if LJ in irrep1 and IK in irrep2
                      else if(nss == nsq .and. nsr == nsp .and. nss /= nsp)then
!                       second restriction to prevent double counting, J<=L in KL IJ
                        if(k == i)then
! old from DMRG interface noccendj = l
!test for 3131 or 2121
                          jstart   = l
                          noccendj = nr_orb_per_sym(nsq)
!end of test for 3131 or 2121
                        else ! otherwise all J are needed
                          noccendj = nr_orb_per_sym(nsq)
                        end if
                      else
                        noccendj = nr_orb_per_sym(nsq)
                      end if

!                     do j = jstart, itmp
                      do j = jstart, noccendj


!                       write (lupri,*) 'i,j,k,l            ',i,j,k,l

                        if(nsp == nsq .and. nsr == nss .and. nsp == nsr)then
                          if(nevpt_host_program == 'dalton ')then
                            if(k == i .and. l == i .and. j < i)then
                              cycle
                            end if
                            if(k == i .and. j < l)then
                               cycle
                            end if
                          else
                            if(k == i .and. l == i .and. j < i) cycle
                            if(k == i .and. j < l) cycle
                          end if
                        end if

                        indn   = indn   + 1
                        indtot = indtot + 1

                        ii = ireost(i+iosh(nsp))
                        jj = ireost(j+iosh(nsq))
                        kk = ireost(k+iosh(nsr))
                        ll = ireost(l+iosh(nss))

                        if(jj > ii)then
                          itamp = ii
                          ii    = jj
                          jj    = itamp
                        end if

                        if(ll > kk)then
                          itamp = kk
                          kk    = ll
                          ll    = itamp
                        end if

                        !> swap indices if k,l > i,j
                        if(kk > ii)then
                          itamp  = kk
                          itamp2 = ll
                          kk     = ii
                          ll     = jj
                          ii     = itamp
                          jj     = itamp2
                        end if

                        !> determine internal ij,kl
                        IJ = nevpt_ijkl%num(ii)+jj
                        KL = nevpt_ijkl%num(kk)+ll

#ifdef IJKL_DEBUG
                        write (lupri,*) 'i,j,k,l',indtot,ii,jj,kk,ll,batch(indn)
                        call flush(lupri)
#endif
                        if (ij < kl)then ! swap indices
                          itamp = ij
                          ij    = kl
                          kl    = itamp
                        endif

                        if(nevpt_ijkl%indic(KL) /= nevpt_ijkl%indic(IJ)) then
                          print * ,'Attenzione errore in ext2int2nev2int'
                          print * ,'KLS,IJS',nevpt_ijkl%indic(KL),nevpt_ijkl%indic(IJ)
                        endif

                        twoint(nevpt_ijkl%num(MAX(nevpt_ijkl%jndic(IJ),nevpt_ijkl%jndic(KL))) +     &
                                              MIN(nevpt_ijkl%jndic(IJ),nevpt_ijkl%jndic(KL))  +     &
                                             NDEB(nevpt_ijkl%indic(KL))                     ) = batch(indn)
                      end do
                    end do
                  end do
                end do
              end do

              deallocate(batch)

            end do
          end do
        end do
      end do


      end subroutine ext2int2nev2int
!-------------------------------------------------------------------------------

      subroutine ext1int2nev1int(                     &
                                 h1_out,              &
                                 nnashx,              &
                                 nr_orb_per_sym,      &
                                 nr_sym,              &
                                 nr_actorb_tot,       &
                                 ireost,              &
                                 lupri)
!
!     purpose: reorder integrals from external (Molcas) format to internal
!              nevpt format.
!
!              1-electron integrals: on input  symmetry-reduced list
!                                    on output triangular list
!
!     -----------------------------------------------------------------
      implicit none
      integer :: nnashx
      integer :: nr_sym
      integer :: nr_orb_per_sym(*)
      integer :: nr_actorb_tot
      integer :: ireost(*)
      integer :: lupri
      real*8  :: h1_out(nnashx)
!     -----------------------------------------------------------------
      integer :: isym, jsym, ijsym
      integer :: offset_ij, offset_external, isorb, jsorb
      integer :: nr_act_isym, nr_act_jsym, orb_tmp
      integer :: i_index, j_index
      integer :: loop_counter_i, isorb_tmp
      character (len=200) :: error_message
      real*8, allocatable  :: h1_in(:)
!     -----------------------------------------------------------------

      allocate(h1_in(nnashx)); h1_in = 0
      datadim(1)   = nnashx; datadim_bound = 1
      call hdf5_get_data(file_id(2),"FockMO",datadim,h1_in)

#ifdef IJKL_DEBUG
      write(lupri,*) ' distributing 1e-ints from h1_in:'
      do isym = 1, nnashx
        write(lupri,*) ' h1_in(',isym,') = ',h1_in(isym)
      end do
      write(lupri,*) ' ireost                          '
      do isym = 1, nr_actorb_tot
        write(lupri,*) ' ireost(',isym,') = ',ireost(isym)
      end do
      write(lupri,*) ' itsym                          '
      do isym = 1, nr_actorb_tot
        write(lupri,*) ' itsym(',isym,') = ',itsym(isym)
      end do
#endif

!     1-electron integrals/density matrix elements
!     --------------------------------------------
      orb_tmp         = 0
      isorb_tmp       = 0
      offset_ij       = 0
      offset_external = 0
      ijsym           = 1

      !> insert check for symmetry of 1-electron operator - quit if not totally symmetric for the time being
      if(ijsym .gt. 1)then
        error_message = '*** error in integral resorting for nevpt! 1-e int-operator is not totally symmetric.'
        print *, error_message
        stop 'quit in ext1int2nev1int'
      end if

      do isym = 1, nr_sym

        jsym        = its(isym,ijsym)

#ifdef IJKL_DEBUG
        write(lupri,*) ' isym, jsym',isym, jsym
#endif

        nr_act_jsym = 0

        !> lower triangle only!
        if(jsym >= isym)then

          nr_act_jsym = nr_orb_per_sym(jsym)

          do jsorb = 1, nr_act_jsym
            do isorb = 1, jsorb

#ifdef IJKL_DEBUG
              print *, 'itsym   is ...                          ',itsym(orb_tmp+isorb)
              print *, 'return_ is ...',&
                        its(itsym(ireost(orb_tmp+isorb)),ijsym)
#endif
              !> restrict to orbitals with the correct jsym
              if(its(itsym(ireost(orb_tmp+isorb)),ijsym) /= jsym) cycle

              !> calculate offset on external symmetry packed triangular 1e-array (lower triangular matrix)
              offset_external  =  offset_external + 1

              !> calculate offset on internal triangular 1e-array (lower triangular matrix)
              i_index   = ireost(orb_tmp+isorb)
              j_index   = (ireost(orb_tmp+jsorb)*(ireost(orb_tmp+jsorb)-1))/2
              offset_ij = i_index + j_index

              !> NEVPT: type-symmetry ordering <==> MOLCAS (or other programs) symmetry-type
!                       array ireost takes care ot that...
              h1_out(offset_ij) = h1_in(offset_external)

#ifdef NEVPT_DEBUG
       write(lupri,*)'ij, ext,jsorb,isorb,h1', offset_ij, offset_external,jsorb,isorb,&
                                               h1_in(offset_external), h1_out(offset_ij)
#endif

            end do
          end do
        end if ! jsym >= isym

        !> total # of active orbitals for each jsym
        orb_tmp = orb_tmp + nr_act_jsym
      end do

#ifdef IJKL_DEBUG
      write(lupri,*) ' final 1e-ints from h1_out:'
      do isym = 1, nnashx
        write(lupri,*) ' h1_out(',isym,') = ',h1_out(isym)
      end do
#endif

      deallocate(h1_in)

      end subroutine ext1int2nev1int
!-------------------------------------------------------------------------------
      real*8 function ai(i,j,k,l)

        integer, intent(in) :: i, j, k, l
        if (Do_Cholesky) then
          ai = ai_cholesky(i,j,k,l)
        else
          ai = ai_twoint(i,j,k,l)
        end if
      end function ai

      real*8 function ai_cholesky(i,j,k,l)

        integer, value :: i, j, k, l !indices
!         integer :: ij, kl ! Cholesky vec. indices, size of Cholesky set
!         integer :: ii, tmp
        real*8, external    :: ddot

        ai_cholesky = 0.0D0
        ! calculate indices
!         if (i < j) then
!           tmp = i
!           i = j
!           j = tmp
!         end if
!         if (k < l) then
!           tmp = k
!           k = l
!           l = tmp
!         end if
!
!         ij = (i-1)*i/2 + j
!         kl = (k-1)*k/2 + l

        ! The implementation below allocates the arrays for Cholesky vectors on every call
        ! This might be a performance bottleneck: if so, change the implementation to use a shared scratch
        ! but make sure it's threadsafe (if we're going to be parallel)
!         allocate(vec_ij(nchovec),vec_kl(nchovec))
        ! this turned out to be too resource intensive

!         call hdf5_read_cholesky(dset_id,dspace_id,ij,vec_ij)
!         call hdf5_read_cholesky(dset_id,dspace_id,kl,vec_kl)

        ! Instead, we now read the Cholesky vectors directly from the nevpt_ijkl.cholesky_array which has been filled at initialisation

        ! a_ijkl = (sum_{L=1,nchovec}(vec_ij^L . vec_kl^L))

!          ai_cholesky = dot_product(nevpt_ijkl%cholesky_array(:,ij), &
!          & nevpt_ijkl%cholesky_array(:,kl))

         ai_cholesky = ddot(nchovec,nevpt_ijkl%cholesky_array(1,indice(i,j)),1, &
       &                nevpt_ijkl%cholesky_array(1,indice(k,l)),1)
#ifdef IJKL_DEBUG
        write (*,*) "a(",i,j,k,l,")=", ai_cholesky, "ij=",indice(i,j),"kl=",indice(k,l)
#endif
!         if(allocated(vec_ij)) deallocate(vec_ij)
!         if(allocated(vec_kl)) deallocate(vec_kl)
      contains
        real*8 function indice(a,b)
          integer, intent(in) :: a,b
          indice=max(a,b)*(max(a,b)-1)/2+min(a,b)
        end function indice
      end function ai_cholesky

      real*8 function ai_twoint(i,j,k,l)
!
!     REPERAGE DE L'INTEGRALE BIELECTRONIQUE IJKL DANS LE TABLEAU twoint
!
      integer, intent(in) :: i, j, k, l
      logical :: zsig
      integer :: ij, kl, ijs, kls, ijkls, nij, nadi, ktyp, nkl, ijkl, nnjk

      IF (I.GE.J) GO TO 10
      IJ=nevpt_ijkl%num(J)+I
      GO TO 15
   10 IJ=nevpt_ijkl%num(I)+J
   15 IF (K.GE.L) GO TO 20
      KL=nevpt_ijkl%num(L)+K
      GO TO 25
   20 KL=nevpt_ijkl%num(K)+L
   25 IJS=nevpt_ijkl%indic(IJ)
      KLS=nevpt_ijkl%indic(KL)
      IF (IJS.LE.KLS) GO TO 30
      IJKLS=nevpt_ijkl%num(IJS)+KLS
      GO TO 35
   30 IJKLS=nevpt_ijkl%num(KLS)+IJS
      NIJ=IJ
      IJ=KL
      KL=NIJ
      IJS=KLS
      KLS=nevpt_ijkl%indic(KL)
   35 CONTINUE
      IF (nevpt_ijkl%NAD(IJKLS)) 40,300,45
   40 ZSIG=.TRUE.
      NADI=-nevpt_ijkl%NAD(IJKLS)
      GO TO 50
   45 ZSIG=.FALSE.
      NADI=nevpt_ijkl%NAD(IJKLS)
   50 KTYP=nevpt_ijkl%KT(IJKLS)
      GO TO (60,70,75,80,85,90,95,100,105,62),KTYP
   62 NIJ=nevpt_ijkl%lndic(IJ)
      NKL=nevpt_ijkl%lndic(KL)
      GO TO 63
   60 NIJ=nevpt_ijkl%jndic(IJ)
      NKL=nevpt_ijkl%jndic(KL)
   63 IF (NIJ.LE.NKL) GO TO 65
      IJKL=nevpt_ijkl%num(NIJ)+NKL
      GO TO 200
   65 IJKL=nevpt_ijkl%num(NKL)+NIJ
      GO TO 200
   70 GO TO 300
   75 GO TO 300
   80 GO TO 300
   85 GO TO 300
   90 GO TO 300
   95 GO TO 300
  100 GO TO 300
  105 GO TO 300
  200 NNJK=nevpt_ijkl%NDEB(NADI)+IJKL

      !> pick the integral from the full list
                ai_twoint = nevpt_ijkl%twoint(nnjk)
      IF (ZSIG) ai_twoint = -ai_twoint
#ifdef IJKL_DEBUG
      write (*,*) "a(",i,j,k,l,")=", ai_twoint, "ij=",ij,"kl=",kl
#endif
      RETURN

  300 ai_twoint=0.D0

      END FUNCTION ai_twoint
!------------------------------------------------------------------------------

#ifdef unused
      subroutine itijklt(norb,nsym)
!
!     LECTURE DES INTEGRALES BIELECTRONIQUES
!
      integer, PARAMETER :: NBUF=4096
      integer :: norb, nsym
      real*8  :: XXm(NBUF)
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
      integer :: itamp, indn, i, j, k, l, ij, kl, idum, idum2, ichange

      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)

      k = 1
!--renzo modif to read Daniel's code for norb >255
      do

        if(k == 0) exit

        call filesplit('READ',50,16*nbuf,ichange,idum2)

        read(50)ind,jnd,knd,lnd,xxm !renzo new

        do indn=1,nbuf
          i=ind(indn)
          j=jnd(indn)
          k=knd(indn)
          l=lnd(indn)

#ifdef IJKL_DEBUG
          write (6,*) 'i,j,k,l',indn,i,j,k,l,XXm(indn)
#endif

          if (i.gt.norb.or.j.gt.norb.or.k.gt.norb.or.l.gt.norb) then
            call flush(6)
            stop 'one of the integral indices is out-of-bounds.'
          endif
          if(k == 0) exit
!--renzo modif end
          IJ = nevpt_ijkl%num(i)+J
          KL = nevpt_ijkl%num(K)+L

          if (ij < kl)then ! swap indices
            itamp = ij
            ij    = kl
            kl    = itamp
          endif

          if(nevpt_ijkl%indic(KL) /= nevpt_ijkl%indic(IJ)) then
            print * ,'Attenzione errore in ITIJKLt'
            print * ,'KLS,IJS',nevpt_ijkl%indic(KL),nevpt_ijkl%indic(IJ)
          endif

          nevpt_ijkl%twoint(nevpt_ijkl%num(MAX(nevpt_ijkl%jndic(IJ),nevpt_ijkl%jndic(KL))) +            &
                                           MIN(nevpt_ijkl%jndic(IJ),nevpt_ijkl%jndic(KL))  +            &
                                          nevpt_ijkl%NDEB(nevpt_ijkl%indic(KL))                     ) = XXm(indn)
!         write (6,7) indn,xxm(indn),i,j,k,l
        end do
      end do
      end subroutine itijklt
#endif
!-------------------------------------------------------------------------------

end module ijkl_utils
