module koop_matrices

#ifdef _OPENMP_
  use omp_lib
#endif
  use eexx_functions
  use info_symmetry
  use hdf5_utils

  implicit none

  public koop_matrices_ctl

  contains

  subroutine koop_matrices_ctl(metat,nact,ncore,atwo,daaa,taa,dal,f,aj)

    integer, intent(in)   :: nact, ncore, metat
    real*8 , intent(in)   :: daaa(metat,nact,nact,nact,nact,nact,nact)
    real*8 , intent(in)   :: taa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: dal(metat,nact,nact)
    real*8 , intent(in)   :: atwo(nact,nact,nact,nact)
    real*8 , intent(in)   :: aj(*)
    real*8 , intent(inout):: f(*)

    call koopE_driver(atwo,taa,dal,nact,ncore,metat,f,aj)

    call koop_matrices_driver(atwo,daaa,taa,dal,nact,ncore,metat,f,aj)


  end subroutine koop_matrices_ctl
! ----------------------------------------------------------------------

  subroutine koopE_driver(atwo,taa,dal,nact,ncore,metat,f,aj)

    integer, intent(in)   :: nact, ncore, metat
    real*8 , intent(in)   :: taa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: dal(metat,nact,nact)
    real*8 , intent(in)   :: atwo(nact,nact,nact,nact)
    real*8 , intent(in)   :: aj(*)
    real*8 , intent(inout):: f(*)

    real*8 , allocatable :: coopipa(:,:,:), coopeaa(:,:,:)

    allocate(coopipa(metat,nact,nact),coopeaa(metat,nact,nact))
    coopipa = 0
    coopeaa = 0

    !print*,' Chiamata di koop'
    call koopE(atwo,coopipa,coopeaa,taa,dal,nact,ncore,metat,f,aj)
    !> put data to file
    datadim(1)   = metat; datadim(2:3) = nact; datadim_bound = 3
    call hdf5_put_data_dp(file_id(1), "1-EAKO" , datadim, coopeaa)
    call hdf5_put_data_dp(file_id(1), "1-IPKO" , datadim, coopipa)
    call hdf5_put_data_dp(file_id(1), "1-RDM"  , datadim, dal)

    deallocate(coopipa,coopeaa)
  end subroutine koopE_driver
! ----------------------------------------------------------------------

  subroutine koop_matrices_driver(atwo,daaa,taa,dal,nact,ncore,metat,f,aj)

    integer, intent(in)   :: nact, ncore, metat
    real*8 , intent(in)   :: daaa(metat,nact,nact,nact,nact,nact,nact)
    real*8 , intent(in)   :: taa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: dal(metat,nact,nact)
    real*8 , intent(in)   :: atwo(nact,nact,nact,nact)
    real*8 , intent(in)   :: aj(*)
    real*8 , intent(inout):: f(*)

    real*8 , allocatable  :: koopaa(:,:,:,:,:)
    integer               :: i, j, k, l, istate

    allocate(koopaa(metat,nact,nact,nact,nact))
    koopaa = 0

    !print*,'costruzione matrici ro2 di buca'
    do i=1,nact
       do j=1,nact
          do k=1,nact
             do l=1,nact
                do istate=1,metat
                   koopaa(istate,i,j,k,l)=ro2t(i,j,k,l,taa,dal,nact,istate,metat)
                enddo
             enddo
          enddo
       enddo
    enddo

    !> put data to file
    datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5
    call hdf5_put_data_dp(file_id(1), "2-HRDM" , datadim, koopaa)

    !print*,' Chiamata di koopman2'
    call koop2E_driver(atwo,daaa,taa,dal,koopaa,nact,ncore,metat,f,aj)

    !print*,' Chiamata di koopman0pE'
    call koopman0pE_driver(atwo,daaa,taa,dal,koopaa,nact,ncore,metat,f,aj)

    deallocate(koopaa)
  end subroutine koop_matrices_driver

! ----------------------------------------------------------------------

  subroutine koop2E_driver(atwo,daaa,taa,dal,koopaa,nact,ncore,metat,f,aj)

    integer, intent(in)   :: nact, ncore, metat
    real*8 , intent(in)   :: daaa(metat,nact,nact,nact,nact,nact,nact)
    real*8 , intent(in)   :: taa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: dal(metat,nact,nact)
    real*8 , intent(in)   :: atwo(nact,nact,nact,nact)
    real*8 , intent(inout):: koopaa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: aj(*)
    real*8 , intent(inout):: f(*)

    real*8 , allocatable  :: koopeaa(:,:,:,:,:)

    allocate(koopeaa(metat,nact,nact,nact,nact))

    call koop2E(atwo,daaa,taa,dal,koopaa,koopeaa,nact,ncore,metat,f,aj)
    !> put data to file
    datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5
    call hdf5_put_data_dp(file_id(1), "2-EAKO" , datadim, koopeaa)
    call hdf5_put_data_dp(file_id(1), "2-IPKO" , datadim, koopaa)
    call hdf5_put_data_dp(file_id(1), "2-RDM"  , datadim, taa)

    deallocate(koopeaa)

  end subroutine koop2E_driver

! ----------------------------------------------------------------------

  subroutine koopman0pE_driver(atwo,daaa,taa,dal,koopaa,nact,ncore,metat,f,aj)

    integer, intent(in)   :: nact, ncore, metat
    real*8 , intent(in)   :: daaa(metat,nact,nact,nact,nact,nact,nact)
    real*8 , intent(in)   :: taa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: dal(metat,nact,nact)
    real*8 , intent(in)   :: atwo(nact,nact,nact,nact)
    real*8 , intent(inout):: koopaa(metat,nact,nact,nact,nact)
    real*8 , intent(in)   :: aj(*)
    real*8 , intent(inout):: f(*)

    real*8 , allocatable  :: koopbb(:,:,:,:,:), koopab(:,:,:)

    allocate(koopbb(metat,nact,nact,nact,nact),koopab(metat,nact,nact))
    koopbb = 0
    koopab = 0

    call koopman0pE(atwo,daaa,taa,dal,koopaa,koopbb,koopab,nact,ncore,metat,f,aj)
    call fdiff_F(koopab,nact,metat)
    !> put data to file
    datadim(1)   = metat; datadim(2:5) = nact; datadim_bound = 5
    call hdf5_put_data_dp(file_id(1), "gk0pa" , datadim, koopaa)
    call hdf5_put_data_dp(file_id(1), "gk0pd" , datadim, koopbb)
    datadim_bound = 3
    call hdf5_put_data_dp(file_id(1), "gk0pf" , datadim, koopab)

    deallocate(koopbb,koopab)

  end subroutine koopman0pE_driver

! ----------------------------------------------------------------------

  subroutine koopE(atwo,koopipa,koopeaa,taa,dal,nact,ncore,metat,f,aj)

      real*8 , intent(in)    :: aj(*)
      real*8 , intent(in)    :: taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      real*8 , intent(in)    :: atwo(nact,nact,nact,nact)
      integer, intent(in)    :: nact,ncore,metat
!---F e' la matrice dei monoelettronici!
      real*8 , intent(inout) :: f(*)
      real*8 , intent(out)   :: koopipa(metat,nact,nact),koopeaa(metat,nact,nact)
      !> local variables
      real*8 , allocatable   :: dumm(:)
      integer                :: a,b,ac,bc,c,ap,apc,cc,d,dc
      integer                :: indice,i,j
      integer                :: ie, ind, irind, istate, jab
      real*8                 :: dum, dum2
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
!
!     Calcolo le matrici di Koopmans
!     Caso -1: matrici di Koopmans per la ionizzazione attiva
!
!     La parte relativa alle singole differenze attive
!     e alle doppie differenze, una attiva e una inattiva
!     viene calcolata mediante rho_ab
!
!----building Dyall's h eff----------
      allocate(dumm(metat))
      do a=1,nact
!cele
       ac=a+ncore
         do b=a,nact
!cele
          bc=b+ncore
            jab=indice(ac,bc)
            f(jab)=aj(jab)
         enddo
      enddo
      do a=1,nact
         ac=a+ncore
         do ap=1,nact
            apc=ap+ncore
            do istate=1,metat
               koopipa(istate,ap,a)=0.d0
               koopeaa(istate,ap,a)=0.d0
            enddo
!            print*,'indici a e ap',a,ap
            do c=1,nact
               cc=c+ncore
               irind=indice(ac,cc)
               dum=f(irind)
               do istate=1,metat
                  koopipa(istate,ap,a)=koopipa(istate,ap,a)-dum*dal(istate,ap,c)
               enddo
               do d=1,nact
                  do ie=1,nact
                     do istate=1,metat
                        koopipa(istate,ap,a)=koopipa(istate,ap,a)-atwo(c,ie,a,d)*taa(istate,ap,c,d,ie)
                     enddo
                  enddo
               enddo
            enddo
            do istate=1,metat
               koopipa(istate,ap,a)=2.d0*koopipa(istate,ap,a)
            enddo
         enddo
      enddo
      do ap=1,nact
         do a=1,nact
            ind=indice(ap+ncore,a+ncore)
            do istate=1,metat
               dumm(istate)=2.d0*f(ind)
            enddo
            do d=1,nact
               do ie=1,nact
                  dum2=2.d0*atwo(d,ie,ap,a)-atwo(ap,ie,d,a)
                  do istate=1,metat
                     dumm(istate)=dumm(istate)+dum2*dal(istate,d,ie)
                  enddo
               enddo
            enddo
            do istate=1,metat
               koopeaa(istate,ap,a)=koopipa(istate,ap,a)+2.d0*dumm(istate)
            enddo
         enddo
      enddo
      deallocate(dumm)
   end subroutine koopE
!----------------------------------------------------------------------

   subroutine koop2E(atwo,daaa,taa,dal,koopaa,koopeaa,nact,ncore,metat,f,aj)

      real*8 , intent(in)    :: aj(*)
      real*8 , intent(in)    :: daaa(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in)    :: taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      real*8 , intent(in)    :: atwo(nact,nact,nact,nact)
      integer, intent(in)    :: nact,ncore,metat
!---F e' la matrice dei monoelettronici!
      real*8 , intent(inout) :: f(*)
      real*8 , intent(out)   :: koopaa(metat,nact,nact,nact,nact),koopeaa(metat,nact,nact,nact,nact)
      !> local variables
      integer                :: a,b,c,d,ac,bc,cc,dc,ap,bp,apc,bpc
      integer                :: indice,i,j
      integer                :: ie, ief, isa, isab, isb, isc, iscd, isd, istate, jab, jad, jbc, jbd, jac
      real*8                 :: dum, dint, dint2, dum1, dum2

      !> statement function
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      !> building heff': in AJ is heff
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
      do a=1,nact
         ac=a+ncore
         isa=itsym(ac)
         do b=1,nact
            bc=b+ncore
            isb=itsym(bc)
            isab=its(isa,isb)
            do ap=1,nact
               apc=ap+ncore
               isd=itsym(apc)
               do bp=1,nact
                  bpc=bp+ncore
                  isc=itsym(bpc)
                  iscd=its(isc,isd)
                  do istate=1,metat
                     koopaa(istate,ap,bp,a,b)=0.d0
                  enddo
                  if (its(isab,iscd).ne.1) goto 1
                  do d=1,nact
                     dc=d+ncore
                     jad=indice(ac,dc)
                     jbd=indice(bc,dc)
                     do istate=1,metat
                        koopaa(istate,ap,bp,a,b)=koopaa(istate,ap,bp,a,b)-f(jad)*taa(istate,ap,bp,d,b)&
                                                                         -f(jbd)*taa(istate,ap,bp,a,d)
                     enddo
                     do ie=1,nact
                        do ief=1,nact
                           dint=atwo(d,ie,a,ief)
                           dint2=atwo(d,ie,b,ief)
                           do istate=1,metat
                              dum=daaa(istate,ap,bp,d,ief,b,ie)
                              if(d.eq.ief)dum=dum+0.5d0*taa(istate,ap,bp,ie,b)
                              if(d.eq.b)dum=dum+0.5d0*taa(istate,ap,bp,ief,ie)
                              koopaa(istate,ap,bp,a,b)=koopaa(istate,ap,bp,a,b)-dint*dum
                              dum=daaa(istate,ap,bp,d,a,ief,ie)
                              if(d.eq.a)dum=dum+0.5d0*taa(istate,ap,bp,ie,ief)
                              if(d.eq.ief)dum=dum+0.5d0*taa(istate,ap,bp,a,ie)
                              koopaa(istate,ap,bp,a,b)=koopaa(istate,ap,bp,a,b)-dint2*dum
                           enddo
                        enddo
                     enddo
                  enddo
 1             enddo
            enddo
         enddo
      enddo
      do a=1,nact
         ac=a+ncore
         isa=itsym(ac)
         do b=1,nact
            bc=b+ncore
            isb=itsym(bc)
            isab=its(isa,isb)
            do ap=1,nact
               apc=ap+ncore
               isd=itsym(apc)
               do bp=1,nact
                  bpc=bp+ncore
                  isc=itsym(bpc)
                  iscd=its(isc,isd)
                  do istate=1,metat
                     koopeaa(istate,ap,bp,a,b)=0.d0
                  enddo
                  if (its(isab,iscd).ne.1) goto 2
                  do c=1,nact
                     cc=c+ncore
                     jac=indice(ac,cc)
                     jbc=indice(bc,cc)
                     do istate=1,metat
                        dum1=ro2t(ap,bp,c,b,taa,dal,nact,istate,metat)
                        dum2=ro2t(ap,bp,a,c,taa,dal,nact,istate,metat)
                        koopeaa(istate,ap,bp,a,b)=koopeaa(istate,ap,bp,a,b)+f(jac)*dum1+f(jbc)*dum2
                     enddo
                     do d=1,nact
                        do ie=1,nact
                           dint=atwo(c,ie,d,a)
                           do istate=1,metat
                              dum=-ro3t(ap,bp,ie,d,b,c,daaa,taa,dal,nact,istate,metat)
                              if(c.eq.ie)dum=dum+2.d0*ro2t(ap,bp,d,b,taa,dal,nact,istate,metat)
                              if(d.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,c,b,taa,dal,nact,istate,metat)
                              if(b.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,d,c,taa,dal,nact,istate,metat)
                              koopeaa(istate,ap,bp,a,b)=koopeaa(istate,ap,bp,a,b)+dint*dum
                           enddo
                           dint=atwo(c,ie,d,b)
                           do istate=1,metat
                              dum=-ro3t(ap,bp,ie,a,d,c,daaa,taa,dal,nact,istate,metat)
                              if(c.eq.ie)dum=dum+2.d0*ro2t(ap,bp,a,d,taa,dal,nact,istate,metat)
                              if(a.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,c,d,taa,dal,nact,istate,metat)
                              if(d.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,a,c,taa,dal,nact,istate,metat)
                              koopeaa(istate,ap,bp,a,b)=koopeaa(istate,ap,bp,a,b)+dint*dum
                           enddo
                        enddo
                     enddo
                  enddo
 2             enddo
            enddo
         enddo
      enddo
   end subroutine koop2E
!-------------------------------------------------------------------

   subroutine koopman0pE(atwo,daaa,taa,dal,koopaa,koopbb,koopab,nact,ncore,metat,f,aj)

      real*8 , intent(in)    :: aj(*)
      real*8 , intent(in)    :: daaa(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in)    :: taa(metat,nact,nact,nact,nact),dal(metat,nact,nact)
      real*8 , intent(in)    :: atwo(nact,nact,nact,nact)
      integer, intent(in)    :: nact,ncore,metat
!---F e' la matrice dei monoelettronici!
      real*8 , intent(inout) :: f(*)
      real*8 , intent(out)   :: koopaa(metat,nact,nact,nact,nact)
      real*8 , intent(out)   :: koopbb(metat,nact,nact,nact,nact), koopab(metat,nact,nact)
      !> local variables

      integer                :: a,b,c,d,ac,bc,cc,dc,ap,bp,apc,bpc
      integer                :: ibc, ica, ie, iec, istate, jab
      integer                :: indice,i,j
      real*8                 :: dum, dum1, dum2

      !> statement function
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      !> build Dyall Hamiltonian in f
      do a=1,nact
         ac=a+ncore
         do b=a,nact
            bc=b+ncore
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nact
               cc=c+ncore
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-0.5d0*dum
            enddo
         enddo
      enddo
      do a=1,nact
         ac=a+ncore
         do b=1,nact
            bc=b+ncore
            do ap=1,nact
               apc=ap+ncore
               do bp=1,nact
                  bpc=bp+ncore
                  do istate=1,metat
                     koopaa(istate,bp,ap,a,b)=0.d0
                  enddo
                  do c=1,nact
                     cc=c+ncore
                     ica=indice(cc,ac)
                     ibc=indice(bc,cc)
                     do istate=1,metat
                        koopaa(istate,bp,ap,a,b)=koopaa(istate,bp,ap,a,b&
                             )+2.d0*(f(ica)*ee2(bp,ap,c,b,taa,dal,nact  &
                             ,istate,metat)-f(ibc)*ee2(bp,ap,a,c        &
                             ,taa,dal,nact,istate,metat))
                     enddo
                     do d=1,nact
                        dc=d+ncore
                        do ie=1,nact
                           iec=ie+ncore
                           do istate=1,metat
                              dum1=eee2(bp,ap,c,ie,d,b,daaa             &
                                   ,taa,dal,nact,istate,metat)          &
                                   +eee2(bp,ap,d,b,c,ie,daaa,taa,dal    &
                                   ,nact,istate,metat)                   
                              dum2=eee2(bp,ap,a,ie,c,d,daaa,taa,dal,nact&
                                   ,istate,metat)+eee2(bp,ap,c,d        &
                                   ,a,ie,daaa,taa,dal,nact,istate,metat)
                              koopaa(istate,bp,ap,a,b)=koopaa(istate,bp &
                                   ,ap,a,b)+atwo(c,ie,d,a)*dum1-atwo(b  &
                                   ,ie,c,d)*dum2
                           enddo
                        enddo
                     enddo
                  enddo
                  do istate=1,metat
                     koopbb(istate,bp,ap,a,b)=0.d0
                  enddo
                  do c=1,nact
                     cc=c+ncore
                     ica=indice(cc,ac)
                     ibc=indice(bc,cc)
                     do istate=1,metat
                        dum1=-ee2(a,ap,bp,c,taa,dal,nact,istate,metat)
                        if(a.eq.ap)dum1=dum1+2.0d0*dal(istate,bp,c)
                        if(bp.eq.ap)dum1=dum1+dal(istate,a,c)
                        dum2=-ee2(c,ap,bp,b,taa,dal,nact,istate,metat)
                        if(ap.eq.c)dum2=dum2+2.0d0*dal(istate,bp,b)
                        if(bp.eq.ap)dum2=dum2+dal(istate,c,b)
                        koopbb(istate,bp,ap,a,b)=koopbb(istate,bp,ap,a,b)-f(ibc)*dum1+f(ica)*dum2
                     enddo
                     do d=1,nact
                        dc=d+ncore
                        do ie=1,nact
                           iec=ie+ncore
                           do istate=1,metat
                              dum1=-eee2(c,ie,a,ap,bp,d,daaa,taa,dal     &
                                   ,nact,istate,metat)-eee2(a,ap         &
                                   ,bp,d,c,ie,daaa,taa,dal,nact,istate   &
                                   ,metat)
                              if(ap.eq.a)dum1=dum1+2.0d0*(ee2(c,ie,bp,d    &
                                   ,taa,dal,nact,istate,metat)           &
                                   +ee2(bp,d,c,ie,taa,dal,nact,istate    &
                                   ,metat))
                              if(bp.eq.ap)dum1=dum1+ee2(c,ie,a,d,taa,dal &
                                   ,nact,istate,metat)+ee2(a,d,c         &
                                   ,ie,taa,dal,nact,istate,metat)
                              if(ap.eq.c)then
                                 dum1=dum1-ee2(a,ie,bp,d,taa,dal,nact,istate,metat)
                                 if(ie.eq.a)dum1=dum1+2.0d0*dal(istate,bp,d)
                              endif
                              if(bp.eq.ie)then
                                 dum1=dum1+ee2(a,ap,c,d,taa,dal,nact,istate,metat)
                                 if(ap.eq.a)dum1=dum1-2.0d0*dal(istate,c,d)
                              endif
                              dum2=-eee2(c,ie,d,ap,bp,b,daaa             &
                                   ,taa,dal,nact,istate,metat)           &
                                   -eee2(d,ap,bp,b,c,ie,daaa,taa,dal     &
                                   ,nact,istate,metat)
                              if(ap.eq.d)dum2=dum2+2.0d0*(ee2(c,ie,bp,b    &
                                   ,taa,dal,nact,istate,metat)           &
                                   +ee2(bp,b,c,ie,taa,dal,nact,istate    &
                                   ,metat))
                              if(bp.eq.ap)dum2=dum2+ee2(c,ie,d,b,taa,    &
                                   dal,nact,istate,metat)+ee2(d          &
                                   ,b,c,ie,taa,dal,nact,istate,metat)
                              if(ap.eq.c)then
                                 dum2=dum2-ee2(d,ie,bp,b,taa,dal,nact,istate,metat)
                                 if(d.eq.ie)dum2=dum2+2.0d0*dal(istate,bp,b)
                              endif
                              if(bp.eq.ie)then
                                 dum2=dum2+ee2(d,ap,c,b,taa,dal,nact,istate,metat)
                                 if(d.eq.ap)dum2=dum2-2.0d0*dal(istate,c,b)
                              endif
                              koopbb(istate,bp,ap,a,b)=koopbb(istate,bp,ap,a,b)-0.5d0*atwo(c,ie,b,d)*dum1&
                                                      +0.5d0*atwo(c,ie,d,a)*dum2
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
#ifdef DEBUG_DMRG_NEVPT
            print *, 'blubb heavy DEBUG for state 1: koopab'
#endif
            do istate=1,metat
               koopab(istate,a,b)=0.d0
            enddo
            do c=1,nact
               cc=c+ncore
               ica=indice(cc,ac)
               ibc=indice(bc,cc)
               do istate=1,metat
                  dum1=dal(istate,c,b)
                  dum2=dal(istate,a,c)
                  koopab(istate,a,b)=koopab(istate,a,b)+2.0d0*(f(ica)*dum1-f(ibc)*dum2)
#ifdef DEBUG_DMRG_NEVPT
                  if(abs(koopab(istate,a,b)) > 1.0d-16 .and. istate == 1)&
                  print '(i1,1i2,3x,d19.12,a,d19.12,2x,d19.12,a,2i2)',&
                  a,b,koopab(istate,a,b), ' dum val ',dum1,dum2, &
' aux idx ',ica,ibc
#endif
               enddo
               do d=1,nact
                  dc=d+ncore
                  do ie=1,nact
                     iec=ie+ncore
                     do istate=1,metat
                        dum1=ee2(c,ie,d,b,taa,dal,nact,istate,metat)+ee2(d,b,c,ie,taa,dal,nact,istate,metat)
                        dum2=ee2(a,ie,c,d,taa,dal,nact,istate,metat)+ee2(c,d,a,ie,taa,dal,nact,istate,metat)
                        koopab(istate,a,b)=koopab(istate,a,b)+atwo(c,ie,d,a)*dum1-atwo(b,ie,c,d)*dum2
#ifdef DEBUG_DMRG_NEVPT
                        if(abs(koopab(istate,a,b)) > 1.0d-16 .and. istate == 1)&
                        print '(i1,1i2,3x,d19.12,a,d19.12,2x,d19.12,a,3i2)',&
                        a,b,koopab(istate,a,b), ' dum val ',dum1,dum2, &
' aux idx ',c,d,ie
#endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

#ifdef DEBUG_DMRG_NEVPT
      print *, 'blubb DEBUG koopab'
      do istate=1,metat
      print *, 'koopab for state ',istate
      do a=1,nact
        do b=1,nact
          if(abs(koopab(istate,a,b)) > 1.0d-16)&
          print '(i1,1i2,3x,d19.12)',&
          a,b,koopab(istate,a,b)
        enddo
      enddo

      enddo ! istate
      call flush(6)
#endif
   end subroutine koopman0pE
!-------------------------------------------------------------------

      subroutine fdiff_F(koopab,nact,metat)

      real*8 , intent(in) :: koopab(metat,nact,nact)
      integer, intent(in) :: nact, metat

      real*8 , parameter  :: threshold = 1.0d-7
      integer             :: istate, i, j
      integer             :: index_ldiff(2)
      real*8              :: bigdiff, val

      index_ldiff(1:2) = 0
      bigdiff=0.d0
      do istate=1,metat
        do i=1,nact
           do j=1,nact
             val=abs(koopab(istate,i,j))
             if(val .gt.bigdiff)then
               bigdiff=val
               index_ldiff(1) = i; index_ldiff(2) = j
             end if
           enddo
        enddo
        if(abs(bigdiff) > threshold)then
          write(6,'(/a,f14.8,a,i3)') ' maximum value of F matrix is ',bigdiff, ' for state ',istate
          write(6,'(a,2i4/)') ' matrix indices of F(i,j): ',index_ldiff(1:2)
        end if
      enddo

      end subroutine fdiff_F
!------------------------------------------------------------------------------

end module koop_matrices
