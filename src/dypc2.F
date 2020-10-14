C-----------------------------------------------------------------------
C
C     Programma per il calcolo della correzione perturbativa
C     dell'energia al secondo ordine secondo la formulazione
C     ``n-electron valence state perturbation theory'' nella versione
C     parzialmente  contratta.
C
C     SOLO H DYALL, STRONGLY E PARTIALLY CONTRACTED
C
C     LAVORA CON FUNZIONI IN PERTURBAZIONE AUTOSTATI DI S2
C
C     Citazioni per NEVPT
C
C     Primo articolo NEVPT
C     - C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J.-P. Malrieu,
C       J. Chem. Phys. 114(23), 10252-10264, (2001).
C
C     Introduzione operatori di Koopmans
C     - C. Angeli, R. Cimiraglia, J-P Malrieu,
C       Chem. Phys. Lett. 350(3-4), 297-305, (2001).
C
C     Formulazione "spin-free"
C     - C. Angeli, R. Cimiraglia, J-P Malrieu,
C       J. Chem. Phys. 117(20), 9138-9153, (2002).
C
C     QD-NEVPT2
C     - C. Angeli, S. Borini, M. Cestari, R. Cimiraglia,
C       J. Chem. Phys. 121(9), 4043-4049, (2004).
C
C     NEVPT3
C     - C. Angeli,  B. Bories, A. Cavallini, R. Cimiraglia,
C       J. Chem. Phys., 124(5), 054108(1-8), (2006).
C
C     Short review of NEVPT
C     - C. Angeli,  M. Pastore, R. Cimiraglia,
C       Theor. Chem. Acc., 117(5-6), 743-754, (2007).
C
C
C-----------------------------------------------------------------------
C
C     ENTREE DES DONNEES
C
C     CARTE 1
C     TITRE   (80 CHARACTERES)
C
C
C     NAMELIST FILES
C     ZVERBOSE T intermediate results are printed
C     ZRO      T perturbed density matrix is calculated (Partially
C     contracted)
C     Three matrices are printed in file 27 (RHOFILE)
C     RHO(i,r),RHO(d,r),RHO(i,d)
C
C     NAMELIST ICINP
C     ZTOUL  = Gli integrali sono letti da MOTRA+MOLCAS+MOLCOST secondo
C     il formalismo di Tolosa. Default = .false.
C     FROM  = if FROM='DALTON' koopmans matrices are scaled with factor
C     ZTHIRD: if T the data needed for the thrd order are stored
C                  in file32 (default)
C             if F file32 is not generated.
C     ZPRT  : if T verbose print in output
C             if F standard (concise) output (default)
C     NCORE : number of doubly inactive orbitals
C     NACT  : number of active orbitals
C     IGELO  = INDICI DEGLI ORBITALI DOPPIAMENTE OCCUPATI
C              CHE SARANNO MANTENUTI INATTIVI NELLA PARTE PERTURBATIVA
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     *     NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     *     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     *     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     *     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      character*80 label
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      character*1280 namel
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      common/debug/zdebug
      common /fil/ file12,file13
      character*64 projf,file12,file13,dir
      COMMON /CPU/ TCPU,PCPU
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM(10)
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
      common/thresho/thresh
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     *     e2mpp(8),psimpp(8),igelo,zverbose,zro
      allocatable eps(:)
      INTEGER*2 NEP,TROUP,PARTP
      integer*1 iocc
      allocatable iocc(:,:)
c     integer*2 icomp,iconf,isegno
c     ele
      common /big/ x,ijkl2,ijkl4,ijkl3,rijkl
      pointer x
      dimension x(:)
      pointer ijkl2,ijkl3,ijkl4
      pointer rijkl
      dimension ijkl2(:),ijkl4(:),ijkl3(:),rijkl(:)
      integer*2 ijkl2
      character*3 ijkl3
      character*6 from
      real*4 rijkl
      common /toul/ ztoul
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
      allocatable gkp1a(:,:),gkm1a(:,:)
      allocatable gkp2aa(:,:,:,:)
      allocatable gkm2aa(:,:,:,:)
      allocatable gk0pa(:,:,:,:),gk0pd(:,:,:,:),gk0pf(:,:)
      allocatable ro3(:,:,:,:,:,:),amat(:,:,:,:,:,:),bmat(:,:,:,:),
     *     cmat(:,:,:,:),dmat(:,:)
      allocatable atmat(:,:,:,:,:,:),btmat(:,:,:,:),
     *     ctmat(:,:,:,:),dtmat(:,:)
      allocatable rho(:,:),w(:),work(:)
C     
C     Density matrices SPIN LESS
C     
C     da  = one-particle density matrix
C     daa = two-particle density matrix
C     daat= two-hole density matrix
C     
      allocatable da(:,:)
      allocatable daa(:,:,:,:)
      allocatable daat(:,:,:,:)
      allocatable fockc(:,:),fockv(:,:)
      DIMENSION NOS(20)
      common /deephole/ ndhb,ndha,zdh
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c     
      print*,'Version compiled with intel IFC compiler'
      zfac=.true.
      ntotdet=0
      ntotconf=0
      ntotcap=0
      zverbose=.false.
      zro=.false.
      from='      '
      write (6,*) '   '
      CALL CPUT(-1)
c     WRITE (6,'(/1X,130(''*''))')
c     WRITE (6,'(/1X,15A8)') ('  CIPPI ',I=1,15)
c     WRITE (6,'(/1X,130(''*''))')
      print '(/,a)',repeat('*',80)
      print '(/,a)',repeat(' NEVPT2 ',10)
      print '(/,a)',repeat('*',80)
      CALL OFILES(zbin,zverbose,zro)
      print '(/,a)',repeat('*',80)
c     WRITE (6,'(/1X,130(''*''))')
      ZTR=.TRUE.
      ZFA=.FALSE.
C     
C     LECTURE DES INTEGRALES
C     
      do i=1,8
       e2mpp(i)=0.d0
       psimpp(i)=0.d0
      enddo
      READ (25) NSYM,NORB,NOC,NDEGEN
      rewind 25
      allocate(itsym(norb+norb+1))
      READ (25) NSYM,NORB,NOC,NDEGEN,(ITSYM(I),I=1,NORB),
     *     nao,ngelo,ngelv
      READ (25) ((ITS(I,J),I=1,NSYM),J=1,NSYM)
      read (4) NORB,NOCB,NOCA,METAT,NCF,II,I
      rewind 4
      allocate(ne(ncf))
      allocate(nd(ncf))
      iditrpa=ii
      allocate(trou(ii))
      allocate(part(ii))
      allocate(e(metat))
      allocate(emp(metat))
      allocate(een(metat))
      allocate(e2en(metat))
      allocate(e2mp(metat))
      allocate(psimp(metat))
      allocate(psien(metat))
      allocate(fmpb(norb+norb))
      allocate(zgel(norb+norb+1))
      allocate(igelo(norb+1))
      allocate(eps(norb))
      allocate(zact(norb+1))
      allocate(num(norb*(norb+1)/2))
      allocate(indic(norb*(norb+1)/2))
      allocate(jndic(norb*(norb+1)/2))
      allocate(lndic(norb*(norb+1)/2))
      allocate(onel(norb*(norb+1)/2))
      read (4) NORB,NOCB,NOCA,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),
     *     (TROU(J),PART(J),J=1,II),(C,J=1,I),(E(J),J=1,METAT),
     *     (EMP(J),J=1,METAT),(FMPB(J),J=1,norb*2),
     *     NNROT,(brdmp,J=1,NNROT),MMETAT,(ZHEFF,J=1,MMETAT),
     *     zhandy,zthree,FACTOR,escf
c     * zhandy,zthree,FACTOR,escf,zenb,(een(j),j=1,metat)

c     Lettura ztoul
c     

      rewind 3
      call cerca(3,'ICINP',itrovo)
      if (itrovo.eq.0) then
         write (6,*) '   *** NAMELIST &ICINP NOT FOUND ***'
         stop 12
      endif

      ztoul=.false.
      ncore=-1
      nact=0
      zmotra=.false.
      zthird=.FALSE.
      zprt=.false.
      call lecnam('&ICINP',namel,3)
      call namell(namel,ztoul,'ZTOUL',1)
      call namell(namel,zthird,'ZTHIRD',1)
      call namell(namel,zprt,'ZPRT',1)
      call nameli(namel,ncore,'NCORE',4)
      call nameli(namel,nact,'NACT',4)
      call namela(namel,from,'FROM')
      jd1=norb
      CALL REIJKL(factor,norb,nsym,itsym,its,nos,jd1)
      finv=1.d0/factor
      allocate(ijkl4(ni4))

      if (ztoul) then
         CALL ITIJKLt(factor,norb,nsym,nos,nijkl,
     *        ijkl4,ijkl2,rijkl,ijkl3)
      else
         CALL ITIJKL(factor,norb,nsym,nos,nijkl,
     *        ijkl4,ijkl2,rijkl,ijkl3)
      endif

      nter=norb*(norb+1)/2
      if (ztoul) then
c     read(26)
c     do
c     read(26) label
c     if(index(label,'==ONEINT MATRIX (TOTAL)').ne.0) exit
c     read(26)
c     read(26)
c     enddo

c     read(26) ! norbttMono,nsymMono,isymMono(1:nsymMono)
c     read(26) onel(1:nter)
         read (26) (onel(l),l=1,nter)
      else
         call filesplit('READ',50,0,ichange,idum2)
         read (50)
         call filesplit('READ',50,8*nter,ichange,idum2)
         read (50) (onel(i),i=1,nter)
      endif
      close(unit=50)
      do i=1,nter
       onel(i)=onel(i)*factor
      enddo

      zcas=.true.
      do i=1,norb+1
       zact(i)=.false.
      enddo
c      print*,'The zero order wave funtion is of CAS type'

      if(nact.eq.0)then
         ii=nd(ncf)+ne(ncf)
         do i=1,II
          if (trou(i).gt.norb)then
             zact(trou(i)-norb)=.true.
          else
             zact(trou(i))=.true.
          endif
          if (part(i).gt.norb)then
             zact(part(i)-norb)=.true.
          else
             zact(part(i))=.true.
          endif
         enddo
      else
         if (ncore.eq.-1) then
            write (6,*)
     $           'Error: NCORE must be set if NACT is set in input'
            stop
         endif
         do i=1,nact
          zact(i+ncore)=.true.
         enddo
      endif

c      write (6,*)
c      write (6,*)
c      do i=1,norb
c       if(zact(i)) write (6,'('' Orbital '',I4,'' active'')') i
c       if(.not.zact(i))write (6,'('' Orbital '',I4,'' inactive'')') i
c       call flush(6)
c      enddo

      call deforb(zact,norb,nocb,ncore,nact)
      if(nact.gt.0)then
         allocate(da(1:nact,1:nact))
         allocate(gkp1a(1:nact,1:nact))
         read(31)
         read(31)((da(k,l),k=1,nact),l=1,nact)
         read(31)((gkp1a(k,l),k=1,nact),l=1,nact)
         if(from.eq.'DALTON')gkp1a=gkp1a*factor
      endif
      
c     
c     
c     Termine V^0 (doppie inattive)
c     
c---  calcolo degli eps secondo la formula:
c--   eps(i)=h(i,i)+Sum(j=1,ncore)(2*J(i,j)-K(i,j))+
c--   +Sum(a,b=1,nact)(<ia|ib>-0.5*<ia|bi>)*da(b,a)
c--   analogo per eps(r)
      call flush(6)
      do i=1,ncore
       indi=indice(i,i)
       eps(i)=onel(indi)
       do j=1,ncore
        eps(i)=eps(i)+2.d0*ai(i,i,j,j)-ai(i,j,j,i)
       enddo
       do ia=1,nact
        do ib=1,nact
         eps(i)=eps(i)+(ai(i,i,ia+ncore,ib+ncore)-0.5d0*ai(i,ib
     $        +ncore,ia+ncore,i))*da(ib,ia)
        enddo
       enddo
      enddo
      allocate(fockc(ncore,ncore))
      do i=1,ncore
       fockc(i,i)=eps(i)
       do j=1,i-1
        indi=indice(i,j)
        fockc(i,j)=onel(indi)
        do k=1,ncore
         fockc(i,j)=fockc(i,j)+2.d0*ai(i,j,k,k)-ai(i,k,k,j)
        enddo
        do ia=1,nact
         do ib=1,nact
          fockc(i,j)=fockc(i,j)+(ai(i,j,ia+ncore,ib+ncore)-0.5d0*ai(i,ib
     $         +ncore,ia+ncore,j))*da(ib,ia)
         enddo
        enddo
        fockc(j,i)=fockc(i,j)
       enddo
      enddo
c     print*,'core Fock matrix'
c     call matout(ncore,ncore,fockc,ncore,factor)
      fmax=0.d0
      do i=1,ncore
       do j=1,i-1
        if(abs(fockc(i,j)).gt.fmax) fmax=abs(fockc(i,j))
       enddo
      enddo
      write (6,*)
      write (6,*)
      print*,'Maximum of off-diagonal core Fock matrix is ',fmax/factor
      deallocate(fockc)
      nvirt=norb-ncore-nact
      do ir=ncore+nact+1,norb
       indi=indice(ir,ir)
       eps(ir)=onel(indi)
       do j=1,ncore
        eps(ir)=eps(ir)+2.d0*ai(ir,ir,j,j)-ai(ir,j,j,ir)
       enddo
       do ia=1,nact
        do ib=1,nact
         eps(ir)=eps(ir)+(ai(ir,ir,ia+ncore,ib+ncore)-0.5d0*ai(ir
     $        ,ib+ncore,ia+ncore,ir))*da(ib,ia)
        enddo
       enddo
      enddo
      allocate(fockv(nvirt,nvirt))
      do ir=ncore+nact+1,norb
       irv=ir-ncore-nact
       fockv(irv,irv)=eps(ir)
       do is=ncore+nact+1,ir-1
        isv=is-ncore-nact
        indi=indice(ir,is)
        fockv(irv,isv)=onel(indi)
        do j=1,ncore
         fockv(irv,isv)=fockv(irv,isv)+2.d0*ai(ir,is,j,j)-ai(ir,j,j,is)
        enddo
        do ia=1,nact
         do ib=1,nact
          fockv(irv,isv)=fockv(irv,isv)+(ai(ir,is,ia+ncore,ib+ncore)-
     $         0.5d0*ai(ir,ib+ncore,ia+ncore,is))*da(ib,ia)
         enddo
        enddo
        fockv(isv,irv)=fockv(irv,isv)
       enddo
      enddo
c     print*,'virtual Fock matrix'
c     call matout(nvirt,nvirt,fockv,nvirt,factor)
      fmax=0.d0
      do i=1,nvirt
       do j=1,i-1
        if(abs(fockv(i,j)).gt.fmax) fmax=abs(fockv(i,j))
       enddo
      enddo
      print*,'Maximum of off-diagonal virtual Fock matrix is ',fmax
     $     /factor
      deallocate(fockv)

      write (6,*)
      write (6,*)
      do i=1,ncore
       write (6,'('' Orbital '',I4,'' inactive: energy = '',f15.8)')
     *   i,eps(i)/factor
       fmpb(i)=eps(i)
       fmpb(i+norb)=eps(i)
      enddo

      do i=ncore+1,ncore+nact
      write (6,'('' Orbital '',I4,'' active'')') i
      enddo
      
      do i=ncore+nact+1,norb
       write (6,'('' Orbital '',I4,'' virtual : energy = '',f15.8)')
     *   i,eps(i)/factor
       fmpb(i)=eps(i)
       fmpb(i+norb)=eps(i)
      enddo
c     
c     Definizione degli orbitali gelati
c     
      do i=1,norb
       igelo(i)=0
       zgel(i)=.false.
       zgel(i+norb)=.false.
      enddo
      zgel(norb+norb+1)=.false.
      rewind 3
      call cerca(3,'ICINP',itrovo)
      if (itrovo.eq.0) then
         write (6,*) '   *** NAMELIST &ICINP NOT FOUND ***'
         stop 12
      endif

      call lecnam('&ICINP',namel,3)
      call nameli(namel,igelo,'IGELO',4)

c     
c     Controllo che nessun orbitale gelato (igelo(i).ne.0) sia
c     attivo nei det dello spazio S
c     
      ngel=0
      do i=1,norb
c     write (6,'('' IGELO '',I4,'' = '',I2)') i,igelo(i)
       if (igelo(i).ne.0) then
          ngel=ngel+1
          zgel(igelo(i))=.true.
          zgel(igelo(i)+norb)=.true.
       endif
      enddo
      if (ngel.gt.0) then
         do i=1,nd(ncf)+ne(ncf)
          if (zgel(trou(i))) then
             write (6,*) 'Attenzione l''orbitale ',trou(i),
     *            ' non puo` essere gelato (e` attivo in S)'
             write (6,*) ' buca n. ',i
             stop 12
          endif
          if (zgel(part(i))) then
             write (6,*) 'Attenzione l''orbitale ',part(i),
     *            ' non puo` essere gelato (e` attivo in S)'
             write (6,*) ' particella n. ',i
             stop 12
          endif
         enddo
         write (6,271) (igelo(i),i=1,ngel)
      endif
 271  format(/'     FROZEN ORBITALS:'/10(5x,18i4/))
      write (6,*)
      write (6,*)
      call flush(6)

      write (6,*)
      write (6,*) 'Calling V(0)'
      call flush(6)
      call v0(ncore,nact,zthird)
      print*,'Done with V(0)'
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(1)/factor
      print*,'Norm =',psimpp(1)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2mpp(1)/factor
      print*,'Norm =',psimpp(1)
      write (6,*)
      CALL CPUT(1)
      call flush(6)
      if(nact.eq.0)goto 3200
c     
c     Termine V^(+1) 
      write (6,*) 'Calling V(+1)'
      call flush(6)
      call v1k(ncore,nact,da,gkp1a,zthird)
      print*,'Done with V(+1)'
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(2)/factor
      print*,'Norm =',psimpp(2)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(2)/factor
      print*,'Norm =',psienep(2)
      write (6,*)
      deallocate (gkp1a)
      CALL CPUT(1)
      call flush(6)
c     
c     Termine V^(-1) 
c     
c      write (6,*) 'Allocable working area for V(-1) =', 
c     *     nact**2*8,' bytes'
c      call flush(6)
      allocate(gkm1a(1:nact,1:nact))
      read(31)((gkm1a(k,l),k=1,nact),l=1,nact)
      if(from.eq.'DALTON')gkm1a=gkm1a*factor

      write (6,*) 'Calling V(-1)'
      call flush(6)
      call vm1k(ncore,nact,da,gkm1a,zthird)
      print*,'Done with V(-1)'
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(3)/factor
      print*,'Norm =',psimpp(3)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(3)/factor
      print*,'Norm =',psienep(3)
      write (6,*)
      deallocate (gkm1a)
      CALL CPUT(1)
      call flush(6)
c     
c     Termine V^(+2) 
c     
      if (zprt)
     * write (6,*) 'Allocable working area for V(+2) =', 
     *     2*nact**4*8,' bytes'
      call flush(6)
      allocate(daat(1:nact,1:nact,1:nact,1:nact))
      allocate(gkp2aa(1:nact,1:nact,1:nact,1:nact))
      do i=1,nact
       do j=1,nact
        read(31)((daat(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
c     write (6,*) 'Spinless two-particle Koopmans matrix'
      do i=1,nact
       do j=1,nact
        read(31)((gkp2aa(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
      if(from.eq.'DALTON')gkp2aa=gkp2aa*factor
      write (6,*) 'Calling V(+2)'
      call flush(6)
      call v2mod(ncore,nact,gkp2aa,daat,zthird)
      print*,'Done with V(+2)'
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(4)/factor
      print*,'Norm =',psimpp(4)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(4)/factor
      print*,'Norm =',psienep(4)
      write (6,*)
      deallocate (gkp2aa)
      deallocate (daat)
      CALL CPUT(1)
      call flush(6)
c     
c     Termine V^(-2) 
c     
      if (zprt)
     * write (6,*) 'Allocable working area for V(-2) =', 
     *     2*nact**4*8,' bytes'
      call flush(6)
      allocate(daa(1:nact,1:nact,1:nact,1:nact))
      allocate(gkm2aa(1:nact,1:nact,1:nact,1:nact))
      do i=1,nact
       do j=1,nact
        read(31)((daa(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
      do i=1,nact
       do j=1,nact
        read(31)((gkm2aa(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
      if(from.eq.'DALTON')gkm2aa=gkm2aa*factor
      write (6,*) 'Calling V(-2)'
c      write (6,*)
      call flush(6)
      call vm2mod(ncore,nact,daa,gkm2aa,zthird)
      print*,'Done with V(-2)'
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(5)/factor
      print*,'Norm =',psimpp(5)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(5)/factor
      print*,'Norm =',psienep(5)
      write (6,*)
      deallocate (gkm2aa)
      CALL CPUT(1)
      call flush(6)
c     
c     Termine V^(0p) 
c     
      if (zprt)
     * write (6,*) 'Allocable working area for V(0)'' =', 
     *     2*nact**4*8+nact**2*8,' bytes'
      call flush(6)
      allocate(gk0pa(1:nact,1:nact,1:nact,1:nact))
      allocate(gk0pd(1:nact,1:nact,1:nact,1:nact))
      allocate(gk0pf(1:nact,1:nact))
      do i=1,nact
       do j=1,nact
        read(31)((gk0pa(i,j,k,l),k=1,nact),l=1,nact)
        read(31)((gk0pd(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
      read(31)((gk0pf(k,l),k=1,nact),l=1,nact)
      if(from.eq.'DALTON')gk0pa=gk0pa*factor
      if(from.eq.'DALTON')gk0pd=gk0pd*factor
      if(from.eq.'DALTON')gk0pf=gk0pf*factor
      write (6,*) 'Calling V(0)'' '
      call flush(6)
      call v0pmod(ncore,nact,da,daa,gk0pa,gk0pd,gk0pf,zthird)
      print*,'Done with V(0'')'
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(8)/factor
      print*,'Norm =',psimpp(8)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(8)/factor
      print*,'Norm =',psienep(8)
      write (6,*)
      call flush(6)
      deallocate (gk0pa)
      deallocate (gk0pd)
      deallocate (gk0pf)
      CALL CPUT(1)
      call flush(6)

c     
c     Termine V^(-1p) 
c     


      allocate(ro3(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(amat(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(bmat(1:nact,1:nact,1:nact,1:nact))
      allocate(cmat(1:nact,1:nact,1:nact,1:nact))
      allocate(dmat(1:nact,1:nact))
c-----lettura matrice spinless a 3 part e matrice amat
      do i=1,nact
       do j=1,nact
        do k=1,nact
         do l=1,nact
          read(31)((ro3(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
          read(31)((amat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
         enddo
        enddo
       enddo
      enddo
      if(from.eq.'DALTON')amat=amat*factor
c---- lettura matrici bmat e cmat, btmat e ctmat
      do i=1,nact
       do j=1,nact
        read(31)((bmat(i,j,k,l),k=1,nact),l=1,nact)
        read(31)((cmat(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
      if(from.eq.'DALTON')bmat=bmat*factor
      if(from.eq.'DALTON')cmat=cmat*factor
c---- lettura matrice dmat e dtmat
      read(31)((dmat(i,j),i=1,nact),j=1,nact)
      if(from.eq.'DALTON')dmat=dmat*factor
      write (6,*) 'Calling V(-1)'' '
      call flush(6)
      call vmupE(ncore,nact,ro3,amat,bmat,cmat,dmat,daa,da,zthird)
      print*,'Done with V(-1)'' '
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(7)/factor
      print*,'Norm =',psimpp(7)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(7)/factor
      print*,'Norm =',psienep(7)
      write (6,*)
      CALL CPUT(1)
      call flush(6)

c-----lettura matrice amat per V+1'
      do i=1,nact
       do j=1,nact
        do k=1,nact
         do l=1,nact
          read(31)((amat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
         enddo
        enddo
       enddo
      enddo
      if(from.eq.'DALTON')amat=amat*factor
c---- lettura matrici bmat e cmat per V+1'
      do i=1,nact
       do j=1,nact
        read(31)((bmat(i,j,k,l),k=1,nact),l=1,nact)
        read(31)((cmat(i,j,k,l),k=1,nact),l=1,nact)
       enddo
      enddo
      if(from.eq.'DALTON')bmat=bmat*factor
      if(from.eq.'DALTON')cmat=cmat*factor
c---- lettura matrice dmat  per V+1'
      read(31)((dmat(i,j),i=1,nact),j=1,nact)
      if(from.eq.'DALTON')dmat=dmat*factor
      write (6,*) 'Calling V(+1)'' '
      call flush(6)
      call vpupE(ncore,nact,ro3,amat,bmat,cmat,dmat,daa,da,zthird)
      print*,'Done with V(+1)'' '
      write (6,*)
      print*,'SC-NEVPT'
      print*,'E2   =',e2mpp(6)/factor
      print*,'Norm =',psimpp(6)
      write (6,*)
      print*,'PC-NEVPT'
      print*,'E2   =',e2enp(6)/factor
      print*,'Norm =',psienep(6)
      write (6,*)
      CALL CPUT(1)
      call flush(6)

 3200 continue
      write (6,*) 
      write (6,*) 
      write (6,'(1X,70(''#''))') 
      write (6,*) 
      write (6,*) ' SUMMARY OF THE VARIOUS TERMS'
      write (6,*) 
c      write (6,*) 'Utilizzando H Dyall'
      print '(10x,a)'
     $     ,'   Norm  SC       Energy SC      Norm  PC       Energy PC '
      write (6,123) '(0)  ',psimpp(1),e2mpp(1)/factor,psienep(1),e2enp(1
     $     )/factor
      write (6,123) '(+1) ',psimpp(2),e2mpp(2)/factor,psienep(2),e2enp(2
     $     )/factor
      write (6,123) '(-1) ',psimpp(3),e2mpp(3)/factor,psienep(3),e2enp(3
     $     )/factor
      write (6,123) '(+2) ',psimpp(4),e2mpp(4)/factor,psienep(4),e2enp(4
     $     )/factor
      write (6,123) '(-2) ',psimpp(5),e2mpp(5)/factor,psienep(5),e2enp(5
     $     )/factor
      write (6,123) '(+1)''',psimpp(6),e2mpp(6)/factor,psienep(6)
     $     ,e2enp(6)/factor
      write (6,123) '(-1)''',psimpp(7),e2mpp(7)/factor,psienep(7)
     $     ,e2enp(7)/factor
      write (6,123) '(0)'' ',psimpp(8),e2mpp(8)/factor,psienep(8)
     $     ,e2enp(8)/factor
      print '(2x,65a1)',('-',i=1,65)
      write (6,123) 'Total',psimp(1),e2mp(1)/factor,psien(1),e2en(1)
     $     /factor
c     123 format (2x,A,4f20.15)
 123  format (2x,A,4f15.10)
      CALL CPUT(1)
      WRITE (6,'(/)')
      WRITE (6,2000)
      WRITE (6,1060)
      WRITE (6,2000)
      WRITE (6,'(/)')
      WRITE (6,1105) ESCF
c     if (zenb) WRITE (6,1070) ' ENB'
c     if (.not.zenb) WRITE (6,1070) 'PC-D'
      WRITE (6,1070) 
      DO 420 M=1,METAT
       E(M)=E(M)*FINV
       E2MP(M)=E2MP(M)*FINV 
       E2EN(M)=E2EN(M)*FINV 
 420   WRITE (6,1100) M,E(M),E2MP(M),E2EN(M)
       WRITE (6,'(//)')
       WRITE (6,'(/)')
c     if (zenb) WRITE (6,1120) ' ENB'
c     if (.not.zenb) WRITE (6,1120) 'PC-D'
       WRITE (6,1120) 
       DO 430 M=1,METAT
        E(M)=E(M)+ESCF
        E2MP(M)=E2MP(M)+E(M)
        E2EN(M)=E2EN(M)+E(M)
 430    WRITE (6,1100) M,E(M),E2MP(M),E2EN(M)
        WRITE (6,'(//)')
        WRITE (6,'(/)')
        if(zro)then
           deallocate(ro3)
           deallocate(amat)
           deallocate(bmat)
           deallocate(cmat)
           deallocate(dmat)
           deallocate(daa)
           allocate(rho(1:norb,1:norb))
           call zeroe(rho,norb**2)
           do i=1,nact
            do j=1,nact
             rho(ncore+i,ncore+j)=da(i,j)
            enddo
           enddo
           do i=1,ncore
            rho(i,i)=2.d0
           enddo
           rewind 27
           read(27)((rho(i,j),i=1,ncore),j=ncore+nact+1,norb)
           do i=1,ncore
            do j=ncore+nact+1,norb
             rho(j,i)=rho(i,j)
            enddo
           enddo
           read(27)((rho(i,j),i=ncore+1,ncore+nact),j=ncore+nact+1,norb)
           do i=ncore+1,ncore+nact
            do j=ncore+nact+1,norb
             rho(j,i)=rho(i,j)
            enddo
           enddo
           read(27)((rho(i,j),i=1,ncore),j=ncore+1,ncore+nact)
           do i=1,ncore
            do j=ncore+1,ncore+nact
             rho(j,i)=rho(i,j)
            enddo
           enddo
           if(zprt)then
              print*,'Density matrix 0+1'
              call matout(norb,norb,rho,norb,1.d0)
           endif
           print*,'calculation of natural occupation numbers'
           allocate(w(norb))
           allocate(work(3*norb))
           call dsyev('V','U',norb,rho,norb,w,work,3*norb,info)
           print*,'core occupation numbers:'
           print '(6f12.8)',(w(i),i=norb-ncore+1,norb)
           print*,'active occupation numbers:'
           print '(6f12.8)',(w(i),i=norb-ncore-nact+1,norb-ncore)
           print*,'virtual occupation numbers:'
           print '(6f12.8)',(w(i),i=1,norb-ncore-nact)
        endif
        CALL CPUT(1)
 1040   FORMAT (i4,f16.5,9x,2f11.6,2x,2f11.6)
 1060   FORMAT (5X,'FINAL RESULTS')
 1070   FORMAT (' STATE',2X,'DIAGONALIZATION        PERTURB. SC-NEVPT',
     *       '        PERTURB. PC-NEVPT')
 1100   FORMAT (I4,F17.8,3F24.8)
 1105   FORMAT (' REFERENCE: ',F16.8,//)
 1120   FORMAT (' STATE',2X,'TOTAL DIAGON.          TOTAL SC-NEVPT',
     *       '          TOTAL PC-NEVPT')
 1225   FORMAT (//5X,'ONT ETE GENERES:',I15,' DETERMINANTS EFFICACES'/
     *       5x,'NEC =',14i12/5x,'DET =',14i12//)
 1400   FORMAT (5X,'NORME DE LA CORRECTION D''ORDRE 1 A LA FONCTION ',
     *       'D''ONDE'/5x,'ETAT            TOTAL   NEC =',i4,14i10)
 1401   format (i8,2x,a4,f12.6,5x,14f10.6)
 2000   FORMAT (5X,13('*'))
        stop
        END
C***********************************************************************
      subroutine read2(a,nact)
      real*8 a(nact,nact)
      read(31)((a(k,l),k=1,nact),l=1,nact)
      return
      end
      subroutine read4(a,nact)
      real*8 a(nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            read(31)((a(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
      return
      end
      subroutine read44(a,b,nact)
      real*8 a(nact,nact,nact,nact)
      real*8 b(nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            read(31)((a(i,j,k,l),k=1,nact),l=1,nact)
            read(31)((b(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
      return
      end
      subroutine read6(a,nact)
      real*8 a(nact,nact,nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  read(31)((a(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine read66(a,b,nact)
      real*8 a(nact,nact,nact,nact,nact,nact)
      real*8 b(nact,nact,nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  read(31)((a(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
                  read(31)((b(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
               enddo
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE OFILES(zbin,zverbose,zro)
      implicit character*64(f)
      logical*1 zbin,zverbose
      character*1280 namel
      character*80 card
      common /fil/ file12,file13
      file03=' '
      file04='CIPPI.FILE04'
      file07=' '
      file09='CIPPI.FILE09'
      file11=' '
      file12=' '
      file13='DIAGO.RESTART'
      file20='CIPPI.FILE20'
      file25='IJKL.FILE25'
      file26=' '
      file31='KOOP.FILE31'
      file50='IJKL.FILE50'
      zbin=.true.
      call lecnam('&FILES',namel,5)
      call namela(namel,file03,'FILE03')
      call namela(namel,file04,'FILE04')
      call namela(namel,file07,'FILE07')
      call namela(namel,file09,'FILE09')
      call namela(namel,file11,'FILE11')
      call namela(namel,file12,'FILE12')
      call namela(namel,file13,'FILE13')
      call namela(namel,file20,'FILE20')
      call namela(namel,file25,'FILE25')
      call namela(namel,file26,'FILE26')
      call namela(namel,file31,'FILE31')
      call namela(namel,file32,'FILE32')
      call namela(namel,file50,'FILE50')
      call namell(namel,zverbose,'ZVERBOSE',1)
      call namell(namel,zro,'ZRO',1)
c      call namell(namel,zbin,'ZBIN',1)
      if (file03.eq.' ') then
      file03=file04
      i=index(file03,'file04')
      if (i.ne.0) file03(i:i+5)='file03'
      i=index(file03,'FILE04')
      if (i.ne.0) file03(i:i+5)='FILE03'
      i=index(file03,'f04')
      if (i.ne.0) file03(i:i+2)='f03'
      i=index(file03,'F04')
      if (i.ne.0) file03(i:i+2)='F03'
      i=index(file03,'file4')
      if (i.ne.0) file03(i:i+4)='file3'
      i=index(file03,'FILE4')
      if (i.ne.0) file03(i:i+4)='FILE3'
      i=index(file03,'f4')
      if (i.ne.0) file03(i:i+1)='f3'
      i=index(file03,'F4')
      if (i.ne.0) file03(i:i+1)='F3'
      endif
      if (file03.eq.file04) file03='CIPPI.FILE03'
      open(3,file=file03,form='FORMATTED',status='UNKNOWN')
      rewind 3
      open(4,file=file04,form='UNFORMATTED',status='UNKNOWN')
      rewind 4

      open(9,file=file09,form='UNFORMATTED',status='UNKNOWN')
      rewind 9
c  apro il file 20 con una subroutine da compilare a parte per la
c  compatibilita` ibm-dec
      open(25,file=file25,form='UNFORMATTED',status='OLD')
      rewind 25
      if (file26.ne.' ') then
      open(26,file=file26,form='UNFORMATTED',status='OLD')
      rewind 26
      endif
      open(31,file=file31,form='UNFORMATTED',status='UNKNOWN')
      rewind 31
      open(32,file=file32,form='UNFORMATTED',status='UNKNOWN')
      rewind 32
      open(50,file=file50,form='UNFORMATTED',status='OLD')
      maxlong=200000000
      call filesplit('OPEN',50,maxlong,idum,idum2)
      rewind 50
      call filesplit('REWIND',50,0,idum,idum2)
      WRITE (6,'(/''    INPUT DATA:''/)')
      write (6,*) ' FILE04 = ',file04
      write (6,*) ' FILE07 = ',file07
      write (6,*) ' FILE11 = ',file11
      write (6,*) ' FILE25 = ',file25
      write (6,*) ' FILE50 = ',file50
      icount=0
   50 read (5,'(a80)',end=99) card
      write (3,'(a80)') card
      icount=icount+1
      if (icount.le.15) write (6,'(1x,a80)') card
      i=index(card,'ICINP')
      if (i.ne.0) icount=0
      goto 50
   99 if (icount.gt.15) WRITE (6,'(8('' ..... etc''))')
      RETURN
      END
C***********************************************************************
      SUBROUTINE MATOUT(N,N2,A,NN,factor)
      real*8 A,factor
c     CHARACTER*3 SYMB(NN)
c     integer symb(*)
      DIMENSION A(NN,*)
C 
   20 FORMAT(//,8X,11(5X,I3,3X))
   30 FORMAT(1X,I4,4X,11(F11.6))
   35 FORMAT(8X,11(5X,I3,3X))
   37 FORMAT(8X,11(4X,A4,3X))
   36 FORMAT(/)
C 
c      DO 40 M=1,N2,11
      DO 40 M=1,N2,6
c      K=M+10
      K=M+5
      IF(K.LE.N2) GOTO 10
      K=N2
   10 WRITE(6,20) (J,J=M,K)
      WRITE(6,36)
c     WRITE(6,35) (ME(I),I=M,K)
c     WRITE(6,37) (SYMB(I),I=M,K)
      WRITE(6,36)
      DO 40 I=1,N
      WRITE(6,30) I,(A(I,J)/factor,J=M,K)
 40   CONTINUE
      call flush(6)
      RETURN
      END
C***********************************************************************
      subroutine v0(ncore,nact,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      integer rs
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
      allocatable c1(:,:),ene(:,:),sqnor(:,:)
c      integer*2 icomp,iconf,isegno
c
c
c     Attenzione il metodo e` state specific!!
c     If Metat.ne.1 STOP!!!!!
c
      if (metat.ne.1) then
      write (6,*) 'Attenzione il metodo NEV-PT SC e` state specific'
      write (6,*) 'Metat deve essere uguale a 1'
      stop
      endif
      nvirt=norb-ncore-nact
      
      if (zthird) then
      allocate(c1(ncore*(ncore+1)/2,nvirt*(nvirt+1)/2))
      allocate(ene(ncore*(ncore+1)/2,nvirt*(nvirt+1)/2))
      allocate(sqnor(ncore*(ncore+1)/2,nvirt*(nvirt+1)/2))
      call zeroe(c1,(ncore*(ncore+1)/2)*(nvirt*(nvirt+1)/2))
      call zeroe(ene,(ncore*(ncore+1)/2)*(nvirt*(nvirt+1)/2))
      call zeroe(sqnor,(ncore*(ncore+1)/2)*(nvirt*(nvirt+1)/2))
      endif

      ij2=0
      do 800 ii=1,ncore
      if (zgel(ii)) goto 800
      isy=itsym(ii)
      do 700 ij=ii,ncore
      if (zgel(ij)) goto 700
      jsy=itsym(ij)
      ijsy=its(isy,jsy)
      ij2=ij2+1
      rs=0
      do 600 ir=ncore+nact+1,norb
      irsy=itsym(ir)
      irijsy=its(ijsy,irsy)

      do 500 is=ir,norb
      issy=itsym(is)
      irsijsy=its(irijsy,issy)
      if(irsijsy.ne.1)goto 500
      rs=rs+1

      zrs=ir.eq.is
      zij=ii.eq.ij
      ab1=ai(ir,ij,is,ii)
      ab2=ai(ir,ii,is,ij)
      if (zrs.and.zij) dint=ab1 !always positive
      if (zrs.neqv.zij) dint=sqrt(2.d0)*abs(ab1)  !I put abs, could be <0
      if (.not.(zrs.or.zij)) dint=2.d0*sqrt(ab1**2+ab2**2-ab1*ab2)
      deno=(fmpb(ir)+fmpb(is)-fmpb(ii)-fmpb(ij))
      co=-dint/deno
      e2mp(1)=e2mp(1)+dint*co
      psimp(1)=psimp(1)+co*co
      e2mpp(1)=e2mpp(1)+dint*co
      psimpp(1)=psimpp(1)+co*co

      if (zthird) then
      c1(ij2,rs)=-abs(dint)/deno
      ene(ij2,rs)=deno
c--dint*finv is a measure of linear dependencies
c      if(dint.ne.0.d0)sqnor(ij2,rs)=1.d0/abs(dint)
      if(dint*finv.gt.1.d-7)sqnor(ij2,rs)=1.d0/abs(dint) !renzo modif
      endif

 500  enddo
  600 enddo

  700 enddo

 800  enddo
      e2enp(1)=e2mpp(1)
      psienep(1)=psimpp(1)
      e2en(1)=e2mp(1)
      psien(1)=psimp(1)

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

      return
      end
c----------------------------------------------------------------
      subroutine v1k(ncore,nact,da,gkp1a,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      real*4 tarray(2)
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension da(nact,nact),gkp1a(nact,nact)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:)
      allocatable c1(:,:,:),ene(:,:,:),sqnor(:,:,:)
c
      nvirt=norb-ncore-nact

      if (zthird) then
      allocate(c1(ncore,ncore,nvirt))
      allocate(ene(ncore,ncore,nvirt))
      allocate(sqnor(ncore,ncore,nvirt))
      call zeroe(c1,ncore*ncore*nvirt)
      call zeroe(ene,ncore*ncore*nvirt)
      call zeroe(sqnor,ncore*ncore*nvirt)
      endif

      do 800 ii=1,ncore
      if (zgel(ii)) goto 800
      isy=itsym(ii)

      do 700 ij=ii,ncore
      if (zgel(ij)) goto 700
      jsy=itsym(ij)
      ijsy=its(isy,jsy)

      do 600 is=ncore+nact+1,norb
      issy=itsym(is)
      isijsy=its(ijsy,issy)

c      write (6,*) 'Eccitazione con gli indici',ii,ij,' -> ',is

      dnormk=0.d0
      epsimk=0.d0
      do ia=1,nact
       iaa=ia+ncore
       do iap=1,nact
        iaap=iap+ncore
	contj=ai(iaa,ii,is,ij)
	contjp=ai(iaap,ii,is,ij)
	if (ij.eq.ii) then
	cint=0.5d0*contj*contjp
	else
	contk=ai(iaa,ij,is,ii)
	contkp=ai(iaap,ij,is,ii)
	cint=contj*contjp+contk*contkp-
     *       0.5d0*contj*contkp-0.5d0*contk*contjp
	endif
	if (ia.eq.iap) then
	ro=2.d0-da(ia,iap)
	else
	ro=-da(ia,iap)
	endif
	dnormk=dnormk+2.d0*ro*cint
	epsimk=epsimk+gkp1a(ia,iap)*cint
      enddo
      enddo

      if (abs(dnormk)*finv**2.lt.1.d-14) goto 600
      epsimk=epsimk/dnormk
      den=fmpb(ii)+fmpb(ij)-fmpb(is)-epsimk
c      write (6,'(A,4F15.10)') 'dn,eps,den',dnormk/factor**2,
c     *            epsimk/factor,den/factor
c      write (6,*) 'E2',e2mpp(2)/factor
      co=dnormk/den
      con=co/den
      e2mp(1)=e2mp(1)+co
      psimp(1)=psimp(1)+con
      e2mpp(2)=e2mpp(2)+co
      psimpp(2)=psimpp(2)+con
      ismeno=is-ncore-nact

      if (zthird) then
      sqnor(ii,ij,ismeno)=sqrt(dnormk)
      c1(ii,ij,ismeno)=sqnor(ii,ij,ismeno)/den
      ene(ii,ij,ismeno)=-den
      sqnor(ii,ij,ismeno)=1.d0/sqnor(ii,ij,ismeno)
      endif

 600  enddo

  700 enddo

  800 enddo

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c----Partially contracted NEV-PT---------------------------------
      tini=etime(tarray)
      allocate(s(1:nact,1:nact))
      allocate(coef(1:nact,1:nact))
      allocate(w(nact))
      allocate(work(3*nact))
      call fill2(s,da,nact)
      call fill(coef,gkp1a,nact**2)
      call halve(coef,nact**2)
c      call dsygv(1,'V','U',nact,coef,nact,s,nact,w,work,3*nact,info)
      call dsyev('V','U',nact,s,nact,w,work,3*nact,info)
c      if(info.ne.0)then
c         print*,'something wrong in diagonalization: sub v1k'
c         stop 'v1k'
c      endif
      thr=1.d-8                 !perhaps too strict
      ndep=0
      do i=1,nact
       if(w(i).lt.0.d0)then
        ndep=ndep+1
        if(zprt) write (6,'(a,i3,a,f12.8)')
     $    'negative eigenvalue: eps(',i,')=',w(i)
        elseif(w(i).lt.thr)then
            ndep=ndep+1
        endif
      enddo
      if (ndep.gt.0.and.zprt) then
      write (6,*)'Linear dependencies analysis:'
      write (6,'(a,i4,a)')' There are ',ndep,' linear dependencies'
      endif
      nnew=nact-ndep
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact,1:nnew))
      call caux0p(caux,s,w,nact,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact,nnew,1.d0,coef,nact,caux,nact,0.d0,s
     $     ,nact)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact,1.d0,caux,nact,s,nact,0.d0
     $     ,coef,nact)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact,w,work,3*nact,info)
c--   back transformation
c--   moltiplico CAUX(nact,nnew) per coef(nnew,nnew)
c--   e metto in S(nact,nnew)
      call dgemm('N','N',nact,nnew,nnew,1.d0,caux,nact,coef,nact,0
     $     .d0,s,nact)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact,nnew)
      deallocate(caux)
      ioffw=1
      if(zprt)then
         print*,'eigenvalues are:'
         do i=1,nact
            print '(a,i2,a,f12.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
      endif
      call e2p1cont(e2,psi2,coef,s,w,da,fmpb,nact,ncore,norb,nnew,zgel)
c      print*,'Partially contracted contribution: class (+1):'
c      print '(a,f20.14)','E2=',e2/factor
c      print '(a,f20.14)','psi2=',psi2
      e2enp(2)=e2
      psienep(2)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
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
      end
c-----------------------------------------------------------
      subroutine vm1k(ncore,nact,da,gkm1a,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension da(nact,nact),gkm1a(nact,nact)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:)
      allocatable c1(:,:,:),ene(:,:,:),sqnor(:,:,:)
c
      nvirt=norb-ncore-nact

      if (zthird) then
      allocate(c1(nvirt,nvirt,ncore))
      allocate(ene(nvirt,nvirt,ncore))
      allocate(sqnor(nvirt,nvirt,ncore))
      call zeroe(c1,nvirt*nvirt*ncore)
      call zeroe(ene,nvirt*nvirt*ncore)
      call zeroe(sqnor,nvirt*nvirt*ncore)
      endif

      do 800 ii=1,ncore
      if (zgel(ii)) goto 800
      isy=itsym(ii)

      do 700 ir=ncore+nact+1,norb 
      irsy=itsym(ir)
      irisy=its(isy,irsy)

      do 600 is=ir,norb
      issy=itsym(is)
      irsisy=its(irisy,issy)

      dnormk=0.d0
      epsimk=0.d0
      do ia=1,nact
       iaa=ia+ncore
       do iap=1,nact
        iaap=iap+ncore
	contj=ai(ir,ii,is,iaa)
	contjp=ai(ir,ii,is,iaap)
	if (ir.eq.is) then
	cint=.5d0*contj*contjp
	else
	contk=ai(ir,iaa,is,ii)
	contkp=ai(ir,iaap,is,ii)
        cint=contj*contjp+contk*contkp-
     *       0.5d0*contj*contkp-0.5d0*contk*contjp
	endif
	dnormk=dnormk+2.d0*da(ia,iap)*cint
	epsimk=epsimk+gkm1a(ia,iap)*cint
       enddo
      enddo

      if (abs(dnormk)*finv**2.lt.1.d-14) goto 600

      epsimk=-epsimk/dnormk
      den=epsimk+fmpb(ii)-fmpb(ir)-fmpb(is)
      co=dnormk/den
      con=co/den
      irmeno=ir-ncore-nact
      ismeno=is-ncore-nact

      if (zthird) then
      sqnor(irmeno,ismeno,ii)=sqrt(dnormk)
      c1(irmeno,ismeno,ii)=sqnor(irmeno,ismeno,ii)/den
      ene(irmeno,ismeno,ii)=-den
      sqnor(irmeno,ismeno,ii)=1.d0/sqnor(irmeno,ismeno,ii)
      endif

c--renzo debug
      if(abs(con).gt.0.1d0)then
         print '(a,i3,a,i3,a,i3,a,3f12.6)','abnormal: eps(',ii,'), eps('
     $        ,ir,') eps(',is,')=',fmpb(ii)/factor,fmpb(ir)/factor
     $        ,fmpb(is)/factor
         print '(a,f12.6)','epsimk=',epsimk/factor
         print '(a,f12.6,a,f12.6)','numeratore=',dnormk/factor
     $        ,' denom.=',den/factor
      endif
c----------
      e2mp(1)=e2mp(1)+co
      psimp(1)=psimp(1)+con
      e2mpp(3)=e2mpp(3)+co
      psimpp(3)=psimpp(3)+con



  600 enddo

  700 enddo

  800 enddo
      
      if(zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c----Partially contracted NEV-PT---------------------------------
      tini=etime(tarray)
      allocate(s(1:nact,1:nact))
      allocate(coef(1:nact,1:nact))
      allocate(w(nact))
      allocate(work(3*nact))
      call fill(s,da,nact**2)
      call fill(coef,gkm1a,nact**2)
      call halve(coef,nact**2)
      call dsyev('V','U',nact,s,nact,w,work,3*nact,info)
c      call dsygv(1,'V','U',nact,coef,nact,s,nact,w,work,3*nact,info)
c      if(info.ne.0)then
c         print*,'something wrong in diagonalization: sub vm1k'
c         stop 'vm1k'
c      endif
      ioffw=1
      thr=1.d-8                 !perhaps too strict
      ndep=0
      do i=1,nact
       if(w(ioffw+i-1).lt.0.d0)then
        ndep=ndep+1
        if (zprt)
     $  write (6,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
        elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
        endif
      enddo
      if (ndep.gt.0.and.zprt) then
      write (6,*)'Linear dependencies analysis:'
      write (6,'(a,i4,a)')' There are ',ndep,' linear dependencies'
      endif
      nnew=nact-ndep
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact,1:nnew))
      call caux0p(caux,s,w,nact,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact,nnew,1.d0,coef,nact,caux,nact,0.d0,s
     $     ,nact)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact,1.d0,caux,nact,s,nact,0.d0
     $     ,coef,nact)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact,w,work,3*nact,info)
c--   back transformation
c--   moltiplico CAUX(nact,nnew) per coef(nnew,nnew)
c--   e metto in S(nact,nnew)
      call dgemm('N','N',nact,nnew,nnew,1.d0,caux,nact,coef,nact,0
     $     .d0,s,nact)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact,nnew)
      deallocate(caux)
      if(zprt)then
         print*,'eigenvalues are:'
         do i=1,nact
            print '(a,i2,a,f12.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
      endif
      call e2m1cont(e2,psi2,coef,s,w,da,fmpb,nact,ncore,norb,nnew,zgel)
c      print*,'Partially contracted contribution: class (-1):'
c      print '(a,f20.14)','E2=',e2/factor
c      print '(a,f20.14)','psi2=',psi2
      e2enp(3)=e2
      psienep(3)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
c-------------------------------------------------------------
      subroutine fill(a,b,n)
      real*8 a(*),b(*)
      do i=1,n
         a(i)=b(i)
      enddo
      return
      end
c----------------------------------------------
      subroutine halve(a,n)
      real*8 a(*)
      do i=1,n
         a(i)=0.5d0*a(i)
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2m1cont(e2,psi2,coef,s,w,da,eps,nact,ncore,norb,nnew,
     $                    zgel)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(*),eps(*),coef(nact,nact),s(nact,nact),da(nact,nact)
      dimension zgel(*)
      integer a,ap
c---costruzione matrice S
      do a=1,nact
         do mu=1,nnew
            s(a,mu)=0.d0
            do ap=1,nact
               s(a,mu)=s(a,mu)+da(ap,a)*coef(ap,mu)
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do i=1,ncore
        if (zgel(i)) goto 1
         do ir=ncore+nact+1,norb
            do is=ir,norb
               do mu=1,nnew
                  rsmui=0.d0
                  rsmuip=0.d0
                  do a=1,nact
                     if(ir.ne.is)then
                        aint=ai(ir,i,is,a+ncore)+ai(is,i,ir,a+ncore)
                        bint=ai(ir,i,is,a+ncore)-ai(is,i,ir,a+ncore)
                        rsmui=rsmui+aint*s(a,mu)
                        rsmuip=rsmuip+bint*s(a,mu)
                     else
                        aint=ai(ir,i,ir,a+ncore)
                        rsmui=rsmui+aint*s(a,mu)
                     endif
                  enddo
                  deno=w(mu)+eps(ir)+eps(is)-eps(i)
                  if(ir.ne.is)then
                     contri=(0.5d0*rsmui**2+1.5d0*rsmuip**2)/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  else
                     contri=rsmui**2/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  endif
               enddo
            enddo
         enddo
 1    enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2m2cont(e2,psi2,coef,s,w,daa,eps,nact,ncore,norb,nnew
     $     ,zgel)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(*),eps(*),coef(nact**2,nact**2),s(nact**2,nact**2)
     $     ,daa(nact,nact,nact,nact)
      dimension zgel(*)
      integer a,ap,b,bp,ab,abp
c--caso r#s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=1,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+daa(ap,bp,a,b)*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb
         do is=ir+1,norb
            do mu=1,nnew
               rsmu=0.d0
               ab=0
               do a=1,nact
                  do b=1,nact
                     ab=ab+1
                     rsmu=rsmu+ai(ir,b+ncore,is,a+ncore)*s(ab,mu)
                  enddo
               enddo
               deno=w(mu)+eps(ir)+eps(is)
               contri=rsmu**2/deno
               e2=e2-contri
               psi2=psi2+contri/deno
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2m2contp(e2,psi2,coef,s,w,daa,eps,nact,ncore,norb,nnew
     $     ,zgel)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(*),eps(*),coef(nact*(nact+1)/2,nact*(nact+1)/2),s(nact
     $     *(nact+1)/2,nact*(nact+1)/2),daa(nact,nact,nact,nact)
      dimension zgel(*)
      integer a,ap,b,bp,ab,abp
c--caso r=s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=a,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=ap,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+0.5d0*(daa(ap,bp,a,b)+daa(ap
     $                    ,bp,b,a))*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb
         do mu=1,nnew
            rsmu=0.d0
            ab=0
            do a=1,nact
               ab=ab+1
                  rsmu=rsmu+ai(ir,a+ncore,ir,a+ncore)*s(ab,mu)
               do b=a+1,nact
                  ab=ab+1
                  rsmu=rsmu+ai(ir,b+ncore,ir,a+ncore)*s(ab,mu)*2.d0
               enddo
            enddo
            deno=w(mu)+eps(ir)+eps(ir)
            contri=rsmu**2/deno
            e2=e2-contri
            psi2=psi2+contri/deno
         enddo
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2p2cont(e2,psi2,coef,s,w,daa,eps,nact,ncore,norb,nnew
     $     ,zgel)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(*),eps(*),coef(nact**2,nact**2),s(nact**2,nact**2)
     $     ,daa(nact,nact,nact,nact),zgel(*)
      integer a,ap,b,bp,ab,abp
c--caso r#s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=1,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+daa(ap,bp,a,b)*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do i=1,ncore
        if (zgel(i)) goto 1
         do j=i+1,ncore
           if (zgel(j)) goto 2
            do mu=1,nnew
               aijmu=0.d0
               ab=0
               do a=1,nact
                  do b=1,nact
                     ab=ab+1
                     aijmu=aijmu+ai(i,b+ncore,j,a+ncore)*s(ab,mu)
                  enddo
               enddo
               deno=w(mu)-eps(i)-eps(j)
               contri=aijmu**2/deno
               e2=e2-contri
               psi2=psi2+contri/deno
            enddo
 2       enddo
 1    enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2p2contp(e2,psi2,coef,s,w,daa,eps,nact,ncore,norb,nnew
     $     ,zgel)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(*),eps(*),coef(nact*(nact+1)/2,nact*(nact+1)/2),s(nact
     $     *(nact+1)/2,nact*(nact+1)/2),daa(nact,nact,nact,nact),zgel(*)
      integer a,ap,b,bp,ab,abp
c--caso r=s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=a,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=ap,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+0.5d0*(daa(ap,bp,a,b)+daa(ap
     $                    ,bp,b,a))*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do i=1,ncore
        if (zgel(i)) goto 1
         do mu=1,nnew
            aijmu=0.d0
            ab=0
            do a=1,nact
               ab=ab+1
               aijmu=aijmu+ai(i,a+ncore,i,a+ncore)*s(ab,mu)
               do b=a+1,nact
                  ab=ab+1
                  aijmu=aijmu+ai(i,b+ncore,i,a+ncore)*s(ab,mu)*2.d0
               enddo
            enddo
            deno=w(mu)-eps(i)-eps(i)
            contri=aijmu**2/deno
            e2=e2-contri
            psi2=psi2+contri/deno
         enddo
 1    enddo
      return
      end
c----------------------------------------------------------------
      subroutine e2p1cont(e2,psi2,coef,s,w,da,eps,nact,ncore,norb,nnew,
     $                    zgel)
      implicit real*8(a-h,o-y),logical*1(z)
      dimension w(*),eps(*),coef(nact,nact),s(nact,nact),da(nact,nact)
     $     ,zgel(*)
      integer a,ap
      logical*1 zcond
c---costruzione matrice S
      do a=1,nact
         do mu=1,nnew
            s(a,mu)=0.d0
            do ap=1,nact
               if(ap.ne.a)then
                  ro=-da(ap,a)
               else
                  ro=2.d0-da(ap,a)
               endif
               s(a,mu)=s(a,mu)+ro*coef(ap,mu)
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb
         do i=1,ncore
	  if (zgel(i)) goto 2
            do j=i,ncore
	      if (zgel(j)) goto 3
               do mu=1,nnew
                  rjimu=0.d0
                  rjimup=0.d0
                  do a=1,nact
                     if(i.ne.j)then
                        aint=ai(ir,j,a+ncore,i)+ai(ir,i,a+ncore,j)
                        bint=ai(ir,j,a+ncore,i)-ai(ir,i,a+ncore,j)
                        rjimu=rjimu+aint*s(a,mu)
                        rjimup=rjimup+bint*s(a,mu)
                     else
                        aint=ai(ir,i,a+ncore,i)
                        rjimu=rjimu+aint*s(a,mu)
                     endif
                  enddo
                  deno=w(mu)+eps(ir)-eps(i)-eps(j)
                  if(i.ne.j)then
                     contri=(0.5d0*rjimu**2+1.5d0*rjimup**2)/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  else
                     contri=rjimu**2/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  endif
               enddo
 3          enddo
 2       enddo
      enddo
      return
      end
c------------------------------------------------------
      subroutine v2mod(ncore,nact,gkp2aa,daat,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension gkp2aa(nact,nact,nact,nact)
      dimension daat(nact,nact,nact,nact)
      allocatable coef(:,:),s(:,:),w(:),work(:),caux(:,:)
      allocatable c1(:,:),ene(:,:),sqnor(:,:)

      if (zthird) then
      allocate(c1(ncore,ncore))
      allocate(ene(ncore,ncore))
      allocate(sqnor(ncore,ncore))
      call zeroe(c1,ncore*ncore)
      call zeroe(ene,ncore*ncore)
      call zeroe(sqnor,ncore*ncore)
      endif

c
c     t ne u
c

      do 1800 ir=1,ncore
      if (zgel(ir)) goto 1800
      irsy=itsym(ir)

      do 1700 is=ir,ncore
      if (zgel(is)) goto 1700
      issy=itsym(is)
      irssy=its(irsy,issy)
c      write (6,*) 'Eccitazione con gli indici',ir,is,' -> '

c
c     Si inizia il ciclo sugli attivi 
c
      
      dn=0.d0
      dk=0.d0

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
                  dum=daat(icpa,idpa,ica,ida)
                  dumk=gkp2aa(icpa,idpa,ica,ida)
                  dint1=ai(is,icp,ir,idp)
                  dint2=ai(is,ic,ir,id)
		  cint=dint1*dint2
		  if (ir.eq.is) then
		  cint=cint*0.25d0
		  dum=dum+daat(icpa,idpa,ida,ica)
		  dumk=dumk+gkp2aa(idpa,icpa,ica,ida)
		  endif
                  dn=dn+cint*dum
                  dk=dk+cint*dumk
 1120          enddo
            enddo
 1121    enddo
      enddo

      if(abs(dn)*finv**2.gt.1.d-14)then
         epsimk=dk/dn                     
         den=fmpb(ir)+fmpb(is)-epsimk

         if (zthird) then
         sqnor(ir,is)=sqrt(dn)
         ene(ir,is)=-den
         c1(ir,is)=sqnor(ir,is)/den
         sqnor(ir,is)=1.d0/sqnor(ir,is)
         endif

         co=dn/den
         con=co/den
         e2mp(1)=e2mp(1)+co
         psimp(1)=psimp(1)+con
         e2mpp(4)=e2mpp(4)+co
         psimpp(4)=psimpp(4)+con
c      write (6,'(A,4F15.10)') 'dn,eps,den',dn/factor**2,
c     *            epsimk/factor,den/factor
c      write (6,*) 'E2',e2mpp(4)/factor
      endif
      
 1700 enddo
 1800 enddo

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c---Partially contracted NEV-PT-----------------------------------
      tini=etime(tarray)
      allocate(s(1:nact**2,1:nact**2))
      allocate(coef(1:nact**2,1:nact**2))
      allocate(w(nact**2))
      allocate(work(3*nact**2))
      call fillcm2(s,daat,nact)
      call fillcm2(coef,gkp2aa,nact)
      ioffw=1
      if(zprt)then
         call provasymm(daat,nact,discrep,discrep2)
         print '(a,f12.8,a,f12.8)','max. err. in herm. of daat=',discrep
     $        ,' in perm.=',discrep2
         call provasymm(gkp2aa,nact,discrep,discrep2)
         print '(a,f12.8,a,f12.8)','max. err. in herm. of K=',discrep
     $        /factor,' in perm.=',discrep2/factor
      endif
c      call matout(nact**2,nact**2,s,nact**2,1.d0)
c      call dsygv(1,'V','U',nact**2,coef,nact**2,s,nact**2,w,work,3*nact
c     $     **2,info)
      nnew=nact**2
      nact3=nact**2
c      info=1
c      if(info.ne.0)then
c         print*,'something wrong in diagonalization: sub v2mod'
c         stop 'v2mod'
c         print*,'trying to eliminate linear dependencies...'
c         print*,'Mind: this is not a normal situation'
c         call fillcm2(s,daat,nact)
c         call fillcm2(coef,gkp2aa,nact)
      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i+ioffw-1)
         enddo
      endif
      if (zprt) print*,'Analisi delle dipendenze lineari:'
      thr=1.d-8                 !perhaps too strict
      ndep=0
      do i=1,nact3
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
      if(nnew.eq.0)then  !nothing to do
         e2=0.d0
         psi2=0.d0
         goto 1
      endif
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
         call flush(6)
      endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef,nact3,0
     $     .d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
      deallocate(caux)
c      endif
      ioffw=1
c      if(zverbose)then
c         print*,'info is ',info
c         print*,'caso i.ne.j'
c         print*,'eigenvalues are:'
c         do i=1,nnew
c            print '(a,i2,a,f12.8)','eps(',i,')=',w(i+ioffw-1)/factor
c         enddo
c      endif
      call e2p2cont(e2,psi2,coef,s,w,daat,fmpb,nact,ncore,norb,nnew,zgel
     $     )
      if(zprt)then
         print*,'Partially contracted contribution: class (-2): i.ne.j'
         print '(a,f20.14)','E2=',e2/factor
         print '(a,f20.14)','psi2=',psi2
      endif
 1    continue
      e2enp(4)=e2
      psienep(4)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      nact2=nact*(nact+1)/2
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(1:nact2,1:nact2))
      allocate(w(1:nact2))
      allocate(work(1:3*nact2))
      call fillcm2p(s,daat,nact)
      call fillcm2p(coef,gkp2aa,nact)
c      call dsygv(1,'V','U',nact2,coef,nact2,s,nact2,w,work,3*nact2
c     $     ,info)
      nact3=nact2
      nnew=nact2
c      info=1
c      if(info.ne.0)then
c         print*,'something wrong in diagonalization: sub v2mod'
c         stop 'v2mod'
c         print*,'trying to eliminate linear dependencies...'
c         print*,'Mind: this is not a normal situation'
c      call fillcm2p(s,daat,nact)
c      call fillcm2p(coef,gkp2aa,nact)
      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i+ioffw-1)
         enddo
      endif
      if (zprt) print*,'Analisi delle dipendenze lineari:'
      thr=1.d-8                 !perhaps too strict
      ndep=0
      do i=1,nact3
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
      if(nnew.eq.0)then  !nothing to do
         e2=0.d0
         psi2=0.d0
         goto 2
      endif
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
         call flush(6)
      endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef,nact3,0
     $     .d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
      deallocate(caux)
c      endif
c      if(zverbose)then
c         print*,'info is ',info
c         print*,'caso i.eq.j'
c         print*,'eigenvalues are:'
c         do i=1,nnew
c            print '(a,i2,a,f12.8)','eps(',i,')=',w(i+ioffw-1)/factor
c         enddo
c      endif
      call e2p2contp(e2,psi2,coef,s,w,daat,fmpb,nact,ncore,norb,nnew
     $     ,zgel)
      if(zprt)then
         print*,'Partially contracted contribution: class (-1): i.eq.j'
         print '(a,f20.14)','E2=',e2/factor
         print '(a,f20.14)','psi2=',psi2
      endif
 2    continue
      e2enp(4)=e2enp(4)+e2
      psienep(4)=psi2+psienep(4)
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
c      print*,'Partially contracted contribution: class (+2): '
c      print '(a,f20.14)','E2=',e2enp(4)/factor
c      print '(a,f20.14)','psi2=',psienep(4)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
c----------------------------------------------------------------
      subroutine vm2mod(ncore,nact,daa,gkm2aa,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension daa(nact,nact,nact,nact)
      dimension gkm2aa(nact,nact,nact,nact)
      allocatable coef(:,:),s(:,:),w(:),work(:),caux(:,:)
      allocatable c1(:,:),ene(:,:),sqnor(:,:)
c
c

      nvirt=norb-ncore-nact

      if (zthird) then
      allocate(c1(nvirt,nvirt))
      allocate(ene(nvirt,nvirt))
      allocate(sqnor(nvirt,nvirt))
      call zeroe(c1,nvirt*nvirt)
      call zeroe(ene,nvirt*nvirt)
      call zeroe(sqnor,nvirt*nvirt)
      endif

      do 1800 ir=ncore+nact+1,norb
      irsy=itsym(ir)

      do 1700 is=ir,norb
      issy=itsym(is)
      irssy=its(irsy,issy)
c      write (6,*) 'Eccitazione con gli indici',' -> ',ir,is

c
c     Si inizia il ciclo sugli attivi 
c
      
      dn=0.d0
      dk=0.d0

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
                  dum=daa(icpa,idpa,ica,ida)
                  dumk=gkm2aa(icpa,idpa,ica,ida)
                  dint1=ai(is,icp,ir,idp)
                  dint2=ai(is,ic,ir,id)
		  cint=dint1*dint2
		  if (ir.eq.is) then
		  cint=cint*0.25d0
		  dum=dum+daa(icpa,idpa,ida,ica)
		  dumk=dumk+gkm2aa(idpa,icpa,ica,ida)
		  endif
                  dn=dn+dum*cint
                  dk=dk+dumk*cint 
 1120          enddo
            enddo
 1121    enddo
      enddo

      if(abs(dn)*finv**2.gt.1.d-14)then
         epsimk=-dk/dn                     
         den=epsimk-fmpb(ir)-fmpb(is)
         irm=ir-ncore-nact
         ism=is-ncore-nact

           if(zthird) then
           sqnor(irm,ism)=sqrt(dn)
           c1(irm,ism)=sqnor(irm,ism)/den
           ene(irm,ism)=-den
           sqnor(irm,ism)=1.d0/sqnor(irm,ism)
           endif

c      write (6,'(A,4F15.10)') 'dn,eps,den',dn/factor**2,
c     *            epsimk/factor,den/factor
c      write (6,*) 'E2',e2mpp(5)/factor
         co=dn/den
         con=co/den
         e2mp(1)=e2mp(1)+co
         psimp(1)=psimp(1)+con
         e2mpp(5)=e2mpp(5)+co
         psimpp(5)=psimpp(5)+con
      endif
      
 1700 enddo
 1800 enddo

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c---Partially contracted NEV-PT-----------------------------------
      tini=etime(tarray)
      allocate(s(1:nact**2,1:nact**2))
      allocate(coef(1:nact**2,1:nact**2))
      allocate(w(nact**2))
      allocate(work(3*nact**2))
      call fillcm2(s,daa,nact)
      call fillcm2(coef,gkm2aa,nact)
      ioffw=1
      if(zprt)then
         call provasymm(daa,nact,discrep,discrep2)
         print '(a,f12.8,a,f12.8)','max. err. in herm. of daa=',discrep
     $        ,' in perm.=',discrep2
         call provasymm(gkm2aa,nact,discrep,discrep2)
         print '(a,f12.8,a,f12.8)','max. err. in herm. of K=',discrep
     $        /factor,' in perm.=',discrep2/factor
      endif
c      call matout(nact**2,nact**2,s,nact**2,1.d0)
c      call flush(6)
c      call dsygv(1,'V','U',nact**2,coef,nact**2,s,nact**2,w,work,3*nact
c     $     **2,info)
      nnew=nact**2
      nact3=nact**2
c      if(info.ne.0)then
c      print*,'something wrong in diagonalization: sub vm2mod'
c      print*,'trying to eliminate linear dependencies...'
c      print*,'Mind: this is not a normal situation'
c      call fillcm2(s,daa,nact)
c      call fillcm2(coef,gkm2aa,nact)
      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i+ioffw-1)
         enddo
      endif
      if (zprt) print*,'Analisi delle dipendenze lineari:'
      thr=1.d-8                 !perhaps too strict
      ndep=0
      do i=1,nact3
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
         call flush(6)
      endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef,nact3,0
     $     .d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
c     stop 'vm2mod'
      deallocate(caux)
c      endif
c      if(zverbose)then
c         print*,'info is ',info
c         print*,'caso r.ne.s'
c         print*,'eigenvalues are:'
c         do i=1,nnew
c            print '(a,i2,a,f12.8)','eps(',i,')=',w(i+ioffw-1)/factor
c         enddo
c      endif
      call e2m2cont(e2,psi2,coef,s,w,daa,fmpb,nact,ncore,norb,nnew,zgel)
      if(zprt)then
         print*,'Partially contracted contribution: class (-2): r.ne.s'
         print '(a,f20.14)','E2=',e2/factor
         print '(a,f20.14)','psi2=',psi2
      endif
      e2enp(5)=e2
      psienep(5)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      nact2=nact*(nact+1)/2
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(1:nact2,1:nact2))
      allocate(w(1:nact2))
      allocate(work(1:3*nact2))
      call fillcm2p(s,daa,nact)
      call fillcm2p(coef,gkm2aa,nact)
c      call dsygv(1,'V','U',nact2,coef,nact2,s,nact2,w,work,3*nact2
c     $     ,info)
      nact3=nact2
      nnew=nact2
c      if(info.ne.0)then
c     print*,'something wrong in diagonalization: sub vm2mod'
c     stop 'vm2mod'
c      print*,'trying to eliminate linear dependencies...'
c      print*,'Mind: this is not a normal situation'
c      call fillcm2p(s,daa,nact)
c      call fillcm2p(coef,gkm2aa,nact)
      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i+ioffw-1)
         enddo
      endif
      if (zprt) print*,'Analisi delle dipendenze lineari:'
      thr=1.d-8                 !perhaps too strict
      ndep=0
      do i=1,nact3
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      if(zprt)then
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
         call flush(6)
      endif
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef,nact3,0.d0
     $     ,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
c     stop 'vm2mod'
      deallocate(caux)
c      endif
      if(zprt)then
         print*,'info is ',info
         print*,'caso r.eq.s'
         print*,'eigenvalues are:'
         do i=1,nnew
            print '(a,i2,a,f12.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
      endif
      call e2m2contp(e2,psi2,coef,s,w,daa,fmpb,nact,ncore,norb,nnew,zgel
     $     )
      if(zprt)then
         print*,'Partially contracted contribution: class (-2): r.eq.s'
         print '(a,f20.14)','E2=',e2/factor
         print '(a,f20.14)','psi2=',psi2
      endif
      e2enp(5)=e2+e2enp(5)
      psienep(5)=psienep(5)+psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
c      print*,'Partially contracted contribution: class (-2): '
c      print '(a,f20.14)','E2=',e2enp(5)/factor
c      print '(a,f20.14)','psi2=',psienep(5)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
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
      end
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
      end
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
      end
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
      end
c----------------------------------------------------------------
      subroutine fillcm2(s,daa,nact)
      real*8 s(nact**2,nact**2),daa(nact,nact,nact,nact)
c---caso r#s
      ij=0
      do i=1,nact
         do j=1,nact
            ij=ij+1
            kl=0
            do k=1,nact
               do l=1,nact
                  kl=kl+1
                  s(ij,kl)=daa(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------
      subroutine fillcm2p(s,daa,nact)
      real*8 s(nact*(nact+1)/2,nact*(nact+1)/2),daa(nact,nact,nact,nact)
c---caso r=s
      ij=0
      do i=1,nact
         do j=i,nact
            ij=ij+1
            kl=0
            do k=1,nact
               do l=k,nact
                  kl=kl+1
                  s(ij,kl)=daa(i,j,k,l)+daa(i,j,l,k)
               enddo
            enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine fillc0p(s,coef,ga,gd,daa,da,nact)
      real*8 s(2*nact**2,2*nact**2),daa(nact,nact,nact,nact)
      real*8 coef(2*nact**2,2*nact**2),da(nact,nact)
      real*8 ga(nact,nact,nact,nact),gd(nact,nact,nact,nact)
      real*8 ee2
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
                  coef(ab,abp)=ga(b,a,ap,bp)
                  s(ab,abp)=2.d0*ee2(bp,ap,a,b,daa,da,nact)
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
                  coef(ab,abp)=-0.5d0*ga(b,a,ap,bp)
                  s(ab,abp)=-ee2(bp,ap,a,b,daa,da,nact)
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
                  coef(ab,abp)=-0.5d0*ga(b,a,ap,bp)
                  s(ab,abp)=-ee2(bp,ap,a,b,daa,da,nact)
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
                  coef(ab,abp)=gd(b,a,ap,bp)
                  s(ab,abp)=-daa(bp,a,b,ap)
                  if(a.eq.ap)s(ab,abp)=s(ab,abp)+2.d0*da(bp,b)
               enddo
            enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine fillcm1p(s,coef,amat,ro3,daa,da,nact,nact3)
      real*8 s(nact3,nact3),daa(nact,nact,nact,nact)
      real*8 coef(nact3,nact3),da(nact,nact)
      real*8 ro3(nact,nact,nact,nact,nact,nact)
      real*8 amat(nact,nact,nact,nact,nact,nact)
      real*8 eee2
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
                        coef(abc,abcp)=amat(ap,bp,cp,a,b,c)
                        s(abc,abcp)=eee2(cp,ap,bp,b,a,c,ro3,daa,da,nact)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine fillcp1p(s,coef,amat,ro3,daa,da,nact,nact3)
      real*8 s(nact3,nact3),daa(nact,nact,nact,nact)
      real*8 coef(nact3,nact3),da(nact,nact)
      real*8 ro3(nact,nact,nact,nact,nact,nact)
      real*8 amat(nact,nact,nact,nact,nact,nact)
      integer a,b,c,ap,bp,cp,abc,abcp
      real*8 eee2t
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
                        coef(abc,abcp)=amat(ap,bp,cp,a,b,c)
                        s(abc,abcp)=eee2t(cp,ap,bp,b,a,c,ro3,daa,da,nact
     $                       )
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c---------------------------------------------------------------
      subroutine v0pmod(ncore,nact,da,daa,gk0pa,gk0pd,gk0pf,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension gk0pa(nact,nact,nact,nact),gk0pd(nact,nact,nact,nact),
     * gk0pf(nact,nact)
      dimension daa(nact,nact,nact,nact)
      dimension da(nact,nact)
      allocatable heffd(:)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:),sp(:,:),spp(:
     $     ,:)
      allocatable c1(:,:),ene(:,:),sqnor(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      allocate(heffd(norb*(norb+1)/2))
      nvirt=norb-ncore-nact

      if (zthird) then
      allocate(c1(ncore,nvirt))
      allocate(ene(ncore,nvirt))
      allocate(sqnor(ncore,nvirt))
      call zeroe(c1,ncore*nvirt)
      call zeroe(ene,ncore*nvirt)
      call zeroe(sqnor,ncore*nvirt)
      endif

c----building Dyall's h eff----------
      do ik=1,ncore
         do ir=ncore+nact+1,norb
            jkr=indice(ik,ir)
            heffd(jkr)=onel(jkr)
            do j=1,ncore
               heffd(jkr)=heffd(jkr)+2.d0*ai(ik,ir,j,j)-ai(ik,j,j,ir)
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

c
c     Si inizia il ciclo sugli attivi 
c
      
      dn=0.d0
      dkoo=0.d0
c
c

      do ia=ncore+1,ncore+nact
         iaa=ia-ncore
         ias=itsym(ia)
         do ib=ncore+1,ncore+nact
            iba=ib-ncore
            ibs=itsym(ib)
	    iabs=its(ias,ibs)
	    if (its(iabs,iirs).ne.1) goto 1121
            do iap=ncore+1,ncore+nact
               iapa=iap-ncore
               iaps=itsym(iap)
               do ibp=ncore+1,ncore+nact
                  ibpa=ibp-ncore
                  ibps=itsym(ibp)
	          iabps=its(iaps,ibps)
	          if (its(iabps,iirs).ne.1) goto 1120
c		  if (its(iabs,iabps).ne.1) goto 1120


		  djp=ai(ir,ii,iap,ibp)
		  dj=ai(ir,ii,ia,ib)
		  dkp=ai(ir,ibp,iap,ii)
		  dk=ai(ir,ib,ia,ii)

		  dint1=2.d0*djp*dj-djp*dk-dkp*dj
		  dint2=dkp*dk

                  dum1=ee2(ibpa,iapa,iaa,iba,daa,da,nact)
		  dum2=ee2tr(ibpa,iba,iapa,iaa,daa,da,nact)
		  if(iaa.eq.iba) dum2=dum2+da(ibpa,iapa)

               dn=dn+dint1*dum1+dint2*dum2

               dumk1=gk0pa(ibpa,iapa,iaa,iba)
               dumk2=gk0pd(ibpa,iapa,iaa,iba)
	       dint1=djp*dj-0.5d0*djp*dk-0.5d0*dj*dkp
	       dint2=dkp*dk
	       dkoo=dkoo+dint1*dumk1+dint2*dumk2

 1120       enddo
         enddo
 1121 enddo
      enddo

      dummy=0.d0
      dummyk=0.d0
      do ia=ncore+1,ncore+nact
         iaa=ia-ncore
         do ib=ncore+1,ncore+nact
            iba=ib-ncore
	    dint=2.d0*ai(ir,ii,ia,ib)-ai(ir,ib,ia,ii)
	    dum=da(iba,iaa)
	    dummy=dummy+dum*dint
	    dint=ai(ir,ii,ia,ib)-0.5d0*ai(ir,ib,ia,ii)
	    dumk=gk0pf(iaa,iba)
	    dummyk=dummyk+dumk*dint
         enddo
      enddo

      iir=indice(ir,ii)
      dn=dn+2.d0*heffd(iir)*dummy+2.d0*heffd(iir)**2
      dkoo=dkoo+heffd(iir)*dummyk
      
c      if(dn.ne.0.d0)then
      if(abs(dn)*finv**2.gt.1.d-14)then
         epsimk=dkoo/dn                     
         den=fmpb(ii)-fmpb(ir)-epsimk
         irm=ir-ncore-nact

         if (zthird) then
         sqnor(ii,irm)=sqrt(dn)
         c1(ii,irm)=sqnor(ii,irm)/den
         ene(ii,irm)=-den
         sqnor(ii,irm)=1.d0/sqnor(ii,irm)
         endif

         co=dn/den
         con=co/den
         e2mp(1)=e2mp(1)+co
         psimp(1)=psimp(1)+con
         e2mpp(8)=e2mpp(8)+co
         psimpp(8)=psimpp(8)+con
      endif
      
 1700 enddo
 1800 enddo

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c---Partially contracted NEV-PT      
      tini=etime(tarray)
      nact2=2*nact**2
      allocate(s(1:nact2,1:nact2))
      allocate(coef(1:nact2,1:nact2))
      allocate(w(nact2))
      allocate(work(3*nact2))
      call fillc0p(s,coef,gk0pa,gk0pd,daa,da,nact)
      if(zprt)then
         call provasymm2(gk0pa,nact,discrep)
         print '(a,f12.8)','max. err. in herm. of gk0pa=',discrep/factor
         call provasymm2(gk0pd,nact,discrep)
         print '(a,f12.8)','max. err. in herm. of gk0pd=',discrep/factor
      endif
      call dsyev('V','U',nact2,s,nact2,w,work,3*nact2,info)
      if(info.ne.0)then
         print*,'something wrong in diagonalization: sub v0pmod'
         stop 'v0pmod'
      endif
      ioffw=1
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact2
            print '(i2,f14.8)',i,w(i+ioffw-1)
         enddo
      endif
c      print*,'prime due colonne autovettori:'
c      ij=0
c      do i=1,nact
c         do j=1,nact
c            ij=ij+1
c            print '(a,2i1,f14.8,f14.8)','ri',i,j,s(ij,1),s(ij,2)
c         enddo
c      enddo
c      do i=1,nact
c         do j=1,nact
c            ij=ij+1
c            print '(a,2i1,f14.8,f14.8)','ri',i,j,s(ij,1),s(ij,2)
c         enddo
c      enddo
      if (zprt) print*,'Analisi delle dipendenze lineari:'
      ndep=0
      ioffw=1
      thr=1.d-8
      do i=1,nact2
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i3,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact2-ndep
c--costruzione della matrice di trasformazione
      allocate (caux(1:nact2,1:nnew))
      call caux0p(caux,s,w,nact2,nnew,thr)
c      jnew=0
c      do j=1,nact2
c         if(abs(w(j)).ge.1.d-8)then
c            jnew=jnew+1
c            val=sqrt(w(j))
c            do i=1,nact2
c               caux(i,jnew)=s(i,j)/val
c            enddo
c         endif
c      enddo
c      if(zverbose)print*,'check: final jnew is',jnew
c---trasformazione della matrice K
c---moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact2,nnew,1.d0,coef,nact2,caux,nact2,0.d0,s
     $     ,nact2)
c---moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact2,1.d0,caux,nact2,s,nact2,0.d0
     $     ,coef,nact2)
c--diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact2,w,work,3*nact2,info)
      if(info.ne.0)then
         print*,'something wrong in diagonalization: sub v0pmod'
         stop 'v0pmod'
      endif
      if(zprt)then
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i2,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
      endif
      call flush(6)
c--back transformation
c--moltiplico CAUX(nact2,nnew) per coef(nnew,nnew)
c-- e metto in S(nact2,nnew)
      call dgemm('N','N',nact2,nnew,nnew,1.d0,caux,nact2,coef,nact2,0.d0
     $     ,s,nact2)
c      print*,'nuovi coefficienti'
c      call matout(nact2,nnew,s,nact2,1.d0)
c--metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact2,nnew)
c      do mu=1,nnew
c         do i=1,nact2
c            coef(i,mu)=s(i,mu)
c         enddo
c      enddo
c--inizio il calcolo perturbativo
      deallocate(s)
      deallocate(caux)
      allocate(s(1:nact**2,1:nnew))
      allocate(sp(1:nact**2,1:nnew))
      allocate(spp(1:nact**2,1:nnew))
      call e20pcont(e2,psi2,heffd,coef,s,sp,spp,w,daa,da,fmpb,nact,ncore
     $     ,norb,nnew,zro,zgel)
c      nn=nact**2
c      call zeroe(s,nn*nnew)
c      call zeroe(sp,nn*nnew)
c      call zeroe(spp,nn*nnew)
c      iab=0
c      do ia=1,nact
c         do ib=1,nact
c            iab=iab+1
c            iabp=0
c            do iap=1,nact
c               do ibp=1,nact
c                  iabp=iabp+1
c                  val=ee2(ibp,iap,ia,ib,daa,da,nact)
c                  val2=-daa(ibp,ia,ib,iap)
c                  if(iap.eq.ia)val2=val2+2.d0*da(ibp,ib)
c                  do mu=1,nnew
c                     s(iab,mu)=s(iab,mu)+val*coef(iabp,mu)
c                     sp(iab,mu)=sp(iab,mu)+val*coef(iabp+nn,mu)
c                     spp(iab,mu)=spp(iab,mu)+val2*coef(iabp+nn,mu)
c                  enddo
c               enddo
c            enddo
c         enddo
c      enddo
c      e2=0.d0
c      psi2=0.d0
c      do ir=ncore+nact+1,norb  !loop sui virtuali
c         do i=1,ncore          !loop sui core
c            iir=indice(ir,i)
c            hir=heffd(iir)
c            do mu=1,nnew
c               rimu=0.d0
c               rimup=0.d0
c               iab=0
c               do ia=1,nact
c                  do ib=1,nact
c                     iab=iab+1
c                     raib=ai(ir,i,ia+ncore,ib+ncore)
c                     rabi=ai(ir,ib+ncore,ia+ncore,i)
c                     rimu=rimu+(2.d0*raib-rabi)*s(iab,mu)+2.d0*da(ib,ia)
c     $                    *hir*coef(iab,mu)
c                     rimup=rimup-raib*sp(iab,mu)+rabi*spp(iab,mu)-hir
c     $                    *da(ib,ia)*coef(iab+nn,mu)
c                  enddo
c               enddo
c               val=(rimu+rimup)**2
c               deno=w(mu)+fmpb(ir)-fmpb(i)
c               cont=val/deno
c               e2=e2-cont
c               psi2=psi2+cont/deno
c            enddo
c         enddo
c      enddo
c      print*,'Partially contracted contribution: class (0''): '
c      print '(a,f20.14)','E2=',e2/factor
c      print '(a,f20.14)','psi2=',psi2
      e2enp(8)=e2
      psienep(8)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      deallocate(s)
      deallocate(sp)
      deallocate(spp)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      deallocate(heffd)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
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
      end
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
      end
c----------------------------------------------------------------
      subroutine e20pcont(e2,psi2,heffd,coef,s,sp,spp,w,daa,da,fmpb,nact
     $     ,ncore,norb,nnew,zro,zgel)
      implicit real*8 (a-h,o-z)
      logical*1 zro,zgel
      allocatable rhoir(:,:),vaux(:)
      dimension coef(2*nact**2,2*nact**2),s(nact**2,nnew),
     $     sp(nact**2,nnew),spp(nact**2,nnew),w(*),fmpb(*),
     $     daa(nact,nact,nact,nact),da(nact,nact),heffd(*),zgel(*)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nn=nact**2
      nvirt=norb-ncore-nact
      if(zro)then
         open(27,file='RHOFILE',form='UNFORMATTED',status='UNKNOWN')
         rewind 27
         allocate(rhoir(1:ncore,1:nvirt))
         allocate(vaux(1:nnew))
         call zeroe(rhoir,ncore*nvirt)
         call zeroe(vaux,nnew)
      endif
      call zeroe(s,nn*nnew)
      call zeroe(sp,nn*nnew)
      call zeroe(spp,nn*nnew)
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
                  val=ee2(ibp,iap,ia,ib,daa,da,nact)
                  val2=-daa(ibp,ia,ib,iap)
                  if(iap.eq.ia)val2=val2+2.d0*da(ibp,ib)
                  do mu=1,nnew
                     s(iab,mu)=s(iab,mu)+val*coef(iabp,mu)
                     sp(iab,mu)=sp(iab,mu)+val*coef(iabp+nn,mu)
                     spp(iab,mu)=spp(iab,mu)+val2*coef(iabp+nn,mu)
                     if(zro.and.iabp.eq.1)vaux(mu)=vaux(mu)+da(ia,ib)*
     $                    (2.d0*coef(iab,mu)-coef(iab+nn,mu))
                  enddo
               enddo
            enddo
         enddo
      enddo
      e2=0.d0
      psi2=0.d0
      do i=1,ncore              !loop sui core
         if (zgel(i)) goto 1
         do ir=ncore+nact+1,norb !loop sui virtuali
            iir=indice(ir,i)
            hir=heffd(iir)
            do mu=1,nnew
               rimu=0.d0
               rimup=0.d0
               iab=0
               do ia=1,nact
                  do ib=1,nact
                     iab=iab+1
                     raib=ai(ir,i,ia+ncore,ib+ncore)
                     rabi=ai(ir,ib+ncore,ia+ncore,i)
                     rimu=rimu+(2.d0*raib-rabi)*s(iab,mu)+2.d0*da(ib,ia)
     $                    *hir*coef(iab,mu)
                     rimup=rimup-raib*sp(iab,mu)+rabi*spp(iab,mu)-hir
     $                    *da(ib,ia)*coef(iab+nn,mu)
                  enddo
               enddo
               co=rimu+rimup
               deno=w(mu)+fmpb(ir)-fmpb(i)
               cirmu=co/deno
               cont=co*cirmu
               e2=e2-cont
               psi2=psi2+cirmu**2
               if(zro)then
                  ivirt=ir-ncore-nact
                  rhoir(i,ivirt)=rhoir(i,ivirt)-cirmu*vaux(mu)
               endif
            enddo
c            print*,'i,ivirt,rho',i,ivirt,rhoir(i,ivirt)
         enddo
 1    enddo
      if(zro)write(27)((rhoir(i,ivirt),i=1,ncore),ivirt=1,nvirt)
      if(zro)deallocate(rhoir)
      return
      end
c----------------------------------------------------------------
      subroutine e21pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb,nact
     $     ,ncore,norb,nnew,zro,zgel)
      implicit real*8 (a-h,o-z)
      logical*1 zro,zgel
      dimension coef(nact**3,nact**3),s(nnew,nact**3),
     $     sp(nnew,nact**3),w(*),fmpb(*),ro3(nact,nact,nact,nact,nact
     $     ,nact),daa(nact,nact,nact,nact),da(nact,nact),heffd(*),zgel(
     $     *)
      real*4 tarray(2)
      logical*1 zcond
      allocatable rhoid(:,:),caux(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      t3=etime(tarray)
      nact3=nact**3
      if(zro)then
         allocate(rhoid(1:ncore,1:nact))
         allocate(caux(1:nact,1:nnew))
         call zeroe(rhoid,nact*ncore)
         call zeroe(caux,nact*nnew)
      endif
      call zeroe(s,nact3*nnew)
      call zeroe(sp,nact3*nnew)
      iabc=0
      do ia=1,nact
         do ib=1,nact
            do ic=1,nact
               zcond=ib.eq.1.and.ic.eq.1
               iabc=iabc+1
               iabcp=0
               do iap=1,nact
                  do ibp=1,nact
                     do icp=1,nact
                        iabcp=iabcp+1
                        val=eee2t(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
                        if(zcond)val2=ee2tr(icp,iap,ibp,ia,daa,da,nact)
                        do mu=1,nnew
                           s(mu,iabc)=s(mu,iabc)+val*coef(mu,iabcp)
                           if(zcond)sp(mu,ia)=sp(mu,ia)+val2*coef(mu
     $                          ,iabcp)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c      if(zro)then
c         iabc=0
c         do ia=1,nact
c            do ib=1,nact
c               do ic=1,nact
c                  iabc=iabc+1
c                  do id=1,nact
c                     val=ee2tl(id,ib,ia,ic,daa,da,nact)
c                     do mu=1,nnew
c                        caux(id,mu)=caux(id,mu)+coef(iabc,mu)*val
c                     enddo
c                  enddo
c               enddo
c            enddo
c         enddo
c      endif
      if(zro)then
         do id=1,nact
            do mu=1,nnew
               iabc=0
               do ia=1,nact
                  do ib=1,nact
                     do ic=1,nact
                        iabc=iabc+1
                        val=ee2tl(id,ib,ia,ic,daa,da,nact)
                        caux(id,mu)=caux(id,mu)+coef(mu,iabc)*val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c-------------------------------------------------------
      t4=etime(tarray)
c      print '(a,f8.2,a)','time to build S',t4-t3,' sec.'
      e2=0.d0
      psi2=0.d0
      do i=1,ncore   !loop sugli orbitali di core
        if (zgel(i)) goto 1 
         do mu=1,nnew
            aimu=0.d0
            iabc=0
            do ia=1,nact
               ira=indice(i,ia+ncore)
               hia=heffd(ira)
               aimu=aimu+hia*sp(mu,ia)
               do ib=1,nact
                  do ic=1,nact
                     iabc=iabc+1
                     baic=ai(ib+ncore,i,ia+ncore,ic+ncore)
                     aimu=aimu+baic*s(mu,iabc)
                  enddo
               enddo
            enddo
            deno=w(mu)-fmpb(i)
            co=aimu/deno
            cont=aimu*co
            e2=e2-cont
            psi2=psi2+co**2
            if(zro)then
               do id=1,nact
                  rhoid(i,id)=rhoid(i,id)-co*caux(id,mu)
               enddo
            endif
         enddo
c         print*,'i,rho',i,(rhoid(i,id),id=1,nact)
 1    enddo
      t5=etime(tarray)
      if(zro)write(27)((rhoid(i,id),i=1,ncore),id=1,nact)
      if(zro)deallocate(rhoid)
c      print '(a,f8.2,a)','time for pert. calcn. ',t5-t4,' sec.'
      return
      end
c------------------------------------------------------------------
      subroutine e2m1pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb
     $     ,nact,ncore,norb,nnew,zro,zgel)
      implicit real*8 (a-h,o-z)
      logical*1 zro,zgel
      dimension coef(nact**3,nact**3),s(nnew,nact**3),
     $     sp(nnew,nact**3),w(*),fmpb(*),ro3(nact,nact,nact,nact,nact
     $     ,nact),daa(nact,nact,nact,nact),da(nact,nact),heffd(*),zgel(
     $     *)
      allocatable rhodr(:,:),caux(:,:)
      real*4 tarray(2)
      logical*1 zcond
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      t3=etime(tarray)
      nact3=nact**3
      nvirt=norb-ncore-nact
      if(zro)then
         allocate(rhodr(1:nact,1:nvirt))
         allocate(caux(1:nact,1:nnew))
         call zeroe(rhodr,nact*nvirt)
         call zeroe(caux,nact*nnew)
      endif
      call zeroe(s,nact3*nnew)
      call zeroe(sp,nact3*nnew)
      iabc=0
      do ia=1,nact
         do ib=1,nact
            do ic=1,nact
               zcond=ib.eq.1.and.ic.eq.1
               iabc=iabc+1
               iabcp=0
               do iap=1,nact
                  do ibp=1,nact
                     do icp=1,nact
                        iabcp=iabcp+1
                        val=eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
                        if(zcond)val2=ee2(icp,iap,ibp,ia,daa,da,nact)
                        do mu=1,nnew
                           s(mu,iabc)=s(mu,iabc)+val*coef(mu,iabcp)
                           if(zcond)sp(mu,ia)=sp(mu,ia)+val2*coef(mu
     $                          ,iabcp)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c      if(zro)then
c         iabc=0
c         do ia=1,nact
c            do ib=1,nact
c               do ic=1,nact
c                  iabc=iabc+1
c                  do id=1,nact
c                     val=ee2(id,ib,ia,ic,daa,da,nact)
c                     do mu=1,nnew
c                        caux(id,mu)=caux(id,mu)+coef(iabc,mu)*val
c                     enddo
c                  enddo
c               enddo
c            enddo
c         enddo
c      endif
      if(zro)then
         do id=1,nact
            do mu=1,nnew
               iabc=0
               do ia=1,nact
                  do ib=1,nact
                     do ic=1,nact
                        iabc=iabc+1
                        val=ee2(id,ib,ia,ic,daa,da,nact)
                        caux(id,mu)=caux(id,mu)+coef(mu,iabc)*val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
c-------------------------------------------------------
      t4=etime(tarray)
c      print '(a,f8.2,a)','time to build S',t4-t3,' sec.'
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb   !loop sui virtuali
         do mu=1,nnew
            rmu=0.d0
            iabc=0
            do ia=1,nact
               ira=indice(ir,ia+ncore)
               hra=heffd(ira)
               rmu=rmu+hra*sp(mu,ia)
               do ib=1,nact
                  do ic=1,nact
                     iabc=iabc+1
                     rabc=ai(ir,ib+ncore,ia+ncore,ic+ncore)
                     rmu=rmu+rabc*s(mu,iabc)
                  enddo
               enddo
            enddo
            deno=w(mu)+fmpb(ir)
            co=rmu/deno
            cont=rmu*co
            e2=e2-cont
            psi2=psi2+co**2
            if(zro)then
               ivirt=ir-nact-ncore
               do id=1,nact
                  rhodr(id,ivirt)=rhodr(id,ivirt)-co*caux(id,mu)
               enddo
            endif
         enddo
c         print*,'ir,rho',ir,(rhodr(id,ivirt),id=1,nact)
      enddo
      t5=etime(tarray)
      if(zro)write(27)((rhodr(id,ivirt),id=1,nact),ivirt=1,nvirt)
      if(zro)deallocate(rhodr)
c      print '(a,f8.2,a)','time for pert. calcn. ',t5-t4,' sec.'
      return
      end
c----------------------------------------------------------------
      subroutine vpupE(ncore,nact,ro3,amat,bmat,cmat,dmat,daa,da,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension nosym(8)
      dimension ro3(nact,nact,nact,nact,nact,nact),
     *          amat(nact,nact,nact,nact,nact,nact),
     *          bmat(nact,nact,nact,nact),
     *          cmat(nact,nact,nact,nact),
     *          daa(nact,nact,nact,nact),
     *          dmat(nact,nact),
     *          da(nact,nact)
      allocatable heffd(:)
      allocatable s(:,:),coef(:,:),sp(:,:),caux(:,:),w(:),work(:)
      allocatable c1(:),ene(:),sqnor(:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      allocate(heffd(norb*(norb+1)/2))

      if (zthird) then
      allocate(c1(ncore))
      allocate(ene(ncore))
      allocate(sqnor(ncore))
      call zeroe(c1,ncore)
      call zeroe(ene,ncore)
      call zeroe(sqnor,ncore)
      endif

c----building Dyall's h eff----------
      do ik=ncore+1,ncore+nact
         do ir=1,ncore
            jkr=indice(ik,ir)
            heffd(jkr)=0.d0
            heffd(jkr)=onel(jkr)
            do j=1,ncore
               heffd(jkr)=heffd(jkr)+2.d0*ai(ik,ir,j,j)-ai(ik,j,j,ir)
            enddo
         enddo
      enddo


      do 800 ii=1,ncore
         if (zgel(ii)) goto 800
         irsy=itsym(ir)

c      write (6,*) 
c      write (6,*) 'Eccitazione con gli indici',ii,' -> '
c      call flush(6)

      dumm=0.d0
      dumk=0.d0

      do ia=1,nact
       iac=ia+ncore
       iia=indice(ii,iac)
       do iap=1,nact
        iapc=iap+ncore
        iiap=indice(ii,iapc)
	rodum=-da(iap,ia)
	if (iap.eq.ia) rodum=rodum+2.d0
        dumm=dumm+heffd(iiap)*rodum*heffd(iia)
        dumk=dumk+heffd(iiap)*dmat(iap,ia)*heffd(iia)
        do ibp=1,nact
         ibpc=ibp+ncore
         do icp=1,nact
          icpc=icp+ncore
          dint=ai(ibpc,ii,iapc,icpc)
          dumm=dumm+2.d0*dint*ee2tr(icp,iap,ibp,ia,daa,da,nact)
     *             *heffd(iia)
          dumk=dumk+dint*bmat(iap,ibp,icp,ia)*heffd(iia)
          dumk=dumk+dint*cmat(ia,iap,ibp,icp)*heffd(iia)
          do ib=1,nact
           ibc=ib+ncore
           do ic=1,nact
            icc=ic+ncore
            dint2=ai(ibc,ii,iac,icc)
            dumm=dumm+dint*eee2t(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
     *           *dint2
            dumk=dumk+dint*amat(iap,ibp,icp,ia,ib,ic)*dint2
            enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      dn=dumm
      dk=dumk
      if(abs(dn)*finv**2.gt.1.d-14)then
         epsimk=dk/dn                     
         den=fmpb(ii)-epsimk
         co=dn/den
         con=co/den

         if (zthird) then
         sqnor(ii)=sqrt(dn)
         c1(ii)=sqnor(ii)/den
         ene(ii)=-den
         sqnor(ii)=1.d0/sqnor(ii)
         endif

         e2mp(1)=e2mp(1)+co
         psimp(1)=psimp(1)+con
         e2mpp(6)=e2mpp(6)+co
         psimpp(6)=psimpp(6)+con
      endif

  800 enddo

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      write (32) e2mp(1)        !second order energy
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c--Partially contracted NEV-PT class (+1')
      tini=etime(tarray)
      nact3=nact**3
      allocate(s(1:nact3,1:nact3))
      allocate(coef(1:nact3,1:nact3))
      allocate(w(nact3))
      allocate(work(10))
      call fillcp1p(s,coef,amat,ro3,daa,da,nact,nact3)
      if(zprt)then
         call provasymm3(amat,nact,discrep)
         print '(a,f12.8)','max. err. in herm. of amat=',discrep/factor
         call provasymm4(s,nact3,discrep)
         print '(a,f12.8)','max. err. in herm. of S=',discrep/factor
      endif
      call dsyev('V','U',nact3,s,nact3,w,work,-1,info)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork),stat=info)
      call dsyev('V','U',nact3,s,nact3,w,work,lwork,info)
      if(info.ne.0)then
         print*,'something wrong in diagonalization: sub vpupe'
         stop 'vpupe'
      endif
      ioffw=1
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i+ioffw-1)
         enddo
      t=etime(tarray)
      print '(a,f8.2,a)','time for diagonalization of S',t-tini,' sec.'
      print*,'Analisi delle dipendenze lineari:'
      endif
      thr=1.d-8  !perhaps too strict
c      thr=1.d-5  !perhaps too large
      ndep=0
         ioffw=1
      do i=1,nact3
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
c--costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w,nact3,nnew,thr)
c---trasformazione della matrice K
c---moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
      t2=etime(tarray)
      if (zprt)print '(a,f8.2,a)','time for matrix transf.',t2-t,' sec.'
c--diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact3,w,work,-1,info)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork),stat=info)
      call dsyev('V','U',nnew,coef,nact3,w,work,lwork,info)
      if(info.ne.0)then
         print*,'something wrong in diagonalization: sub vpupe'
         stop 'vpupe'
      endif
      t3=etime(tarray)
      if(zprt)then
      print '(a,f8.2,a)','time for diagonalization of H ',t3-t2,' sec.'
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
      endif
c--back transformation
c--moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c-- e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef,nact3,0.d0
     $     ,s,nact3)
c      call matout(nact3,nnew,s,nact3,1.d0)
c--metto in coef i coefficienti trasformati (for clarity sake)
      do mu=1,nnew
         do i=1,nact3
            coef(mu,i)=s(i,mu)
         enddo
      enddo
c--inizio il calcolo perturbativo
c--costruiamo le matrici ausiliarie S(abc,mu) e S'(a,mu)
      deallocate(s)
      deallocate(caux)
      allocate(s(1:nnew,1:nact3))
      allocate(sp(1:nnew,1:nact3))
      call e21pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb,nact,ncore
     $     ,norb,nnew,zro,zgel)
c      print*,'Partially contracted contribution: class (+1''): '
c      print '(a,f20.14)','E2=',e2/factor
c      print '(a,f20.14)','psi2=',psi2
      e2enp(6)=e2
      psienep(6)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      deallocate(s)
      deallocate(sp)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      deallocate(heffd)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
c----------------------------------------------------------------
      subroutine vmupE(ncore,nact,ro3,amat,bmat,cmat,dmat,daa,da,zthird)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      pointer e,een,emp,e2en,e2mp,fmpb,psimp,psien
      dimension e(:),een(:),emp(:),e2en(:),e2mp(:),fmpb(:),psimp(:)
     $     ,psien(:)
      COMMON /VEC/ E,EEN,EMP,E2EN,E2MP,FMPB,psimp,psien,NTOT
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,itsym
      dimension nd(:),ne(:),trou(:),part(:),itsym(:)
      COMMON /DET/ ND,NE,TROU,PART,ITSYM,ITS(20,20)
      pointer onel,igelo
      dimension onel(:),igelo(:)
      COMMON /PERTURBER/onel,e2enp(8),psienep(8),
     * e2mpp(8),psimpp(8),igelo,zverbose,zro
      real*4 tarray(2)
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 NEP,TROUP,PARTP
      common/debug/zdebug
      common /deephole/ ndhb,ndha,zdh
      pointer zact
      dimension zact(:)
      common/detinf/ntotdet,ntotconf,ntotcap,zact
c      integer*2 icomp,iconf,isegno
      dimension nosym(8)
      dimension ro3(nact,nact,nact,nact,nact,nact),
     *          amat(nact,nact,nact,nact,nact,nact),
     *          bmat(nact,nact,nact,nact),
     *          cmat(nact,nact,nact,nact),
     *          daa(nact,nact,nact,nact),
     *          dmat(nact,nact),
     *          da(nact,nact)
      allocatable heffd(:)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:),sp(:,:)
      allocatable c1(:),ene(:),sqnor(:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      allocate(heffd(norb*(norb+1)/2))
      nvirt=norb-ncore-nact

      if (zthird) then
      allocate(c1(nvirt))
      allocate(ene(nvirt))
      allocate(sqnor(nvirt))
      call zeroe(c1,nvirt)
      call zeroe(ene,nvirt)
      call zeroe(sqnor,nvirt)
      endif

c----building Dyall's h eff----------
      do ik=ncore+1,ncore+nact
         do ir=ncore+nact+1,norb
            jkr=indice(ik,ir)
            heffd(jkr)=0.d0
            heffd(jkr)=onel(jkr)
            do j=1,ncore
               heffd(jkr)=heffd(jkr)+2.d0*ai(ik,ir,j,j)-ai(ik,j,j,ir)
            enddo
	    do j=ncore+1,ncore+nact
                heffd(jkr)=heffd(jkr)-ai(ir,j,j,ik)
	    enddo
         enddo
      enddo


      do 800 ir=ncore+nact+1,norb
      irsy=itsym(ir)

c      write (6,*) 
c      write (6,*) 'Eccitazione con gli indici',' -> ',ir
c      call flush(6)

      dumm=0.d0
      dumk=0.d0

      do ia=1,nact
       iac=ia+ncore
       ira=indice(ir,iac)
       do iap=1,nact
        iapc=iap+ncore
        irap=indice(ir,iapc)
        dumm=dumm+heffd(irap)*da(iap,ia)*heffd(ira)
        dumk=dumk+heffd(irap)*dmat(iap,ia)*heffd(ira)
        do ibp=1,nact
         ibpc=ibp+ncore
         do icp=1,nact
          icpc=icp+ncore
          dint=ai(ir,ibpc,iapc,icpc)
          dumm=dumm+2.d0*dint*ee2(icp,iap,ibp,ia,daa,da,nact)
     *             *heffd(ira)
          dumk=dumk+dint*bmat(iap,ibp,icp,ia)*heffd(ira)
          dumk=dumk+dint*cmat(ia,iap,ibp,icp)*heffd(ira)
          do ib=1,nact
           ibc=ib+ncore
           do ic=1,nact
            icc=ic+ncore
            dint2=ai(ir,ibc,iac,icc)
            dumm=dumm+dint*eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
     *           *dint2
            dumk=dumk+dint*amat(iap,ibp,icp,ia,ib,ic)*dint2
            enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      dn=dumm
      dk=dumk
      if(abs(dn)*finv**2.gt.1.d-14)then
         epsimk=-dk/dn                     
         den=epsimk-fmpb(ir)
         co=dn/den
         con=co/den
         irm=ir-ncore-nact

         if (zthird) then
         sqnor(irm)=sqrt(dn)
         c1(irm)=sqnor(irm)/den
         ene(irm)=-den
         sqnor(irm)=1.d0/sqnor(irm)
         endif

         e2mp(1)=e2mp(1)+co
         psimp(1)=psimp(1)+con
         e2mpp(7)=e2mpp(7)+co
         psimpp(7)=psimpp(7)+con
      endif

  800 enddo

      if (zthird) then
      write (32) c1
      write (32) ene
      write (32) sqnor
      deallocate(c1)
      deallocate(ene)
      deallocate(sqnor)
      endif

c--Partially contracted NEV-PT class (-1')
      tini=etime(tarray)
      nact3=nact**3
      allocate(s(1:nact3,1:nact3))
      allocate(coef(1:nact3,1:nact3))
      allocate(w(nact3))
      allocate(work(10))
      call fillcm1p(s,coef,amat,ro3,daa,da,nact,nact3)
      if(zprt)then
         call provasymm3(amat,nact,discrep)
         print '(a,f12.8)','max. err. in herm. of amat=',discrep/factor
         call provasymm4(s,nact3,discrep)
         print '(a,f12.8)','max. err. in herm. of S=',discrep/factor
      endif
      call dsyev('V','U',nact3,s,nact3,w,work,-1,info)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork),stat=info)
      call dsyev('V','U',nact3,s,nact3,w,work,lwork,info)
      if(info.ne.0)then
         print*,'something wrong in diagonalization: sub vmupe'
         stop 'vmupe'
      endif
      ioffw=1
      if(zprt)then
         print*,'matrice metrica diagonalizzata con info=',info
         print*,'autovalori matrice metrica:'
         do i=1,nact3
            print '(i4,f14.8)',i,w(i+ioffw-1)
         enddo
      t=etime(tarray)
      print '(a,f8.2,a)','time for diagonalization of S',t-tini,' sec.'
      print*,'Analisi delle dipendenze lineari:'
      endif
      thr=1.d-8  !perhaps too strict
c      thr=1.d-5  !perhaps too large
      ndep=0
         ioffw=1
      do i=1,nact3
         if(w(ioffw+i-1).lt.0.d0)then
            ndep=ndep+1
            if (zprt)
     $      print '(a,i3,a,f12.8)','negative eigenvalue: eps(',i,')='
     $           ,w(ioffw+i-1)
         elseif(w(ioffw+i-1).lt.thr)then
            ndep=ndep+1
         endif
      enddo
      if (zprt) print '(a,i4,a)','Ci sono',ndep,' dipendenze lineari'
      nnew=nact3-ndep
c--costruzione della matrice di trasformazione
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,s,w,nact3,nnew,thr)
      if(zverbose)print*,'check: final jnew is',jnew
c---trasformazione della matrice K
c---moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
      t2=etime(tarray)
      if (zprt)print '(a,f8.2,a)','time for matrix transf.',t2-t,' sec.'
c--diagonalizzo la matrice coef(hamiltoniana trasformata)
      call dsyev('V','U',nnew,coef,nact3,w,work,-1,info)
      lwork=work(1)
c      print*,'best dimension for lwork is ',lwork
      deallocate(work)
      allocate(work(lwork),stat=info)
      call dsyev('V','U',nnew,coef,nact3,w,work,lwork,info)
      if(info.ne.0)then
         print*,'something wrong in diagonalization: sub vmupe'
         stop 'vmupe'
      endif
      t3=etime(tarray)
      if(zprt)then
      print '(a,f8.2,a)','time for diagonalization of H ',t3-t2,' sec.'
      call flush(6)
         print*,'diagonalization got info=',info
         print*,'eigenvalues:'
         do i=1,nnew
            print '(a,i4,a,f15.8)','eps(',i,')=',w(i+ioffw-1)/factor
         enddo
         call flush(6)
      endif
c--back transformation
c--moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c-- e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,coef,nact3,0.d0
     $     ,s,nact3)
c      print*,'nuovi coefficienti'
c      call matout(nact3,nnew,s,nact3,1.d0)
c--metto in coef i coefficienti trasformati (for clarity sake)
c      call copia(coef,s,nact3,nnew)
      do mu=1,nnew
         do i=1,nact3
            coef(mu,i)=s(i,mu)
         enddo
      enddo
c--inizio il calcolo perturbativo
c--costruiamo le matrici ausiliarie S(abc,mu) e S'(a,mu)
      deallocate(s)
      deallocate(caux)
      allocate(s(1:nnew,1:nact3),stat=info)
      if(info.ne.0)then
         print*,'cannot allocate s'
         call flush(6)
         stop
      endif
      allocate(sp(1:nnew,1:nact3),stat=info)
      if(info.ne.0)then
         print*,'cannot allocate s'
         call flush(6)
         stop
      endif
      call e2m1pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb,nact
     $     ,ncore,norb,nnew,zro,zgel)
c      print*,'Partially contracted contribution: class (-1''): '
c      print '(a,f20.14)','E2=',e2/factor
c      print '(a,f20.14)','psi2=',psi2
c      call flush(6)
      e2enp(7)=e2
      psienep(7)=psi2
      e2en(1)=e2en(1)+e2
      psien(1)=psien(1)+psi2
      deallocate(s)
      deallocate(sp)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      deallocate(heffd)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t-tini,' sec.'
      return
      end
C***********************************************************************
      subroutine deforb(zact,norb,nocb,ncore,nact)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      dimension zact(*)

      ncore=0
      nact=0
      zfin=.false.
      do i=1,nocb
      if (.not.zact(i)) then
       if (zfin) then
       write (6,*) 'Errore: orbitali non ordinati correttamente'
       stop
      endif
      ncore=ncore+1
      else
      zfin=.true. 
      nact=nact+1
      endif
      enddo

      zfin=.false.
      do i=nocb+1,norb
      if (zact(i)) then
       if (zfin) then
       write (6,*) 'Errore: orbitali non ordinati correttamente'
       stop
      endif
      nact=nact+1
      else
      zfin=.true. 
      endif
      enddo

c      write (6,*) 'ncore =',ncore
c      write (6,*) 'nact =',nact

      return
      end
c-------------------------------------------------------------------
      function eee2(a,b,c,d,e,f,daaa,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
c--the same as eee but with the spinless density matrices as input
      eee2=daaa(a,c,e,b,d,f)
      if(e.eq.d)eee2=eee2+taa(a,c,b,f)
      if(e.eq.b)eee2=eee2+taa(a,c,f,d)
      if(b.eq.c)eee2=eee2+taa(a,e,d,f)
      if(b.eq.c.and.e.eq.d)eee2=eee2+dal(a,f)
      return
      end
c-------------------------------------------------------------------
      function ee2(a,b,c,d,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      ee2=taa(a,c,b,d)
      if(b.eq.c)ee2=ee2+dal(a,d)
      return
      end
c-------------------------------------------------------------------
      function eee2t(a,b,c,d,e,f,daaa,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
      eee2t=-eee2(a,b,d,c,e,f,daaa,taa,dal,nact)
      if(c.eq.d)eee2t=eee2t+2.d0*ee2(a,b,e,f,taa,dal,nact)
      return
      end
c-------------------------------------------------------------------
      function ee2tr(a,b,c,d,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      ee2tr=-ee2(a,b,d,c,taa,dal,nact)
      if(c.eq.d)ee2tr=ee2tr+2.d0*dal(a,b)
      return
      end
c-------------------------------------------------------------------
      function ee2tl(a,b,c,d,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      ee2tl=-ee2(b,a,c,d,taa,dal,nact)
      if(a.eq.b)ee2tl=ee2tl+2.d0*dal(c,d)
      return
      end
c-------------------------------------------------------------------
      SUBROUTINE REIJKL(factor,norb,nsym,itsym,its,nos,jd1)
C 
C     PREPARATION POUR LE STOCKAGE DES INTEGRALES BIELECTONIQUES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      dimension itsym(*)
      DIMENSION NOS(20),ITS(20,20)
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      common/caratte/cara1,cara2
      character*1 cara1,cara2
      dimension norbsym(8)
      common /toul/ ztoul
c      cara1=char(-1)
      cara2=char(0)
c      print*,'REIJKL compiled with intel IFC compiler'
      rewind 25
      READ (25) NSYM,NORB,NOC,NDEGEN,(ITSYM(I),I=1,NORB)
      READ (25) ((ITS(I,J),I=1,NSYM),J=1,NSYM)
      IJ=0
      NUM(1)=0
      DO 1 I=2,JD1
      DO 1 J=1,I
      IJ=IJ+1
    1 NUM(IJ+1)=NUM(IJ)+IJ
      DO 801 I=1,NSYM
  801 NOS(I)=0
      DO 802 I=1,NORB
      IS=ITSYM(I)
      DO 803 J=1,I
      JS=ITSYM(J)
      IJS=ITS(IS,JS)
      NOS(IJS)=NOS(IJS)+1
      IJ=NUM(I)+J
      INDIC(IJ)=IJS
  803 JNDIC(IJ)=NOS(IJS)
  802 CONTINUE
      NDEB(1)=0
      NSM1=NSYM-1
      DO 804 K=1,NSM1
      KP1=K+1
      if (kp1.gt.500) then
	 write (6,*) ' errore, KP1 =',kp1
	 stop 12
      endif
  804 NDEB(KP1)=NDEB(K)+(NOS(K)*(NOS(K)+1))/2

      dnijkl=1.d0*NDEB(NSYM)+(1.d0*NOS(NSYM)*(NOS(NSYM)+1))/2
      if (dnijkl.gt.2147483646.d0) then
      write (6,*) 'Too many integrals:'
      write (6,*) 'the dimension of IJKL is larger than 2147483646)'
      stop
      endif
c      NIJKL=NDEB(NSYM)+(NOS(NSYM)*(NOS(NSYM)+1))/2
      nijkl=nint(dnijkl)

      NJP=0
      DO 820 IS=1,NSYM
      DO 820 JS=1,IS
      NJP=NJP+1
      if (njp.gt.2000) then
	 write (6,*) ' errore, NJP =',njp
	 stop 12
      endif
      IF(JS.EQ.IS)GOTO 821
      NAD(NJP)=0
      KT(NJP)=0
      GOTO 820
  821 NAD(NJP)=IS
      KT(NJP)=1
  820 CONTINUE
      NSYT=NJP
      ni2=1
      nr4=1
      ni3=1
      ni4=nijkl
      WRITE (6,74) NIJKL,4
      if (zfac) then
c_toul
         if (ztoul) then
          call intmaxt(bigint,norb,nsym,nos)
         else
          call intmax(bigint,norb,nsym,nos)
         endif
c_toul
         print'(''     MAXIMUM INTEGRAL IS '',D10.4)',bigint
         factor=2147483646.d0/bigint
      endif
      print'(''     FACTOR WILL BE '',D10.4)',factor
 74   FORMAT (//5X,'DIMENSION OF THE IJKL TABLE:',I10,' *',i1,' BYTES')
      RETURN
      END
C***********************************************************************
      SUBROUTINE intmax(bigint,norb,nsym,nos)
C 
C     trova l'integrale piu' grande
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      PARAMETER(MBUF=3255,NBUF=4096)
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      bigint=0.d0
c      WRITE (6,1000)
c 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
      if (.not.zmotra) then
      call filesplit('READ',50,8,ichange,idum2)
      READ (50) FACDUM
      endif
 20   continue
      call filesplit('READ',50,12*mbuf,ichange,idum2)
      READ(50)XX,IX
      do n=1,mbuf
         lab=ix(N)
         if (lab.eq.0) return
         A8=abs(XX(N))
         if (a8.gt.bigint) bigint=a8
      enddo
      goto 20
      END
C***********************************************************************
      SUBROUTINE intmaxt(bigint,norb,nsym,nos)
C
C     trova l'integrale piu' grande
C
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      PARAMETER(MBUF=3255,NBUF=4096)
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
      bigint=0.d0
c      WRITE (6,1000)
c 1000 FORMAT (5X,'LECTURE DES INTEGRALES tolosa')
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
 20   continue
      call filesplit('READ',50,16*nbuf,ichange,idum2)
      read(50)ind,jnd,knd,lnd,xxm
      do n=1,nbuf
       lab=knd(n)
       if (lab.eq.0) return
       a8=abs(xxm(n))
       if (a8.gt.bigint) bigint=a8
      enddo
      goto 20
      END
C***********************************************************************
      SUBROUTINE ITIJKL(factor,norb,nsym,nos,nijkl,ijkl4,
     * ijkl2,rijkl,ijkl3)
C 
C     LECTURE DES INTEGRALES BIELECTRONIQUES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      PARAMETER(MBUF=3255,NBUF=4096)
c     character*2000 namel
      common /hst/ host
      character*8 host
      dimension ijkl4(*)
      integer*2 ijkl2(*)
      character*3 ijkl3(*)
      character*4 cara
      character*1 pivec(4),pivec2(4)
      equivalence (itamp,pivec2(1))
      equivalence (jtamp,pivec(1))
      equivalence (idummy,cara)
      real*4 rijkl(*)
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      idummy=0
      WRITE (6,1000)
 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      call getenv('HOSTTYPE',host)
      write (6,'(a,a)') '     HOSTTYPE=',host
      if (.not.zmotra) write (6,*) '    INTEGRALI DA HONDO + IJKLPI'
      if (zmotra) write (6,*) '    INTEGRALI DA MOLCAS + MOTRA'
c      do 500 k=1,ni2
c 500     ijkl2(k)=0
      do 501 k=1,ni4
 501     ijkl4(k)=0
c      do k=1,ni3
c         ijkl3(k)=cara(2:4)
c      enddo
      khan=0
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
      call filesplit('READ',50,8,ichange,idum2)
      READ (50) FACDUM
      WRITE (6,10) FACTOR
 10   FORMAT (5X,'INTEGRAL PREFACTOR =',D15.2)
      NNT=NUM(NORB)+NORB
      NTT=NDEB(NSYM)+(NOS(NSYM)*(NOS(NSYM)+1))/2
      M=0
      N=1
      call filesplit('READ',50,12*mbuf,ichange,idum2)
      READ(50)XX,IX
      DO 820 I=1,NORB
      II=NUM(I)
      DO 820 J=1,I
      IJ=II+J
      IJS=INDIC(IJ)
      IJP=JNDIC(IJ)
      DO 830 K=1,I
      LMAX=K
      IF(K.EQ.I)LMAX=J
      KK=NUM(K)
      DO 830 L=1,LMAX
      KL=KK+L
      KLS=INDIC(KL)
      M=M+1
      IF(KLS.NE.IJS) GO TO 830
      KLP=JNDIC(KL)
  831 LAB=IX(N)
      IF(LAB.EQ.0)GOTO 860
      IF(LAB.GT.M)GOTO 830
      N1=MAX0(IJP,KLP)
      N2=MIN0(IJP,KLP)
      NA=NUM(N1)+N2+NDEB(KLS)
      A8=XX(N)*factor
      CALL HANDY(NA,A8,KHAN,ZHANDY,zthree,ijkl4,ijkl3,ijkl2,rijkl)
      N=N+1
      IF(N.LE.MBUF)GOTO 830
      call filesplit('READ',50,12*mbuf,ichange,idum2)
      READ(50)XX,IX
      N=1
  830 CONTINUE
  820 CONTINUE
  860 RETURN
      return
      END
C***********************************************************************
      SUBROUTINE HANDY(I,A8,khan,zhandy,zthree,ijkl4,ijkl3,ijkl2,rijkl)
C 
C     STOCKAGE DES INTEGRALES A LA HANDY (N.C.HANDY, CHEM.PHYS.LETT.
C     74, 280 (1980))
C 
      integer*2 ijkl2(*)
      real*4 rijkl(*)
      character*3 ijkl3(*)
      dimension ijkl4(*)
      logical*1 zhandy,zthree
      character*4 cara
      equivalence (idummy,cara)
      REAL*8 A8
      if (zhandy) then
         A=A8
         if(.not.zthree)then
            IF (ABS(A).GT.20000.) GOTO 10
            IJKL2(I)=A
            RETURN
 10         KHAN=KHAN+1
            IF (KHAN.GT.25534) then
 20            WRITE (6,30)
 30       FORMAT (//5X,'TOO MANY INTEGRALS > 20000, DECREASE FACTOR'//)
               STOP 12
            endif
            RIJKL(KHAN)=A
            IF (KHAN.GT.12767) GOTO 15
            IJKL2(I)=KHAN+20000
            RETURN
 15         IJKL2(I)=-KHAN-7233
            RETURN
         else
            idummy=a
            if(idummy.gt.8388607)then
               print*,' Un integrale troppo grande ',idummy
               print*,' Il massimo si ha con 8388607'
               print*,' Riduci FACTOR'
               stop 'riduci factor'
            endif
            ijkl3(i)=cara(2:4)
         endif
      else
C
C     Aumento di 1 il max intero 2147483647 per sicurezza
c
         if (abs(a8).gt.2147483648.d0) then
            write (6,40)
            STOP 12
         endif
 40      format (//5X,'INTEGRAL OVERFLOW, DECREASE FACTOR'//)
         ijkl4(i)=nint(a8)
         return
      endif
      END
C***********************************************************************
      FUNCTION AI(I,J,K,L)
C
C     REPERAGE DE L'INTEGRALE BIELECTRONIQUE IJKL DANS LE TABLEAU IJKL2
C
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM1,LOCM2,LOCM3,LOCM4,LOCM(6)
      integer pick3
      common /big/ x,ijkl2,ijkl4,ijkl3,rijkl
      pointer x
      dimension x(:)
      pointer ijkl2,ijkl3,ijkl4
      pointer rijkl
      dimension ijkl2(:),ijkl4(:),ijkl3(:),rijkl(:)
      integer*2 ijkl2
      character*3,ijkl3
      real*4 rijkl
      IF (I.GE.J) GO TO 10
      IJ=NUM(J)+I
      GO TO 15
   10 IJ=NUM(I)+J
   15 IF (K.GE.L) GO TO 20
      KL=NUM(L)+K
      GO TO 25
   20 KL=NUM(K)+L
   25 IJS=INDIC(IJ)
      KLS=INDIC(KL)
      IF (IJS.LE.KLS) GO TO 30
      IJKLS=NUM(IJS)+KLS
      GO TO 35
   30 IJKLS=NUM(KLS)+IJS
      NIJ=IJ
      IJ=KL
      KL=NIJ
      IJS=KLS
      KLS=INDIC(KL)
   35 CONTINUE
      IF (NAD(IJKLS)) 40,300,45
   40 ZSIG=.TRUE.
      NADI=-NAD(IJKLS)
      GO TO 50
   45 ZSIG=.FALSE.
      NADI=NAD(IJKLS)
   50 KTYP=KT(IJKLS)
      GO TO (60,70,75,80,85,90,95,100,105,62),KTYP
   62 NIJ=LNDIC(IJ)
      NKL=LNDIC(KL)
      GO TO 63
   60 NIJ=JNDIC(IJ)
      NKL=JNDIC(KL)
   63 IF (NIJ.LE.NKL) GO TO 65
      IJKL=NUM(NIJ)+NKL
      GO TO 200
   65 IJKL=NUM(NKL)+NIJ
      GO TO 200
   70 NIJ=JNDIC(IJ)
      NKL=JNDIC(KL)
      IJKL=NBO(KLS)*(NIJ-1)+NKL
      GO TO 200
   75 NIJ=LNDIC(IJ)
      NKL=JNDIC(KL)
      IJKL=NBO(KLS)*(NIJ-1)+NKL
      GO TO 200
   80 NIJ=JNDIC(IJ)
      NKL=LNDIC(KL)
      IJKL=NBO(KLS)*(NIJ-1)+NKL
      GO TO 200
   85 NIJ=LNDIC(IJ)
      NKL=LNDIC(KL)
      IJKL=NBO(KLS)*(NIJ-1)+NKL
      GO TO 200
   90 NIJ=JNDIC(IJ)
      NKL=JNDIC(KL)
      IJKL=NBO(IJS)*(NKL-1)+NIJ
      GO TO 200
   95 NIJ=LNDIC(IJ)
      NKL=JNDIC(KL)
      IJKL=NBO(IJS)*(NKL-1)+NIJ
      GO TO 200
  100 NIJ=JNDIC(IJ)
      NKL=LNDIC(KL)
      IJKL=NBO(IJS)*(NKL-1)+NIJ
      GO TO 200
  105 NIJ=LNDIC(IJ)
      NKL=LNDIC(KL)
      IJKL=NBO(IJS)*(NKL-1)+NIJ
  200 NNJK=NDEB(NADI)+IJKL
      AI=PICK4(nnjk,ijkl4)
      IF (ZSIG) AI=-AI
      RETURN
  300 AI=0.D0
      RETURN
      END
C***********************************************************************
      FUNCTION PICK2(NNJK,IJKL2,RIJKL)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      integer*2 ijkl2(1)
      real*4 rijkl(1)
      LAI=IJKL2(NNJK)
      IF (LAI.GT.20000) then
      PICK2=RIJKL(LAI-20000)
      ELSEIF (LAI.LT.-20000) then
      PICK2=RIJKL(-LAI-7233)
      ELSE
      PICK2=DBLE(LAI)
      endif
      return
      END
C***********************************************************************
      FUNCTION PICK4(NNJK,IJKL4)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      dimension ijkl4(1)
      PICK4=ijkl4(nnjk)
      return
      END
C***********************************************************************
      integer FUNCTION PICK3(NNJK,IJKL3)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      character*3 ijkl3
      dimension ijkl3(1)
      common/caratte/c1,c2
      character*1 c1,c2
      character*4 cara
c     integer*2 i2
      equivalence (idummy,cara)
      equivalence (z2,cara(2:2))
c c1=char(-1),c2=char(0)      
      cara(2:4)=ijkl3(nnjk)
      if (z2) then
         cara(1:1)=c2
      else
         cara(1:1)=c1
      endif
      pick3=idummy
      return
      END
c----------------------------------------------------
      subroutine zeroe(a,n)
      real*8 a(*)
      do i=1,n
         a(i)=0.d0
      enddo
      return
      end
C***********************************************************************
      SUBROUTINE ITIJKLt(factor,norb,nsym,nos,nijkl,ijkl4,
     * ijkl2,rijkl,ijkl3)
C 
C     LECTURE DES INTEGRALES BIELECTRONIQUES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      PARAMETER(NBUF=4096)
c     character*2000 namel
      common /hst/ host
      character*8 host
      dimension ijkl4(*)
      integer*2 ijkl2(*)
      character*3 ijkl3(*)
      character*4 cara
      character*1 pivec(4),pivec2(4)
      equivalence (itamp,pivec2(1))
      equivalence (jtamp,pivec(1))
      equivalence (idummy,cara)
      real*4 rijkl(*)
      COMMON /HND/ ZHANDY,ZTHREE,zmotra,zfac
      pointer num,indic,jndic,lndic
      dimension num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,NUM,NDEB(500),NBO(99),NAD(2000)
     $     ,KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XXm(NBUF),IXm(NBUF)
c_toul
c      character*1 chin1,chin2,chin3,chin4,chin5
c      dimension  chin1(nbuf),chin2(nbuf),chin3(nbuf),chin4(nbuf)
c     $     ,chin5(nbuf)
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
c_toul
      idummy=0
c      WRITE (6,1000)
c 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
c      call getenv('HOSTTYPE',host)
c      write (6,'(a,a)') '     HOSTTYPE=',host
c      write (6,*) '    INTEGRALI LETTI CON IL FORMALISMO DI TOLOSA'
c      do 500 k=1,ni2
c 500     ijkl2(k)=0
      do 501 k=1,ni4
 501     ijkl4(k)=0
c      do k=1,ni3
c         ijkl3(k)=cara(2:4)
c      enddo
      khan=0
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
      WRITE (6,10) FACTOR
   10 FORMAT (5X,'INTEGRAL PREFACTOR =',D15.2)
      NNT=NUM(NORB)+NORB
      NTT=NDEB(NSYM)+(NOS(NSYM)*(NOS(NSYM)+1))/2
      M=0
      N=1
c--renzo modif to read Daniel's code for norb >255
   5  continue
c      if(norb.le.255)then
c         READ(50)chin1,chin2,chin3,chin4,xxm
c      else
c         READ(50)chin1,chin2,chin3,chin4,xxm,chin5
c      endif
      call filesplit('READ',50,16*nbuf,ichange,idum2)
      read(50)ind,jnd,knd,lnd,xxm !renzo new
      do indn=1,nbuf
       i=ind(indn)
       j=jnd(indn)
       k=knd(indn)
       l=lnd(indn)
       if(k.eq.0)return
c--renzo modif end
      II=NUM(i)
      IJ=II+J
      KK=NUM(K)
      KL=KK+L
      if (ij.lt.kl) then
      itamp=ij
      ij=kl
      kl=itamp
      endif
c      write (6,*) 'ij,kl',ij,kl
      IJS=INDIC(IJ)
      IJP=JNDIC(IJ)
      KLS=INDIC(KL)
c      write (6,*) 'inic',(indic(icel),icel=1,10)
c      write (6,*) 'nperm',(nperm(icel),icel=1,7)
      IF(KLS.NE.IJS) then
      print * ,'Attenzione errore in ITIJKL'
      print * ,'KLS,IJS',KLS,IJS
      endif
      KLP=JNDIC(KL)
      N1=MAX0(IJP,KLP)
      N2=MIN0(IJP,KLP)
      NA=NUM(N1)+N2+NDEB(KLS)
      A8=XXm(indn)*factor
c      write (6,7) indn,xxm(indn),i,j,k,l,a8
      CALL HANDY(NA,A8,KHAN,ZHANDY,zthree,ijkl4,ijkl3,ijkl2,rijkl)
c     CALL HANDY(NA,A8,KHAN,ZHANDY)
      enddo
      goto 5
    7 format (i4,F12.8,4I4,F20.2)
   12 NTER=NORB*(NORB+1)/2
      nfind=4096-indn
   11 if (nfind.ge.nter) return
      call filesplit('READ',50,0,ichange,idum2)
      read (50)
      nfind=nfind+4096
      goto 11  
      END
