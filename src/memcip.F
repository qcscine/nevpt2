#ifndef INTEL
      SUBROUTINE MEMCIP(ACTION)
#else
      SUBROUTINE MEMCIP(ACTION,doit)
#endif      
C
C     DYNAMIC MEMORY ALLOCATION, PROGRAM CIPPI
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Y),LOGICAL*1(Z)
      LOGICAL NOREL
      CHARACTER*5 ACTION
      character*240 namel
      include 'parameter.cippi'
      PARAMETER (ID20=110,ID132=ID13*2)
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM(10)
#ifndef INTEL
      COMMON /BIG/ X(1)
#else
      common /big/ x,ijkl2,ijkl4,ijkl3,rijkl
      pointer x
      dimension x(:)
      pointer ijkl2,ijkl3,ijkl4
      pointer rijkl
      dimension ijkl2(:),ijkl4(:),ijkl3(:),rijkl(:)
      integer*2 ijkl2
      character*3,ijkl3
      real*4 rijkl
      logical doit
#endif
      COMMON /INT/ NI4,NI3,NI2,NR4
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel(id3),ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      integer*2 ietats
      COMMON /DOR/ NSAB,NREFMO,NRMO2,NRFM,ZIDEN
C--------- HARDWARE SWITCH ------------
*DEC     EXTERNAL LOC
c     external loc
*AIX     NOREL=.FALSE.
*DEC     NOREL=.TRUE.
*HPUX    NOREL=.FALSE.
*LINUX   NOREL=.FALSE.
c     NOREL=.true.
      NOREL=.false.
C--------------------------------------
C
C     Release memory, if possible
C
      if (action.eq.'RELSE') then
#ifndef INTEL 
      if (.not.norel) call memrel(locmem)
#else  
      norel=.false.
      deallocate (x)
      if(doit)then
      deallocate (ijkl2)
      deallocate (ijkl3)
      deallocate (ijkl4)
      deallocate (rijkl)
      endif
#endif      
      return
      endif
C
C     Compute the amount of memory needed (IWANT*8 bytes)
C     and compare with memory available (MEMSIZ*8 bytes)
C     If ACTION.eq.TOTAL (first call) the maximum amount
C     of memory requested for integrals and other items
C     is computed.  This amount can be increased by setting
C     IWANT in the input (namelist &MEMRY).
C     The following calls of MEMCIP compute the memory
C     requirements for single steps of the calculation.
C     In systems where the release of memory is not possible
C     (NOREL=.TRUE.), the first call determines once for all
C     the allocated memory (MEMSIZ) and the following calls
C     only check that IWANT<MEMSIZ.
C     In systems where the allocated memory can be released
C     (NOREL=.FALSE), every call of MEMCIP allocates a
C     different amount of memory.
C
      iwant=0
      if (action.eq.'TOTAL') then
      memsiz=0
      iwant=10000
      rewind 3
      call cerca(3,'MEMRY',itrovo)
      if (itrovo.ne.0) then
      call lecnam('&MEMRY',namel,3)
      call nameli(namel,IWANT,'IWANT',4)
      call namell(namel,NOREL,'NOREL',4)
      call namell(namel,ZPRT,'ZPRT',1)
      endif
      endif
      locm(1)=0
C
C     ... for REIJKL + HMAT
C
      if (action.eq.'IJKLH'.or.action.eq.'TOTAL') then
#ifndef INTEL         
      locm(2)=locm(1)+nwrd(ni4*4)
      locm(3)=locm(2)+nwrd(ni3*3)
      locm(4)=locm(3)+nwrd(ni2*2)
      locm(5)=locm(4)+nwrd(nr4*4)
#else      
      if(doit)then
         allocate(ijkl2(ni2))
         allocate(ijkl3(ni3))
         allocate(ijkl4(ni4))
         allocate(rijkl(nr4))
      endif
      locm(2)=locm(1)+1
      locm(3)=locm(2)
      locm(4)=locm(3)
      locm(5)=locm(4)
#endif
      locm(6)=locm(5)+nwrd(5000*2)
      locm(7)=locm(6)+5000
      locm(8)=locm(7)+id4
      locm(9)=locm(8)+nwrd(id4)
      nwtot  =locm(9)+nwrd(id4)
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     ... for DAVID
C
      if (action.eq.'DAVID') then
      ndb=maxm*ncf
      ndab=ndb
      if (ndab.lt.id4) ndab=id4
      if (zvec.and.projf.ne.' '.and.ndab.lt.id8/2) ndab=id8/2+1
      locm(2)=locm(1)+ndb
      nwtot  =locm(2)+ndab
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     ... for ATOM
C
      if (action.eq.'ATOML'.or.action.eq.'TOTAL') then
      nwtot  =locm(1)+id13*id1
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     ... for DIABAT
C
      if (action.eq.'DIABA') then
      locm(2)=locm(1)+nrefmo*id13
      locm(3)=locm(2)+nsab*nsab
      locm(4)=locm(3)+nsab
      locm(5)=locm(4)+nsab
      locm(6)=locm(5)+nwrd(nsab*4)
      locm(7)=locm(6)+id13*nrfm
      locm(8)=locm(7)+nwrd(id20*nrmo2)
      nwtot  =locm(8)+nwrd(nrefmo*4)
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     ... for PREPER
C
      if (action.eq.'PREPR'.or.action.eq.'TOTAL') then
      locm(2)=locm(1)+nwrd(id5*20)
      locm(3)=locm(2)+nwrd(id6*20*2)
      nwtot  =locm(3)+nwrd(id5*10*2)
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     ... for MKF04
C
      if (action.eq.'MKF04'.or.action.eq.'TOTAL') then
      locm(2)=locm(1)+nwrd(id4*4)
      nwtot  =locm(2)+nwrd(id4*2)
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     ... for REIJKL (perturbation step)
C
      if (action.eq.'IJKLP'.or.action.eq.'TOTAL') then
#ifndef INTEL         
      locm(2)=locm(1)+nwrd(ni4*4)
      locm(3)=locm(2)+nwrd(ni3*3)
      locm(4)=locm(3)+nwrd(ni2*2)
      nwtot  =locm(4)+nwrd(nr4*4)
#else
      if(doit)then
      allocate(ijkl2(ni2))
      allocate(ijkl3(ni3))
      allocate(ijkl4(ni4))
      allocate(rijkl(nr4))
      endif
      locm(2)=locm(1)+1
      locm(3)=locm(2)
      locm(4)=locm(2)
      nwtot  =locm(4)
#endif      
      if (iwant.lt.nwtot) iwant=nwtot
      endif
C
C     First call, systems where the release of allocated memory
C     is not allowed,         or
C     Following calls, systems where the release of allocated memory
C     is allowed:         allocate IWANT*8 bytes of memory
C
C     First call, systems where the release of allocated memory
C     is allowed:         no action
C
C     Following calls, systems where the release of allocated memory
C     is not allowed:     compare IWANT and MEMSIZ
C
      if ((norel.and.action.eq.'TOTAL').or.
     *    (.not.norel.and.action.ne.'TOTAL')) then
#ifndef INTEL         
      LOCMEM = MEMGET(IWANT)
         IF (LOCMEM.EQ.0) THEN
         WRITE(6,50) iwant
         STOP 12
      ENDIF
#else
      allocate(x(iwant))
#endif      
      MEMSIZ = IWANT
      elseif (.not.norel.and.action.eq.'TOTAL') then
	 if (zprt) write (6,60) iwant
      elseif (norel.and.action.ne.'TOTAL') then
	 if (iwant.gt.memsiz) then
	 write (6,70) iwant,memsiz
	 stop 12
	 endif
      endif
   50 format (//'  *********************************************'/
     *          '  SYSTEM HAS NOT ENOUGH MEMORY AVAILABLE'/
     *          '  REQUIRED MEMORY: IWANT =',i10,' *8 BYTES'/
     *          '  *********************************************'/)
   60 format (//'  PRELIMINARY EVALUATION OF DYNAMIC MEMORY REQUIRED'/
     *          '  FOR THIS RUN: IWANT =',i10,' *8 BYTES'/)
   70 format (//'  *********************************************'/
     *          '  NOT ENOUGH MEMORY ALLOCATED, INCREASE IWANT'/
     *          '  REQUIRED  MEMORY: ',i10,' *8 BYTES'/
     *          '  ALLOCATED MEMORY: ',i10,' *8 BYTES'/
     *          '  *********************************************'/)
C
C         COMPUTE THE OFFSET -LOFFS- FROM X(1) TO THE BEGINNING
C         OF THE DYNAMIC POOL (E.G. BYTE TO WORD CONVERSION).
C         -LOFFS- WILL BE ZERO FOR STATIC IMPLEMENTATIONS.
C
C--------- HARDWARE SWITCH ------------
*AIX    LOCX = LADDRS(X)
*DEC    LOCX = LOC(X)
*HPUX   LOCX = LOC(X)
*LINUX  LOCX = LOC(X)
#ifndef INTEL      
      LOCX = laddrs(x)
#endif      
c     LOCX = loc(x)
C--------------------------------------
#ifndef INTEL      
      LOFFS = LOCMEM - LOCX
      LOFFS = (LOFFS+7)/8 + 1
#else
      loffs=1
#endif      
      do i=1,10
      locm(i)=locm(i)+loffs
      enddo
      RETURN
      END
C***********************************************************************
      SUBROUTINE MEMDIG
C
C     DYNAMIC MEMORY ALLOCATION, PROGRAMS CIPDIAGR AND CIPNESBET
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Y),LOGICAL*1(Z)
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM(10)
#ifndef INTEL
      COMMON /BIG/ X(1)
#else
#define CIP
#ifdef NES      
      common /big/ ijkl2,ijkl4,ijkl3,rijkl
      pointer ijkl2,ijkl3,ijkl4
      pointer rijkl
      dimension ijkl2(:),ijkl3(:),ijkl4(:),rijkl(:)
      integer*2 ijkl2
      real*4 rijkl
      character*3 ijkl3
      print*,'You shouln''t be calling memdig'
      stop 'you should not call memdig'
#endif
#ifdef CIP    
      common /big/ px
      pointer x
      dimension x(:)
#endif      
#endif
      COMMON /INT/ NI4,NI3,NI2,NR4
C--------- HARDWARE SWITCH ------------
*DEC    EXTERNAL LOC
c     external loc
C--------------------------------------
C
C     Compute the amount of memory needed (IWANT*8 bytes)
C     to store integrals.
C
      locm(1)=0
      locm(2)=locm(1)+nwrd(ni4*4)
      locm(3)=locm(2)+nwrd(ni3*3)
      locm(4)=locm(3)+nwrd(ni2*2)
      iwant  =locm(4)+nwrd(nr4*4)
C
C     Allocate memory
C
#ifndef INTEL      
      LOCMEM = MEMGET(IWANT)
         IF (LOCMEM.EQ.0) THEN
         WRITE(6,50) iwant
         STOP 12
      ENDIF
#else
      allocate(x(iwant))
#endif      
      MEMSIZ = IWANT
   50 format (//'  *********************************************'/
     *          '  SYSTEM HAS NOT ENOUGH MEMORY AVAILABLE'/
     *          '  REQUIRED MEMORY: IWANT =',i10,' *8 BYTES'/
     *          '  *********************************************'/)
C
C         COMPUTE THE OFFSET -LOFFS- FROM X(1) TO THE BEGINNING
C         OF THE DYNAMIC POOL (E.G. BYTE TO WORD CONVERSION).
C         -LOFFS- WILL BE ZERO FOR STATIC IMPLEMENTATIONS.
C
C--------- HARDWARE SWITCH ------------
*AIX    LOCX = LADDRS(X)
*DEC    LOCX = LOC(X)
*HPUX   LOCX = LOC(X)
*LINUX  LOCX = LOC(X)
#ifndef INTEL      
      LOCX = laddrs(x)
#endif      
c      LOCX = loc(x)
C--------------------------------------
#ifndef INTEL      
      LOFFS = LOCMEM - LOCX
      LOFFS = (LOFFS+7)/8 + 1
#else
      loffs=1
#endif      
      do i=1,10
      locm(i)=locm(i)+loffs
      enddo
      RETURN
      END
C***********************************************************************
      SUBROUTINE MEMROS
C
C     DYNAMIC MEMORY ALLOCATION, PROGRAM ROSPI
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Y),LOGICAL*1(Z)
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM(10)
#ifndef INTEL
      COMMON /BIG/ X(1)
#else
      common /big/ x
      pointer x
      dimension x(:)
#endif
      COMMON /INT/ NI4,NI3,NI2,NR4
C--------- HARDWARE SWITCH ------------
*DEC    EXTERNAL LOC
c     external loc
C--------------------------------------
C
C     Compute the amount of memory needed (IWANT*8 bytes)
C     to store integrals.
C
      locm(1)=0
      locm(2)=locm(1)+nwrd(ni4*4)
      locm(3)=locm(2)+nwrd(ni3*3)
      locm(4)=locm(3)+nwrd(ni2*2)
      iwant  =locm(4)+nwrd(nr4*4)
C
C     Allocate memory
C
#ifndef INTEL      
      LOCMEM = MEMGET(IWANT)
      IF (LOCMEM.EQ.0) THEN
         WRITE(6,50) iwant
         STOP 12
      ENDIF
#else   
      allocate(x(iwant))
#endif      
      MEMSIZ = IWANT
   50 format (//'  *********************************************'/
     *          '  SYSTEM HAS NOT ENOUGH MEMORY AVAILABLE'/
     *          '  REQUIRED MEMORY: IWANT =',i10,' *8 BYTES'/
     *          '  *********************************************'/)
C
C         COMPUTE THE OFFSET -LOFFS- FROM X(1) TO THE BEGINNING
C         OF THE DYNAMIC POOL (E.G. BYTE TO WORD CONVERSION).
C         -LOFFS- WILL BE ZERO FOR STATIC IMPLEMENTATIONS.
C
C--------- HARDWARE SWITCH ------------
*AIX    LOCX = LADDRS(X)
*DEC    LOCX = LOC(X)
*HPUX   LOCX = LOC(X)
*LINUX  LOCX = LOC(X)
#ifndef INTEL      
      LOCX = laddrs(x)
#endif      
c     LOCX = loc(x)
C--------------------------------------
#ifndef INTEL      
      LOFFS = LOCMEM - LOCX
      LOFFS = (LOFFS+7)/8 + 1
#else
      loffs=1
#endif      
      do i=1,10
      locm(i)=locm(i)+loffs
      enddo
      RETURN
      END
C***********************************************************************
      FUNCTION NWRD(NBYTE)
      NWRD=(NBYTE+7)/8
      RETURN
      END
