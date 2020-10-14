      implicit real*8 (a-h,o-y),logical*1 (z)
      allocatable coopipa(:,:),coopeaa(:,:)
      allocatable dal(:,:),koopaa(:,:,:,:),koopeaa(:,:,:,:)
      allocatable koopbb(:,:,:,:),koopab(:,:)
      allocatable daaa(:,:,:,:,:,:)
      allocatable taa(:,:,:,:)
      allocatable amat(:,:,:,:,:,:),bmat(:,:,:,:),cmat(:,:,:,:),dmat(:,:
     $     )
      allocatable atmat(:,:,:,:,:,:),btmat(:,:,:,:),ctmat(:,:,:,:)
     $     ,dtmat(:,:)
      allocatable ro4(:,:,:,:,:)
      allocatable atwo(:,:,:,:)
      real*8 koopaa,koopbb,koopab,koopeaa,koopebb,koopeab
c      dimension dal(idim,idim),dbe(idim,idim)
c      pointer ne,nd,ispin,iorb,itsym
      allocatable ne(:),nd(:),ispin(:),iorb(:),itsym(:)
      allocatable c(:)
c      dimension c(:),ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON/COE/C
      COMMON /CIP/ 
     *ITS(20,20),IGELS(2),NCFG,
     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
c      pointer trou,part
      allocatable trou(:),part(:)
c      COMMON /BUPA/TROU,PART,NREQ
      common /actspace/mact,nact2,nact3,nact4
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     * ,normind(18,18,18,18)
      real*4 tarray(2)
      external index
      logical z
c      pointer f,aj,ak,num,indic,jndic,lndic
      allocatable f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * NDEB(500),NBO(99),
     * NAD(2000),KT(2000)
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      character*80 file04,file50,file25,file31,file26
      allocatable iocc(:,:)
      integer*1 iocc
      integer*2 icomp,iconf,isegno
      common /toul/ ztoul
      allocatable icomp(:),iconf(:),isegno(:),zcapos(:)
      NAMELIST /LEGGI/ NCORE,NACT,THR,NSPIN,file04,file50,file26,ztoul,
     * file25,file31,zro4,ZORDER,zverbose
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c      read '(a)',file04
c      read '(a)',file50
c      read '(a)',file25
c      read '(a)',file31
      write (6,*)                      
      write (6,*)                      
      print '(/,a)',repeat('*',80)
      print '(/,a)',repeat(' KOOPRO ',10)
      print '(/,a)',repeat('*',80)
      write (6,*)                      
      write (6,*)                      
      ncore=0
      nact=0
      thr=1.d-10
      nspin=1   ! nspin is 2*S+1
      zro4=.false. !doesn't write ro4 by default
      ztoul=.false.
      file04=' '
      file50=' '
      file25=' '
      file31=' '
      file26=' '
      ZORDER=.true.
      zverbose=.false.
      READ (5,nml=LEGGI)
      if(.not.zorder) then
      write (6,*) 'The determinants are not correctly ordered:'
      write (6,*) 'each determinant is considered explicitly.'
      endif
      if (zverbose) print '(a,l1)',' ztoul is ',ztoul
      if(file04.eq.' '.or.file50.eq.' '.or.file25.eq.' '.or.file31.eq
     $     .' ')then
         print*,' file04, file50, file25, file31 are mandatory input'
         stop 'missing files'
      endif
      if(ztoul.and.file26.eq.' ')then
         print*,' file26 is mandatory when ztoul eq true'
         stop 'missing file26'
      endif
      open(4,file=file04,form='UNFORMATTED',status='OLD')
      open(50,file=file50,form='UNFORMATTED',status='OLD')
      maxlong=200000000
      call filesplit('OPEN',50,maxlong,idum,idum2)
      open(25,file=file25,form='UNFORMATTED',status='OLD')
      open(31,file=file31,form='UNFORMATTED',status='UNKNOWN')
      if(ztoul) open(26,file=file26,form='UNFORMATTED',status='OLD')
cele      READ (4) NORB,NOCB,NOCA,IDUM,NCF,II
      READ (4) NORB,NOCA,NOCB,METAT,NCF,II,I
cele      rewind 4
      allocate(ne(ncf+1000))
      allocate(nd(ncf+1000))
      allocate(icomp(ncf+1000))
      allocate(iconf(ncf+1000))
      allocate(isegno(ncf+1000))
      allocate(zcapos(ncf+1000))
      allocate(trou(ii+1000))
      allocate(part(ii+1000))
      iditrpa=ii
cele      READ (4) NORB,NOCA,NOCB,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),
cele     *     (TROU(J),PART(J),J=1,II)
      allocate(c(i+1000))
      allocate(itsym(norb+norb+1))
      allocate(ispin(norb+norb+1))
      allocate(iorb(norb+norb+1))
      allocate(f(norb*(norb+1)/2))
      allocate(aj(norb*(norb+1)/2))
      allocate(ak(norb*(norb+1)/2))
      allocate(num(norb*(norb+1)/2))
      allocate(indic(norb*(norb+1)/2))
      allocate(jndic(norb*(norb+1)/2))
      allocate(lndic(norb*(norb+1)/2))
      write (6,*) 
      write (6,*) 'FILE04 informations:'
      write (6,*) 'NORB = ',norb
      write (6,*) 'NOCA = ',noca
      write (6,*) 'NOCB = ',nocb
      write (6,*) 'METAT= ',metat
      write (6,*) 'NCF  = ',ncf
      write (6,*) 'II   = ',ii 
      write (6,*) 'I    = ',i  
      call flush(6)
C      rewind 4
C      READ (4) NORB,NOCA,NOCB,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),
C     *     (TROU(J),PART(J),J=1,II),(C(J),J=1,I),
C     *     (E,J=1,METAT),(EMP,J=1,METAT),(EPS,J=1,NORB),
C     *     (EPS,J=1,NORB),NNROT,(dummy,J=1,NNROT),
C     *     MMETAT,(ZHEFF,J=1,MMETAT),zhandy,zthree,FACTOR
      rewind 4
      READ (4) NORB,NOCA,NOCB,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),
     *     (TROU(J),PART(J),J=1,II),(C(J),J=1,I),
     *     (E,J=1,METAT),(EMP,J=1,METAT),(EPS,J=1,NORB),
     *     (EPS,J=1,NORB),NNROT,(dummy,J=1,NNROT),
     *     MMETAT,(ZHEFF,J=1,MMETAT),zhandy,zthree,FACTOR,escf
c--renzo debug begin
c      print*,'determinanti appena letti'
c--renzo debug end
c      call printtp(ne,nd,trou,part,ncf,norb,c)

      zion=.false.
      if (ne(1).ne.0.and.part(ne(1)).gt.2*norb) zion=.true.
c      zion=.false.
      if(zion) write(6,*) ' Odd number of electrons'                
c      print*,'leggo ncore e nact'
c      read*,ncore,nact,thr
      write (6,*) 
      write (6,'(a,I4)') ' Number of doubly occup. MOs      =',ncore
      write (6,'(a,I4)') ' Number of active MOs             =',nact
      write (6,'(a,D10.3)') ' Threshold to discard dets coeff. =',thr
      write (6,'(a,f15.8)')  ' Energy of the Fermi vacuum       = ',escf
      write (6,*)
      call flush(6)
      mact=nact
      nact2=nact**2
      nact3=nact**3
      nact4=nact**4
      nele=(noca-ncore)*2
      if(zion)nele=nele-1
      write (6,*) 'Number of active electrons  =',nele
      if(nele.eq.0)then
         print*,' There are no active electrons; stopping...'
         stop
      endif
      ndim2=nact*(nact+1)/2
      ndim3=ndim2*(ndim2+1)/2
      do i=1,nact
      lindi(i)=(i-1)*i*(i+1)*(i+2)/24
      lindj(i)=(i-1)*i*(i+1)/6
      lindk(i)=i*(i-1)/2
      enddo
      do i=1,nact
      do j=1,nact
      do k=1,nact
      do l=1,nact
      call getnorm(i,j,k,l,norm)
      normind(i,j,k,l)=norm
      enddo
      enddo
      enddo
      enddo
c      if(nact.gt.12)then
c         print*,'sorry, nact > 12 not yet implemented'
c         stop 'nact greater than 12'
c      endif
      nwords=nact*(nact+1)*(nact+2)*(nact+3)/24
      allocate(daaa(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(dal(1:nact,1:nact))
      allocate(taa(1:nact,1:nact,1:nact,1:nact))
      allocate(amat(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(atwo(1:nact,1:nact,1:nact,1:nact))
      allocate(ro4(1:nwords,1:nact,1:nact,1:nact,1:nact))
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
      call readint(atwo,factor,escf,nact,ncore,indic,jndic,lndic,num,f
     $     ,aj,ne,nd,ispin,iorb,itsym)
c      print*,'after readint: factor=',factor,' nact=',nact,' ncore='
c     $     ,ncore
      write (6,*) 
      write (6,*) 'Calling ro4 (four particle density matrix).'
      call flush(6)
      call zeroe(ro4,nwords*nact**4)
      call zeroe(daaa,nact**6)
      call zeroe(taa,nact**4)
      call zeroe(dal,nact**2)
      nocb=min(noca,nocb)
      noca=nocb
      call esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,nocb,norb
     $     ,nspin)
c      print*,'fatto'
      ncapos=0
      inext=1
      do i=1,ncf
         if(i.eq.inext)then
            zcapos(i)=.true.
            if (zorder) then
             inext=i+icomp(i)
            else
             inext=i+1
cele 16-03-2010
             icomp(i)=1
            endif
         else
            zcapos(i)=.false.
         endif
         if(zcapos(i))ncapos=ncapos+1
      enddo
      write (6,*) 
      write (6,'(a,i6,a)') ' There are ',ncapos,' "progenitors".'
      write (6,*) 
      call flush(6)
      allocate(iocc(1:norb,1:ncapos))
c      print*,'Allocata la memoria per iocc'
c      call flush(6)
      icapos=0
      do i=1,ncf
         if(zcapos(i))then
            icapos=icapos+1
            call giveocc(iocc(1,icapos),i,nd,ne,trou,part,nocb,norb)
         endif
      enddo
      if(nele.ge.4)then
         call rofour(c,ro4,ncore,nact,thr,icomp,iconf,iocc,zcapos,
     $        ncapos,trou,part,ne,nd,ispin,iorb,itsym,zorder,zverbose)
         call bro3(ro4,daaa,nact,nele)
         call bro2(daaa,taa,nact,nele)
         call bro1(taa,dal,nact,nele)
      elseif(nele.eq.3)then
         call ro3(c,daaa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,itsym
     $        )
         call bro2(daaa,taa,nact,nele)
         call bro1(taa,dal,nact,nele)
      elseif(nele.eq.2)then
         call ro2(c,taa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,itsym)
         call bro1(taa,dal,nact,nele)
      elseif(nele.eq.1)then
         print*,'just one active electron: are you kidding?'
         call ro1(c,dal,ncore,nact,thr,trou,part,ne,nd,ispin,iorb,itsym)
      endif
c      print * ,'one el den matr'
c      call matout(nact,nact,dal,nact)
      t=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t,' sec.'
      t2=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t2-t,' sec.'
      n=nact
      nstop=(n+n*n)/2
 9    format(6f12.6)
c 100  format(//20x,' Matrice alpha-alpha diagonale'//)
c 101  format(//20x,' Matrice alpha-beta diagonale'//)
c 910  format(i3,9f12.7,5(/3x,9f12.7))
c      write (6,*) 'Bulding the Koopmans matrices for double ionization' 
c---calcolo Koopmans esteso doppie ionizzazioni---------------------
      allocate(coopipa(1:nact,1:nact))
      allocate(coopeaa(1:nact,1:nact))
      allocate(koopaa(1:nact,1:nact,1:nact,1:nact))
      allocate(koopbb(1:nact,1:nact,1:nact,1:nact))
      allocate(koopab(1:nact,1:nact))
      allocate(koopeaa(1:nact,1:nact,1:nact,1:nact))
c      print*,'Chiamata di koop'
      call koopE(atwo,coopipa,coopeaa,taa,dal,factor,nact,ncore,f,aj)
c      print * ,'one el koo matr'
c      do i=1,nact
c      do j=1,nact
c      coopipa(i,j)=coopipa(i,j)/factor
c      enddo
c      enddo
c      call matout(nact,nact,coopipa,nact)
c      do i=1,nact
c      do j=1,nact
c      coopipa(i,j)=coopipa(i,j)*factor
c      enddo
c      enddo
      call wr1(dal,coopeaa,coopipa,nact)
      call wr2(koopaa,taa,dal,nact)
c      print*,'done with wr1 and wr2'
c      call flush(6)
c      write(31) nact
c      write(31)((dal(k,l),k=1,nact),l=1,nact)
c      write(31)((coopeaa(k,l),k=1,nact),l=1,nact)
c      write(31)((coopipa(k,l),k=1,nact),l=1,nact)
c      print*,'costruzione matrici ro2 di buca'
c      do i=1,nact
c         do j=1,nact
c            do k=1,nact
c               do l=1,nact
c                  adum=ro2t(i,j,k,l,taa,dal,nact)
c                  koopaa(i,j,k,l)=adum
c               enddo
c            enddo
c         enddo
c      enddo
c      print*,'scrittura matrici densita di buca su file 31'
c      do i=1,nact
c         do j=1,nact
c            write(31)((koopaa(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c      print*,'Chiamata di koopman2'
c      call flush(6)
c      call koopman2(atwo,daaa,dbbb,daab,dbba,taa,tbb,tab,dal,dbe,koopaa
c     $     ,koopbb,koopab,koopeaa,koopebb,koopeab,factor,nact,ncore)
      call koop2E(atwo,daaa,taa,dal,koopaa,koopeaa,factor,nact,ncore
     $     ,f,aj,itsym,its)
      if (zverbose) then
      t3=etime(tarray)
      print '(a,f10.2,a)','koopman2 has taken ',t3-t2,' seconds'
      call flush(6) 
      endif
      call wr3(koopeaa,taa,koopaa,nact)
c      print*,'done with wr3'
c      call flush(6)
c      do i=1,nact
c         do j=1,nact
c            write(31)((koopeaa(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c-------scrittura matrici densita a due particelle
c      do i=1,nact
c         do j=1,nact
c            write(31)((taa(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c-------scrittura matrici koopman doppia ionizzazione
c      do i=1,nact
c         do j=1,nact
c            write(31)((koopaa(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c      print '(a,f8.2,a)','elapsed time ',t3-t2,' sec.'
c      print*,'Chiamata di koopman0pE'
      call koopman0pE(atwo,daaa,taa,dal,koopaa,koopbb,koopab,factor,nact
     $     ,ncore,f,aj,itsym,its)
      write (6,*) 
      write (6,*) 'Checking  the elements of the F matrix'
      fmax=0.d0
      do i=1,nact
         do j=1,nact
            fprov=abs(koopab(i,j))
            if(fprov.gt.fmax)fmax=fprov
         enddo
      enddo
      write (6,*) 'Maximum value of F matrix is',fmax/factor
c      t4=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t4-t3,' sec.'
      call wr4(koopaa,koopbb,koopab,nact)
c      print*,'done with wr4'
c      call flush(6)
c-----scrittura matrici A e D-----------------------------------
c      do i=1,nact
c         do j=1,nact
c            write(31)((koopaa(i,j,k,l),k=1,nact),l=1,nact)
c            write(31)((koopbb(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c      write(31)((koopab(k,l),k=1,nact),l=1,nact)
c      t5=etime(tarray)
c      print '(a,f8.2,a)','elapsed time ',t5-t4,' sec.'
      allocate(atmat(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(bmat(1:nact,1:nact,1:nact,1:nact))
      allocate(cmat(1:nact,1:nact,1:nact,1:nact))
      allocate(dmat(1:nact,1:nact))
      allocate(btmat(1:nact,1:nact,1:nact,1:nact))
      allocate(ctmat(1:nact,1:nact,1:nact,1:nact))
      allocate(dtmat(1:nact,1:nact))
      call zeroe(amat,nact**6)
      call zeroe(atmat,nact**6)
c      print*,'finishing up construction of matrix AMAT'
      call bamat(amat,atmat,atwo,ro4,daaa,taa,dal,nact,ncore,f,aj,itsym
     $     ,its,norb)
      t6=etime(tarray)
c      print '(a,f10.2,a)','bamat has taken ',t6-t5,' seconds'
c      print*,'Building matrix BMAT'
      call bbmat(atwo,bmat,btmat,daaa,taa,dal,nact,ncore,f,aj)
      t7=etime(tarray)
c      print '(a,f10.2,a)','bbmat has taken ',t7-t6,' seconds'
c      print*,'Building matrix CMAT'
      call bcmat(atwo,cmat,ctmat,daaa,taa,dal,nact,ncore,f,aj)
      t8=etime(tarray)
c      print*,'finding differences in bmat and cmat'
      call fdiff(bmat,cmat,btmat,ctmat,factor,nact)
c      print*,'done with fdiff'
c      call flush(6)
c      do i =1,nact
c         do j=1,nact
c            do ii=1,nact
c               do jj=1,nact
c                  val=bmat(i,j,ii,jj)/factor
c                  val2=cmat(jj,i,j,ii)/factor
c                  if(abs(val-val2).gt.1.d-8)then
c                     print '(a,4i3,a,f12.8,a,f12.8)','indici',i,j,ii,jj,
c     $                    'bmat=',val,' cmat=',val2
c                  endif
c               enddo
c            enddo
c         enddo
c      enddo
c      print*,'finding differences in btmat and ctmat'
c      do i =1,nact
c         do j=1,nact
c            do ii=1,nact
c               do jj=1,nact
c                  val=btmat(i,j,ii,jj)/factor
c                  val2=ctmat(jj,i,j,ii)/factor
c                  if(abs(val-val2).gt.1.d-8)then
c                     print '(a,4i3,a,f12.8,a,f12.8)','indici',i,j,ii,jj,
c     $                    'btmat=',val,' ctmat=',val2
c                  endif
c               enddo
c            enddo
c         enddo
c      enddo
c--------
c      print '(a,f10.2,a)','bcmat has taken ',t8-t7,' seconds'
c      print*,'Building matrix DMAT'
      call bdmat(atwo,dmat,dtmat,taa,dal,nact,ncore,f,aj)
      t9=etime(tarray)
c      print '(a,f10.2,a)','bdmat has taken ',t9-t8,' seconds'
c      print '(a,f8.2,a)','Total elapsed time ',t9,' sec.'
c      call flush(6)
      call wr5(daaa,amat,atmat,bmat,btmat,cmat,ctmat,dmat,dtmat,nact)
c      print*,'done with wr5'
c      call flush(6)
      if(zro4)call wr6(ro4,nact,nwords)
c-----scrittura matrice spinless a 3 part e matrice amat
c      do i=1,nact
c         do j=1,nact
c            do k=1,nact
c               do l=1,nact
c                  write(31)((daaa(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
c                  write(31)((amat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
c               enddo
c            enddo
c         enddo
c      enddo
c----scrittura matrici bmat e cmat
c      do i=1,nact
c         do j=1,nact
c            write(31)((bmat(i,j,k,l),k=1,nact),l=1,nact)
c            write(31)((cmat(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c----scrittura matrice dmat
c      write(31)((dmat(i,j),i=1,nact),j=1,nact)
c-----scrittura matrice atmat
c      do i=1,nact
c         do j=1,nact
c            do k=1,nact
c               do l=1,nact
c                  write(31)((atmat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
c               enddo
c            enddo
c         enddo
c      enddo
c----scrittura matrici btmat e ctmat
c      do i=1,nact
c         do j=1,nact
c            write(31)((btmat(i,j,k,l),k=1,nact),l=1,nact)
c            write(31)((ctmat(i,j,k,l),k=1,nact),l=1,nact)
c         enddo
c      enddo
c----scrittura matrice dtmat
c      write(31)((dtmat(i,j),i=1,nact),j=1,nact)
c----that's all folks
      write (6,*)
      write (6,*)
      end
c-------------------------------------------------------------------
      subroutine rofour(c,amat,ncore,nact,thr,icomp,iconf,iocc,zcapos,
     $     ncaptot,trou,part,ne,nd,ispin,iorb,itsym,zorder,zverbose)
      implicit real*8 (a-h,o-y),logical*1 (z)
c      include 'rotwo.par'
      dimension amat(nwords,nact,nact,nact,nact)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     *     ,normind(18,18,18,18)
c      pointer ne,nd,ispin,iorb,itsym
      dimension c(*),zcapos(*)
      integer*2 icomp(*),iconf(*)
      integer*1 iocc(norb,ncaptot),flip
      allocatable flip(:)
c      COMMON/COE/C
      COMMON /CIP/ 
     *ITS(20,20),IGELS(2),NCFG,
     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      dimension ne(*),nd(*),ispin(*),iorb(*),itsym(*)
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      integer*8 nhits,ijcal
c      pointer trou,part
c      dimension trou(:),part(:)
      dimension trou(*),part(*)
c      COMMON /BUPA/TROU,PART,NREQ
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      allocatable iocca(:),ioccb(:)
c      dimension iocca(id1),ioccb(id1)
      integer op,pp,a,b,op1,op2,op3,pp1,pp2,pp3
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * F,AJ,AK,NUM,NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer segno,segno0,segno1,segno2,segno3,segno4,segno5
      nocc=noca-ncore
      allocate(iocca(norb+1))
      allocate(ioccb(norb+1))
      allocate(flip(ncf+1000))
      do i=1,ncf
       flip(i)=0
      enddo
c      print*,'nocc,noca',nocc,noca
c      print*,'subr. ro3'
c      print*,'nact=',nact,' nocc=',nocc,' ncf=',ncf
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
      if (zorder) call esame(flip,ne,nd,trou,part,ispin,iorb,ncf)
c      print*,'spin flips'
c      print '(30i2)',(flip(i),i=1,ncf)
      do  mcap=1,ncaptot
         minit=mfin+1
         mfin=minit+icomp(minit)-1
         if (mod(mcap,1000).eq.0)print *,'  mcap=',mcap
c         if (mod(mcap,100).eq.0)print *,'nhits,nhitsd',nhits,nhitsd
         call flush(6)
         nec=ne(minit)
         ndm=nd(minit)
         nfin=0
         do  ncap=1,mcap      !do 2
            ninit=nfin+1
            nfin=ninit+icomp(ninit)-1
            nedif=abs(ne(minit)-ne(ninit))
            znocalc=nedif.gt.4
            if(znocalc)nhitsd=nhitsd+1
            if(.not.znocalc)znocalc=ndiff(iocc(ncore+1,mcap)
     $           ,iocc(ncore+1,ncap),nact).gt.8
            if(znocalc)nhits=nhits+1
            if(znocalc)goto 2
            do m=minit,mfin
               nec=ne(m)
               ndm=nd(m)
               call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc
     $              ,nec,ndm,trou,part,iorb,ispin) 
               if(abs(c(m)).lt.thr) goto 3
               if(flip(m).eq.2)goto 3
               if(mcap.eq.ncap)then
                  nbeg=minit
                  nend=m-1
c                  ci2=cmval*c(m)
c                  val=ci2
                  call case0
               else
                  nbeg=ninit
                  nend=nfin
               endif
               do n=nbeg,nend
                  if(abs(c(n)).lt.thr) goto 4
                  ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,4,zseg)
                  ijcal=ijcal+1
                  if(ij.eq.0)nbdif=5
                  if(nbdif.gt.4)goto 4
                  if(zseg)then
                     segno=-1
                  else
                     segno=1
                  endif
                  if(nbdif.eq.1)then
                     nsings=nsings+1
                     n1=iorb(nv(1))
                     ns1=ispin(nv(1))
                     n2=iorb(nw(1))
                     ns2=ispin(nw(1))
c     goto 10
                     call case10
                  elseif(nbdif.eq.2)then
                     ndoubs=ndoubs+1
                     call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorb(nv(1))
                     ns1=ispin(nv(1))
                     n2=iorb(nv(2))
                     ns2=ispin(nv(2))
                     n3=iorb(nw(1))
                     ns3=ispin(nw(1))
                     n4=iorb(nw(2))
                     ns4=ispin(nw(2))
c     goto 20
                     call case20
                  elseif(nbdif.eq.3)then
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
c     goto 30
                     call case30
                  else          !nbdif.eq.4
c---  renzo debug
                     nfours=nfours+1
c---- 
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
c                     goto 40
                  call case40
                  endif
 4             enddo
 3          enddo
 2       enddo
 1    enddo
c---renzo debug
      if (zverbose) then
      print*,'number of single    differences=',nsings
      print*,'number of double    differences=',ndoubs
      print*,'number of triple    differences=',ntrips
      print*,'number of quadruple differences=',nfours
      print*,'number of jdifgen calls ',ijcal
      endif
c---contributi diagonali
      return
c     una differenza
c 10         continue
      contains
      subroutine case0
c      print*,'case 0 m=',m
c---  contributo a daaaa
            if(flip(m).eq.1)then
               valcm=2.d0*c(m)
            else
               valcm=c(m)
            endif
            cij=valcm*c(m)
            val=cij
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ic=iocca(j)
                  do k=j+1,nocca
                     ie=iocca(k)
                     do l=k+1,nocca
                        ig=iocca(l)
                        call ordsame(ia,ic,ie,ig,ia,ic,ie,ig,val,amat
     $                       ,nact,m,m)
                     enddo
                  enddo
               enddo
            enddo
c---  contributo a dbbbb
            do i=1,noccb
               ia=ioccb(i)
               do j=i+1,noccb
                  ic=ioccb(j)
                  do k=j+1,noccb
                     ie=ioccb(k)
                     do l=k+1,noccb
                        ig=ioccb(l)
                        call ordsame(ia,ic,ie,ig,ia,ic,ie,ig,val,amat
     $                       ,nact,m,m)
                     enddo
                  enddo
               enddo
            enddo
c---  contributi a daaab,daaba,dabaa,dbaaa
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ic=iocca(j)
                  do k=j+1,nocca
                     ie=iocca(k)
                     do l=1,noccb
                        ig=ioccb(l)
                        call ord31(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact
     $                       ,m,m)
                     enddo
                  enddo
               enddo
            enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
            do i=1,noccb
               ia=ioccb(i)
               do j=i+1,noccb
                  ic=ioccb(j)
                  do k=j+1,noccb
                     ie=ioccb(k)
                     do l=1,nocca
                        ig=iocca(l)
                        call ord31(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact
     $                       ,m,m)
                     enddo
                  enddo
               enddo
            enddo
c---  contributo a daabb,dbbaa,dabab,dbaba,dabba,dbaab
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ic=iocca(j)
                  do k=1,noccb
                     ie=ioccb(k)
                     do l=k+1,noccb
                        ig=ioccb(l)
                        call ord22(ia,ic,ie,ig,ia,ic,ie,ig,val,amat,nact
     $                       ,m,m)
                     enddo
                  enddo
               enddo
            enddo
            return
            end subroutine case0
      subroutine case10
            n1=n1-ncore
            n2=n2-ncore
            if(flip(m).eq.1.and.n.ne.m-1)then
               valcm=2.d0*c(m)
            else
               valcm=c(m)
            endif
            cij=valcm*c(n)*segno
            val=cij
            if(ns1.eq.1.and.ns2.eq.1)then
c---a sinistra il piu' grande
               if(n1.ge.n2)then
                  nl=n1
                  nr=n2
               else
                  nl=n2
                  nr=n1
               endif
c---  contributo a daaaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ic=iocca(j)
                     do k=j+1,nocca
                        ie=iocca(k)
                        if(ia.ne.n1.and.ic.ne.n1.and.ie.ne.n1.and.
     $                       ia.ne.n2.and.ic.ne.n2.and.ie.ne.n2)then
c			write (6,*) 'calling ordsame'
                           call ordsame(ia,ic,ie,nl,ia,ic,ie,nr,val
     $                          ,amat,nact,m,n)
                        endif
                     enddo
                  enddo
               enddo
c---  contributo a dbbba
               do i=1,noccb
                  ia=ioccb(i)
                  do j=i+1,noccb
                     ic=ioccb(j)
                     do k=j+1,noccb
                        ie=ioccb(k)
c			write (6,*) 'calling ord31'
                           call ord31(ia,ic,ie,nl,ia,ic,ie,nr,val
     $                          ,amat,nact,m,n)
                     enddo
                  enddo
               enddo
c---  contributi a daaab,daaba,dabaa,dbaaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ic=iocca(j)
                     if(ia.ne.n1.and.ic.ne.n1.and.ia.ne.n2.and.ic.ne.n2
     $                    )then
                        do k=1,noccb
                           ie=ioccb(k)
c			write (6,*) 'calling ord31 aaab'
c			write (6,*) ia,ic,ie,nl,nr
                           call ord31(ia,ic,nl,ie,ia,ic,nr,ie,val,amat
     $                          ,nact,m,n)
                        enddo
                     endif
                  enddo
               enddo
c---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
               do i=1,nocca
                  ia=iocca(i)
                  do j=1,noccb
                     ic=ioccb(j)
                     do k=j+1,noccb
                        ie=ioccb(k)
                        if(ia.ne.n1.and.ia.ne.n2)then
c			write (6,*) 'calling ord22'
                           call ord22(ia,nl,ic,ie,ia,nr,ic,ie,val,amat
     $                          ,nact,m,n)
                        endif
                     enddo
                  enddo
               enddo
            elseif(ns1.eq.0.and.ns2.eq.0)then
c---a sinistra il piu' grande
               if(n1.ge.n2)then
                  nl=n1
                  nr=n2
               else
                  nl=n2
                  nr=n1
               endif
c---  contributo a dbbbb
               do i=1,noccb
                  ia=ioccb(i)
                  do j=i+1,noccb
                     ic=ioccb(j)
                     do k=j+1,noccb
                        ie=ioccb(k)
                        if(ia.ne.n1.and.ic.ne.n1.and.ie.ne.n1.and.
     $                       ia.ne.n2.and.ic.ne.n2.and.ie.ne.n2)then
c			write (6,*) 'calling ordsame'
                           call ordsame(ia,ic,ie,nl,ia,ic,ie,nr,val,amat
     $                          ,nact,m,n)
                        endif
                     enddo
                  enddo
               enddo
c---  contributo a daaab
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ic=iocca(j)
                     do k=j+1,nocca
                        ie=iocca(k)
c			write (6,*) 'calling ord31'
                           call ord31(ia,ic,ie,nl,ia,ic,ie,nr,val,amat
     $                          ,nact,m,n)
                     enddo
                  enddo
               enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
               do i=1,noccb
                  ia=ioccb(i)
                  do j=i+1,noccb
                     ic=ioccb(j)
                     if(ia.ne.n1.and.ic.ne.n1.and.ia.ne.n2.and.ic.ne.n2
     $                    )then
                        do k=1,nocca
                           ie=iocca(k)
c			write (6,*) 'calling ord31'
                           call ord31(ia,ic,nl,ie,ia,ic,nr,ie,val,amat
     $                          ,nact,m,n)
                        enddo
                     endif
                  enddo
               enddo
c---  contributi a dbbaa,daabb,dbaba,dabab,dbaab,dabba
               do i=1,noccb
                  ia=ioccb(i)
                  do j=1,nocca
                     ic=iocca(j)
                     do k=j+1,nocca
                        ie=iocca(k)
                        if(ia.ne.n1.and.ia.ne.n2)then
c			write (6,*) 'calling ord22'
                           call ord22(ic,ie,ia,nl,ic,ie,ia,nr,val,amat
     $                          ,nact,m,n)
                        endif
                     enddo
                  enddo
               enddo
            endif
c            goto 2
            return
            end subroutine case10
c     due differenze
c 20         continue
c            contains 
            subroutine case20
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            if(flip(m).eq.1.and.n.ne.m-1)then
               valcm=2.d0*c(m)
            else
               valcm=c(m)
            endif
            cij=valcm*c(n)*segno
            val=cij
c--- n1l,n2l a sinistra 
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
c---  contributo a daaaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ic=iocca(j)
                     if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4
     $                    .and.ic.ne.n1.and.ic.ne.n2.and.ic.ne.n3.and.ic
     $                    .ne.n4)then
                        call ordsame(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val
     $                       ,amat,nact,m,n)
                     endif
                  enddo
               enddo
c---  contributo a daaab,daaba,dabaa,dbaaa
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4
     $                 )then
                     do j=1,noccb
                        ic=ioccb(j)
                        call ord31(ia,n1l,n2l,ic,ia,n1r,n2r,ic,val,amat
     $                       ,nact,m,n)
                     enddo
                  endif
               enddo
c---  contributo a dbbaa,daabb,dabab,dbaba,dabba,dbaab
               do i=1,noccb
                  ia=ioccb(i)
                  do j=i+1,noccb
                     ic=ioccb(j)
                     call ord22(n1l,n2l,ia,ic,n1r,n2r,ia,ic,val,amat
     $                    ,nact,m,n)
                  enddo
               enddo
            elseif(ns1.eq.0.and.ns2.eq.0)then
c---  contributo a dbbbb
               do i=1,noccb
                  ia=ioccb(i)
                  do j=i+1,noccb
                     ic=ioccb(j)
                     if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4
     $                    .and.ic.ne.n1.and.ic.ne.n2.and.ic.ne.n3.and.ic
     $                    .ne.n4)then
                        call ordsame(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val
     $                       ,amat,nact,m,n)
                     endif
                  enddo
               enddo
c---  contributo a dbbba,dbbab,dbabb,dabbb
               do i=1,noccb
                  ia=ioccb(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4
     $                 )then
                     do j=1,nocca
                        ic=iocca(j)
                        call ord31(ia,n1l,n2l,ic,ia,n1r,n2r,ic,val,amat
     $                       ,nact,m,n)
                     enddo
                  endif
               enddo
c---  contributo a daabb,dbbaa,dbaba,dabab,dbaab,dabba
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ic=iocca(j)
                     call ord22(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val,amat
     $                    ,nact,m,n)
                  enddo
               enddo
            elseif(ns1.eq.1.and.ns2.eq.0)then
c---  contributi a daaab,daaba,dabaa,dbaaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ic=iocca(j)
                     if(ia.ne.n1.and.ia.ne.n3.and.ic.ne.n1.and.ic.ne.n3
     $                    )then
c---- daaab
                        call ord31(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val,amat
     $                       ,nact,m,n)
                     endif
                  enddo
               enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
               do i=1,noccb
                  ia=ioccb(i)
                  do j=i+1,noccb
                     ic=ioccb(j)
                     if(ia.ne.n2.and.ia.ne.n4.and.ic.ne.n2.and.ic.ne.n4
     $                    )then
c---- dbbba
                        call ord31(ia,ic,n2l,n1l,ia,ic,n2r,n1r,val,amat
     $                       ,nact,m,n)
                     endif
                  enddo
               enddo
c---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
               do i=1,nocca
                  ia=iocca(i)
                  do j=1,noccb
                     ic=ioccb(j)
                     if(ia.ne.n1.and.ia.ne.n3.and.ic.ne.n2.and.ic.ne.n4
     $                    )then
c---- daabb
                        call ord22(ia,n1l,ic,n2l,ia,n1r,ic,n2r,val,amat
     $                       ,nact,m,n)
                     endif
                  enddo
               enddo    
            endif
c            goto 2
            return
            end subroutine case20
c     tre differenze
c 30         continue
c            contains 
            subroutine case30
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            n5=n5-ncore
            n6=n6-ncore
            if(flip(m).eq.1.and.n.ne.m-1)then
               valcm=2.d0*c(m)
            else
               valcm=c(m)
            endif
            cij=valcm*c(n)*segno
            val=cij
c--- n1l,n2l,n3l a sinistra 
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
c     casi aaa, bbb, aab, bba
            if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
c---  contributo daaaa
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4.and
     $                 .ia.ne.n5.and.ia.ne.n6)then
                     call ordsame(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r,val,amat
     $                    ,nact,m,n)
                  endif
               enddo
c---  contributi a daaab,daaba,dabaa,dbaaa
               do i=1,noccb
                  ia=ioccb(i)
                  call ord31(n1l,n2l,n3l,ia,n1r,n2r,n3r,ia,val,amat,nact
     $                 ,m,n)
               enddo
            elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
c---  contributo dbbbb
               do i=1,noccb
                  ia=ioccb(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4.and
     $                 .ia.ne.n5.and.ia.ne.n6)then
                     call ordsame(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r,val,amat
     $                    ,nact,m,n)
                  endif
               enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
               do i=1,nocca
                  ia=iocca(i)
                  call ord31(n1l,n2l,n3l,ia,n1r,n2r,n3r,ia,val,amat,nact
     $                 ,m,n)
               enddo
            elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
c---  contributi a daaab,daaba,dabaa,dbaaa
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n4.and.ia.ne.n5
     $                 )then
c---- daaab
                     call ord31(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r,val,amat
     $                    ,nact,m,n)
                  endif
               enddo
c--- contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
               do i=1,noccb
                  ia=ioccb(i)
                  if(ia.ne.n3.and.ia.ne.n6)then
c---- daabb
                     call ord22(n1l,n2l,ia,n3l,n1r,n2r,ia,n3r,val,amat
     $                    ,nact,m,n)
                  endif
               enddo
            elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
c---  contributi a dbbba,dbbab,dbabb,dabbb
               do i=1,noccb
                  ia=ioccb(i)
                  if(ia.ne.n2.and.ia.ne.n3.and.ia.ne.n5.and.ia.ne.n6
     $                 )then
c---- dbbba
                     call ord31(ia,n2l,n3l,n1l,ia,n2r,n3r,n1r,val,amat
     $                    ,nact,m,n)
                  endif
               enddo
c--- contributi a dbbaa,daabb,dbaba,dabab,dbaab,dabba
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n4)then
c---- dbbaa
                     call ord22(n1l,ia,n2l,n3l,n1r,ia,n2r,n3r,val,amat
     $                    ,nact,m,n)
                  endif
               enddo
            endif
c            goto 2
            return
            end subroutine case30
c     quattro differenze
c 40         continue
c            contains 
            subroutine case40
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            n5=n5-ncore
            n6=n6-ncore
            n7=n7-ncore
            n8=n8-ncore
            if(flip(m).eq.1.and.n.ne.m-1)then
               valcm=2.d0*c(m)
            else
               valcm=c(m)
            endif
            cij=valcm*c(n)*segno
            val=cij
c--- n1l,n2l,n3l a sinistra 
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
c     contributo a daaaa
               call ordsame(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat
     $              ,nact,m,n)
            elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0.and.ns4.eq.0)then
c     contributo a dbbbb
               call ordsame(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat
     $              ,nact,m,n)
            elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1.and.ns4.eq.0)then
c     contributo a daaab
               call ord31(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat
     $              ,nact,m,n)
            elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0.and.ns4.eq.0)then
c     contributo a dbbba
               call ord31(n2l,n3l,n4l,n1l,n2r,n3r,n4r,n1r,val,amat
     $              ,nact,m,n)
            elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0.and.ns4.eq.0)then
c     contributo a daabb
               call ord22(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val,amat
     $              ,nact,m,n)
            endif
            return
            end subroutine case40
      end
c--------------------------------------------------------------------
      subroutine ro3(c,daaa,ncore,nact,thr,trou,part,ne,nd,ispin,iorb
     $     ,itsym)
      implicit real*8 (a-h,o-y),logical*1 (z)
c      include 'rotwo.par'
      allocatable dbbb(:,:,:,:,:,:),daab(:,:,:,:,:,:),dbba(:,:,:,:,:,:)
      dimension daaa(nact,nact,nact,nact,nact,nact)
c      pointer ne,nd,ispin,iorb,itsym
      dimension c(*)
c      COMMON/COE/C
      COMMON /CIP/ 
     *ITS(20,20),IGELS(2),NCFG,
     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*),itsym(*)
      integer segno
c      pointer trou,part
c      dimension trou(:),part(:)
      dimension trou(*),part(*)
c      COMMON /BUPA/TROU,PART,NREQ
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      allocatable iocca(:),ioccb(:)
c      dimension iocca(id1),ioccb(id1)
      allocate(iocca(norb))
      allocate(ioccb(norb))
      nocc=noca-ncore
c      print*,'nocc,noca',nocc,noca
c      print*,'subr. ro3'
c      print*,'nact=',nact,' nocc=',nocc,' ncf=',ncf
      allocate(dbbb(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(daab(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(dbba(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
c---azzeramento
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  do m=1,nact
                     do n=1,nact
                        daaa(i,j,k,l,m,n)=0.d0
                        dbbb(i,j,k,l,m,n)=0.d0
                        daab(i,j,k,l,m,n)=0.d0
                        dbba(i,j,k,l,m,n)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      do 1 m=1,ncf
         if (mod(m,100).eq.0)print *,'m=',m
         call flush(6)
         if(abs(c(m)).lt.thr)goto 1
         nec=ne(m)
         ndm=nd(m)
         call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec,ndm
     $        ,trou,part,iorb,ispin) 
         do 2 n=1,m
            if(abs(c(n)).lt.thr)goto 2
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
            ci2=c(m)*c(m)
c     costruzione di daaa
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do k=j+1,nocca
                     ka=iocca(k)
                     daaa(ia,ja,ka,ia,ja,ka)=daaa(ia,ja,ka,ia,ja,ka)
     $                    +ci2
                     call perm3(daaa,ia,ja,ka,ia,ja,ka,nact)
                  enddo
               enddo
            enddo
c     costruzione di dbbb            
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do k=j+1,noccb
                     kb=ioccb(k)
                     dbbb(ib,jb,kb,ib,jb,kb)=dbbb(ib,jb,kb,ib,jb,kb)
     $                    +ci2
                     call perm3(dbbb,ib,jb,kb,ib,jb,kb,nact)
                  enddo
               enddo
            enddo
c     costruzione di daab
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  do k=1,noccb
                     kb=ioccb(k)
                     daab(ia,ja,kb,ia,ja,kb)=daab(ia,ja,kb,ia,ja,kb)
     $                    +ci2
                     call perm2(daab,ia,ja,kb,ia,ja,kb,nact)
                  enddo
               enddo
            enddo
c     costruzione di dbba            
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  do k=1,nocca
                     ka=iocca(k)
                     dbba(ib,jb,ka,ib,jb,ka)=dbba(ib,jb,ka,ib,jb,ka)
     $                    +ci2
                     call perm2(dbba,ib,jb,ka,ib,jb,ka,nact)
                  enddo
               enddo
            enddo
            goto 2
c     una differenza
 10         continue
            n1=n1-ncore
            n2=n2-ncore
            cij=c(m)*c(n)*segno
            if(ns1.eq.1.and.ns2.eq.1)then
c---contributo a daaa
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ja=iocca(j)
                     if(ia.ne.n1.and.ja.ne.n1.and.ia.ne.n2.and.ja.ne
     $                    .n2)then
                        daaa(n1,ia,ja,n2,ia,ja)=daaa(n1,ia,ja,n2,ia
     $                       ,ja)+cij
                        daaa(n2,ia,ja,n1,ia,ja)=daaa(n2,ia,ja,n1,ia
     $                       ,ja)+cij
                        call perm3(daaa,n1,ia,ja,n2,ia,ja,nact)
                        call perm3(daaa,n2,ia,ja,n1,ia,ja,nact)
                     endif
                  enddo
               enddo
c----contributo a dbba
               do i=1,noccb
                  ib=ioccb(i)
                  do j=i+1,noccb
                     jb=ioccb(j)
                     dbba(ib,jb,n1,ib,jb,n2)=dbba(ib,jb,n1,ib,jb,n2)
     $                    +cij
                     call perm2(dbba,ib,jb,n1,ib,jb,n2,nact)
                  enddo
               enddo
c----contributo a daab
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2)then
                     do j=1,noccb
                        jb=ioccb(j)
                        daab(n1,ia,jb,n2,ia,jb)=daab(n1,ia,jb,n2,ia
     $                       ,jb)+cij
                        call perm2(daab,n1,ia,jb,n2,ia,jb,nact)
                     enddo
                  endif
               enddo
            elseif(ns1.eq.0.and.ns2.eq.0)then
c----contributo a dbbb
               do i=1,noccb
                  ib=ioccb(i)
                  do j=i+1,noccb
                     jb=ioccb(j)
                     if(ib.ne.n1.and.jb.ne.n1.and.ib.ne.n2.and.jb.ne
     $                    .n2)then
                        dbbb(n1,ib,jb,n2,ib,jb)=dbbb(n1,ib,jb,n2,ib
     $                       ,jb)+cij
                        dbbb(n2,ib,jb,n1,ib,jb)=dbbb(n2,ib,jb,n1,ib
     $                       ,jb)+cij
                        call perm3(dbbb,n1,ib,jb,n2,ib,jb,nact)
                        call perm3(dbbb,n2,ib,jb,n1,ib,jb,nact)
                     endif
                  enddo
               enddo
c---contrinuto a daab
               do i=1,nocca
                  ia=iocca(i)
                  do j=i+1,nocca
                     ja=iocca(j)
                     daab(ia,ja,n1,ia,ja,n2)=daab(ia,ja,n1,ia,ja,n2)
     $                    +cij
                     call perm2(daab,ia,ja,n1,ia,ja,n2,nact)
                  enddo
               enddo
c----contributo a dbba
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n1.and.ib.ne.n2)then
                     do j=1,nocca
                        ja=iocca(j)
                        dbba(n1,ib,ja,n2,ib,ja)=dbba(n1,ib,ja,n2,ib
     $                       ,ja)+cij
                        call perm2(dbba,n1,ib,ja,n2,ib,ja,nact)
                     enddo
                  endif
               enddo
            endif
            goto 2
c     due differenze
 20         continue
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            cij=c(m)*c(n)*segno
            if(ns1.eq.1.and.ns2.eq.1)then
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia.ne.n4
     $                 )then
                     daaa(n1,n2,ia,n3,n4,ia)=daaa(n1,n2,ia,n3,n4,ia)
     $                    +cij
                     daaa(n3,n4,ia,n1,n2,ia)=daaa(n3,n4,ia,n1,n2,ia)
     $                    +cij
                     call perm3(daaa,n1,n2,ia,n3,n4,ia,nact)
                     call perm3(daaa,n3,n4,ia,n1,n2,ia,nact)
                  endif
               enddo
               do i=1,noccb
                  ib=ioccb(i)
                  daab(n1,n2,ib,n3,n4,ib)=daab(n1,n2,ib,n3,n4,ib)+cij
                  call perm2(daab,n1,n2,ib,n3,n4,ib,nact)
               enddo

            elseif(ns1.eq.0.and.ns2.eq.0)then
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n1.and.ib.ne.n2.and.ib.ne.n3.and.ib.ne.n4
     $                 )then
                     dbbb(n1,n2,ib,n3,n4,ib)=dbbb(n1,n2,ib,n3,n4,ib)
     $                    +cij
                     dbbb(n3,n4,ib,n1,n2,ib)=dbbb(n3,n4,ib,n1,n2,ib)
     $                    +cij
                     call perm3(dbbb,n1,n2,ib,n3,n4,ib,nact)
                     call perm3(dbbb,n3,n4,ib,n1,n2,ib,nact)
                  endif
               enddo
               do i=1,nocca
                  ia=iocca(i)
                  dbba(n1,n2,ia,n3,n4,ia)=dbba(n1,n2,ia,n3,n4,ia)+cij
                  call perm2(dbba,n1,n2,ia,n3,n4,ia,nact)
               enddo
            elseif(ns1.eq.1.and.ns2.eq.0)then
               do i=1,nocca
                  ia=iocca(i)
                  if(ia.ne.n1.and.ia.ne.n3)daab(n1,ia,n2,n3,ia,n4)
     $                 =daab(n1,ia,n2,n3,ia,n4)+cij
                  call perm2(daab,n1,ia,n2,n3,ia,n4,nact)
               enddo
               do i=1,noccb
                  ib=ioccb(i)
                  if(ib.ne.n2.and.ib.ne.n4)dbba(n2,ib,n1,n4,ib,n3)
     $                 =dbba(n2,ib,n1,n4,ib,n3)+cij
                  call perm2(dbba,n2,ib,n1,n4,ib,n3,nact)
               enddo
            endif
            goto 2
c     tre differenze
 30         continue
            n1=n1-ncore
            n2=n2-ncore
            n3=n3-ncore
            n4=n4-ncore
            n5=n5-ncore
            n6=n6-ncore
            cij=c(m)*c(n)*segno
c     casi aaa, bbb, aab, bba
            if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
               daaa(n1,n2,n3,n4,n5,n6)=daaa(n1,n2,n3,n4,n5,n6)+cij
               daaa(n4,n5,n6,n1,n2,n3)=daaa(n4,n5,n6,n1,n2,n3)+cij
               call perm3(daaa,n1,n2,n3,n4,n5,n6,nact)
               call perm3(daaa,n4,n5,n6,n1,n2,n3,nact)
            elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
               dbbb(n1,n2,n3,n4,n5,n6)=dbbb(n1,n2,n3,n4,n5,n6)+cij
               dbbb(n4,n5,n6,n1,n2,n3)=dbbb(n4,n5,n6,n1,n2,n3)+cij
               call perm3(dbbb,n1,n2,n3,n4,n5,n6,nact)
               call perm3(dbbb,n4,n5,n6,n1,n2,n3,nact)
            elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
               daab(n1,n2,n3,n4,n5,n6)=daab(n1,n2,n3,n4,n5,n6)+cij
               call perm2(daab,n1,n2,n3,n4,n5,n6,nact)
            elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
               dbba(n2,n3,n1,n5,n6,n4)=dbba(n2,n3,n1,n5,n6,n4)+cij
               call perm2(dbba,n2,n3,n1,n5,n6,n4,nact)
            endif
 2       continue
 1    continue
c--costruzione matrice ro3 spinless
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do ip=1,nact
                  do jp=1,nact
                     do kp=1,nact
                        daaa(i,j,k,ip,jp,kp)=daaa(i,j,k,ip,jp,kp)+
     $                       dbbb(i,j,k,ip,jp,kp)+daab(i,j,k,ip,jp,kp)
     $                       +dbba(i,j,k,ip,jp,kp)+daab(i,k,j,ip,kp,jp)
     $                       +daab(j,k,i,jp,kp,ip)+dbba(i,k,j,ip,kp,jp)
     $                       +dbba(k,j,i,kp,jp,ip)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(dbbb)
      deallocate(daab)
      deallocate(dbba)
      return
      end
c-------------------------------------------------------
      subroutine ro2(c,d2alal,ncore,nact,thr,trou,part,ne,nd,ispin,iorb
     $     ,itsym)
      implicit real*8 (a-h,o-y),logical*1 (z)
c      include 'rotwo.par'
      dimension d2alal(nact,nact,nact,nact)
      allocatable d2bebe(:,:,:,:),d2albe(:,:,:,:)
c      pointer ne,nd,ispin,iorb,itsym
      dimension c(*)
c      COMMON/COE/C
      COMMON /CIP/ 
     *ITS(20,20),IGELS(2),NCFG,
     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*),itsym(*)
c      pointer trou,part
c      dimension trou(:),part(:)
      dimension trou(*),part(*)
c      COMMON /BUPA/TROU,PART,NREQ
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      logical z1,z2,z3,z4,z5
      integer segno
      allocatable iocca(:),ioccb(:)
c      dimension iocca(id1),ioccb(id1)
      external index
      parameter (pt5=0.5d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c common e dimensionamenti mettere dopo
      allocate(d2bebe(1:nact,1:nact,1:nact,1:nact))
      allocate(d2albe(1:nact,1:nact,1:nact,1:nact))
      allocate(iocca(norb))
      allocate(ioccb(norb))
      n=nact
      nocc=noca-ncore
c      print*,'nocc,noca',nocc,noca
      nbig=(n*n+n)/2
      nstop=nbig
      nbig=(nbig*nbig+nbig)/2
      do i=1,n
         do j=1,n
            do k=1,n
               do l=1,n
                  d2alal(i,j,k,l)=0.d0
                  d2bebe(i,j,k,l)=0.d0
                  d2albe(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo
c      print*,'subr. ro2'
c      print*,'n=',n,' nocc=',nocc,' ncf=',ncf
      do 1 m=1,ncf
         if (mod(m,100).eq.0)print *,'m=',m
         call flush(6)
         if(abs(c(m)).lt.thr)goto 1
         ndi=nd(m)
         nec=ne(m)
         call giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec
     $        ,ndi,trou,part,iorb,ispin)
         do 2 n=1,m
            if(abs(c(n)).lt.thr)goto 2
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
c                  print '(a,4i3,l2)','due diff: ',n1,n3,n2,n4,zseg
                  goto 20
               endif
c            goto(10,20),nbdif
 3          continue
c     nessuna differenza
            ci2=c(m)*c(m)
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  d2alal(ia,ja,ia,ja)=d2alal(ia,ja,ia,ja)+ci2
                  call permbis(d2alal,ia,ja,ia,ja,nact)
               enddo
            enddo
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  d2bebe(ib,jb,ib,jb)=d2bebe(ib,jb,ib,jb)+ci2
                  call permbis(d2bebe,ib,jb,ib,jb,nact)
               enddo
            enddo
            do i=1,nocca
               ia=iocca(i)
               do j=1,noccb
                  jb=ioccb(j)
                  d2albe(ia,jb,ia,jb)=d2albe(ia,jb,ia,jb)+ci2
               enddo
            enddo
               goto 2
 10            continue
c     una differenza
               n1=n1-ncore
               n2=n2-ncore
               cij=c(m)*c(n)*segno
               if(ns1.eq.1.and.ns2.eq.1)then
                  do i=1,nocca
                     ia=iocca(i)
                     if(ia.ne.n1.and.ia.ne.n2)then
                        d2alal(n1,ia,n2,ia)=d2alal(n1,ia,n2,ia)+cij
                        call permbis(d2alal,n1,ia,n2,ia,nact)
                     endif
                  enddo
                  do i=1,noccb
                     ib=ioccb(i)
                     d2albe(n1,ib,n2,ib)=d2albe(n1,ib,n2,ib)+cij
                     d2albe(n2,ib,n1,ib)=d2albe(n2,ib,n1,ib)+cij
                  enddo
               else
                  do i=1,noccb
                     ib=ioccb(i)
                     if(ib.ne.n1.and.ib.ne.n2)then
                        d2bebe(n1,ib,n2,ib)=d2bebe(n1,ib,n2,ib)+cij
                        call permbis(d2bebe,n1,ib,n2,ib,nact)
                     endif
                  enddo
                  do i=1,nocca
                     ia=iocca(i)
                     d2albe(ia,n1,ia,n2)=d2albe(ia,n1,ia,n2)+cij
                     d2albe(ia,n2,ia,n1)=d2albe(ia,n2,ia,n1)+cij
                  enddo
               endif
               goto 2
c     due differenze
 20            continue
               n1=n1-ncore
               n2=n2-ncore
               n3=n3-ncore
               n4=n4-ncore
               cij=c(m)*c(n)*segno
c               print '(a,i1,a,i1,a,i1,a,i1)','Due diff.: ns1=',ns1,
c     $              ' ns2=',ns2,' ns3=',ns3,' ns4=',ns4
               if(ns1.eq.1.and.ns2.eq.1)then
                  d2alal(n1,n2,n3,n4)=d2alal(n1,n2,n3,n4)+cij
                  call permbis(d2alal,n1,n2,n3,n4,nact)
               elseif(ns1.eq.0.and.ns2.eq.0)then
                  d2bebe(n1,n2,n3,n4)=d2bebe(n1,n2,n3,n4)+cij
                  call permbis(d2bebe,n1,n2,n3,n4,nact)
               elseif(ns1.eq.1.and.ns2.eq.0)then
                  d2albe(n1,n2,n3,n4)=d2albe(n1,n2,n3,n4)+cij
                  d2albe(n3,n4,n1,n2)=d2albe(n3,n4,n1,n2)+cij
               endif
 2          continue
 1       continue
c--costruzione ro2 spinless
         do i=1,nact
            do j=1,nact
               do ip=1,nact
                  do jp=1,nact
                     d2alal(i,j,ip,jp)=d2alal(i,j,ip,jp)+
     $                    d2albe(i,j,ip,jp)+d2albe(j,i,jp,ip)
     $                    +d2bebe(i,j,ip,jp)
                  enddo
               enddo
            enddo
         enddo
         deallocate(d2bebe)
         deallocate(d2albe)
         return
         end
c-----------------------------------------------------
      subroutine ro1(c,d1,ncore,nact,thr,trou,part,ne,nd,ispin,iorb
     $     ,itsym)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d1(nact,nact)
      integer segno
c      pointer ne,nd,ispin,iorb,itsym
      dimension c(*)
c      COMMON/COE/C
      COMMON /CIP/ 
     *ITS(20,20),IGELS(2),NCFG,
     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispin(*),iorb(*),itsym(*)
c      pointer trou,part
c      dimension trou(:),part(:)
      dimension trou(*),part(*)
c      COMMON /BUPA/TROU,PART,NREQ
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      n=nact
      nocc=noca-ncore
      nn=(n*n+n)/2
c      print*,'subr. ro1'
c      print*,'n=',n,' nocc=',nocc,' ncf=',ncf
      nac=nocc
      do j=1,nac
         d1(j,j)=2.d0
      enddo
      do 1 i=1,ncf
         if(abs(c(i)).lt.thr)goto 1
         do 2 j=1,i
            if(abs(c(j)).lt.thr)goto 2
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
            d1(n1,n2)=d1(n1,n2)+c(i)*c(j)*segno
            d1(n2,n1)=d1(n1,n2)
            goto 2
 3          nec=ne(i)
            if(nec.eq.0)goto 2
            cmi2=c(i)*c(i)
            ndi=nd(i)
            do 30 l=1,nec
               lt=trou(l+ndi)
               lp=part(l+ndi)
               lt=iorb(lt)-ncore
               lp=iorb(lp)-ncore
               d1(lt,lt)=d1(lt,lt)-cmi2
               d1(lp,lp)=d1(lp,lp)+cmi2
 30         continue
 2       continue
 1    continue
      return
      end
c-------------------------------------------------
      subroutine ord2(nv,nw,zseg,segno)
c--prima gli alfa poi i beta
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
      return
      end
c----------------------------------------------------------------
      subroutine ord3(nv,nw,zseg,segno)
c--prima gli alfa poi i beta
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
      return
      end
c----------------------------------------------------------------
      subroutine ord4(nv,nw,zseg,segno)
c--prima gli alfa poi i beta
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
      return
      end
c----------------------------------------------------------------
      subroutine giveocc2(m,iocca,ioccb,nact,ncore,nocca,noccb,nocc,nec
     $     ,ndm,trou,part,iorb,ispin) 
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension iocca(*),ioccb(*),trou(*),part(*),iorb(*),ispin(*)
      integer*2 trou,part,iorb,ispin
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
         its=ispin(it)
         ips=ispin(ip)
         it=iorb(it)-ncore
         ip=iorb(ip)-ncore
         if(its.eq.1)then
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
      return
      end
c--------------------------------------------------------------------
      subroutine perm3(d,i,j,k,l,m,n,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(nact,nact,nact,nact,nact,nact)
c--block 1
      a=d(i,j,k,l,m,n)
      d(j,i,k,l,m,n)=-a
      d(i,k,j,l,m,n)=-a
      d(k,j,i,l,m,n)=-a
      d(k,i,j,l,m,n)=a
      d(j,k,i,l,m,n)=a
c--block 2
      d(i,j,k,m,l,n)=-a
      d(j,i,k,m,l,n)=a
      d(i,k,j,m,l,n)=a
      d(k,j,i,m,l,n)=a
      d(k,i,j,m,l,n)=-a
      d(j,k,i,m,l,n)=-a
c--block 3
      d(i,j,k,l,n,m)=-a
      d(j,i,k,l,n,m)=a
      d(i,k,j,l,n,m)=a
      d(k,j,i,l,n,m)=a
      d(k,i,j,l,n,m)=-a
      d(j,k,i,l,n,m)=-a
c--block 4
      d(i,j,k,n,m,l)=-a
      d(j,i,k,n,m,l)=a
      d(i,k,j,n,m,l)=a
      d(k,j,i,n,m,l)=a
      d(k,i,j,n,m,l)=-a
      d(j,k,i,n,m,l)=-a
c--block 5
      d(i,j,k,n,l,m)=a
      d(j,i,k,n,l,m)=-a
      d(i,k,j,n,l,m)=-a
      d(k,j,i,n,l,m)=-a
      d(k,i,j,n,l,m)=a
      d(j,k,i,n,l,m)=a
c--block 6
      d(i,j,k,m,n,l)=a
      d(j,i,k,m,n,l)=-a
      d(i,k,j,m,n,l)=-a
      d(k,j,i,m,n,l)=-a
      d(k,i,j,m,n,l)=a
      d(j,k,i,m,n,l)=a
      return
      end
c-----------------------------------------------------------
      subroutine perm2(d,i,j,k,l,m,n,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(nact,nact,nact,nact,nact,nact)
c--block 1
      a=d(i,j,k,l,m,n)
      d(j,i,k,l,m,n)=-a
      d(l,m,n,i,j,k)=a
      d(l,m,n,j,i,k)=-a
c--block 2
      d(i,j,k,m,l,n)=-a
      d(j,i,k,m,l,n)=a
      d(m,l,n,i,j,k)=-a
      d(m,l,n,j,i,k)=a
      return
      end
c--------------------------------------------------------------------
      subroutine permbis(d,i,j,k,l,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension d(nact,nact,nact,nact)
c--block 1
      a=d(i,j,k,l)
      d(j,i,k,l)=-a
      d(k,l,i,j)=a
      d(k,l,j,i)=-a
c--block 2
      d(i,j,l,k)=-a
      d(j,i,l,k)=a
      d(l,k,i,j)=-a
      d(l,k,j,i)=a
      return
      end
c--------------------------------------------------------------------
      subroutine permut123(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
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
c---permutazioni dispari
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
      return
      end
c----------------------------------------------------------
      subroutine permut124(i,j,l,k,ip,jp,lp,kp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
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
c---permutazioni dispari
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
      return
      end
c--------------------------------------------------------------------
      subroutine permut134(i,l,j,k,ip,lp,jp,kp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
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
c---permutazioni dispari
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
      return
      end
c--------------------------------------------------------------------
      subroutine permut234(l,i,j,k,lp,ip,jp,kp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
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
c---permutazioni dispari
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
      return
      end
c--------------------------------------------------------------------
      subroutine permut1234(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
c--permutazioni pari
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
c---permutazioni dispari
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
      return
      end
c--------------------------------------------------------------------
      subroutine permut1324(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
c--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
c         ip=k
c         jp=i
c         kp=l
c         lp=j
         ip=k
         jp=l
         kp=i
         lp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.3)then
c         ip=k
c         jp=i
c         kp=j
c         lp=l
         ip=k
         jp=j
         kp=i
         lp=l
         segno=-1
      elseif(ivolte.eq.4)then
c         ip=i
c         jp=k
c         kp=l
c         lp=j
         ip=i
         jp=l
         kp=k
         lp=j
         segno=-1
      endif
      return
      end
c--------------------------------------------------------------------
      subroutine permut1423(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
c--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
c         ip=l
c         jp=i
c         kp=k
c         lp=j
         ip=l
         jp=k
         kp=j
         lp=i
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.3)then
c         ip=l
c         jp=i
c         kp=j
c         lp=k
         ip=l
         jp=j
         kp=k
         lp=i
         segno=-1
      elseif(ivolte.eq.4)then
c         ip=i
c         jp=l
c         kp=k
c         lp=j
         ip=i
         jp=k
         kp=j
         lp=l
         segno=-1
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine permut4(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer segno
      if(ivolte.eq.1)then
c--permutazioni pari      
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
c---permutazioni dispari
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
      return
      end
c---------------------------------------------------------------------
      FUNCTION Jdifgen(ne,nd,trou,part,II,JJ,nv1,nv2,nma,zsig)
      implicit real*8(a-h,o-y),logical*1(z)
C     
C     SSP DE RECHERCHE DES ELEMENTS D'INTERACTION
C     ENTRE DEUX DETERMINANTS, NUMEROTES II ET JJ
C     
c      include 'rotwo.par'
      parameter (nmax=4)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      pointer trou,part
c      dimension trou(:),part(:)
c      COMMON /BUPA/TROU,PART,NREQ
c      INTEGER*2 ISPIN,IORB
      dimension ne(*),nd(*),trou(*),part(*)
      INTEGER*2 NE,TROU,PART
      DIMENSION ITA(99),ITB(99),NV1(nmax),NV2(nmax),NATU1(nmax)
     $     ,NATU2(nmax)
c      EQUIVALENCE(NV11,NV1(1)),(NV12,NV1(2)),(NV21,NV2(1)),(NV22,NV2(2))
      common/diff/n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      KR=II
      KL=JJ
C     
C     LE DETERMINANT - 1 - EST CHOISI COMME LE
C     DETERMINANT AYANT LE PLUS D'O.M. EXCITEES
C     
      zshift=.false.
      IF (NE(KR)-NE(KL)) 5,7,7
    5 KKKK=KR
      KR=KL
      KL=KKKK
      zshift=.true.
    7 NE1=NE(KR)
C     TEST SUR LE FONDAMENTAL
      IF (NE1.EQ.0) GOTO 44
      NE2=NE(KL)
C     NBDIF : EST LE NOMBRE DE SPIN-ORBITALES
C     DIFFERENTES
      NBDIF=NE1-NE2
      IF (NBDIF.GT.nma) GOTO 24
C     CONSTRUCTION POUR LE DETERMINANT -1- ET -2-
C     DE IT ET NS, QUI CONTIENNENT LES NUMEROS
C     DES ORBITALES ET LEURS SPINS
C     
      ND1=ND(KR)
      ND2=ND(KL)
      DO 15 J=1,NE1
         ITA(J)=TROU(J+ND1)
         ITA(J+NE1)=PART(J+ND1)
 15   CONTINUE
C     
C     ZSIG :  EST LA PARITE DU  NOMBRE DE CROISEMENTS
C     DANS LE DIAGRAMME D INTERACTION
      ZSIG=.FALSE.
C     SI NE2=0 LE DETERMINANT -2- EST LE FONDAMENTAL
      IF (NE2.EQ.0) GOTO 125
 121  CONTINUE
      DO 25 J= 1,NE2
         ITB(J)=TROU(J+ND2)
 25   ITB(J+NE2)=PART(J+ND2)
    4 DO 8 K=1,2
         K1=(K-1)*NE1
         K2=(K-1)*NE2
         DO 10I=1,NE2
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
C     
C     NAR : EXCITATION SUPLEMENTAIRE DU DETERMINANT -1-
C     PAR RAPPORT AU DETERMINANT -2-
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
C     CALCUL DE L ELEMENT DE MATRICE
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
         END
c------------------------------------------------------------------
      integer function index(i,j,k,l)
      if(i.ge.k)then
         ik=(i*i-i)/2+k
      else
         ik=(k*k-k)/2+i
      endif
      if(j.ge.l)then
         jl=(j*j-j)/2+l
      else
         jl=(l*l-l)/2+j
      endif
      if(ik.ge.jl)then
         index=(ik*ik-ik)/2+jl
      else
         index=(jl*jl-jl)/2+ik
      endif
      return
      end
c------------------------------------------------------------------
      subroutine accum(amat,ip,jp,kp,lp,mp,np,op,pp,val,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
      dimension amat(nwords,nact,nact,nact,nact)
      integer op,pp
      ind=lindice(ip,jp,kp,lp)
      amat(ind,mp,np,op,pp)=amat(ind,mp,np,op,pp)+val
c      write (6,'(''accum'',8I3,I6,2f15.8)') 
c     * ip,jp,kp,lp,mp,np,op,pp,ind,val,amat(ind,mp,np,op,pp)
      return
      end
c------------------------------------------------------------------
      function lindiceold(i,j,k,l,m,n,o,p)
      integer o,p,op,pp
      logical zi,zj,zjk,zk,zl,zij,zijk
      dimension ivec(12)
      data ivec/1,48,639,4371,20327,73321,220741,580353,
     $     1373052,2985159,6055496,11593084/
      ind=ivec(i-1)
      do jp=1,j
         zj=jp.eq.j
         if(zj)then
            klim=k
         else
            klim=jp
         endif
         do kp=1,klim
            zk=kp.eq.k
            zjk=zj.and.zk
            if(zjk)then
               llim=l
            else
               llim=kp
            endif
            do lp=1,llim
c                  zl=lp.eq.l
               measure=i+jp+kp+lp
               do mp=1,i
                  do np=1,i
                     do op=1,i
                        do pp=1,i
                           measurep=mp+np+op+pp
                           if(measurep.le.measure)ind=ind+1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      lindiceold=ind
      return
      end
c------------------------------------------------------------------
      function lindice(i,j,k,l)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
      lindice=lindi(i)+lindj(j)+lindk(k)+l
      return
      end
c-----------------------------------------------------------------
      SUBROUTINE MATOUT(N,N2,A,NN)
      DOUBLE PRECISION A
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
      WRITE(6,30) I,(A(I,J),J=M,K)
 40   CONTINUE
      call flush(6)
      RETURN
      END
c-------------------------------------------------------------------
      subroutine readint(atwo,factor,escf,nact,ncore,indic,jndic,lndic
     $     ,num,f,aj,ne,nd,ispin,iorb,itsym)
      implicit real*8 (a-h,o-y),logical*1 (z)
      allocatable ijkl4(:)
      character*80 label
      dimension ne(*),nd(*),ispin(*),iorb(*),itsym(*)
      COMMON /CIP/ 
     *ITS(20,20),IGELS(2),NCFG,
     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
      INTEGER*2 ISPIN,IORB
      INTEGER*2 NE,TROU,PART
      dimension f(*),aj(*)
c      pointer trou,part
c      dimension trou(:),part(:)
c      dimension trou(*),part(*)
c      COMMON /BUPA/TROU,PART,NREQ
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * NDEB(500),NBO(99),
     * NAD(2000),KT(2000)
      dimension INDIC(*),JNDIC(*),LNDIC(*),num(*)
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      dimension nos(20)
      integer a,b,ac,bc
      dimension atwo(nact,nact,nact,nact)
      common /toul/ ztoul
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c      jd1=id1
      jd1=norb
      call flush(6)
      CALL REIJKL(factor,norb,nsym,itsym,its,nos,jd1,indic,jndic,lndic
     $     ,num)
      finv=1.d0/factor
c      print*,'    escf=',escf
c      allocate(ijkl2(1))
c      allocate(ijkl3(1))
c      allocate(rijkl(1))
c      allocate(x(1))
      allocate(ijkl4(ni4))
      if (ztoul) then
         CALL ITIJKLt(factor,norb,nsym,nos,
     * ijkl4,ijkl2,rijkl,ijkl3,indic,jndic,lndic,num)
      else
         CALL ITIJKL(factor,norb,nsym,nos,
     *        ijkl4,ijkl2,rijkl,ijkl3,indic,jndic,lndic,num)
      endif
      NTER=NORB*(NORB+1)/2
      if(.not.ztoul)then
         call filesplit('READ',50,8*NTER,ichange,idum2)
         READ (50) (F(L),L=1,NTER) !matrice di Fock
         call filesplit('READ',50,8*NTER,ichange,idum2)
         READ (50) (AJ(L),L=1,NTER) !monoelettronici puri
      else
c      read(26)
c      do
c         read(26) label
c         if(index(label,'==ONEINT MATRIX (TOTAL)').ne.0) exit
c         read(26)
c         read(26)
c      enddo

c      read(26) ! norbttMono,nsymMono,isymMono(1:nsymMono)
c      read(26) aj(1:nter)
         READ (26) (AJ(L),L=1,NTER) !monoelettronici puri
      endif
c----building Dyall's h eff in aj----------
      do a=1,nact
         ac=a+ncore
         do b=a,nact
            bc=b+ncore
            jab=indice(ac,bc)
            aj(jab)=aj(jab)*factor
            do j=1,ncore
               aj(jab)=aj(jab)+2.d0*ai(ac,bc,j,j,indic,jndic,lndic,num
     $              ,ijkl4)-ai(ac,j,j,bc,indic,jndic,lndic,num,ijkl4)
            enddo
         enddo
      enddo
c---building matrix atwo
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  atwo(i,j,k,l)=ai(i+ncore,j+ncore,k+ncore,l+ncore,indic
     $                 ,jndic,lndic,num,ijkl4)
               enddo
            enddo
         enddo
      enddo
c------deallocating memory
c      deallocate(ijkl2)
c      deallocate(ijkl3)
c      deallocate(rijkl)
c      deallocate(x)
      deallocate(ijkl4)
      return
      end
c-------------------------------------------------------------------
      subroutine koop2E(atwo,daaa,taa,dal,koopaa,koopeaa,factor,nact
     $     ,ncore,f,aj,itsym,its)
      implicit real*8 (a-h,o-y),logical*1 (z)
c---F e' la matrice dei monoelettronici!
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
      dimension atwo(nact,nact,nact,nact)
      dimension koopaa(nact,nact,nact,nact),koopeaa(nact,nact,nact,nact)
      real*8 koopaa,koopbb,koopab,koopeaa,koopebb,koopeab
c      pointer ne,nd,ispin,iorb,itsym
c      dimension ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      INTEGER*2 ISPIN,IORB
c      INTEGER*2 NE,TROU,PART
      dimension itsym(*),its(20,20)
      dimension f(*),aj(*)
c      pointer trou,part
c      dimension trou(:),part(:)
c      dimension trou(*),part(*)
c      COMMON /BUPA/TROU,PART,NREQ
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
      dimension nos(20)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,apc,bpc
      dimension w(100),work(1000)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nocc=noca
c      print*,'nocc=',nocc,' ncore=',ncore,' nact=',nact
c      call flush(6)
c--building heff'
c--in AJ is heff
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
c-------------------
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
                  koopaa(ap,bp,a,b)=0.d0
		  if (its(isab,iscd).ne.1) goto 1
                  do d=1,nact
                     dc=d+ncore
                     jad=indice(ac,dc)
                     jbd=indice(bc,dc)
                     koopaa(ap,bp,a,b)=koopaa(ap,bp,a,b)-f(jad)*taa(ap
     $                    ,bp,d,b)-f(jbd)*taa(ap,bp,a,d)
c                     write (6,*) 'a,b,ap,bp,d',a,b,ap,bp,d
c		     write (6,*) 'k',koopaa(ap,bp,a,b)/factor
c     $                     ,f(jad),taa(ap
c     $                    ,bp,d,b),f(jbd),taa(ap,bp,a,d)
                     do ie=1,nact
                        do ief=1,nact
                           dint=atwo(d,ie,a,ief)
                           dum=daaa(ap,bp,d,ief,b,ie)
                           if(d.eq.ief)dum=dum+0.5d0*taa(ap,bp,ie,b)
                           if(d.eq.b)dum=dum+0.5d0*taa(ap,bp,ief,ie)
                           koopaa(ap,bp,a,b)=koopaa(ap,bp,a,b)-dint*dum
c                     write (6,*) 'a,b,ap,bp,d,e,f',a,b,ap,bp,d,ie,ief
c           write (6,*) 'k',koopaa(ap,bp,a,b)/factor,f(jad),dint,dum
                           dint=atwo(d,ie,b,ief)
                           dum=daaa(ap,bp,d,a,ief,ie)
                           if(d.eq.a)dum=dum+0.5d0*taa(ap,bp,ie,ief)
                           if(d.eq.ief)dum=dum+0.5d0*taa(ap,bp,a,ie)
                           koopaa(ap,bp,a,b)=koopaa(ap,bp,a,b)-dint*dum
c                     write (6,*) 'a,b,ap,bp,d,e,f',a,b,ap,bp,d,ie,ief
c           write (6,*) 'k',koopaa(ap,bp,a,b)/factor,f(jad),dint,dum
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
                  koopeaa(ap,bp,a,b)=0.d0
		  if (its(isab,iscd).ne.1) goto 2
                  do c=1,nact
                     cc=c+ncore
                     jac=indice(ac,cc)
                     jbc=indice(bc,cc)
                     dum1=ro2t(ap,bp,c,b,taa,dal,nact)
                     dum2=ro2t(ap,bp,a,c,taa,dal,nact)
                     koopeaa(ap,bp,a,b)=koopeaa(ap,bp,a,b)+f(jac)*dum1
     $                    +f(jbc)*dum2
                     do d=1,nact
                        do ie=1,nact
                           dint=atwo(c,ie,d,a)
                           dum=-ro3t(ap,bp,ie,d,b,c,daaa,taa,dal,nact)
                           if(c.eq.ie)dum=dum+2.d0*ro2t(ap,bp,d,b,taa
     $                          ,dal,nact)
                           if(d.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,c,b,taa
     $                          ,dal,nact)
                           if(b.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,d,c,taa
     $                          ,dal,nact)
                           koopeaa(ap,bp,a,b)=koopeaa(ap,bp,a,b)+dint
     $                          *dum
                           dint=atwo(c,ie,d,b)
                           dum=-ro3t(ap,bp,ie,a,d,c,daaa,taa,dal,nact)
                           if(c.eq.ie)dum=dum+2.d0*ro2t(ap,bp,a,d,taa
     $                          ,dal,nact)
                           if(a.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,c,d,taa
     $                          ,dal,nact)
                           if(d.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,a,c,taa
     $                          ,dal,nact)
                           koopeaa(ap,bp,a,b)=koopeaa(ap,bp,a,b)+dint
     $                          *dum
                        enddo
                     enddo
                  enddo
 2             enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine koopman0pE(atwo,daaa,taa,dal,koopaa,koopbb,koopab
     $     ,factor,nact,ncore,f,aj,itsym,its)
      implicit real*8 (a-h,o-y),logical*1 (z)
c---F e' la matrice dei monoelettronici!
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),atwo(nact,nact,nact,nact),dal(nact,nact)
      dimension koopaa(nact,nact,nact,nact),koopbb(nact,nact,nact,nact)
     $     ,koopab(nact,nact)
      real*8 koopaa,koopbb,koopab
      dimension itsym(*),its(20,20)
c      pointer ne,nd,ispin,iorb,itsym
c      dimension ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      INTEGER*2 ISPIN,IORB
c      INTEGER*2 NE,TROU,PART
c      pointer trou,part
c      dimension trou(:),part(:)
c      dimension trou(*),part(*)
      dimension f(*),aj(*)
c      COMMON /BUPA/TROU,PART,NREQ
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,apc,bpc
      dimension w(100),work(1000)
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nocc=noca
c      print*,'nocc=',nocc,' ncore=',ncore,' nact=',nact
c      call flush(6)
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
                  koopaa(bp,ap,a,b)=0.d0
                  do c=1,nact
                     cc=c+ncore
                     ica=indice(cc,ac)
                     ibc=indice(bc,cc)
                     koopaa(bp,ap,a,b)=koopaa(bp,ap,a,b)+2.d0*(f(ica)
     $                    *ee2(bp,ap,c,b,taa,dal,nact)-f(ibc)
     $                    *ee2(bp,ap,a,c,taa,dal,nact))
                     do d=1,nact
                        dc=d+ncore
                        do ie=1,nact
                           iec=ie+ncore
                           dum1=eee2(bp,ap,c,ie,d,b,daaa
     $                          ,taa,dal,nact)+eee2(bp,ap,d,b
     $                          ,c,ie,daaa,taa,dal,nact)
                           dum2=eee2(bp,ap,a,ie,c,d,daaa,taa,dal,nact)
     $                          +eee2(bp,ap,c,d,a,ie,daaa,taa,dal,nact)
                           koopaa(bp,ap,a,b)=koopaa(bp,ap,a,b)+atwo(c
     $                          ,ie,d,a)*dum1-atwo(b,ie,c,d)*dum2
                        enddo
                     enddo
                  enddo
                  koopbb(bp,ap,a,b)=0.d0
                  do c=1,nact
                     cc=c+ncore
                     ica=indice(cc,ac)
                     ibc=indice(bc,cc)
                     dum1=-ee2(a,ap,bp,c,taa,dal,nact)
                     if(a.eq.ap)dum1=dum1+two*dal(bp,c)
                     if(bp.eq.ap)dum1=dum1+dal(a,c)
                     dum2=-ee2(c,ap,bp,b,taa,dal,nact)
                     if(ap.eq.c)dum2=dum2+two*dal(bp,b)
                     if(bp.eq.ap)dum2=dum2+dal(c,b)
                     koopbb(bp,ap,a,b)=koopbb(bp,ap,a,b)-f(ibc)*dum1
     $                    +f(ica)*dum2
                     do d=1,nact
                        dc=d+ncore
                        do ie=1,nact
                           iec=ie+ncore
                           dum1=-eee2(c,ie,a,ap,bp,d,daaa,taa,dal,nact)
     $                          -eee2(a,ap,bp,d,c,ie,daaa,taa,dal,nact)
                           if(ap.eq.a)dum1=dum1+two*(ee2(c,ie,bp,d,taa
     $                          ,dal,nact)+ee2(bp,d,c,ie,taa
     $                          ,dal,nact))
                           if(bp.eq.ap)dum1=dum1+ee2(c,ie,a,d,taa,dal
     $                          ,nact)+ee2(a,d,c,ie,taa,dal,nact)
                           if(ap.eq.c)then
                              dum1=dum1-ee2(a,ie,bp,d,taa,dal,nact)
                              if(ie.eq.a)dum1=dum1+two*dal(bp,d)
                           endif
                           if(bp.eq.ie)then
                              dum1=dum1+ee2(a,ap,c,d,taa,dal,nact)
                              if(ap.eq.a)dum1=dum1-two*dal(c,d)
                           endif
                           dum2=-eee2(c,ie,d,ap,bp,b,daaa
     $                          ,taa,dal,nact)-eee2(d,ap,bp,b
     $                          ,c,ie,daaa,taa,dal,nact)
                           if(ap.eq.d)dum2=dum2+two*(ee2(c,ie,bp,b,taa
     $                          ,dal,nact)+ee2(bp,b,c,ie,taa
     $                          ,dal,nact))
                           if(bp.eq.ap)dum2=dum2+ee2(c,ie,d,b,taa,
     $                          dal,nact)+ee2(d,b,c,ie,taa,dal,nact)
                           if(ap.eq.c)then
                              dum2=dum2-ee2(d,ie,bp,b,taa,dal,nact)
                              if(d.eq.ie)dum2=dum2+two*dal(bp,b)
                           endif
                           if(bp.eq.ie)then
                              dum2=dum2+ee2(d,ap,c,b,taa,dal,nact)
                              if(d.eq.ap)dum2=dum2-two*dal(c,b)
                           endif
                           koopbb(bp,ap,a,b)=koopbb(bp,ap,a,b)-0.5d0*
     $                          atwo(c,ie,b,d)*dum1+0.5d0*atwo(c,ie,d
     $                          ,a)*dum2
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            koopab(a,b)=0.d0
            do c=1,nact
               cc=c+ncore
               ica=indice(cc,ac)
               ibc=indice(bc,cc)
               dum1=dal(c,b)
               dum2=dal(a,c)
               koopab(a,b)=koopab(a,b)+two*(f(ica)*dum1
     $              -f(ibc)*dum2)
               do d=1,nact
                  dc=d+ncore
                  do ie=1,nact
                     iec=ie+ncore
                     dum1=ee2(c,ie,d,b,taa,dal,nact)+ee2(d,b,c
     $                    ,ie,taa,dal,nact)
                     dum2=ee2(a,ie,c,d,taa,dal,nact)+ee2(c,d,a
     $                    ,ie,taa,dal,nact)
                     koopab(a,b)=koopab(a,b)+atwo(c,ie,d
     $                    ,a)*dum1-atwo(b,ie,c,d)*dum2
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bamat(amat,atmat,atwo,ro4,daaa,taa,dal,nact,ncore,f,aj
     $     ,itsym,its,norb)
      implicit real*8 (a-h,o-y),logical*1 (z)
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     *     ,normind(18,18,18,18)
      dimension ro4(nwords,nact,nact,nact,nact)
      dimension amat(nact,nact,nact,nact,nact,nact)
      dimension atmat(nact,nact,nact,nact,nact,nact)
      dimension atwo(nact,nact,nact,nact)
      dimension itsym(*),its(20,20)
c      pointer ne,nd,ispin,iorb,itsym
c      dimension ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      INTEGER*2 ISPIN,IORB
c      INTEGER*2 NE,TROU,PART
c      pointer trou,part
c      dimension trou(:),part(:)
c      COMMON /BUPA/TROU,PART,NREQ
      dimension f(*),aj(*)
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
c      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      allocatable indice(:,:)
c      dimension indice(id1,id1)
      allocate(indice(norb,norb))
      nocc=noca
c      print*,'nocc=',nocc,' ncore=',ncore,' nact=',nact
c      call flush(6)
c-------------------
      do i=1,norb
       do j=1,norb
       indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
       enddo
      enddo
c--building heff'
c--in AJ is heff
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
c-------------------
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
			goto 1
			endif
                        do d=1,nact
                           dc=d+ncore
			   isd=itsym(dc)
                           ida=indice(dc,ac)
                           icd=indice(cc,dc)
                           ibd=indice(bc,dc)
                           amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                          ,c)+f(ida)*eee2(cp,ap,bp,b,d,c,daaa,taa
     $                          ,dal,nact)-f(icd)*eee2(cp,ap,bp,b,a,d
     $                          ,daaa,taa,dal,nact)-f(ibd)*eee2(cp,ap,bp
     $                          ,d,a,c,daaa,taa,dal,nact)
                           atmat(ap,bp,cp,a,b,c)=atmat(ap,bp,cp,a,b
     $                          ,c)+f(ida)*eee2t(cp,ap,bp,b,d,c,daaa,taa
     $                          ,dal,nact)-f(icd)*eee2t(cp,ap,bp,b,a,d
     $                          ,daaa,taa,dal,nact)+f(ibd)*eee2t(cp,ap
     $                          ,bp,d,a,c,daaa,taa,dal,nact)
                           do ie=1,nact
                              iec=ie+ncore
			      ise=itsym(iec)
                              dint1=atwo(c,ie,d,a)
                              dint2=atwo(d,ie,ie,a)
                              dint3=atwo(c,d,d,ie)
                              dint4=atwo(b,d,d,ie)
                              dum1=eee2(cp,ap,bp,b,d,ie,daaa,taa,dal
     $                             ,nact)
                              dum2=eee2(cp,ap,bp,b,d,c,daaa,taa,dal
     $                             ,nact)
                              dum3=eee2(cp,ap,bp,b,a,ie,daaa,taa,dal
     $                             ,nact)
                              dum4=eee2(cp,ap,bp,ie,a,c,daaa,taa,dal
     $                             ,nact)
c                              dum1t=eee2t(cp,ap,bp,b,d,ie,daaa,taa,dal
c     $                             ,nact)
c                              dum2t=eee2t(cp,ap,bp,b,d,c,daaa,taa,dal
c     $                             ,nact)
c                              dum3t=eee2t(cp,ap,bp,b,a,ie,daaa,taa,dal
c     $                             ,nact)
                              dum1t=-dum1
                              dum2t=-dum2
                              dum3t=-dum3
                              if(bp.eq.b)then
                                 dum1t=dum1t+2.d0*ee2(cp,ap,d
     $                                ,ie,taa,dal,nact)
                                 dum2t=dum2t+2.d0*ee2(cp,ap,d,c
     $                                ,taa,dal,nact)
                                 dum3t=dum3t+2.d0*ee2(cp,ap,a,ie
     $                                ,taa,dal,nact)
                              endif
                              atmat(ap,b,cp,a,bp,c)=atmat(ap,b,cp,a,bp,c
     $                             )+dum1t*dint1-0.5d0*(dum2t*dint2
     $                             +dum3t*dint3)
                              dum4t=eee2t(cp,ap,bp,ie,a,c,daaa,taa,dal
     $                             ,nact)
                              amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b,c)+
     $                             dum1*dint1-0.5d0*(dum2*dint2+dum3
     $                             *dint3-dum4*dint4)
                              atmat(ap,bp,cp,a,b,c)=atmat(ap,bp,cp,a,b,c
     $                             )+0.5d0*dum4t*dint4
                              do ief=1,nact
                                 iefc=ief+ncore
				 isf=itsym(iefc)
				 isfd=its(isf,isd)
				 isfde=its(isfd,ise)
                                 istot=its(isabcp,isfde)
				 if (its(isfde,isa).eq.1.and.its(istot
     $                                ,isbc).eq.1) then
                                 dum1=eeee(cp,ap,bp,b,d,ief,ie,c,daaa
     $                                   ,taa,dal,ro4,nact)
                                 dum1t=-dum1
                                 if(bp.eq.b)dum1t=dum1t+2.d0*eee2(cp,ap
     $                                ,d,ief,ie,c,daaa,taa,dal,nact)
c                                 dum1t=eeeet(cp,ap,bp,b,d,ief,ie,c,daaa
c     $                                ,taa,dal,ro4,nact)
                                 dint1=atwo(d,ief,ie,a)
                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                                ,c)+dum1*dint1
                                 atmat(ap,b,cp,a,bp,c)=atmat(ap,b,cp,a
     $                                ,bp,c)+dum1t*dint1
				 endif
				 if (its(isfde,isc).eq.1.and.its(istot
     $                                ,isab).eq.1) then
                                 dum2=eeee(cp,ap,bp,b,d,ief,a,ie,daaa
     $                                   ,taa,dal,ro4,nact)
                                 dum2t=-dum2
                                 if(bp.eq.b)dum2t=dum2t+2.d0*eee2(cp,ap
     $                                ,d,ief,a,ie,daaa,taa,dal,nact)
c                                 dum2t=eeeet(cp,ap,bp,b,d,ief,a,ie,daaa
c     $                                ,taa,dal,ro4,nact)
                                 dint2=atwo(d,ief,c,ie)
                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                                ,c)-dum2*dint2
                                 atmat(ap,b,cp,a,bp,c)=atmat(ap,b,cp,a
     $                                ,bp,c)-dum2t*dint2
				 endif
				 if (its(isfde,isb).eq.1.and.its(istot
     $                                ,isac).eq.1) then
                                 dum3=eeee(cp,ap,bp,ie,d,ief,a,c,daaa
     $                                ,taa,dal,ro4,nact)
                                 dum3t=eeeet(cp,ap,bp,ie,d,ief,a,c,daaa
     $                                ,taa,dal,ro4,nact)
                                 dint3=atwo(d,ief,b,ie)
                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                                ,c)-dum3*dint3
                                 atmat(ap,bp,cp,a,b,c)=atmat(ap,bp,cp,a
     $                                ,b,c)+dum3t*dint3
                                 endif
c                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
c     $                                ,c)+dum1*dint1-dum2*dint2-dum3
c     $                                *dint3
                              enddo
                           enddo
                        enddo
   1                 enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bbmat(atwo,bmat,btmat,daaa,taa,dal,nact,ncore,f,aj)
      implicit real*8 (a-h,o-y),logical*1 (z)
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
      dimension bmat(nact,nact,nact,nact)
      dimension btmat(nact,nact,nact,nact)
      dimension atwo(nact,nact,nact,nact)
c      pointer ne,nd,ispin,iorb,itsym
c      dimension ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      INTEGER*2 ISPIN,IORB
c      INTEGER*2 NE,TROU,PART
c      pointer trou,part
c      dimension trou(:),part(:)
c      COMMON /BUPA/TROU,PART,NREQ
      dimension f(*),aj(*)
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nocc=noca
c      print*,'nocc=',nocc,' ncore=',ncore,' nact=',nact
c      call flush(6)
c---building matrix heff''
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
c-------------------------------------
      do ap=1,nact
         apc=ap+ncore
         do bp=1,nact
            bpc=bp+ncore
            do cp=1,nact
               cpc=cp+ncore
               do a=1,nact
                  ac=a+ncore
                  bmat(ap,bp,cp,a)=0.d0
                  btmat(ap,bp,cp,a)=0.d0
                  do c=1,nact
                     cc=c+ncore
                     iac=indice(cc,ac)
                     bmat(ap,bp,cp,a)=bmat(ap,bp,cp,a)-f(iac)
     $                    *ee2(cp,ap,bp,c,taa,dal,nact)
                     btmat(ap,bp,cp,a)=btmat(ap,bp,cp,a)+aj(iac)
     $                    *ee2tr(cp,ap,bp,c,taa,dal,nact)
                     do ie=1,nact
                        iec=ie+ncore
                        do ief=1,nact
                           iefc=ief+ncore
                           dum=eee2(cp,ap,bp,ie,c,ief,daaa,taa
     $                          ,dal,nact)
                           dumt=eee2t(cp,ap,bp,ie,c,ief,daaa,taa
     $                          ,dal,nact)
                           dint=atwo(a,ie,c,ief)
                           bmat(ap,bp,cp,a)=bmat(ap,bp,cp,a)-dum
     $                          *dint
                           btmat(ap,bp,cp,a)=btmat(ap,bp,cp,a)+dumt
     $                          *dint
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bcmat(atwo,cmat,ctmat,daaa,taa,dal,nact,ncore,f,aj)
      implicit real*8 (a-h,o-y),logical*1 (z)
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
      dimension cmat(nact,nact,nact,nact)
      dimension ctmat(nact,nact,nact,nact)
      dimension atwo(nact,nact,nact,nact)
c      pointer ne,nd,ispin,iorb,itsym
c      dimension ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      INTEGER*2 ISPIN,IORB
c      INTEGER*2 NE,TROU,PART
c      pointer trou,part
c      dimension trou(:),part(:)
c      COMMON /BUPA/TROU,PART,NREQ
      dimension f(*),aj(*)
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nocc=noca
c      print*,'nocc=',nocc,' ncore=',ncore,' nact=',nact
c      call flush(6)
c--building heff'
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
c-------------------
      do ap=1,nact
         apc=ap+ncore
         do a=1,nact
            ac=a+ncore
            do b=1,nact
               bc=b+ncore
               do c=1,nact
                  cc=c+ncore
                  cmat(ap,a,b,c)=0.d0
                  ctmat(ap,a,b,c)=0.d0
                  do d=1,nact
                     dc=d+ncore
                     ida=indice(dc,ac)
                     icd=indice(cc,dc)
                     ibd=indice(bc,dc)
                     cmat(ap,a,b,c)=cmat(ap,a,b,c)+f(ida)*ee2(ap,b,d,c
     $                    ,taa,dal,nact)-f(icd)*ee2(ap,b,a,d
     $                    ,taa,dal,nact)-f(ibd)*ee2(ap,d,a,c
     $                    ,taa,dal,nact)
                     ctmat(ap,a,b,c)=ctmat(ap,a,b,c)+f(ida)*ee2t(ap,b,d
     $                    ,c,taa,dal,nact)-f(icd)*ee2t(ap,b,a,d,taa,dal
     $                    ,nact)+f(ibd)*ee2t(ap,d,a,c,taa,dal,nact)
                     do ie=1,nact
                        iec=ie+ncore
                        dum1=ee2(ap,b,d,ie,taa,dal,nact)
                        dum2=ee2(ap,b,d,c,taa,dal,nact)
                        dum3=ee2(ap,b,a,ie,taa,dal,nact)
                        dum4=ee2(ap,ie,a,c,taa,dal,nact)
                        dum1t=ee2t(ap,b,d,ie,taa,dal,nact)
                        dum2t=ee2t(ap,b,d,c,taa,dal,nact)
                        dum3t=ee2t(ap,b,a,ie,taa,dal,nact)
                        dum4t=ee2t(ap,ie,a,c,taa,dal,nact)
                        dint1=atwo(c,ie,d,a)
                        dint2=atwo(d,ie,ie,a)
                        dint3=atwo(c,d,d,ie)
                        dint4=atwo(b,d,d,ie)
                        cmat(ap,a,b,c)=cmat(ap,a,b,c)+dum1*dint1-0.5d0
     $                       *(dum2*dint2+dum3*dint3-dum4*dint4)
                        ctmat(ap,a,b,c)=ctmat(ap,a,b,c)+dum1t*dint1-0
     $                       .5d0*(dum2t*dint2+dum3t*dint3-dum4t*dint4)
                        do ief=1,nact
                           iefc=ief+ncore
                           dum1=eee2(ap,b,d,ief,ie,c,daaa,taa
     $                          ,dal,nact)
                           dum1t=eee2tl(ap,b,d,ief,ie,c,daaa,taa
     $                          ,dal,nact)
                           dint1=atwo(d,ief,ie,a)
                           dum2=eee2(ap,b,d,ief,a,ie,daaa,taa
     $                          ,dal,nact)
                           dum2t=eee2tl(ap,b,d,ief,a,ie,daaa,taa
     $                          ,dal,nact)
                           dint2=atwo(d,ief,c,ie)
                           dum3=eee2(ap,ie,d,ief,a,c,daaa,taa
     $                          ,dal,nact)
                           dum3t=eee2tl(ap,ie,d,ief,a,c,daaa,taa
     $                          ,dal,nact)
                           dint3=atwo(d,ief,b,ie)
                           cmat(ap,a,b,c)=cmat(ap,a,b,c)+dum1*dint1-dum2
     $                          *dint2-dum3*dint3
                           ctmat(ap,a,b,c)=ctmat(ap,a,b,c)+dum1t*dint1
     $                          -dum2t*dint2+dum3t*dint3
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bdmat(atwo,dmat,dtmat,taa,dal,nact,ncore,f,aj)
      implicit real*8 (a-h,o-y),logical*1 (z)
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      dimension dmat(nact,nact)
      dimension dtmat(nact,nact)
      dimension atwo(nact,nact,nact,nact)
c      pointer ne,nd,ispin,iorb,itsym
c      dimension ne(:),nd(:)
c     $     ,ispin(:),iorb(:),itsym(:)
c      COMMON /CIP/ 
c     *NE,ND,ISPIN,
c     *IORB,ITSYM,ITS(20,20),IGELS(2),NCFG,
c     *NORB,NOCA,NOCB,MNORB,NCF,NSYM,ISYM,NTRSY,METAT,METAT1,ZQDPT,
c     *ZPRT,ZKEEP,ZUNI,ZMP,ZWRT,ZSEG
c      INTEGER*2 ISPIN,IORB
c      INTEGER*2 NE,TROU,PART
      dimension f(*),aj(*)
c      pointer trou,part
c      dimension trou(:),part(:)
c      COMMON /BUPA/TROU,PART,NREQ
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nocc=noca
c      print*,'nocc=',nocc,' ncore=',ncore,' nact=',nact
c      call flush(6)
c---building matrix heff''
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
c---------------------------
      do ap=1,nact
         apc=ap+ncore
         do a=1,nact
            ac=a+ncore
            dmat(ap,a)=0.d0
            dtmat(ap,a)=0.d0
            do c=1,nact
               cc=c+ncore
               iac=indice(cc,ac)
               dmat(ap,a)=dmat(ap,a)-f(iac)*dal(ap,c)
               dum=-dal(ap,c)
               if(ap.eq.c)dum=dum+2.d0
               dtmat(ap,a)=dtmat(ap,a)+aj(iac)*dum
               do ie=1,nact
                  iec=ie+ncore
                  do ief=1,nact
                     iefc=ief+ncore
                     dum=ee2(ap,ie,c,ief,taa,dal,nact)
                     dumt=ee2t(ap,ie,c,ief,taa,dal,nact)
                     dint=atwo(a,ie,c,ief)
                     dmat(ap,a)=dmat(ap,a)-dum*dint
                     dtmat(ap,a)=dtmat(ap,a)+dumt*dint
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      function eeee(a,b,c,d,e,f,g,h,d3,d2,d1,ro4,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f,g,h
      dimension d3(nact,nact,nact,nact,nact,nact),d2(nact,nact,nact,nact
     $     ),d1(nact,nact)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     *     ,normind(18,18,18,18)
      dimension ro4(nwords,nact,nact,nact,nact)
c----four-particle density matrix dealt with here
c--a,c,e,g    b,d,f,h
      lindice1(i,j,k,l)=lindi(i)+lindj(j)+lindk(k)+l
      call ord8(a,c,e,g,b,d,f,h,i1,i2,i3,i4,i5,i6,i7,i8)
      ind=lindice1(i1,i2,i3,i4)
      eeee=ro4(ind,i5,i6,i7,i8)
      if(b.eq.c)eeee=eeee+d3(a,e,g,d,f,h)
      if(b.eq.e)eeee=eeee+d3(a,c,g,f,d,h)
      if(b.eq.g)eeee=eeee+d3(a,c,e,h,d,f)
      if(d.eq.e)eeee=eeee+d3(a,g,c,b,h,f)
      if(d.eq.g)eeee=eeee+d3(a,c,e,b,h,f)
      if(f.eq.g)eeee=eeee+d3(a,e,c,b,h,d)
      if(b.eq.c)then
         if(d.eq.e)eeee=eeee+d2(a,g,f,h)
         if(d.eq.g)eeee=eeee+d2(a,e,h,f)
         if(f.eq.g)eeee=eeee+d2(a,e,d,h)
      endif
      if(b.eq.e.and.d.eq.g)eeee=eeee+d2(a,c,f,h)
      if(b.eq.e.and.f.eq.g)eeee=eeee+d2(a,c,h,d)
      if(b.eq.g.and.d.eq.e)eeee=eeee+d2(a,c,h,f)
      if(d.eq.e.and.f.eq.g)eeee=eeee+d2(a,c,b,h)
      if(b.eq.c.and.d.eq.e.and.f.eq.g)eeee=eeee+d1(a,h)
      return
      end
c-------------------------------------------------------------------
      function eeeet(a,b,c,d,e,f,g,h,d3,d2,d1,ro4,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f,g,h
      dimension d3(nact,nact,nact,nact,nact,nact),d2(nact,nact,nact,nact
     $     ),d1(nact,nact)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     *     ,normind(18,18,18,18)
      dimension ro4(nwords,nact,nact,nact,nact)
      eeeet=-eeee(a,b,d,c,e,f,g,h,d3,d2,d1,ro4,nact)
      if(c.eq.d)eeeet=eeeet+2.d0*eee2(a,b,e,f,g,h,d3,d2,d1,nact)
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
      function eee2tl(a,b,c,d,e,f,daaa,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d,e,f
      dimension daaa(nact,nact,nact,nact,nact,nact),taa(nact,nact,nact
     $     ,nact),dal(nact,nact)
      eee2tl=-eee2(b,a,c,d,e,f,daaa,taa,dal,nact)
      if(a.eq.b)eee2tl=eee2tl+2.d0*ee2(c,d,e,f,taa,dal,nact)
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
      function ee2t(a,b,c,d,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      ee2t=-ee2(b,a,c,d,taa,dal,nact)
      if(a.eq.b)ee2t=ee2t+2.d0*dal(c,d)
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
c----------------------------------------------------------------------
      function daaat(i,j,k,l,m,n,daaa,daa,da,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daaa(nact,nact,nact,nact,nact,nact),daa(nact,nact,nact
     $     ,nact),da(nact,nact)
      v=-daaa(i,j,k,l,m,n)
      if(k.eq.n)v=v+daat(i,j,l,m,daa,da,nact)
      if(j.eq.n)v=v-daat(i,k,l,m,daa,da,nact)
      if(i.eq.n)v=v+daat(j,k,l,m,daa,da,nact)
      if(m.eq.k)then
         v=v-daa(n,l,j,i)
         if(i.eq.l)v=v+da(n,j)
         if(j.eq.l)v=v-da(n,i)
      endif
      if(m.eq.j)then
         v=v+daa(n,l,k,i)
         if(k.eq.l)v=v+da(n,i)
         if(i.eq.l)v=v-da(n,k)
      endif
      if(m.eq.i)then
         v=v-daa(n,l,k,j)
         if(k.eq.l)v=v-da(n,j)
         if(j.eq.l)v=v+da(n,k)
      endif
      if(k.eq.l)v=v+daa(i,j,m,n)
      if(l.eq.j)v=v-daa(i,k,m,n)
      if(l.eq.i)v=v+daa(j,k,m,n)
      daaat=v
      return
      end
c----------------------------------------------------------------------
      function dbbbt(i,j,k,l,m,n,dbbb,dbb,db,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension dbbb(nact,nact,nact,nact,nact,nact),dbb(nact,nact,nact
     $     ,nact),db(nact,nact)
      v=-dbbb(i,j,k,l,m,n)
      if(k.eq.n)v=v+dbbt(i,j,l,m,dbb,db,nact)
      if(j.eq.n)v=v-dbbt(i,k,l,m,dbb,db,nact)
      if(i.eq.n)v=v+dbbt(j,k,l,m,dbb,db,nact)
      if(m.eq.k)then
         v=v-dbb(n,l,j,i)
         if(i.eq.l)v=v+db(n,j)
         if(j.eq.l)v=v-db(n,i)
      endif
      if(m.eq.j)then
         v=v+dbb(n,l,k,i)
         if(k.eq.l)v=v+db(n,i)
         if(i.eq.l)v=v-db(n,k)
      endif
      if(m.eq.i)then
         v=v-dbb(n,l,k,j)
         if(k.eq.l)v=v-db(n,j)
         if(j.eq.l)v=v+db(n,k)
      endif
      if(k.eq.l)v=v+dbb(i,j,m,n)
      if(l.eq.j)v=v-dbb(i,k,m,n)
      if(l.eq.i)v=v+dbb(j,k,m,n)
      dbbbt=v
      return
      end
c----------------------------------------------------------------------
      function daabt(i,j,k,l,m,n,daab,daa,dab,da,db,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daab(nact,nact,nact,nact,nact,nact),dab(nact,nact,nact
     $     ,nact),daa(nact,nact,nact,nact),da(nact,nact),db(nact,nact)
      v=-daab(i,j,k,l,m,n)
      if(k.eq.n)v=v+daat(i,j,l,m,daa,da,nact)
      if(m.eq.j)v=v+dab(l,n,i,k)
      if(m.eq.i)v=v-dab(l,n,j,k)
      if(l.eq.j)v=v-dab(i,k,m,n)
      if(l.eq.i)v=v+dab(j,k,m,n)
      if(m.eq.j.and.i.eq.l)v=v-db(n,k)
      if(m.eq.i.and.j.eq.l)v=v+db(n,k)
      daabt=v
      return
      end
c----------------------------------------------------------------------
      function dbbat(i,j,k,l,m,n,dbba,dbb,dab,da,db,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension dbba(nact,nact,nact,nact,nact,nact),dab(nact,nact,nact
     $     ,nact),dbb(nact,nact,nact,nact),da(nact,nact),db(nact,nact)
      dbbat=v
      v=-dbba(i,j,k,l,m,n)
      if(k.eq.n)v=v+dbbt(i,j,l,m,dbb,db,nact)
      if(m.eq.j)v=v+dab(n,l,k,i)
      if(m.eq.i)v=v-dab(n,l,k,j)
      if(l.eq.j)v=v-dab(k,i,n,m)
      if(l.eq.i)v=v+dab(k,j,n,m)
      if(m.eq.j.and.i.eq.l)v=v-da(n,k)
      if(m.eq.i.and.j.eq.l)v=v+da(n,k)
      dbbat=v
      return
      end
c--------------------------------------------------------------------------
      function daat(a,b,c,d,daa,da,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension daa(nact,nact,nact,nact),da(nact,nact)
      v=daa(a,b,c,d)
      if(b.eq.d)v=v+dat(c,a,da,nact)
      if(a.eq.d)v=v-dat(c,b,da,nact)
      if(b.eq.c)v=v+da(a,d)
      if(a.eq.c)v=v-da(b,d)
      daat=v
      return
      end
c--------------------------------------------------------------------------
      function dbbt(a,b,c,d,dbb,db,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension dbb(nact,nact,nact,nact),db(nact,nact)
      v=dbb(a,b,c,d)
      if(b.eq.d)v=v+dbt(c,a,db,nact)
      if(a.eq.d)v=v-dbt(c,b,db,nact)
      if(b.eq.c)v=v+db(a,d)
      if(a.eq.c)v=v-db(b,d)
      dbbt=v
      return
      end
c--------------------------------------------------------------------------
      function dabt(a,b,c,d,dab,da,db,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,d
      dimension dab(nact,nact,nact,nact),da(nact,nact),db(nact,nact)
      v=dab(a,b,c,d)
      if(b.eq.d)v=v+dat(c,a,da,nact)
      if(a.eq.c)v=v-db(b,d)
      dabt=v
      return
      end
c------------------------------------------------------------------------
      function dat(a,b,da,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b
      dimension da(nact,nact)
      v=-da(a,b)
      if(a.eq.b)v=v+1.0d0
      dat=v
      return
      end
c------------------------------------------------------------------------
      function dbt(a,b,db,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b
      dimension db(nact,nact)
      v=-db(a,b)
      if(a.eq.b)v=v+1.0d0
      dbt=v
      return
      end
c------------------------------------------------------------------------
      subroutine koopE(atwo,koopipa,koopeaa,taa,dal,factor,nact,ncore,f
     $     ,aj)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      dimension koopipa(nact,nact),koopeaa(nact,nact)
      dimension atwo(nact,nact,nact,nact)
      real*8 koopipa,koopipb,koopeaa,koopeab
      dimension f(*),aj(*)
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
c      dimension nos(20)
      integer a,b,ac,bc,c,ap,apc,cc,d,dc
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
C
C     Calcolo le matrici di Koopmans
C     Caso -1: matrici di Koopmans per la ionizzazione attiva
C
C     La parte relativa alle singole differenze attive
C     e alle doppie differenze, una attiva e una inattiva
C     viene calcolata mediante rho_ab
C
c----building Dyall's h eff----------
      do a=1,nact
cele
       ac=a+ncore
         do b=a,nact
cele
          bc=b+ncore
            jab=indice(ac,bc)
c            aj(jab)=aj(jab)*factor
            f(jab)=aj(jab)
         enddo
      enddo
      do a=1,nact
         ac=a+ncore
         do ap=1,nact
            apc=ap+ncore
            koopipa(ap,a)=0.d0
c            print*,'indici a e ap',a,ap
            do c=1,nact
               cc=c+ncore
               irind=indice(ac,cc)
               dum=f(irind)
               koopipa(ap,a)=koopipa(ap,a)-dum*dal(ap,c)
c               print*,'indice c=',c
c               print*,'1: ',dum/factor,dal(ap,c),koopipa(ap,a)/factor
               do d=1,nact
                  do ie=1,nact
                     koopipa(ap,a)=koopipa(ap,a)-atwo(c,ie,a,d)*taa(ap,c
     $                    ,d,ie)
c                     print*,'indici d, ie',d,ie
c                     print*,'2: ',atwo(c,ie,a,d)/factor,taa(ap,c,d,ie)
c     $                    ,koopipa(ap,a)/factor
                  enddo
               enddo
            enddo
            koopipa(ap,a)=2.d0*koopipa(ap,a)
         enddo
      enddo
      do ap=1,nact
         do a=1,nact
            ind=indice(ap+ncore,a+ncore)
            dum=2.d0*f(ind)
            do d=1,nact
               do ie=1,nact
                  dum=dum+(2.d0*atwo(d,ie,ap,a)-atwo(ap,ie,d,a))*dal(d
     $                 ,ie)
               enddo
            enddo
            koopeaa(ap,a)=koopipa(ap,a)+2.d0*dum
         enddo
      enddo
      return
      end
c-----------------------------------------------------
      subroutine zeroe(a,n)
      real*8 a(*)
      do i=1,n
         a(i)=0.d0
      enddo
      return
      end
c------------------------------------------------------
      subroutine ord31(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nact,
     * m,n)
      implicit none
      integer segno,segno2,op,pp,op2,pp2,m,n,mp2,np2,mp,np
      integer ia,ic,ie,ig,ia2,ic2,ie2,ig2,nact,ivolte,jvolte
      integer i1,i2,i3,i4,i5,i6,i7,i8,norm1,norm2,nwords,normind
      integer lindi,lindj,lindk
      real*8 amat(nwords,nact,nact,nact,nact),val,vals
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     * ,normind(18,18,18,18)

c      call getnorm(ia,ic,ie,ig,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia,ic,ie,ig)
      norm2=normind(ia2,ic2,ie2,ig2)
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
c      write (6,'(''ord31 1'',8I3)') mp,np,op,pp
        do jvolte=1,6   
        call permut123(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !aaba
      call permut124(i1,i2,i4,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
c      write (6,'(''ord31 2'',8I3)') mp,np,op,pp
        do jvolte=1,6   
        call permut124(i5,i6,i8,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !abaa
      call permut134(i1,i4,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
c      write (6,'(''ord31 3'',8I3)') mp,np,op,pp
        do jvolte=1,6   
        call permut134(i5,i8,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !baaa
      call permut234(i4,i1,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
c      write (6,'(''ord31 4'',8I3)') mp,np,op,pp
        do jvolte=1,6   
        call permut234(i8,i5,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
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
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !aaba
      call permut124(i1,i2,i4,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut124(i5,i6,i8,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !abaa
      call permut134(i1,i4,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut134(i5,i8,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,6  !baaa
      call permut234(i4,i1,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut234(i8,i5,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      return
      end
c-----------------------------------------------------------
      subroutine ord22(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nact,
     * m,n)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     * ,normind(18,18,18,18)
      real*8 amat(nwords,nact,nact,nact,nact),val,vals
      integer segno,segno2,op,pp,op2,pp2,m,n

c      write (6,'(''ord22 in'',8I3)') ia,ic,ie,ig,ia2,ic2,ie2,ig2     
c      call getnorm(ia,ic,ie,ig,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
c      write(6,*) 'norme',norm1,norm2
      norm1=normind(ia,ic,ie,ig)
      norm2=normind(ia2,ic2,ie2,ig2)
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
c      write (6,'(''ord22 in'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8

      do ivolte=1,4  !aabb
      call permut1234(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 1'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abab
      call permut1324(i1,i3,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 2'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i5,i7,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abba
      call permut1423(i1,i3,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 3'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i5,i7,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baab
      call permut1423(i3,i1,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 4'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i7,i5,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baba
      call permut1324(i3,i1,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 5'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i7,i5,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !bbaa
      call permut1234(i3,i4,i1,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 6'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i7,i8,i5,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
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

c      write (6,'(''ord22 in'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8

      do ivolte=1,4  !aabb
      call permut1234(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 1'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abab
      call permut1324(i1,i3,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 2'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i5,i7,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !abba
      call permut1423(i1,i3,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 3'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i5,i7,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baab
      call permut1423(i3,i1,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 4'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i7,i5,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !baba
      call permut1324(i3,i1,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 5'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i7,i5,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      do ivolte=1,4  !bbaa
      call permut1234(i3,i4,i1,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 6'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i7,i8,i5,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
        enddo
      endif
      enddo

      return
      end
c------------------------------------------------------
      subroutine ordsame(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nact,
     * m,n)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     * ,normind(18,18,18,18)
      real*8 amat(nwords,nact,nact,nact,nact),val,vals
      integer segno,segno2,op,pp,op2,pp2,n,m

c      call getnorm(ia,ic,ie,ig,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia,ic,ie,ig)
      norm2=normind(ia2,ic2,ie2,ig2)
c      write (6,*) 'calling ordsame'
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
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
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
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nact)
	enddo
      endif
      enddo


      return
      end
c-----------------------------------------------------
      subroutine ord9(ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2,
     * i1,i2,i3,i4,i5,i6,i7,i8)
c---  ordina in maniera decrescente ia1, ic1, ie1, ig1 in i1 i2 i3 i4
c---  ed fa la stessa perm anche su ia2, ic2, ie2, ig2 in i5 i6 i7 i8
      dimension nv1(4),nv2(4)
      common /actspace/nact,nact2,nact3,nact4
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     * ,normind(18,18,18,18)
c      write (6,'(''ord8 in'',8I3)') ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2

c      call getnorm(ia1,ic1,ie1,ig1,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia1,ic1,ie1,ig1)
      norm2=normind(ia2,ic2,ie2,ig2)


      if (norm1.ge.norm2) then
      nv1(1)=ia1
      nv1(2)=ic1
      nv1(3)=ie1
      nv1(4)=ig1
      nv2(1)=ia2
      nv2(2)=ic2
      nv2(3)=ie2
      nv2(4)=ig2
      else
      nv1(1)=ia2
      nv1(2)=ic2
      nv1(3)=ie2
      nv1(4)=ig2
      nv2(1)=ia1
      nv2(2)=ic1
      nv2(3)=ie1
      nv2(4)=ig1
      endif
      do i=1,4
         imax=nv1(i)
        do j=i+1,4
            if(nv1(j).gt.imax)then
               n=imax
               imax=nv1(j)
               nv1(j)=n
	       n=nv2(i)
	       nv2(i)=nv2(j)
	       nv2(j)=n
            endif
            nv1(i)=imax
         enddo
      enddo
      i1=nv1(1)
      i2=nv1(2)
      i3=nv1(3)
      i4=nv1(4)
      i5=nv2(1)
      i6=nv2(2)
      i7=nv2(3)
      i8=nv2(4)
c      write (6,'(''ord8 ou'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8         
      return
      end
c-----------------------------------------------------
      subroutine ord8(ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2,
     * i1,i2,i3,i4,i5,i6,i7,i8)
c---  ordina in maniera decrescente ia1, ic1, ie1, ig1 in i1 i2 i3 i4
c---  ed fa la stessa perm anche su ia2, ic2, ie2, ig2 in i5 i6 i7 i8
      dimension nv1(4),nv2(4)
      common /actspace/nact,nact2,nact3,nact4
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     * ,normind(18,18,18,18)
c      write (6,'(''ord8 in'',8I3)') ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2

c      call getnorm(ia1,ic1,ie1,ig1,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia1,ic1,ie1,ig1)
      norm2=normind(ia2,ic2,ie2,ig2)


      if (norm1.ge.norm2) then
      nv1(1)=ia1
      nv1(2)=ic1
      nv1(3)=ie1
      nv1(4)=ig1
      nv2(1)=ia2
      nv2(2)=ic2
      nv2(3)=ie2
      nv2(4)=ig2
      else
      nv1(1)=ia2
      nv1(2)=ic2
      nv1(3)=ie2
      nv1(4)=ig2
      nv2(1)=ia1
      nv2(2)=ic1
      nv2(3)=ie1
      nv2(4)=ig1
      endif
      do i=1,4
         imax=nv1(i)
	 jmax=i
         do j=i+1,4
            if(nv1(j).gt.imax)then
               imax=nv1(j)
               jmax=j  
            endif
         enddo
         if (jmax.ne.i) then
         ndum=nv1(i)
         nv1(i)=nv1(jmax)
         nv1(jmax)=ndum
         ndum=nv2(i)
         nv2(i)=nv2(jmax)
         nv2(jmax)=ndum
         endif
c      write (6,'(''ord8 ou'',i6,3x,8I3)') i,nv1(1),nv1(2),nv1(3),nv1(4),
c     *                                      nv2(1),nv2(2),nv2(3),nv2(4)
      enddo
      i1=nv1(1)
      i2=nv1(2)
      i3=nv1(3)
      i4=nv1(4)
      i5=nv2(1)
      i6=nv2(2)
      i7=nv2(3)
      i8=nv2(4)
c      write (6,'(''ord8 ou'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8         
      return
      end
c-----------------------------------------------------
      subroutine getnorm(ia,ic,ie,ig,norm)
      implicit none
      dimension nv(4)
      integer ia,ic,ie,ig,norm,nv,i,j,jmax,imax,ndum

c      norm=ia+ic+ie+ig
c      return

      nv(1)=ia
      nv(2)=ic
      nv(3)=ie
      nv(4)=ig
      

c      write (6,'(''norm in'',8I3)') ia,ic,ie,ig

      do i=1,4
       imax=nv(i)
       jmax=i
        do j=i+1,4
	if(nv(j).gt.imax) then
	imax=nv(j)
	jmax=j
	endif
	enddo
	if(jmax.ne.i) then
	ndum=nv(i)
	nv(i)=nv(jmax)
	nv(jmax)=ndum
	endif
      enddo
c      write (6,'(''norm ou'',8I3)') nv(1),nv(2),nv(3),nv(4)

      norm=(nv(1)-1)*nv(1)*(nv(1)+1)*(nv(1)+2)/24+
     *     (nv(2)-1)*nv(2)*(nv(2)+1)/6+
     *     nv(3)*(nv(3)-1)/2+nv(4)

      return
      end
c-------------------------------------------------------------------
      function ro3t(a,b,c,ap,bp,cp,d3,d2,d1,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,ap,bp,cp
      dimension d3(nact,nact,nact,nact,nact,nact),d2(nact,nact,nact,nact
     $     ),d1(nact,nact)
      ro3t=-d3(ap,bp,cp,a,b,c)
      if(ap.eq.a)ro3t=ro3t+2.d0*d2(bp,cp,b,c)
      if(bp.eq.b)ro3t=ro3t+2.d0*d2(ap,cp,a,c)
      if(cp.eq.c)ro3t=ro3t+2.d0*d2(ap,bp,a,b)
      if(bp.eq.c)ro3t=ro3t-d2(ap,cp,a,b)
      if(ap.eq.c)ro3t=ro3t-d2(cp,bp,a,b)
      if(cp.eq.b)ro3t=ro3t-d2(bp,ap,c,a)
      if(ap.eq.b)ro3t=ro3t-d2(bp,cp,a,c)
      if(cp.eq.a)ro3t=ro3t-d2(ap,bp,c,b)
      if(bp.eq.a)ro3t=ro3t-d2(ap,cp,b,c)
      if(cp.eq.c.and.bp.eq.b)ro3t=ro3t-4.d0*d1(ap,a)
      if(cp.eq.b.and.bp.eq.c)ro3t=ro3t+2.d0*d1(ap,a)
      if(cp.eq.c.and.ap.eq.a)ro3t=ro3t-4.d0*d1(bp,b)
      if(cp.eq.a.and.ap.eq.c)ro3t=ro3t+2.d0*d1(bp,b)
      if(ap.eq.a.and.bp.eq.b)ro3t=ro3t-4.d0*d1(cp,c)
      if(ap.eq.b.and.bp.eq.a)ro3t=ro3t+2.d0*d1(cp,c)
      if(cp.eq.c.and.ap.eq.b)ro3t=ro3t+2.d0*d1(bp,a)
      if(cp.eq.b.and.ap.eq.c)ro3t=ro3t-d1(bp,a)
      if(cp.eq.c.and.bp.eq.a)ro3t=ro3t+2.d0*d1(ap,b)
      if(cp.eq.a.and.bp.eq.c)ro3t=ro3t-d1(ap,b)
      if(bp.eq.b.and.ap.eq.c)ro3t=ro3t+2.d0*d1(cp,a)
      if(bp.eq.c.and.ap.eq.b)ro3t=ro3t-d1(cp,a)
      if(bp.eq.b.and.cp.eq.a)ro3t=ro3t+2.d0*d1(ap,c)
      if(bp.eq.a.and.cp.eq.b)ro3t=ro3t-d1(ap,c)
      if(ap.eq.a.and.cp.eq.b)ro3t=ro3t+2.d0*d1(bp,c)
      if(ap.eq.b.and.cp.eq.a)ro3t=ro3t-d1(bp,c)
      if(ap.eq.a.and.bp.eq.c)ro3t=ro3t+2.d0*d1(cp,b)
      if(ap.eq.c.and.bp.eq.a)ro3t=ro3t-d1(cp,b)
      if(ap.eq.a.and.bp.eq.b.and.cp.eq.c)ro3t=ro3t+8.d0
      if(ap.eq.b.and.bp.eq.a.and.cp.eq.c)ro3t=ro3t-4.d0
      if(ap.eq.a.and.bp.eq.c.and.cp.eq.b)ro3t=ro3t-4.d0
      if(ap.eq.c.and.bp.eq.b.and.cp.eq.a)ro3t=ro3t-4.d0
      if(ap.eq.b.and.bp.eq.c.and.cp.eq.a)ro3t=ro3t+2.d0
      if(ap.eq.c.and.bp.eq.a.and.cp.eq.b)ro3t=ro3t+2.d0
      return
      end
c-------------------------------------------------------------------
      function ro2t(ap,bp,a,b,d2,d1,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      integer a,b,c,ap,bp,cp
      dimension d2(nact,nact,nact,nact),d1(nact,nact)
      ro2t=d2(a,b,ap,bp)
      if(ap.eq.a)ro2t=ro2t-2.d0*d1(b,bp)
      if(bp.eq.b)ro2t=ro2t-2.d0*d1(a,ap)
      if(ap.eq.b)ro2t=ro2t+d1(a,bp)
      if(bp.eq.a)ro2t=ro2t+d1(b,ap)
      if(ap.eq.a.and.bp.eq.b)ro2t=ro2t+4.d0
      if(ap.eq.b.and.bp.eq.a)ro2t=ro2t-2.d0
      return
      end
c-----------------------------------------------------------------
      subroutine bro3(ro4,daaa,nact,nele)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daaa(nact,nact,nact,nact,nact,nact)
      common/lndc/nwords,lindi(18),lindj(18),lindk(18)
     *     ,normind(18,18,18,18)
      dimension ro4(nwords,nact,nact,nact,nact)
cc---controllo su contrazione della ro4
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do ii=1,nact
                  do jj=1,nact
                     do kk=1,nact
                        dum4=0.d0
                        do ll=1,nact
			call ord8(i,j,k,ll,ii,jj,kk,ll,
     $                            i1,i2,i3,i4,i5,i6,i7,i8)
			ind=lindice(i1,i2,i3,i4)
                        dum4=dum4+ro4(ind,i5,i6,i7,i8)
                        enddo
                        daaa(i,j,k,ii,jj,kk)=dum4/(nele-3)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------
      subroutine bro2(daaa,taa,nact,nele)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension daaa(nact,nact,nact,nact,nact,nact)
      dimension taa(nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            do ii=1,nact
               do jj=1,nact
                  dum4=0.d0
                  do ll=1,nact
                     dum4=dum4+daaa(i,j,ll,ii,jj,ll)
                  enddo
                  taa(i,j,ii,jj)=dum4/(nele-2)
               enddo
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------
      subroutine bro1(taa,dal,nact,nele)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      do i=1,nact
         do j=1,nact
            dum4=0.d0
            do ll=1,nact
               dum4=dum4+taa(i,ll,j,ll)
            enddo
            dal(i,j)=dum4/(nele-1)
         enddo
      enddo
      return
      end
c------------------------------------------------------------------------
      subroutine wr1(dal,coopeaa,coopipa,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension dal(nact,nact),coopeaa(nact,nact)
     $     ,coopipa(nact,nact)
      write(31) nact
      write(31)((dal(k,l),k=1,nact),l=1,nact)
      write(31)((coopeaa(k,l),k=1,nact),l=1,nact)
      write(31)((coopipa(k,l),k=1,nact),l=1,nact)
      return
      end
c------------------------------------------------------------------------
      subroutine wr2(koopaa,taa,dal,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      real*8 dal(nact,nact),koopaa(nact,nact,nact,nact
     $     ),taa(nact,nact,nact,nact)
c      print*,'costruzione matrici ro2 di buca'
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  adum=ro2t(i,j,k,l,taa,dal,nact)
                  koopaa(i,j,k,l)=adum
               enddo
            enddo
         enddo
      enddo
c      print*,'scrittura matrici densita di buca su file 31'
      do i=1,nact
         do j=1,nact
            write(31)((koopaa(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
      return
      end
c------------------------------------------------------------------------
      subroutine wr3(koopeaa,taa,koopaa,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      real*8 koopaa(nact,nact,nact,nact
     $     )
     $     ,taa(nact,nact,nact,nact),koopeaa(nact,nact,nact,nact)
      do i=1,nact
         do j=1,nact
            write(31)((koopeaa(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
c-------scrittura matrici densita a due particelle
      do i=1,nact
         do j=1,nact
            write(31)((taa(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
c-------scrittura matrici koopman doppia ionizzazione
      do i=1,nact
         do j=1,nact
            write(31)((koopaa(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine wr4(koopaa,koopbb,koopab,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      real*8 koopaa(nact,nact,nact,nact
     $     ),koopbb(nact,nact,nact,nact),koopab(nact,nact)
      do i=1,nact
         do j=1,nact
            write(31)((koopaa(i,j,k,l),k=1,nact),l=1,nact)
            write(31)((koopbb(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
      write(31)((koopab(k,l),k=1,nact),l=1,nact)
      return
      end
c------------------------------------------------------------------
      subroutine fdiff(bmat,cmat,btmat,ctmat,factor,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension bmat(nact,nact,nact,nact),cmat(nact,nact,nact,nact)
      dimension btmat(nact,nact,nact,nact),ctmat(nact,nact,nact,nact)
      write (6,*) 
      write (6,*) 'Finding differences in bmat and cmat'
      bigdiff=0.d0
      do i =1,nact
         do j=1,nact
            do ii=1,nact
               do jj=1,nact
                  val=bmat(i,j,ii,jj)/factor
                  val2=cmat(jj,i,j,ii)/factor
c                  if(abs(val-val2).gt.1.d-8)then
c                     print '(a,4i3,a,f12.8,a,f12.8)','indici',i,j,ii,jj,
c     $                    'bmat=',val,' cmat=',val2
c                  endif
                  pdiff=abs(val-val2)
                  if(pdiff.gt.bigdiff)bigdiff=pdiff
               enddo
            enddo
         enddo
      enddo
      write(6,'(a,a,d20.10)')' Largest difference between bmat',
     $     ' and cmat  =',bigdiff
      write (6,*) 
      write (6,*) 'Finding differences in btmat and ctmat'
      bigdiff=0.d0
      do i =1,nact
         do j=1,nact
            do ii=1,nact
               do jj=1,nact
                  val=btmat(i,j,ii,jj)/factor
                  val2=ctmat(jj,i,j,ii)/factor
c                  if(abs(val-val2).gt.1.d-8)then
c                     print '(a,4i3,a,f12.8,a,f12.8)','indici',i,j,ii,jj,
c     $                    'btmat=',val,' ctmat=',val2
c                  endif
                  pdiff=abs(val-val2)
                  if(pdiff.gt.bigdiff)bigdiff=pdiff
               enddo
            enddo
         enddo
      enddo
      write(6,'(a,a,d20.10)')' Largest difference between btmat and',
     $     ' ctmat=',bigdiff
      call flush(6)
      return
      end
c------------------------------------------------------------------------------
      subroutine wr5(daaa,amat,atmat,bmat,btmat,cmat,ctmat,dmat,dtmat
     $     ,nact)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension bmat(nact,nact,nact,nact),cmat(nact,nact,nact,nact)
      dimension amat(nact,nact,nact,nact,nact,nact),dmat(nact,nact)
      dimension btmat(nact,nact,nact,nact),ctmat(nact,nact,nact,nact)
      dimension atmat(nact,nact,nact,nact,nact,nact),dtmat(nact,nact)
      dimension daaa(nact,nact,nact,nact,nact,nact)
c-----scrittura matrice spinless a 3 part e matrice amat
c      print*,'wr5: nact=',nact
c      call flush(6)
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  write(31)((daaa(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
                  write(31)((amat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
               enddo
            enddo
         enddo
      enddo
c      print*,'done one'
c      call flush(6)
c----scrittura matrici bmat e cmat
      do i=1,nact
         do j=1,nact
            write(31)((bmat(i,j,k,l),k=1,nact),l=1,nact)
            write(31)((cmat(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
c      print*,'done two'
c      call flush(6)
c----scrittura matrice dmat
      write(31)((dmat(i,j),i=1,nact),j=1,nact)
c      print*,'done three'
c      call flush(6)
c-----scrittura matrice atmat
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  write(31)((atmat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
               enddo
            enddo
         enddo
      enddo
c      print*,'done four'
c      call flush(6)
c----scrittura matrici btmat e ctmat
      do i=1,nact
         do j=1,nact
            write(31)((btmat(i,j,k,l),k=1,nact),l=1,nact)
            write(31)((ctmat(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
c      print*,'done five'
c      call flush(6)
c----scrittura matrice dtmat
      write(31)((dtmat(i,j),i=1,nact),j=1,nact)
c      print*,'done six'
c      call flush(6)
      return
      end
c-----------------------------------------------------
      subroutine wr6(ro4,nact,nwords)
      implicit real*8 (a-h,o-y),logical*1 (z)
      dimension ro4(nwords,nact,nact,nact,nact)
c-----lettura matrice spinless a 4 part
      print*,'wr6: nact=',nact
      call flush(6)
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  write (31) (ro4(nwo,i,j,k,l),nwo=1,nwords)
               enddo
            enddo
         enddo
      enddo
      return
      end
c------------------------------------------------------------
      SUBROUTINE REIJKL(factor,norb,nsym,itsym,its,nos,jd1,indic,jndic
     $     ,lndic,num)
C 
C     PREPARATION POUR LE STOCKAGE DES INTEGRALES BIELECTONIQUES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      DIMENSION NOS(20),ITSYM(*),ITS(20,20)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * NDEB(500),NBO(99),
     * NAD(2000),KT(2000)
      dimension INDIC(*),JNDIC(*),LNDIC(*),num(*)
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      common/caratte/cara1,cara2
      character*1 cara1,cara2
      dimension norbsym(8)
      common /toul/ ztoul
c      cara1=char(-1)
      cara2=char(0)
cele 08-12-04
      zfac=.true.
cele 08-12-04
#ifdef INTEL
c      print*,'REIJKL compiled with intel ifc compiler'
#endif
c      if (id1.ne.jd1) then
c      write (6,*) '   ATTENTION: LA DIMENSION ID1 N''EST PAS LA MEME'
c      write (6,*) '   DANS REIJKL ET DANS LE PROGRAMME PRINCIPAL'
c      stop 12
c      endif
c      k28=id2*8
c      if (id1*20*2.gt.k28) then
c      write (6,*) '   ERREUR DIMENSIONS DANS REIJKL (EQUIV. F-NOMS)'
c      stop 12
c      endif
      rewind 25
      READ (25) NSYM,NORB,NOC,NDEGEN,(ITSYM(I),I=1,NORB)
      READ (25) ((ITS(I,J),I=1,NSYM),J=1,NSYM)
      IJ=0
      NUM(1)=0
      DO 1 I=2,norb
      DO 1 J=1,I
      IJ=IJ+1
c      write (6,*) 'num',i,j,ij,NUM(IJ)+IJ
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
      if (zhandy) then
         if(.not.zthree)then
            ni2=nijkl
            nr4=25534
            ni4=1
            ni3=1
            WRITE (6,74) NIJKL,2
         else
            ni3=nijkl
            nr4=1
            ni4=1
            ni2=1
            write (6,74) nijkl,3
            if (zfac) then
               if (ztoul) then
                  call intmaxt(bigint,norb,nsym,nos)
               else
                  call intmax(bigint,norb,nsym,nos)
               endif
               print'(''     MAXIMUM INTEGRAL IS '',D10.4)',bigint
               factor=8388607.d0/bigint
            endif
           print'(''     FACTOR WILL BE '',D10.4)',factor
         endif
      else
         ni2=1
         nr4=1
         ni3=1
         ni4=nijkl
         WRITE (6,74) NIJKL,4
         call flush(6)
         if (zfac) then
            if (ztoul) then
               call intmaxt(bigint,norb,nsym,nos)
            else
            call intmax(bigint,norb,nsym,nos)
         endif
            print'(''     MAXIMUM INTEGRAL IS '',D10.4)',bigint
            call flush(6)
            factor=2147483646.d0/bigint
         endif
         print'(''     FACTOR WILL BE '',D10.4)',factor
         call flush(6)
      endif
   74 FORMAT (//5X,'DIMENSION OF THE IJKL TABLE:',I10,' *',i1,' BYTES')
      RETURN
      END
C***********************************************************************
      SUBROUTINE intmax(bigint,norb,nsym,nos)
C 
C     trova l'integrale piu' grande
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      PARAMETER(MBUF=3255,NBUF=4096)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
c      COMMON /INT/ NI4,NI3,NI2,NR4,
c     * NDEB(500),NBO(99),
c     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
c      INTEGER*2 INDIC,JNDIC,LNDIC,NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      bigint=0.d0
c      k22=id2*2
c      if (nbuf.gt.k22) then
c      write (6,*) '   *** ERREUR DIMENSIONS DANS INTMAX ***'
c      stop
c      endif
c      WRITE (6,1000)
c 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
      call filesplit('READ',50,8,ichange,idum2)    
      READ (50) FACDUM
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
c-----------------------------------------------------
      SUBROUTINE intmaxt(bigint,norb,nsym,nos)
C 
C     trova l'integrale piu' grande
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      PARAMETER(MBUF=3255,NBUF=4096)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
      bigint=0.d0
c      WRITE (6,1000)
c 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
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
      SUBROUTINE ITIJKL(factor,norb,nsym,nos,ijkl4,
     * ijkl2,rijkl,ijkl3,indic,jndic,lndic,num)
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
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * NDEB(500),NBO(99),
     * NAD(2000),KT(2000)
      dimension INDIC(*),JNDIC(*),LNDIC(*),num(*)
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      idummy=0
c      k22=id2*2
c      if (mbuf.gt.k22) then
c      write (6,*) '   *** ERREUR DIMENSIONS DANS ITIJKL ***'
c      stop
c      endif
      WRITE (6,1000)
 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
c      call getenv('HOSTTYPE',host)
      host='LINUX'
      write (6,'(a,a)') '    HOSTTYPE=',host
      write (6,*) '    INTEGRALI DA DALTON + IJKLPI'
      do 500 k=1,ni2
 500     ijkl2(k)=0
      do 501 k=1,ni4
 501     ijkl4(k)=0
      do k=1,ni3
         ijkl3(k)=cara(2:4)
      enddo
      khan=0
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
      call filesplit('READ',50,8,ichange,idum2)    
      READ (50) FACDUM
      WRITE (6,10) FACTOR
   10 FORMAT (5X,'INTEGRAL PREFACTOR =',D15.2)
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
c      write (6,*) '0',xx(n),i,j,k,l
      CALL HANDY(NA,A8,KHAN,ZHANDY,zthree,ijkl4,ijkl3,ijkl2,rijkl)
c     CALL HANDY(NA,A8,KHAN,ZHANDY)
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
c-------------------------------------------------------------
      SUBROUTINE ITIJKLt(factor,norb,nsym,nos,ijkl4,
     * ijkl2,rijkl,ijkl3,indic,jndic,lndic,num)
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
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * NDEB(500),NBO(99),
     * NAD(2000),KT(2000)
      dimension INDIC(*),JNDIC(*),LNDIC(*),num(*)
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
c_toul
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
c_toul
      idummy=0
c      WRITE (6,1000)
c 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
c      call getenv('HOSTTYPE',host)
c      host='LINUX'
c      write (6,'(a,a)') '     HOSTTYPE=',host
c      write (6,*) '    INTEGRALI LETTI CON IL FORMALISMO DI TOLOSA'
      do 500 k=1,ni2
 500     ijkl2(k)=0
      do 501 k=1,ni4
 501     ijkl4(k)=0
      do k=1,ni3
         ijkl3(k)=cara(2:4)
      enddo
c
c      NUM(1)=0
c      DO 1 I=2,norb
c      DO 1 J=1,I
c      IJ=IJ+1
c      write (6,*) 'num',i,j,ij,NUM(IJ)+IJ
c    1 NUM(IJ+1)=NUM(IJ)+IJ

      khan=0
      REWIND 50
      call filesplit('REWIND',50,0,idum,idum2)
c      WRITE (6,10) FACTOR
c   10 FORMAT (5X,'INTEGRAL PREFACTOR =',D15.2)
      NNT=NUM(NORB)+NORB
      NTT=NDEB(NSYM)+(NOS(NSYM)*(NOS(NSYM)+1))/2
      M=0
      N=1
c--renzo modif to read Daniel's code for norb >255
   5  continue
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
      FUNCTION AI(I,J,K,L,indic,jndic,lndic,num,ijkl4)
C
C     REPERAGE DE L'INTEGRALE BIELECTRONIQUE IJKL DANS LE TABLEAU IJKL2
C
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      pointer f,aj,ak,num,indic,jndic,lndic
c      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
c      common/fajk/f,aj,ak,num
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * NDEB(500),NBO(99),
     * NAD(2000),KT(2000)
      dimension INDIC(*),JNDIC(*),LNDIC(*),num(*),ijkl4(*)
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM1,LOCM2,LOCM3,LOCM4,LOCM(6)
      integer pick3
c      common /big/ x,ijkl2,ijkl4,ijkl3,rijkl
c      pointer x
c      dimension x(:)
c      pointer ijkl2,ijkl3,ijkl4
c      pointer rijkl
c      dimension ijkl2(:),ijkl4(:),ijkl3(:),rijkl(:)
c      integer*2 ijkl2
c      character*3,ijkl3
c      real*4 rijkl
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
c      if (zhandy) then
c         if(.not.zthree)then
c      AI=PICK2(nnjk,ijkl2,rijkl)
c      else
c         ai=pick3(nnjk,ijkl3)
c         endif
c      else
      AI=PICK4(nnjk,ijkl4)
c      endif
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
c-----------------------------------------------------------
      subroutine giveocc(iocc,m,nd,ne,trou,part,nocc,norb)
      integer*2 ne(*),trou(*),part(*)
      integer*1 iocc(*)
      integer nd(*)
      do i=1,nocc
         iocc(i)=2
      enddo
      do i=nocc+1,norb
         iocc(i)=0
      enddo
      nstart=nd(m)
      do i=1,ne(m)
         ib=trou(nstart+i)
         ip=part(nstart+i)
         if(ib.gt.norb)ib=ib-norb
         if(ip.gt.norb)ip=ip-norb
         iocc(ib)=iocc(ib)-1
	 if (ip.le.norb) then
         iocc(ip)=iocc(ip)+1
	 endif
      enddo
      return
      end
c------------------------------------------------------------------
      integer function ndiff(iocc,jocc,norb)
      integer*1 iocc(*),jocc(*)
      ndiff=0
      do i=1,norb
         ia=iocc(i)-jocc(i)
         if(ia.ne.0)ndiff=ndiff+abs(ia)
         if(ndiff.gt.8)return
      enddo
      return
      end
c--------------------------------------------------------
      subroutine esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,nocc
     $     ,norb,nspin)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      integer*2 ne,trou,part,icomp,iconf,isegno
      dimension nd(*),ne(*),trou(*),part(*),icomp(*),iconf(*),isegno(*)
      equivalence (ip1,ib1),(ip2,ib2)
      allocatable iocc(:)
      integer*1 iocc
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
c      print*,'esclass found ',ncapos,' capostipiti'
c      call flush(6)
      return
      end
c------------------------------------------------------------------------
      SUBROUTINE ESAME(flip,ne,nd,trou,part,ispin,iorb,ncf)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      integer*2 ne,trou,part,ispin,iorb
      integer*1 flip
      dimension flip(*),ne(*),nd(*),trou(*),part(*),ispin(*),iorb(*)
      DIMENSION IT(12),IP(12)
      DO2M=1,NCF
      IF(M.EQ.1)GOTO5
      IF(FLIP(M-1).EQ.2)GOTO22
    5 NEC=NE(M)
      IF(NEC.EQ.0)GOTO3
      NAD=ND(M)
      ISP=0
      DO10I=1,NEC
      IT(I)=TROU(NAD+I)
      IP(I)=PART(NAD+I)
      IF(ISPIN(IT(I)).EQ.0)ISP=ISP-1
      IF(ISPIN(IT(I)).EQ.1)ISP=ISP+1
      IF(ISPIN(IP(I)).EQ.0)ISP=ISP+1
      IF(ISPIN(IP(I)).EQ.1)ISP=ISP-1
   10 CONTINUE
      IF(ISP.NE.0)RETURN
      IF((NEC/2)*2.NE.NEC)GOTO21
      DO20I=1,NEC
      ZREV1=.TRUE.
      ZREV2=.TRUE.
      IIT=IORB(IT(I))
      IIP=IORB(IP(I))
      DO30J=1,NEC
      IF(J.EQ.I)GOTO30
      JJT=IORB(IT(J))
      JJP=IORB(IP(J))
      IF(JJT.EQ.IIT)ZREV1=.FALSE.
      IF(JJP.EQ.IIP)ZREV2=.FALSE.
   30 CONTINUE
      ZREV=ZREV1.OR.ZREV2
      IF(ZREV)GOTO21
   20 CONTINUE
      FLIP(M)=0
      GOTO2
   21 FLIP(M)=2
      GOTO2
   22 FLIP(M)=1
      GOTO2
    3 FLIP(M)=0
    2 CONTINUE
      RETURN
      END
c---------------------------------------------------------------------
      integer function noverk(n,k)
      if(k.eq.0)then
         noverk=1
         return
      elseif(k.eq.1)then
         noverk=n
         return
      else
         num=n
         iden=1
         do i=2,k
            num=num*(n-i+1)
            iden=iden*i
         enddo
         noverk=num/iden
         return
      endif
      end
c-------------------------------------------------------------------
      subroutine printtp(ne,nd,trou,part,ncf,norb,c)
      integer*2 ne(*),trou(*),part(*)
      integer nd(*)
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
      return
      end
