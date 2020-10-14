      SUBROUTINE lecnam(NAME,NAMEL,IUNIT)
      IMPLICIT REAL*8 (A-H,O-X,Z),LOGICAL*1(Y)
C
C***********************************************************************
C   LECTURE DE LA NAMELIST DE NOM NAME. ON RANGE TOUTES LES DONNEES
C   DANS LA VARIABLE NAMEL.
C   LA NAMELIST PEUT CONTENIR UN MAXIMUN DE 500 CARACTERES NON BLANCS
C   (NAME ET &END EXCLUS), SOIT ENVIRON 7 CARTES PLEINES. LES DONNEES
C   DOIVENT SE TROUVER ENTRE LES COLONNES 2 ET 80.
C   IUNIT C'EST L'UNITE D'OU LES DONNES SONT LUS
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C  Occorre fare attenzione al caso in cui due o piu' nomi della
C  namelist siano contenuti l'uno dentro l'altro. Ad esempio se
C  a e ab sono due variabili dello stesso tipo dichiarate in
C  namelist bisogna che il call namel.. per ab preceda il call per a
C  Quindi la regola e': sistemare per prime le variabili i cui nomi
C  sono piu' lunghi !!!!!!!!
C
      PARAMETER (IN=5,IOUT=6)
      CHARACTER*(*) NAME
      character*1 der
      CHARACTER*1 BLANK,CC,COT1,COT2,COTEST
      CHARACTER*80 CARTE
C     CHARACTER*(NLEN) NAMEL
      CHARACTER*(*) NAMEL
      LOGICAL LEND,LFIRST
      common/dder/ider
C
C     DATA BLANK/' '/COT1/''''/,COT2/'"'/
      DATA BLANK/' '/,COT1/''''/
      COT2=CHAR(34)
      yw=.false.
      iecr=0
      nlen=len(namel)
      ierr=0
      if(yw) print *,'name(nom de la namelist) ',name
      if(yw)       PRINT*, 'nlen,longueur du nom ',nlen,len(name)
      if(yw)       call flush(6)
C
 9999 FORMAT(A80)
 9998 FORMAT(5X,'*** WAITING FOR THE NAMELIST ',A7,' I FOUND => ',A80,
     $' ** STOP ***')
 9997 FORMAT(5X,'*** TOO MANY CHARACTERS IN NAMELIST ',A8,' ** STOP **')
C
 46   continue
      READ(IUNIT,9999,end=1000) CARTE
      if(yw)      print*,'carte lue',carte
      call majnam(carte,80)
c     call majnam(name,len(name))
      if(yw)      print*,'carte mise en majuscules ',carte
C     print*,'name',name
      INAM=INDEX(CARTE,NAME)
      IF(INAM.EQ.0) THEN
         if(yw)WRITE(IOUT,9998) NAME,CARTE
c        STOP 200
         goto 46
      ENDIF
C
C                                   MISE A BLANK DE LA NAMELIST
C
      DO 100 I=1,NLEN
      NAMEL(I:I)=BLANK
  100 CONTINUE
      IC=0
      IST=INAM+LEN(NAME)
      itrouv=0
  300 IEND=INDEX(CARTE,'&END') + INDEX(CARTE,'$END')
c--   ajout fin de namelist = '/'
c     if(yw)       print*,'len(carte)',len(carte)
      if(yw)       print*,'dernier caractere de la carte ',der(carte,80)
      if(iend.eq.0.and.der(carte,80).eq.'/') iend=1
      if(yw)      print*,'iend (.ne.0 si fin de namelist trouvee)',iend
      if(yw)      print*,' longueur de la ligne ',ider
c     if(yw)      print*,'carte4',carte
c--   si le signe de fin n'est pas seul sur la ligne :
      if(iend.ne.0) then
         if(yw) then
            if(ider.ge.2)print*,'ider,carte(ider)',
     *      ider,carte(ider:ider)
         endif
         if(der(carte,80).eq.'/') then
            nbcar=0
         else
            nbcar=3
         endif
         if(yw)print*,'nbcar',nbcar
         do k=1,ider-1-nbcar
            if(carte(k:k).ne.' ') itrouv=1
         enddo
         if(yw)print*,'nbcar,itrouv',nbcar,itrouv
         if(itrouv.eq.1) then
            iend=0
            do k=ider-nbcar,min(ider+3,80)
               carte(k:k)=' '
            enddo
         endif
      endif
      if(yw)      print*,'carte6',carte
      LEND=IEND.EQ.0
      LFIRST=.FALSE.
      IF(LEND) THEN
         if(yw)print*,'la carte n''est pas la fin de la namelist'
         call ajvir(carte)
         if(yw)print*,'on ajoute (eventuellement) une virgule ',carte
         IEND=80
      ELSE
         IEND=IEND-1
         if(yw)         print*,'fin namelist1'
         if(iecr.eq.1) then
            print*,' '
            print*,'=============== NAMELIST  ',name,
     *      ' ==( en donnee )============='
            do k=1,ic,60
               print'(50(1x,a60))', namel(k:k+59)
            enddo
            print*,
     *      '=========================='//
     *      '=================================='
            print*,' '
         endif
         return
      ENDIF
      DO 200 I=IST,IEND
      CC=CARTE(I:I)
c     print*,'i,cc',i,cc
      IF(CC.EQ.COT1.OR.CC.EQ.COT2) THEN
         IF(LFIRST.AND.CC.EQ.COTEST) THEN
            LFIRST=.NOT.LFIRST
            COTEST=BLANK
         ELSE IF(.NOT.LFIRST) THEN
            LFIRST=.NOT.LFIRST
            COTEST=CC
         ENDIF
      ENDIF
      IF(CC.NE.BLANK.OR.LFIRST) THEN
         IC=IC+1
         IF(IC.GT.NLEN) THEN
            WRITE(IOUT,9997) NAME
            STOP 201
         ENDIF
         NAMEL(IC:IC)=CC
c        if(yw)       print*,'ic,cc',ic,cc
c        if(yw)       print*,'namel',namel(1:ic)
      ENDIF
  200 CONTINUE
      if(yw) print*,'la ligne est rajoutee a l''info totale '
      if(yw) print*,namel(1:ic)
      IF(.NOT.LEND) RETURN
      IST=1
      if(itrouv.eq.1) then
         if(yw)      print*,'fin namelist2'
         if(iecr.eq.1) then
            print*,' '
            print*,'=============== NAMELIST  ',name,
     *      ' ==( en donnee )============='
            do k=1,ic,60
               print'(50(1x,a60))', namel(k:k+59)
            enddo
            print*,
     *      '=========================='//
     *      '=================================='
            print*,' '
         endif
         return
      endif
      READ(IUNIT,9999) CARTE
      if(yw)       print*,'carte suivante ',carte
      call majnam(carte,80)
      if(yw)       print*,'carte suivante (en majuscules) ',carte
      GOTO 300
C     RETURN
 1000 ierr=1
      if(yw)      print*,'lecnam ierr',ierr
      return
      END
       subroutine ajvir(carte)
      character*80 carte
      if(carte(80:80).ne.','.and.carte(80:80).ne.' ') then
         print*, carte
         print*, 'ligne de 80 caracteres sans vigule a la fin'
         stop
      endif
      do k=80,1,-1
         if(carte(k:k).eq.',') return
         if(carte(k:k).eq.'=') return
         if(carte(k:k).ne.' ') then
            carte(k+1:k+1)=','
            return
         endif
      enddo
      return
      end
      character*1 function der(carte,n)
      character*80 carte
      common/dder/ider
      do ider=n,1,-1
         if(carte(ider:ider).ne.' ') goto 1
      enddo
      return
 1    der=carte(ider:ider)
      return
      end
 
       SUBROUTINE NINPI(IA,IST,IEND,JBUF)
C
C***********************************************************************
C     TRANSCODAGE D'UN INTEGER*4
C
C ===>  REQUESTED ROUTINES: NECERR
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C
      CHARACTER*1 JCHAR(14),ITEMP,IA*(*)
C
      DATA JCHAR/'0','1','2','3','4','5','6','7','8','9','+','&',' ',
     *'-'/
C
      IFACT=1
      LL=LEN(IA)
      LI=MIN0(IST+79,LL)
      JBUF=0
      DO 10 NSTRT=IEND,IST,-1
      ITEMP=IA(NSTRT:NSTRT)
      DO 43 J=1,14
      IF(JCHAR(J).EQ.ITEMP) GOTO 45
   43 CONTINUE
      CALL NECERR(IA(IST:LI),NSTRT,2)
   45 IF(J.GE.11) THEN
         IF(NSTRT.NE.IST) CALL NECERR(IA(IST:LI),NSTRT,2)
         IF(J.GE.14) THEN
            JBUF=-JBUF
            RETURN
         ENDIF
      ELSE
         JBUF=JBUF+(J-1)*IFACT
         IFACT=IFACT*10
      ENDIF
   10 CONTINUE
      RETURN
      END
       SUBROUTINE NINPL(IA,IST,IEND,LBUF)
C
C***********************************************************************
C     TRANSCODAGE D'UN LOGICAL*4
C
C ===>  REQUESTED ROUTINES: NECERR
C***********************************************************************
C
      LOGICAL LBUF
      CHARACTER IA*(*)
      IFACT=1
      LL=LEN(IA)
      LI=MIN0(IST+79,LL)
      IND =INDEX(IA(IST:IEND),'T')
      IF(IND.GT.0)THEN
         LBUF=.TRUE.
         RETURN
      ENDIF
      IND=INDEX(IA(IST:IEND),'F')
      IF(IND.GT.0)THEN
         LBUF=.FALSE.
         RETURN
      ENDIF
      WRITE(6,*)' ERRORE NELLA LETTURA DI UN LOGICO'
      STOP
      END
       SUBROUTINE NAMELD(NAMEL,TBUF,CBUF2,isize)
C
C***********************************************************************
C   LECTURE DES DONNEES DE TYPE DOUBLE PRECISION DANS LE TABLEAU TBUF.
C   CBUF EST LE NOM DU TABLEAU.
C-----------------------------------------------------------------------
C  UN CALL LECNAM(NAME,NAMEL) EST NECESSAIRE AVANT L'APPEL A NAMELD
C
C ===>  REQUESTED ROUTINES: NINPI,NINPF,NECERR
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C     isize e' la lunghezza del reale (4 o 8)
C
      IMPLICIT REAL*8 (A-H,O-X,Z),LOGICAL*1(Y)
      PARAMETER (IN=5,IOUT=6)
      DOUBLE PRECISION TBUFF
      CHARACTER*(*) CBUF2
      character(len=len(cbuf2))cbuf
      CHARACTER*1 BLANK
C     CHARACTER*(NLEN) NAMEL
      CHARACTER*(*) NAMEL
      CHARACTER*1 CH
      LOGICAL LPAR
      real*4 tbuf,tbu(2)
      equivalence (tbu(1),tbuff)
      DIMENSION TBUF(*)
C
      DATA BLANK/' '/
C
      cbuf=cbuf2
      yw=.false.
 9999 FORMAT(5X,'*** ERREUR EN LISANT LA VARIABLE ',A6,' DANS ',A80,
     $' ***   STOP ***')
C
C             RECHERCHE LA VARIABLE ET DELIMITE LA ZONE DE SES DONNEES.
C
      LUN=LEN(NAMEL)
      if(yw)print*,'entree nameld'
!      call maj(cbuf,len(cbuf)) !
	call ucase(cbuf)
      if(yw)       print*,'namel ',namel
      if(yw)       print*,'lun ',lun
      if(yw)       print*,'cbuf ',cbuf
      nlen=lun
      IUNO=1
      JPOS=0
  400 IPOS=INDEX(NAMEL(IUNO:LUN),CBUF)+JPOS
      if(yw)      print*,'ipos',ipos,INDEX(NAMEL(IUNO:LUN),CBUF)
      IF(IPOS.EQ.JPOS)RETURN
      JPOS=IPOS
      IF(IPOS.NE.0)THEN
         LUNG=LEN(CBUF)
         CH=NAMEL(IPOS+LUNG:IPOS+LUNG)
         IF(CH.NE.' '.AND.CH.NE.'('.AND.CH.NE.'=')THEN
            IUNO=IPOS+LUNG
            GOTO400
         ENDIF
      ENDIF
      JPOS=0
      IUNO=1
      IST=INDEX(NAMEL(IPOS:NLEN),'=')+IPOS
      IEG=INDEX(NAMEL(IST:NLEN),'=')
      IF(IEG.EQ.0) THEN
         IEG=NLEN
      ELSE
         IEG=IEG+IST-1
      ENDIF
      IFIN=MIN0(IPOS+80,NLEN)
      DO 100 I=IEG,IST,-1
      IF(NAMEL(I:I).EQ.',') GOTO 110
  100 CONTINUE
      WRITE(IOUT,9999) CBUF,NAMEL(IPOS:IFIN)
      STOP 210
  110 INEXT=I-1
C
C                                    IST= 1ERE COLONNE DES DONNEES,
C                                    INEXT= DERNIERE COLONNE
C
      IPAR1=INDEX(NAMEL(IPOS:IST),'(')
      LPAR=IPAR1.NE.0
      IF(LPAR) THEN
         IPAR1=IPAR1+IPOS-1
         IPAR2=INDEX(NAMEL(IPAR1:IST),')')+IPAR1-2
C
C                                    DECODE L'INDICE DU TABLEAU
C
         CALL NINPI(NAMEL,IPAR1+1,IPAR2,IDEX)
C
C                        LECTURE DE LA VALEUR
C
         CALL NINPF(NAMEL,IST,INEXT,TBUFF)
         if(isize.eq.4)then
            tbuf(idex)=tbuff
         elseif(isize.eq.8)then
            tbuf(2*idex-1)=tbu(1)
            tbuf(2*idex)=tbu(2)
         else
            write(6,*)' il valore',isize,' per un reale e'' intollerabil
     &      e'
            stop
         endif
C
C              MET A BLANC LA ZONE DE LA NAMELIST QUI VIENT D'ETRE LUE.
C
         DO 200 I=IPOS,INEXT+1
         NAMEL(I:I)=BLANK
  200    CONTINUE
      ELSE
C
C                   LE TABLEAU EST DONNE AVEC LES ELEMENTS CONSECUTIFS
C
         IDEX=0
  310    IVIRG=INDEX(NAMEL(IST:INEXT),',')
         IF(IVIRG.EQ.0) THEN
            IVIRG=INEXT
         ELSE
            IVIRG=IVIRG-2+IST
         ENDIF
C
C                                 CHERCHE UN FACTEUR DE REPETITION
C
         IETOIL=INDEX(NAMEL(IST:IVIRG),'*')
         IF(IETOIL.EQ.0) THEN
            IETOIL=IST
            ICOUNT=1
         ELSE
            IETOIL=IETOIL+IST-2
            CALL NINPI(NAMEL,IST,IETOIL,ICOUNT)
            IETOIL=IETOIL+2
         ENDIF
C
C                                   LIT LA VALEUR
C
         CALL NINPF(NAMEL,IETOIL,IVIRG,tbuff)
         DO 300 I=1,ICOUNT
         IDEX=IDEX+1
         if(isize.eq.4)then
            tbuf(idex)=tbuff
         elseif(isize.eq.8)then
            tbuf(2*idex-1)=tbu(1)
            tbuf(2*idex)=tbu(2)
         else
            write(6,*)' il valore',isize,' per un reale e'' intollerabil
     &      e'
            stop
         endif
  300    CONTINUE
         IST=IVIRG+2
         IF(IST.LE.INEXT) GOTO 310
C
C                                MET LA ZONE DE LA NAMELIST A BLANC
C
         DO 320 I=IPOS,INEXT+1
         NAMEL(I:I)=BLANK
  320    CONTINUE
      ENDIF
      GOTO 400
C     RETURN
      END
       SUBROUTINE NAMELL(NAMEL,LBUF,CBUF2,isize)
C***********************************************************************
C   LECTURE DES DONNEES DE TYPE LOGICAL*4 DANS LE TABLEAU IBUF.
C   CBUF EST LE NOM DU TABLEAU.
C-----------------------------------------------------------------------
C  UN CALL LECNAM(NAME,NAMEL) EST NECESSAIRE AVANT L'APPEL A NAMELL
C
C ===> REQUESTED ROUTINES: NINPL
C***********************************************************************
C     isize e' la lunghezza del logico (1 o 4)
C
      IMPLICIT REAL*8 (A-H,O-X,Z),LOGICAL*1(Y)
      PARAMETER (IN=5,IOUT=6)
      CHARACTER*(*) CBUF2
      character(len=len(cbuf2))cbuf
      CHARACTER*1 BLANK
C     CHARACTER*(NLEN) NAMEL
      CHARACTER*(*) NAMEL
      CHARACTER*1 CH
      LOGICAL LPAR,SIGN
      LOGICAL*1 LBUF,lbu(4)
      logical lbuff
      equivalence (lbu(1),lbuff)
      DIMENSION LBUF(*)
C
      DATA BLANK/' '/
C
      cbuf=cbuf2
      yw=.false.
 9999 FORMAT(5X,'*** ERREUR EN LISANT LA VARIABLE ',A6,' DANS ',A80,
     $' ***   STOP ***')
C
C            RECHERCHE LA VARIABLE ET DELIMITE LA ZONE DE SES DONNEES.
C
      LUN=LEN(NAMEL)
!      call maj(cbuf,len(cbuf))
	call ucase(cbuf)
      if(yw)       print*,'namel',namel
      if(yw)       print*,'lun',lun
      if(yw)       print*,'cbuf',cbuf
      nlen=lun
      IUNO=1
      JPOS=0
  400 IPOS=INDEX(NAMEL(IUNO:LUN),CBUF)+JPOS
      if(yw)      print*,'ipos',ipos,INDEX(NAMEL(IUNO:LUN),CBUF)
      IF(IPOS.EQ.JPOS)RETURN
      JPOS=IPOS
      SIGN=.FALSE.
      IF(IPOS.NE.0)THEN
         LUNG=LEN(CBUF)
         IPOSL=IPOS+LUNG
         CH=NAMEL(IPOSL:IPOSL)
         IF(CH.NE.' '.AND.CH.NE.'('.AND.CH.NE.'=')THEN
            IUNO=IPOS+LUNG
            SIGN=.TRUE.
         ENDIF
      ENDIF
      IF(SIGN)GOTO400
      IUNO=1
      JPOS=0
      IST=INDEX(NAMEL(IPOS:NLEN),'=')+IPOS
      IEG=INDEX(NAMEL(IST:NLEN),'=')
      IF(IEG.EQ.0) THEN
         IEG=NLEN
      ELSE
         IEG=IEG+IST-1
      ENDIF
      IFIN=MIN0(IPOS+80,NLEN)
      DO 100 I=IEG,IST,-1
      IF(NAMEL(I:I).EQ.',') GOTO 110
  100 CONTINUE
      WRITE(IOUT,9999) CBUF,NAMEL(IPOS:IFIN)
      STOP 210
  110 INEXT=I-1
C
C                                    IST= 1ERE COLONNE DES DONNEES,
C                                    INEXT= DERNIERE COLONNE
C
      IPAR1=INDEX(NAMEL(IPOS:IST),'(')
      LPAR=IPAR1.NE.0
      IF(LPAR) THEN
         IPAR1=IPAR1+IPOS-1
         IPAR2=INDEX(NAMEL(IPAR1:IST),')')+IPAR1-2
C
C                             DECODE L'INDICE DU TABLEAU
C
         CALL NINPI(NAMEL,IPAR1+1,IPAR2,IDEX)
C
C                             LECTURE DE LA VALEUR
C
         CALL NINPL(NAMEL,IST,INEXT,lbuff)
         if(isize.eq.1)then
            lbuf(idex)=lbuff
         elseif(isize.eq.4)then
            lbuf(4*idex-3)=lbu(1)
            lbuf(4*idex-2)=lbu(2)
            lbuf(4*idex-1)=lbu(3)
            lbuf(4*idex)=lbu(4)
         else
            write(6,*)' il valore',isize,' per un logico e'' intollerabi
     &      le'
            stop
         endif
C
C              MET A BLANC LA ZONE DE LA NAMELIST QUI VIENT D'ETRE LUE.
C
         DO 200 I=IPOS,INEXT+1
         NAMEL(I:I)=BLANK
  200    CONTINUE
      ELSE
C
C                   LE TABLEAU EST DONNE AVEC LES ELEMENTS CONSECUTIFS
C
         IDEX=0
  310    IVIRG=INDEX(NAMEL(IST:INEXT),',')
         IF(IVIRG.EQ.0) THEN
            IVIRG=INEXT
         ELSE
            IVIRG=IVIRG-2+IST
         ENDIF
C
C                           CHERCHE UN FACTEUR DE REPETITION
C
         IETOIL=INDEX(NAMEL(IST:IVIRG),'*')
         IF(IETOIL.EQ.0) THEN
            IETOIL=IST
            ICOUNT=1
         ELSE
            IETOIL=IETOIL+IST-2
            CALL NINPI(NAMEL,IST,IETOIL,ICOUNT)
            IETOIL=IETOIL+2
         ENDIF
C
C                               LIT LA VALEUR
C
         CALL NINPL(NAMEL,IETOIL,IVIRG,lbuff)
         DO 300 I=1,ICOUNT
         IDEX=IDEX+1
         if(isize.eq.1)then
            lbuf(idex)=lbuff
         elseif(isize.eq.4)then
            lbuf(4*idex-3)=lbu(1)
            lbuf(4*idex-2)=lbu(2)
            lbuf(4*idex-1)=lbu(3)
            lbuf(4*idex)=lbu(4)
         else
            write(6,*)' il valore',isize,' per un logico e'' intollerabi
     &      le'
            stop
         endif
  300    CONTINUE
         IST=IVIRG+2
         IF(IST.LE.INEXT) GOTO 310
C
C                            MET LA ZONE DE LA NAMELIST A BLANC
C
         DO 320 I=IPOS,INEXT+1
         NAMEL(I:I)=BLANK
  320    CONTINUE
      ENDIF
      GOTO 400
C     RETURN
      END
       SUBROUTINE NINPF(IA,IST,IEND,BUF)
C
C***********************************************************************
C    TRANSCODAGE D'UN REAL*8
C
C ===>  REQUESTED ROUTINES: NINPI,NECERR
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C
      real*8  ZERO,UN,DIX,PT1,BUF,FACT,EXPO
      PARAMETER (ZERO=0.0D+00,UN=1.0D+00,DIX=10.0D+00,PT1=0.1D+00)
      CHARACTER*1 JCHAR(17),ITEMP,IA*(*)
C
      DATA JCHAR/'0','1','2','3','4','5','6','7','8','9','+','&',' ',
     *'-','.','D','E'/
C
      LIMIT=15
      JEND=IEND
      LL=LEN(IA)
      LI=MIN0(IST+79,LL)
C
C                                 RECHERCHE D'UN EVENTUEL EXPOSANT
C
      IEXP1=INDEX(IA(IST:IEND),JCHAR(16))
      IF(IEXP1.EQ.0) IEXP1=INDEX(IA(IST:IEND),JCHAR(17))
      IF(IEXP1.NE.0) THEN
         CALL NINPI(IA,IEXP1+IST,IEND,JEXP)
         EXPO=DIX**JEXP
         JEND=IEXP1-2+IST
      ELSE
         EXPO=UN
      ENDIF
C
C                                  LECTURE DE LA CARACTERISTIQUE.
C
      IFACT2=0
      FACT=UN
      I=0
      BUF=ZERO
      DO 9 NSTRT=JEND,IST,-1
      I=I+1
      ITEMP=IA(NSTRT:NSTRT)
      DO 3 J=1,LIMIT
      IF(JCHAR(J).EQ.ITEMP) GOTO 5
    3 CONTINUE
      CALL NECERR(IA(IST:LI),NSTRT,1)
    5 IF(J.GE.11) THEN
         IF(J.GT.14) THEN
            IFACT2=I-1
            LIMIT=14
         ELSE
            IF(NSTRT.NE.IST) CALL NECERR(IA(IST:LI),NSTRT,1)
            IF(J.EQ.14) BUF=-BUF
            BUF=(PT1**IFACT2)*BUF*EXPO
            RETURN
         ENDIF
      ELSE
         BUF=BUF+FLOAT(J-1)*FACT
         FACT=FACT*DIX
      ENDIF
    9 CONTINUE
      BUF=(PT1**IFACT2)*BUF*EXPO
      RETURN
      END
       SUBROUTINE NECERR(IA,M,IER)
C
C***********************************************************************
C    SORTIE EN CAS D'ERREUR DE LECTURE
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 BLANK,HYPHEN,STAR
      CHARACTER*6 TYPE(2)
      CHARACTER*(*) IA,IB*80
C
      DATA BLANK/' '/,HYPHEN/'-'/,STAR/'*'/,TYPE/'REEL  ','ENTIER'/
C
 9999 FORMAT(/,3X,' === ERREUR DANS LA LECTURE D''UN ',A6,
     $//,3X,A80,/,3X,A80)
C
      DO 10 I=1,M-1
      IB(I:I)=HYPHEN
   10 CONTINUE
      DO 20 I=M+1,80
      IB(I:I)=BLANK
   20 CONTINUE
      IB(M:M)=STAR
      WRITE(6,9999) TYPE(IER),IA,IB
      STOP 211
      END
       SUBROUTINE NAMELI(NAMEL,IBUF,CBUF2,ISIZE)
C
C***********************************************************************
C   LECTURE DES DONNEES DE TYPE ENTIER DANS LE TABLEAU IBUF.
C   CBUF EST LE NOM DU TABLEAU.
C-----------------------------------------------------------------------
C  UN CALL LECNAM(NAME,NAMEL) EST NECESSAIRE AVANT L'APPEL A NAMELI
C
C ===> REQUESTED ROUTINES: NINPI,NECERR
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C     ISIZE PUO' ESSERE 2 o 4 A SECONDA DELLA LUNGHEZZA DELL'INTERO
C
      IMPLICIT REAL*8 (A-H,O-X,Z),LOGICAL*1(Y)
      PARAMETER (IN=5,IOUT=6)
      CHARACTER*(*) CBUF2
      character(len=len(cbuf2))cbuf
      CHARACTER*1 BLANK
C     CHARACTER*(NLEN) NAMEL
      CHARACTER*(*) NAMEL
      CHARACTER*1 CH
      LOGICAL LPAR
      integer*2 ibuf,ibu(2)
      equivalence (ibu(1),ibuff)
      DIMENSION IBUF(*)
C
      DATA BLANK/' '/
      cbuf=cbuf2
      yw=.false.
C
 9999 FORMAT(5X,'*** ERREUR EN LISANT LA VARIABLE ',A6,' DANS ',A80,
     $' ***   STOP ***')
!      call maj(cbuf,len(cbuf))
	call ucase(cbuf)
C
C            RECHERCHE LA VARIABLE ET DELIMITE LA ZONE DE SES DONNEES.
C
      LUN=LEN(NAMEL)
      if(yw)       print*,'namel ',namel
      if(yw)       print*,'lun ',lun
      if(yw)       print*,'cbuf ',cbuf
      if(yw)       call flush(6)
      nlen=lun
      IUNO=1
      JPOS=0
  400 IPOS=INDEX(NAMEL(IUNO:LUN),CBUF)+JPOS
      if(yw)      print*,'ipos',ipos,jpos
      IF(IPOS.EQ.JPOS)RETURN
      JPOS=IPOS
      IF(IPOS.NE.0)THEN
         LUNG=LEN(CBUF)
         CH=NAMEL(IPOS+LUNG:IPOS+LUNG)
         if(yw)print*,'lung,ch ',lung,ch
         IF(CH.NE.' '.AND.CH.NE.'('.AND.CH.NE.'=')THEN
            IUNO=IPOS+LUNG
            GOTO400
         ENDIF
      ENDIF
      JPOS=0
      IUNO=1
      IST=INDEX(NAMEL(IPOS:NLEN),'=')+IPOS
      IEG=INDEX(NAMEL(IST:NLEN),'=')
      IF(IEG.EQ.0) THEN
         IEG=NLEN
      ELSE
         IEG=IEG+IST-1
      ENDIF
      IFIN=MIN0(IPOS+80,NLEN)
      DO 100 I=IEG,IST,-1
      IF(NAMEL(I:I).EQ.',') GOTO 110
  100 CONTINUE
      WRITE(IOUT,9999) CBUF,NAMEL(IPOS:IFIN)
      STOP 210
  110 INEXT=I-1
C
C                                    IST= 1ERE COLONNE DES DONNEES,
C                                    INEXT= DERNIERE COLONNE
C
      IPAR1=INDEX(NAMEL(IPOS:IST),'(')
      LPAR=IPAR1.NE.0
      IF(LPAR) THEN
         IPAR1=IPAR1+IPOS-1
         IPAR2=INDEX(NAMEL(IPAR1:IST),')')+IPAR1-2
C
C                             DECODE L'INDICE DU TABLEAU
C
         CALL NINPI(NAMEL,IPAR1+1,IPAR2,IDEX)
C
C                             LECTURE DE LA VALEUR
C
         CALL NINPI(NAMEL,IST,INEXT,ibuff)
         if(isize.eq.2)then
            ibuf(idex)=ibuff
         elseif(isize.eq.4)then
            ibuf(2*idex-1)=ibu(1)
            ibuf(2*idex)=ibu(2)
         else
            write(6,*)' il valore',isize,' per un intero e'' intoller
     &      abile'
            stop
         endif
C
C              MET A BLANC LA ZONE DE LA NAMELIST QUI VIENT D'ETRE LUE.
C
         DO 200 I=IPOS,INEXT+1
         NAMEL(I:I)=BLANK
  200    CONTINUE
      ELSE
C
C                   LE TABLEAU EST DONNE AVEC LES ELEMENTS CONSECUTIFS
C
         IDEX=0
  310    IVIRG=INDEX(NAMEL(IST:INEXT),',')
         IF(IVIRG.EQ.0) THEN
            IVIRG=INEXT
         ELSE
            IVIRG=IVIRG-2+IST
         ENDIF
C
C                           CHERCHE UN FACTEUR DE REPETITION
C
         IETOIL=INDEX(NAMEL(IST:IVIRG),'*')
         IF(IETOIL.EQ.0) THEN
            IETOIL=IST
            ICOUNT=1
         ELSE
            IETOIL=IETOIL+IST-2
            CALL NINPI(NAMEL,IST,IETOIL,ICOUNT)
            IETOIL=IETOIL+2
         ENDIF
C
C                               LIT LA VALEUR
C
         CALL NINPI(NAMEL,IETOIL,IVIRG,ibuff)
         DO 300 I=1,ICOUNT
         IDEX=IDEX+1
         if(isize.eq.2)then
            ibuf(idex)=ibuff
         elseif(isize.eq.4)then
            ibuf(2*idex-1)=ibu(1)
            ibuf(2*idex)=ibu(2)
         else
            write(6,*)' il valore',isize,' per un intero e'' intoller
     &      abile'
            stop
         endif
  300    CONTINUE
         IST=IVIRG+2
         IF(IST.LE.INEXT) GOTO 310
C
C                            MET LA ZONE DE LA NAMELIST A BLANC
C
         DO 320 I=IPOS,INEXT+1
         NAMEL(I:I)=BLANK
  320    CONTINUE
      ENDIF
      GOTO 400
C     RETURN
      END
       SUBROUTINE NAMELA(NAMEL,TBUF,CBUF2)
C
C***********************************************************************
C   LECTURE DES DONNEES DE TYPE CHARACTER DANS LE TABLEAU TBUF.
C   CBUF EST LE NOM DU TABLEAU.
C-----------------------------------------------------------------------
C  UN CALL LECNAM(NAME,NAMEL) EST NECESSAIRE AVANT L'APPEL A NAMELA
C
C ===>  REQUESTED ROUTINES: NINPI,NECERR
C-----------------------------------------------------------------------
C  WRITTEN BY J.P. FLAMENT
C             GROUPE DE CHIMIE THEORIQUE
C             ECOLE POLYTECHNIQUE    91128 PALAISEAU CEDEX FRANCE
C***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-X,Z),LOGICAL*1(Y)
      PARAMETER (IN=5,IOUT=6)
      CHARACTER*(*) CBUF2,TBUF
      character(len=len(cbuf2))cbuf
      CHARACTER*1 CC,BLANK,CHRTST
C     CHARACTER*(NLEN) NAMEL
      CHARACTER*(*) NAMEL
      CHARACTER*1 CH
      LOGICAL LPAR
      DIMENSION TBUF(*)
C
      DATA BLANK/' '/
C
 9999 FORMAT(5X,'*** ERREUR EN LISANT LA VARIABLE ',A6,' DANS ',A80,
     $' ***   STOP ***')
      cbuf=cbuf2
      yw=.false.
!      call maj(cbuf,len(cbuf))
	call ucase(cbuf)
C
C             RECHERCHE LA VARIABLE ET DELIMITE LA ZONE DE SES DONNEES.
C
C NELLA PROSSIMA SCHEDA C'ERA MENO IN ORIGINE
      ILEN=LEN(TBUF(1))+1
      LUN=LEN(NAMEL)
      if(yw)       print*,'namel',namel
      if(yw)       print*,'lun',lun
      if(yw)       print*,'cbuf',cbuf
      nlen=lun
      IUNO=1
      JPOS=0
  400 IPOS=INDEX2(NAMEL(IUNO:LUN),CBUF)+JPOS
      IF(IPOS.EQ.JPOS)RETURN
      JPOS=IPOS
      IF(IPOS.NE.0)THEN
         LUNG=LEN(CBUF)
         CH=NAMEL(IPOS+LUNG:IPOS+LUNG)
         IF(CH.NE.' '.AND.CH.NE.'('.AND.CH.NE.'=')THEN
            IUNO=IPOS+LUNG
            GOTO400
         ENDIF
      ENDIF
      JPOS=0
      IUNO=1
      IEG=INDEX(NAMEL(IPOS:NLEN),'=')+IPOS-1
C
C                   CHRTST EST LE 1ER CARACTERE APRES LE SIGNE = .
C                   IL SERA CONSIDERE COMME LE DELIMITEUR DE LA CHAINE.
C
      ICOT1=INDEX(NAMEL(IPOS:NLEN),'''')
      ICOT2=INDEX(NAMEL(IPOS:NLEN),CHAR(34))
      ICOT=MIN0(ICOT1,ICOT2)
      IF(ICOT1.EQ.0) ICOT=ICOT2
      IF(ICOT2.EQ.0) ICOT=ICOT1
      IF(ICOT1.EQ.0.AND.ICOT2.EQ.0) THEN
         WRITE(IOUT,9999) CBUF,NAMEL(IPOS:IFIN)
         STOP 210
      ENDIF
      IST=ICOT+IPOS-1
      CHRTST=NAMEL(IST:IST)
C
C                                  IST= 1ERE COLONNE DES DONNEES,
C                                  INEXT= DERNIERE COLONNE
C
      IPAR1=INDEX(NAMEL(IPOS:IST),'(')
      LPAR=IPAR1.NE.0
      IF(LPAR) THEN
         IPAR1=IPAR1+IPOS-1
         IPAR2=INDEX(NAMEL(IPAR1:IST),')')+IPAR1-2
C
C                                  DECODE L'INDICE DU TABLEAU
C
         CALL NINPI(NAMEL,IPAR1+1,IPAR2,IDEX)
C
C                                  LECTURE DE LA VALEUR
C
         IST=IST+1
         DO 100 I=IST,IST+ILEN
         IF(NAMEL(I:I).EQ.CHRTST) GOTO 110
  100    CONTINUE
         IFIN=MIN0(IPOS+79,NLEN)
         WRITE(IOUT,9999) CBUF,NAMEL(IPOS:IFIN)
         STOP 211
  110    INEXT=I-1
         TBUF(IDEX)=NAMEL(IST:INEXT)
C
C              MET A BLANC LA ZONE DE LA NAMELIST QUI VIENT D'ETRE LUE.
C
         DO 200 I=IPOS,INEXT+2
         NAMEL(I:I)=BLANK
  200    CONTINUE
      ELSE
C
C                  LE TABLEAU EST DONNE AVEC LES ELEMENTS CONSECUTIFS
C
         IDEX=0
C
C                               CHERCHE UN FACTEUR DE REPETITION
C
  310    IST=IST+1
         IETOIL=INDEX(NAMEL(IEG:IST),'*')
         IF(IETOIL.EQ.0) THEN
            IETOIL=IST
            ICOUNT=1
         ELSE
            IETOIL=IETOIL+IEG-2
            CALL NINPI(NAMEL,IEG+1,IETOIL,ICOUNT)
            IETOIL=IETOIL+2
         ENDIF
C
C                                LIT LA VALEUR
C
         DO 300 I=1,ICOUNT
         IDEX=IDEX+1
         JFIRST=IST
         JLAST=IST+ILEN
         DO 220 J=JFIRST,JLAST
         IF(NAMEL(J:J).EQ.CHRTST) GOTO 230
  220    CONTINUE
         IFIN=MIN0(IPOS+79,NLEN)
         WRITE(IOUT,9999) CBUF,NAMEL(IPOS:IFIN)
         STOP 212
  230    INEXT=J-1
         TBUF(IDEX)=NAMEL(IST:INEXT)
  300    CONTINUE
         IEG=INEXT+2
         IST=INEXT+3
         CC=NAMEL(IST:IST)
         IF(CC.EQ.CHRTST) GOTO 310
         ICC=ICHAR(CC)
         IF(ICC.GE.ICHAR('0').AND.ICC.LE.ICHAR('9')) THEN
            DO 330 I=IST,NLEN
            IF(NAMEL(I:I).EQ.CHRTST) GOTO 340
  330       CONTINUE
            GOTO 350
  340       IST=I
            GOTO 310
         ENDIF
C
C                                 MET LA ZONE DE LA NAMELIST A BLANC
C
  350    DO 320 I=IPOS,INEXT+2
         NAMEL(I:I)=BLANK
  320    CONTINUE
      ENDIF
      GOTO 400
      END
       function index2(namel,c)
      logical y
      character*(*)namel,c
      l=len_trim77(namel)
      lc=len_trim77(c)
      y=.false.
      do k=1,l-lc+1
         if(namel(k:k).eq.'''') y=.not.y
         if(y) goto 1
         if(namel(k:k+lc-1).eq.c) then
            index2=k
            return
         endif
 1       continue
      enddo
      index2=0
      return
      end
       function len_trim77(namel)
      character*(*)namel
      l=len(namel)
      len_trim77=0
      do k=l,1,-1
         if(namel(k:k).ne.' ') then
            len_trim77=k
            return
         endif
      enddo
      return
      end
      SUBROUTINE CERCA(IUNIT,STRING,NCARD)
      character*(*) string
      character*16 card
      ncard=0
   10 read (iunit,20,err=99,end=99) card
   20 format (a16)
      ncard=ncard+1
      i=index(card,string)
      if (i.eq.0) goto 10
      j=index(card,'$')
      if (j.eq.0) j=index(card,'&')
      if (j.eq.0) goto 10
      if (i-j.ne.1) goto 10
      backspace iunit
      return
   99 ncard=0
      return
      end
       subroutine majnam(a2,l)
      IMPLICIT REAL*8 (A-H,O-X,Z),LOGICAL*1(Y)
      character*1 a2(l)
      dimension ichang(0:300)
      yw=.false.
      ichang(0)=1
      do k=1,l
         if(a2(k).eq.'''') then
            ichang(k)=-ichang(k-1)
         else
            ichang(k)=ichang(k-1)
         endif
      enddo
      if(yw)       write(6,'('' ichang'',20i2)') (ichang(k),k=1,20)
      do k=1,l
c         if(ichang(k).eq.1) call maj(a2(k),1)
         if(ichang(k).eq.1) call ucase(a2(k))
      enddo
      return
      end
      subroutine ucase(string)
      implicit none
! from lower to upper case
      character(len=*), intent(inout)::string
      integer :: i,length
      length=len(string)
      do i=1,length
       if(lge(string(i:i),'a').and.lle(string(i:i),'z'))then
          string(i:i)=achar(iachar(string(i:i))-32)
       endif
      enddo
      end subroutine
