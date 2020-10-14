c-----------------------------------------------------------------------
C 
C     PROGRAMME CIPSI - TOULOUSE
C     VERSION MODIFIEE "CIPPI" - PISA  JANVIER 1993
C 
C     INTERACTION DE CONFIGURATION PAR DIAGONALISATION DANS UN ESPACE
C     DETERMINANTAL D'ORDRE ZERO, S, ET PERTURBATION AU DEUXIEME ORDRE
C     (SINGLES ET DOUBLES EXCITATIONS A PARTIR DES DETERMINANTS DE S)
C  
C     LAVORA CON CONFIGURAZIONI AUTOSTATI DI S2
c
c     Si possono trattare eccitazioni fino alle dodecuple.
C     Ferrara 2000
C-----------------------------------------------------------------------
C 
C     ENTREE DES DONNEES
C 
C     CARTE 1
C     TITRE   (80 CHARACTERES)
C 
C 
C     NAMELIST FILES
C     ZBIN   = F, file07 is written in formatted form (default);
C              T, file07 is written in binary form.
C     ZCOST  = F(default), if TRUE file14 (file .Info from COST) and 
C                 file28 (file .Mono from COST) are read.
C     ZWCONF = T, the CSF's coefficients are written on the file29
C              F, default.
C 
C     NAMELIST ICINP
C     ZDET   = T, LES DETERMINANTS SONT LUS DU FICHIER 4
C              F, LES DETERMINANTS SONT LUS DU FICHIER 5 (DEFAUT)
C     ISZ    = COMPOSANTE Z DU SPIN TOTAL (DEFAUT, ISZ=0)
C     ZION   = F, NOMBRE PAIR D'ELECTRONS (DEFAUT)
C            = T, NOMBRE IMPAIR D'ELECTRONS; ON PLACE UN ELECTRON DE
C              SPIN +1/2 DANS UNE ORBITALE FICTIVE, NORB+1, NON
C              INTERAGISSANTE AVEC LES AUTRES; S'EMPLOIE OBLIGATOIREMENT
C              AVEC ISZ=1, LE DOUBLET ETANT CONSIDERE COMME UN TRIPLET
C              DEGENERE
C     METAT  = NOMBRE D'ETATS TRAITES EN PERTURBATION (DEFAUT, 1)
C     NROT   = NOMBRE DE VECTEURS PROPRES DETERMINES DANS L'ESPACE S
C              VOIR AUSSI IETATS, NAMELIST DIAGO, ET NREF, NAMELIST
C              DIABA (NROT.GE.METAT; DEFAUT, NROT=METAT)
C     ZPERTU = T, ON FAIT LE CALCUL PERTURBATIF (DEFAUT)
C            = F, PAS DE PERTURBATION
C     ZENB   = F, CALCUL EPSTEIN-NESBET VALEURS PROPRES (DEFAULT)
C            = T, CALCUL EPSTEIN-NESBET BARICENTRIQUE
C     ZMPPUN = T, i dati del file07 sono relativi alla part. MP
C            = F, i dati del file07 sono relativi alla part. EN (default)
C     TEST   = SEUIL D'IMPRESSION DES DETERMINANTS.
C              LE DETERMINANT I PERTURBANT L'ETAT M EST IMPRIME SI:
C              C(M,I) .GE. TEST*MAX(C0(M,K)), OU C0(M,K) SONT LES
C              COEFFICIENTS DES DETERMINANTS DE S POUR L'ETAT M
C              ET C(M,I) EST LE COEFFICIENT PERTURBATIF (ENVP) DE I
C              (DEFAUT, TEST=0.05)
C     TAU    = SEUIL POUR LA PERFORATION DES DETERMINANTS, ANALOGUE
C              A TEST (TAU.LE.TEST, PAR DEFAUT TAU=TEST)
C     ZQD    = T, LES DETERMINANTS I POUR LESQUELS EMPB(I) .LE.
C              MAX(EMPB(M)+DELTA U.A.), OU EMPB SONT LES ENERGIES MPB
C              DE I ET DE L'ETAT M, SONT IMPRIMES ET PERFORES (DEFAUT)
C            = F, PAS D'IMPRESSION DES DETERMINANTS QUASI-DEGENERES
C     DELTA  = SEUIL POUR LA DEFINITION DES DETERMINANTS QUASI-DEGENERES
C              (DEFAUT, DELTA=0.10).
C     ZBRD   = T, CALCUL DE L=HAMILTONIEN EFFECTIF QDPT
C            = F, PAS DE QDPT (DEFAUT)
C     ZPRT   = T, ECRITURE DE RESULTATS INTERMEDIAIRES (INTEGRALES J ET
C              K, TABLEAUX LOGIQUES, ITERATIONS MAXOVL, ETC)
C            = F, OUTPUT NORMAL (DEFAUT)
C     ZPUN   = T, LES DETERMINANTS SELECTIONNES POUR L'IMPRESSION SONT
C              AUSSI PERFORES
C            = F, PAS DE PERFORATION (DEFAUT)
C     ZHMAT  = T, ON CALCULE LA MATRICE D'I.C.  (SOUSPROGRAMME HMAT,
C              DEFAUT)
C            = F, ON LIT LA MATRICE DU FICHIER 20
C     KGEN   = NUMEROS DES DETERMINANTS GENERATEURS (VECTEUR, A
C              UTILISER SEULEMENT AVEC ZDET=F, PAR DEFAUT TOUS LES
C              DETERMINANTS DE S SONT GENERATEURS)
C     ZRAP   = T, ON FAIT L'APPROXIMATION DE NEGLIGER L'INTERACTION
C              ENTRE DETERMINANTS DE S NON GENERATEURS ET DETERMINANTS
C              PERTURBATEURS (CALCUL PLUS RAPIDE)
C            = F, CALCUL RIGUREUX DE L'INTERACTION (DEFAUT)
C     ZDIAB  = T, ROTATION DES VECTEURS D'ORDRE ZERO, PAR LE CRITERE
C              DE L'OVERLAP MAXIMUM AVEC UNE BASE (DIABATIQUE) DE
C              REFERENCE; IMPLIQUE ZBRD=T.
C            = F, CALCUL ADIABATIQUE (DEFAUT)
C     ZSPIN  = T, ON CALCULE LES ELEMENTS DE MATRICE DE S2 (CARRE DU
C              SPIN TOTAL).  DEFAUT, ZSPIN=F.
C     ZSROT  = T, LES VECTEURS PROPRES DE H DANS L'ESPACE S SONT MODIFIES
C              DE FACON QUE S2 SOIT DIAGONALISE; IMPLIQUE ZBRD=T.
C              DEFAUT ZSROT=F.
C     ZATOM  = T, ON CALCULE LES ELEMENTS DE MATRICE DE L2 (CARRE DU
C              MOMENT ANGULAIRE TOTAL) ET DE LZ2 (CARRE DE LA
C              PROJECTION DE L SUR L'AXE Z). POUR ATOMES SEULEMENT,
C              VOIR LE NAMELIST ATOMO. DEFAUT ZATOM=F.
C     ZDIALZ = T, LES VECTEURS PROPRES DE H DANS L'ESPACE S SONT MODIFIES
C              DE FACON QUE LZ2 SOIT DIAGONALISE PAR BLOCS. CHAQUE
C              BLOC EST DEFINI PAR UN SOUS-ESPACE DE VECTEURS DEGENERES
C              EN ENERGIE ET L2 (DEFAUT, ZDIALZ=F).
C     TDEGEN = SEUIL UTILISE DANS LES SSP DEGEN ET ATOM POUR DEFINIR LA
C              DEGENERATION DES VECTEURS PROPRES (DEFAUT, 1.D-7).
C     SHIFT(I),I=1,NSHIFT, DEPLACEMENTS DES ENERGIES MP ET EN DES
C              DETERMINANTS PERTURBATEURS, GENERANT DES SOMMES
C              PERTURBATIVES MODIFIEES DANS LES DENOMINATEURS.
C     NSHIFT = NOMBRE D'ALTERATIONS DES DENOMINATEURS ENVISAGEES
C              (DEFAUT, NSHIFT=0).
C     ZCAS   = The zero order wf is of CAS type (false by default)
C              can speed up the controls considerably !
C     ZTHIRD = printing of coefficients for successive third order 
C              calculation
C     ZMP3   = Moeller Plesset third order (default FALSE)
C     ZEN3   = Epstein Nesbet third order  (default TRUE if ZTHIRD.EQ.T)
C              ZMP3 and ZEN3 are mutually exclusive (ZEN3 prevails)
C     ZHANDY = T, INTEGRALES STOCKEES I*2 
C              F, INTEGRALES STOCKEES I*4. (DEFAUT)
C     FACTOR = FACTEUR MULTIPLICATIF DES INTEGRALES, QUI PERMET DE LES
C              STOCKER DANS LE TABLEAU IJKL2 (INTEGER*2).
C              SI  ABS(A*FACTOR).GT.20000, L'INTEGRALE A EST STOCKEE
C              DANS LE TABLEAU RIJKL (REAL*4), DE DIMENSION 25534.
C              SI ZHANDY=F, ON UTILISE LE TABLEAU IJKL4 (INTEGER*4)
C              (DEFAUT, FACTOR=1.D6 SI ZHANDY=T, 1.D9 SI ZHANDY=F).
C     FACTOR is deprecated. Use ZFAC instead
C     ZFAC   = T, FACTOR e` calcolato automaticamente sulla base del 
C                 massimo integrale
C              F, FACTOR ha il valore di default o quello dato in 
C                 input (se FACTOR viene dichiarato) (DEFAULT)
C     IGEL   = INDICI DEGLI ORBITALI CHE SARANNO MANTENUTI INATTIVI
C              (DOPPI OCCUPARI O VUOTI) NELLA PARTE PERTURBATIVA
C              DEFAULT: NESSUN ORBITALE GELATO
C     ZDAV   = F, the EN denominators are calculated using the energy of
C                 the perturbation deterninants (default);
C              T, the EN denominators are calculated using the average
C                 procedure od Davidson (average over all the determinants
C                 with a given orbital occupation).
C     DIR    = NOM de la directory qui contient les fichiers des
C              matrices T (des templates)
C     ZDH    = T, tutti i determinanti perturbatori hanno l'orbitale
C              NDHB ``bucato'': utilizzato nel calcolo degli IP di core
C              F, default
C     NDHB   = Orbitale spaziale di core ionizzato (ZDH)
C
C     ZDDH   = T, tutti i determinanti dello spazio S e perturbatori hanno gli orbitali
C              NDHB1 NDHB2 ``bucati'': utilizzato nel calcolo degli IP di core
C              nel caso di sistemi simmetrici
C              F, default
C     NDHB1   = Orbitale spaziale di core ionizzato (ZDDH)
C     NDHB2   = Orbitale spaziale di core ionizzato (ZDDH)
C
C     THRESH = Soglia per la scrittura dell'occupazione orbitalica.
C              Solo i coefficienti > thresh verranno considerati
C              Default = 0.0
C
C     ZTOUL  = Gli integrali sono letti da MOTRA+MOLCAS+MOLCSD secondo
C              il formalismo di Tolosa. Default = .false.
C 
C     NAMELIST DIAGO
C     SFO    = SEUIL DE CONVERGENCE DE LA METHODE DE DAVIDSON, POUR
C              LA FONCTION D'ONDE  (DEFAUT, 1.D-5)
C              SFO EST UN VECTEUR (DIM. 5).  ON REPETE LA RECHERCHE        
C              DES VECTEURS PROPRES JUSQU'A 5 FOIS, AVEC DES SEUILS
C              DECROISSANTS, EN DONNANT SFO(2), SFO(3) ETC.   CHAQUE
C              FOIS ON PREND COMME VECTEURS D'ESSAI CEUX DE LA TOURNEE
C              PRECEDENTE.
C     SEN    = SEUIL POUR L'ENERGIE  (DEFAUT, 0.)
C     NITER  = NOMBRE MAXIMAL D'ITERATIONS  (DEFAUT, 50)
C     MAXM   = DIMENSION MAXIMALE DU SOUS ESPACE DIAGONALISE
C              (MAXM.LE.ID14, PAR DEFAUT MAXM=12)
C     MAXGUE = NUMERO DI DETERMINANTI USATI PER PREPARARE I VETTORI DI PROVA
C              (VEDI ZVEC). DEFAULT = 200, MASSIMO=600 (IL VALORE MASSIMO
C              NON DEVE ECCEDERE (ID7/2)**.5
C     IOP    = 1, ON EVITE L'ACCUMULATION DES ERREURS LORS DE LA
C              REDUCTION DE LA DIMENSION M DU SOUS ESPACE DIAGONALISE
C              0, CALCUL NORMAL (PLUS RAPIDE)
C     ZVEC   = T, LES VECTEURS D'ESSAI SE TROUVENT DANS LE FICHIER 4
C              (INCOMPATIBLE AVEC ZDET=F) OU BIEN DANS LE FICHIER 14
C              (VOIR PROJF).
C              F, ON PREPARE LES VECTEURS D'ESSAI EN DIAGONALISANT UNE
C              PETITE SOUS MATRICE, DE DIMENSION AU PLUS ID12=700,ID13=150
C              (DEFAUT).
C     PROJF  = FICHIER 14, D'OU L'ON PEUT LIRE DES VECTEURS D'ESSAI POUR
C              L'ALGORITHME DE DAVIDSON.  CA DOIT ETRE LE FICHIER 4
C              D'UN CALCUL PRECEDENT AVEC UN SOUS-ESPACE DE L'ESPACE S
C              COURANT.  UTILISE SEULEMENT SI ZVEC=T (DEFAUT, PROFJ=' ', 
C              ON LIT LES VECTEURS D'ESSAI DU FICHIER 4, PAS DE
C              PROJECTION SUR UN SOUS-ESPACE).
C     ZSS    = T, ON SELECTIONNE LES VECTEURS D'ESSAI QUI APPARTIENNENT
C              A UNE SYMETRIE DONNEE
C              F, PAS DE SELECTION PAR SYMETRIE (DEFAUT)
C     CSEL,ISEL = VECTEUR REPRESENTATIF DE LA SYMETRIE VOULUE (SEULS
C              LES VECTEURS D'ESSAI QUI ONT UN RECOUVREMENT SIGNIFICATIF
C              AVEC CSEL SERONT RETENUS). CSEL=COEFFICIENTS, ISEL=
C              NUMEROS DES DETERMINANTS. DEFAUT, CSEL=0., ISEL=0
C     ZFOL   = T, LE VECTEUR PROPRE AUQUEL ON APPLIQUE L'ALGORITHME DE
C              DAVIDSON EST DETERMINE PAR SES COEFFICIENTS (A UTILISER
C              PREFERABLEMENT DANS LE CAS OU DES VECTEURS DE SYMETRIE
C              DIFFERENTE SONT CALCULES A LA FOIS)
C            = F, LE VECTEUR EST DETERMINE PAR L'ORDRE ENERGETIQUE
C              PROVISOIRE DE L'ITERATION COURANTE (DEFAUT)
C     ZWRT   = T, DES RESULTATS INTERMEDIAIRES SONT IMPRIMES (DEFAUT, F)
C     SECR   = SEUIL POUR L'IMPRESSION DES VECTEURS PROPRES (SEULS LES
C              COEFFICIENTS PLUS GRANDS QUE SECR SERONT IMPRIMES). PAR
C              DEFAUT, SECR=0.05
C     IETATS(I),I=1,NROT = NUMEROS DES ETATS QUE L'ON VEUT GARDER APRES
C              LA DIAGONALISATION.  ON DETERMINE IETATS(NROT) VECTEURS
C              PROPRES (ATTENTION AUX DIMENSIONS MAXIMALES DANS DAVID),
C              ET ON EN CHOISIT NROT SELON L'ORDRE ENERGETIQUE.
C              SI ZDIALZ=T (ZSROT=T), LA SELECTION SUIVRA LES MODIFICATIONS
C              DES VECTEURS QUE L'ON FAIT LORS DU CALCUL DE L2 (S2).
C              PAR DEFAUT, ON GARDE LES PREMIERS NROT ETATS.
C     ZREST  = T, RESTART DURANTE LA PROCEDURA DI DIAGONALIZZAZIONE DI
C              UN CALCOLO PRECEDENTE. I VETTORI SONO LETTI DAL FILE13
C     ZSPINSEL=T, SELEZIONE SUGLI AUTOVALORI DI S2 IN MODO DA AVERE SOLO
C              STATI DI SPIN VOLUTO
C     SPINMAX,SPINMIN=ESTREMI DELL'AUTOVALORE DI S2 PER CUI UNO STATO
C              VERRA' SELEZIONATO
C 
C     CARTES DEFINISSANT LES DETERMINANTS DE S
C     UNE CARTE CORRESPOND A UN DETERMINANT OU A UN GROUPE DE
C     DETERMINANTS. CHAQUE CARTE COMMENCE PAR UN CHARACTERE EN
C     COLONNE 1 :
C     F  PAS D'EXCITATION  (DOUBLE OCCUPATION; L'ENERGIE DE REFERENCE
C                           EST TOUJOURS CALCULEE PAR RAPPORT A CE
C                           DETERMINANT, MEME SI LE CALCUL SCF A ETE
C                           FAIT AVEC UNE OCCUPATION DIFFERENTE)
C     M  DETERMINANT MONOEXCITE
C     D  DETERMINANT BIEXCITE
C     T  TRIEXCITE
C     Q  QUADRIEXCITE
C     P  PENTAEXCITE
C     H  HEXAEXCITE
C     E  EPTAEXCITE
C     O  OTTAEXCITE
C     N  NONAEXCITE
C     X  DECAEXCITE
C     U  UNDECAEXCITE
C     Z  DODECAEXCITE
C     A PARTIR DE COLONNE 3, SEPARES PAR UN SEUL BLANC (OU VIRGULE)
C     LES TROUS, SUIVIS PAR UN OU PLUSIEURS CHOIX DE PARTICULES.
C     TROUS ET PARTICULES EN GENERAL SONT INDIQUES COMME SPINORBITALES
C     (NUMERO DE L'O.M. SUIVI PAR LE SPIN, + OU -). SI LE SPIN EST
C     OMIS, LE PROGRAMME GENERE AUTOMATIQUEMENT TOUTES LES PARTIES DE
C     SPIN COMPATIBLES AVEC ISZ (SEULEMENT POUR MONO, BI ET TRIEXCITES).
C     SI DES TRANSFORMATIONS DE SYMETRIE SPATIALE ONT ETE DEFINIES AU
C     NIVEAU IJKL, ELLES SONT APPLIQUEES A TOUT DETERMINANT ENGENDRE
C     ET IL EST INUTILE D'ENTRER LES DETERMINANTS EQUIVALENTS DANS S.
C     UNE CARTE BLANCHE TERMINE LA SEQUENCE DES DETERMINANTS.
C 
C     NAMELIST DIABA (SEULEMENT SI ZDIAB=T; PARAMETRES POUR LA
C     DIABATIZATION PAR ROTATION DES VECTEURS ADIABATIQUES, SUIVANT
C     LE CRITERE DE L'OVERLAP MAXIMUM AVEC UNE BASE D'ETATS DE
C     REFERENCE)
C     NREF   = NOMBRE D'ETATS DE REFERENCE. NORMALEMENT, NREF=NROT=
C              =METAT. SI NREF=METAT<NROT, ON CHOISIT DE PERTURBER DES
C              ETATS QUASI-DIABATIQUES QUI RASSEMBLENT LE PLUS POSSIBLE
C              AUX REFERENCES, EN LES DECOUPLANT DE NROT-METAT ETATS
C              AVEC LESQUELS ILS SONT MELANGES A L'ORDRE ZERO ADIA-
C              BATIQUE.  SI NREF-1=METAT=NROT (ZCOMB=T), ON UTILISE UNE
C              COMBINAISON LINEAIRE DES DEUX DERNIERES REFERENCES, ET
C              NREF => NREF-1.
C     ZCOMB  = T, LES DEUX DERNIERS ETATS DE REFERENCE SONT UTILISES EN
C              COMBINAISON LINEAIRE, DONT LES COEFFICIENTS SONT DETERM-
C              MINES PAR LE CRITERE DE L'OVERLAP MAXIMUM
C            = F, PAS DE COMBINAISON LINEAIRE (DEFAUT)
C     ZLOWD  = T, LES ETATS DE REFERENCE SONT ORTHOGONALISES A LA
C              LOWDIN
C            = F, PAS D'ORTHOGONALISATION (NORMALISATION SEULEMENT,
C              DEFAUT)
C     ZGRAM  = T, LES ETATS DE REFERENCE SONT ORTHOGONALISES A LA
C              GRAM-SCHMIDT
C            = F, PAS D'ORTHOGONALISATION (NORMALISATION SEULEMENT,
C              DEFAUT)
C     ZRD    = T, LE CRITERIOM DE L'OVERLAP MAXIMUM EST MODIFIE DE
C              FACON QUE, SI LA REFERENCE I N'A PAS DE SUPERPOSITION
C              AVEC AUCUN ETAT ADIABATIQUE, L'ETAT DIABATIQUE I SERA
C              IDENTIFIE AVEC L'ETAT ADIABATIQUE I (L'ORDRE DES REFE-
C              RENCES N'EST PAS SANS IMPORTANCE, DANS CE CAS-CI)
C     NREFMO = NOMBRE D'ORBITALES UTILISES POUR CONSTRUIRE LES ETATS
C              DE REFERENCE
C     NREFCI = NOMBRE DE DETERMINANTS UTILISES POUR CONSTRUIRE LES
C              ETATS DE REFERENCE
C     ZIDEN  = T, LES ORBITALES DE REFERENCE COINCIDENT AVEC LES O.M.
C              UTILISES POUR L'IC
C            = F, ON LIT LES ORBITALES DE REFERENCE (DEFAUT)
C     IREFMO = NUMEROS DES OM UTILISES DANS LES DETERMINANTS DE
C              REFERENCE.  VALIDE SEULEMENT POUR ZIDEN=T.  DEFAULT,
C              IREFMO=1,2,3... NREFMO
C     MODISK = UNITE' DE LECTURE DES O.M. DE REFERENCE, 5 (DEFAUT)
C              OU BIEN 12.
C     MOFORM = FORMAT DE LECTURE DES OM DE REFERENCE (DEFAUT, 5F14.10)
C              ON LIT UN VECTEUR A LA FOIS, (REFMO(I,J),I=1,NAO)).
C            = STAR, LIBRE FORMAT
C            = BIN (SEULEMENT SI MODISK=12), ON LIT LES OM EN BINAIRE,
C              UN SEUL ENREGISTREMENT, ((REFMO(I,J),I=1,NAO,J=1,NREFMO).
C     IWRT   = 0,1,2, OUTPUT DE MINIMALE A ETENDU (DEFAUT, IWRT=1)
C 
C     CARTES DEFINISSANT LES OM DE REFERENCE (VOIR MODISK ET ZIDEN,
C     CARTE PRECEDENTE)
C 
C     CARTES DEFINISSANT LES DETERMINANTS DE REFERENCE
C     UNE CARTE POUR CHAQUE DETERMINANT, INDIQUANT L'OCCUPATION DES
C     SPINORBITALES BETA, PUIS CELLE DES ALPHA (T POUR OCCUPE, F OU
C     BLANK POUR VIDE). LE SIGNE DES DETERMINANTS SUIT LA CONVENTION:
C     TOUS LES ELECTRONS BETA AVANT, TOUS LES ALPHA APRES.
C 
C     DEUX TYPES DE CARTES POUR CHAQUE ETAT DE REFERENCE:
C     CARTE A, LABEL DE L'ETAT, COLONNES 1-16
C     CARTES B, VECTEUR DANS LA BASE DES DETERMINANTS DE REFERENCE,
C     5F14.10
C 
C     NAMELIST ATOMO             (SEULEMENT SI ZATOM=T)
C     PARAMETRES POUR LE CALCUL DES ELEMENTS DE MATRICES DE LX
C     (PROJECTION DU MOMENT ANGULAIRE SUR L'AXE X), DE LY ET DE LZ DANS
C     LA BASE DES O.M. 
C     IP     = NOMBRE QUI SPECIFIE LA POSITION DES FONCTION P DANS LA
C              BASE DES A.O.; POUR CHAQUE SET DE TROIS FONCTIONS P IL 
C              FAUT SPECIFIER SEULEMENT LE NUMERO DE LA PREMIERE
C     ID     = COMME IP, MAIS POUR LES FONCTIONS D
C     ATTOL  = SEUIL UTILISE POUR DEFINIR UNE O.M. COMME S, OU P, OU D
C              (DEFAUT, ATTOL=1.D-6)
C-----------------------------------------------------------------------
C 
C     DIMENSIONS MAXIMALES
C     BASE D'ORBITALES MOLECULAIRES ............... 150 = ID1
C     ESPACE DETERMINANTAL D'ORDRE ZERO ........ 100000 = ID4
C     DETERMINANTS LUS DU FICHIER 5 ............ 100000 = ID5
C     DETERMINANTS GENERATEURS ................. 100000 = ID5
C     MONOEXCITATIONS .......................... 100000 = MIN(ID4,ID8)
C     DIEXCITATIONS ............................ 100000 = MIN(ID4,ID8/2)
C     TRI ...................................... 100000 = MIN(ID4,ID8/3)
C     QUADRI ...................................  75000 = MIN(ID4,ID8/4)
C     PENTA ....................................  60000 = MIN(ID4,ID8/5)
C     HEXA .....................................  50000 = MIN(ID4,ID8/6)
C     EPTA .....................................  42857 = MIN(ID4,ID8/7)
C     OCTA .....................................  37500 = MIN(ID4,ID8/8)
C     NONA .....................................  33333 = MIN(ID4,ID8/9)
C     DECA .....................................  30000 = MIN(ID4,ID8/10)
C     UNDECA....................................  27272 = MIN(ID4,ID8/11)
C     DODECA....................................  25000 = MIN(ID4,ID8/12)
C     ETATS TRAITES ................................ 18 = ID6
C     ETATS * DETERMINANTS ..................... 110001 = ID7
C     COUCHES OUVERTES DANS UN DETERMINANT .......... 9
C     TROUS ET PARTICULES ...................... 650000 = ID8
C     BASE ATOMIQUE (SSP. ATOM ET DIAOVL) ......... 150 = ID13
C     NOMBRE DE VECTEURS D'ESSAI (DAVIDSON) ........ 30 = ID14
C 
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)

      pointer morb(:),mspin(:),isingle(:),iel(:),morbp(:)
      integer*2 morb,morbp,mspin
      common /ric/ ihalf,nec1,morb,isingle,mspin,iel,morbp
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN,itm,itk
      allocatable lgen(:),zef(:,:),itm(:,:),itk(:,:)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      common/third/zthird,zmp3,zen3
      common/debug/zdebug
      common /fil/ file12,file13,file26,zwconf
      COMMON /ALE/SPINMAX,SPINMIN,ZSPINSEL,ZINCRE
      character*64 projf,file12,file13,file26,dir
      INTEGER*2 ietats
      COMMON /CPU/ TCPU,PCPU
      COMMON /MEM/ MEMSIZ,LOFFS,LOCMEM,LOCM(10)
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      common/thresho/thresh
      common/classi/E2ENc(ID6,8),psienec(id6,8),E2MPc(ID6,8),
     * psimpec(id6,8)
      common /toul/ ztoul
      allocatable iocc(:,:)
      integer*1 iocc
      integer*2 icomp,iconf,isegno
#ifdef SINGLET
      common/zmat/ zint(460673)
#endif
#ifdef DOUBLET
      common/zmat/ zint(115009)
#endif
#ifdef TRIPLET
      common/zmat/ zint(336842)
#endif
      common /big/ x,ijkl2,ijkl4,ijkl3,rijkl
      pointer x
      dimension x(:)
      pointer ijkl2,ijkl3,ijkl4
      dimension ijkl2(:),ijkl3(:),ijkl4(:),rijkl(:)
      pointer rijkl
      integer*2 ijkl2
      character*3,ijkl3
      real*4 rijkl
      allocatable indis(:,:),tocom(:)
c      INTEGER*2 INDIS(ID1,20)
      INTEGER*2 INDIS
c      DIMENSION TOCOM(ID1)
      DIMENSION NOS(20)
      external lshift
      common /deephole/ ndhb,zdh,ndhb1,ndhb2,zddh
      character*4 strarg
c
      print*,'Version compiled with intel ifc compiler'
c      if(iargc().eq.0)then
c         print*,'program called with no arguments: id1 assumed as',id1
c         idimorb=id1
c      else
c         call getarg(1,strarg)
c         print '(a,a)','strarg is ',strarg
c         read(strarg,'(i4)')idimorb
c         print*,'id1 will be set equal to',idimorb
c      endif
      ntotdet=0
      ntotconf=0
      ntotcap=0
      write (6,*) '   '
      CALL CPUT(-1)
      WRITE (6,'(/1X,130(''*''))')
      WRITE (6,'(/1X,15A8)') ('  CIPPI ',I=1,15)
      WRITE (6,'(/1X,130(''*''))')
      CALL OFILES(zbin,ndims,ndimtp,zcost)
      WRITE (6,'(/1X,130(''*''))')
      ZTR=.TRUE.
      ZFA=.FALSE.
      if(zcost)call crea25
      rewind 25
      READ (25) NSYM,NORB
      rewind 25
      idimorb=norb
#ifdef DOUBLET
      idimorb=idimorb+2
#endif
      idimorb2=idimorb*(idimorb+1)/2
      idimorb3=idimorb*2+1
c----provisional allocations
c      allocate(zgel(idimorb3))
      allocate(fmpb(idimorb3))
cc      allocate(c(id7))
c      allocate(nd(ndimtp2))
c      allocate(ne(ndimtp2))
c      allocate(trou(ndimtp))
c      allocate(part(ndimtp))
c      allocate(kgen(ndims))
c      allocate(ndeg(ndims))
c      allocate(iorb(idimorb3))
c      allocate(ispin(idimorb3))
c      allocate(itsym(idimorb3))
c      allocate(zsaut(ndims))
      allocate(zoc(idimorb3))
      allocate(zocs(idimorb3))
c      allocate(zfree(idimorb))
      allocate(zheff(id11))
      allocate(f(idimorb2))
      allocate(aj(idimorb2))
      allocate(ak(idimorb2))
      allocate(num(idimorb2))
      allocate(indic(idimorb2))
      allocate(jndic(idimorb2))
      allocate(lndic(idimorb2))
      allocate(indis(idimorb,20))   !!!!!
      allocate(tocom(idimorb))      !!!!!
      allocate(morb(idimorb))       !!!!!
      allocate(mspin(idimorb))      !!!!!
      allocate(isingle(idimorb))    !!!!!
      allocate(iel(idimorb))        !!!!!
      allocate(morbp(idimorb))      !!!!!
C     Chiamo definfo che costruisce alcune matrici utili in
C     seguito
c
      call definfo(idimorb)
C 
C     LECTURE DES DONNEES
C 
      CALL DETER(zdet,TOCOM,INDIS,tgen,zrest,zdav,zdavid2,dir
     $     ,idimorb,iditrpa,ndims)
      write (6,*) 'Zdavid2',zdavid2
      call flush(6)
C 
C     LECTURE DES INTEGRALES
C 
c      deallocate(c)
c      if(allocated(kgen))deallocate(kgen)
c      if(allocated(ndeg))deallocate(ndeg)
c      deallocate(zsaut)
c      allocate(c(NCF*METAT+1000))
c      allocate(cconf(ncf*metat+1000))
      allocate(kgen(ncfg+1000))
      allocate(ndeg(ncfg+1000))
c      allocate(zsaut(ncf+1000))
      allocate(icomp(ncf+1000))
      allocate(iconf(ncf+1000))
      allocate(isegno(ncf+1000))
      allocate(zcapos(ncf+1000))
      allocate(zact(idimorb))
c      jd1=id1
      jd1=idimorb
      CALL REIJKL(factor,norb,nsym,itsym,its,nos,jd1)
      finv=1.d0/factor
      CALL MEMCIP('TOTAL',.false.)
      CALL MEMCIP('IJKLH',.true.)
      if (ztoul) then
      CALL ITIJKLt(factor,norb,nsym,nos,nijkl,
     * ijkl4,ijkl2,rijkl,ijkl3)
      else
      CALL ITIJKL(factor,norb,nsym,nos,nijkl,
     * ijkl4,ijkl2,rijkl,ijkl3)
      endif
      CALL IJKF(1,TOCOM)
C 
C     I.C. VARIATIONNELLE: CALCUL DE LA MATRICE D'I.C.
C 
c      CALL HMAT(X(locm(5)),X(locm(6)),X(locm(7)),X(locm(8)),X(locm(9)))
cele16-11-04
       do i=1,norb !spostato piu` su ...
       zact(i)=.false.
       enddo
cele16-11-04
      if(zcas)then
         print*,'The zero order wave funtion is of CAS type'
cele16-11-04         do i=1,norb !spostato piu` su ...
cele16-11-04            zact(i)=.false.
cele16-11-04         enddo
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
         do i=1,norb
            if(zact(i))print*,'orbitale ',i,' attivo'
            if(.not.zact(i))print*,'orbitale ',i,' inattivo'
         enddo
      endif
c--renzo--modif--begin
      call esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,nocb,norb)
      ncapos=0
      inext=1
      do i=1,ncf
         if(i.eq.inext)then
            zcapos(i)=.true.
            inext=i+icomp(i)
         else
            zcapos(i)=.false.
         endif
         if(zcapos(i))ncapos=ncapos+1
      enddo
      print '(a,i6,a)','Ci sono',ncapos,' capostipiti'
      call flush(6)   !tolgo flush
      allocate(iocc(1:norb,1:ncapos))
      print*,'Allocata la memoria per iocc'
      call flush(6)   !tolgo flush
      icapos=0
      do i=1,ncf
         if(zcapos(i))then
            icapos=icapos+1
         call giveocc(iocc(1,icapos),i,nd,ne,trou,part,nocb,norb)
      endif
      enddo
      call hmat(nconf,dir,iocc,norb)
c--renzo--modif--end      
      CALL CPUT(1)
      CALL MEMCIP('RELSE',.true.)
C 
C     I.C. VARIATIONNELLE: VECTEURS ET VALEURS PROPRES
C 
      IF (METAT.GT.NCF) METAT=NCF
      nvec=ietats(nrot)
c alessandro: luglio 98 
      nvecsave=nvec
 100  continue
      if (ZSPINSEL) then
       zspin=.TRUE.
       do i=1,id6
	ietats(i)=0
       enddo
      endif
c alessandro: luglio 98 
      IF (nvec.GT.NCF) then
      write (6,*)
      write (6,*) ' TROP D''ETATS EN DIAGONALISATION, NVEC,NCF,ID6 =',
     * NVEC,NCF,ID6
      stop 12
      endif
      if(.not.zdavid2)CALL MEMCIP('DAVID',.false.)
c      CALL DAVID(NCF,NVEC,X(locm(1)),X(locm(2)),zrest)
      if(.not.zdavid2)then
      allocate(c(NCF*NVEC+1000))
      allocate(cconf(ncf*NVEC+1000))
         CALL DAVID(ncf,NCONF,NVEC,X(locm(1)),X(locm(2)),zrest,secr
     $        ,iditrpa)
      else
         call david2(ncf,nconf,nvec,secr,zrest)
      endif
      CALL CPUT(1)
      if(.not.zdavid2)CALL MEMCIP('RELSE',.false.)
      ZVBIS=.FALSE.
      CALL WRTVAR(NVEC,ZVBIS,zbin)
c
c     Cele: si definisce un set di generatori nel caso ZPERTU=F e
c     ZDET=T. Servono per costruire le opportune variabili che 
c     seguono (ZEF, ZOCD, ZOCS e ZFREE). Non vogliamo usare i
c     generatori perche` potrebbero essere ben maggiori di ID5.
c
      if (zdet.and..not.zpertu) then
         ncfg=0
         do k=1,ncf
            if (.not.zsaut(k)) then
               ztake=.false.
               do m=1,nrot
                  km=(m-1)*ncf+k
                  if (dabs(c(km)).gt.TGEN) ztake=.true.
               enddo
            endif
               if (ztake) then
                  ncfg=ncfg+1
                  kgen(ncfg)=k
               endif
         enddo
      WRITE (6,240) NCFG
  240 FORMAT (5X,'Definizione di un set di ''generatori'' interni',/
     * 5X,'NOMBRE DE DETERMINANTS GENERATEURS, NCFG =',I6/)
      call flush(6)   !tolgo flush
c      if (ncfg.gt.id5) then
c      write (6,*) 'Errore: NCFG maggiore di id5',NCFG,ID5
c      write (6,*) 'Aumentare TGEN (default .05) o ID5'
c      endif
      endif
C
C     SPIN, MOMENT ANGULAIRE, SYMETRISATION ET DIABATISATION
C
      NNROT=NROT*(NROT+1)/2
      DO  I=1,NROT
         do  j=1,nrot
            brdmp(i,j)=0.d0
            brden(I,j)=0.D0
         enddo
      enddo
      if (zspin)  then
c alessandro: luglio 98 
      call spin(nvec,zvbis)  !useless with cippiS2
       if (zincre) then
	write(6,*)
	write(6,*)'   NON E` STATO RAGGIUNTO IL NUMERO DI STATI',
     *'RICHIESTO;   AUMENTO IL NUMERO DI AUTOVETTORI'
	write(6,*)
	nvec=nvec+1
	goto 100
      endif
      nvec=nvecsave
c alessandro: luglio 98 
      if (zsrot) zvbis=.true.
      endif
      if (zatom) then
      CALL MEMCIP('ATOML',.false.)
      call atom(nvec,zvbis,X(locm(1)))
      CALL MEMCIP('RELSE',.false.)
      endif
      call elim(zvbis)
      CALL DEGEN(zdet,ZVBIS)
      IF (ZDIAB) then
      CALL DIABAT(nconf)
      ZVBIS=.TRUE.
      endif
      DO 80 I=1,NROT
      do 80 j=1,nrot
      brdmp(I,j)=brdmp(I,j)*FACTOR*2.D0
   80 brden(i,j)=brdmp(i,j)
      IF (ZVBIS) CALL WRTVAR(METAT,ZVBIS,zbin)
c      if (zspin.and.zvbis) call spin(nrot,zvbis) !useless with cippiS2
      if (zatom.and.zvbis) then
      CALL MEMCIP('ATOML',.false.)
      call atom(nrot,zvbis,X(locm(1)))
      CALL MEMCIP('RELSE',.false.)
      endif
      DO 90 I=1,METAT
   90 E(I)=E(I)*FACTOR
C 
C     PREPARATION DE L'I.C. PERTURBATIVE
C 
c      CALL MEMCIP('PREPR',.false.)
      allocate(zef(ncfg,20))
      allocate(itm(id6,20))
      allocate(itk(ncfg,10))
      CALL PREPER(zdet,idimorb3,
     * zef,itm,itk)
c      CALL MEMCIP('RELSE',.false.)
      deallocate(zef)
      deallocate(itm)
      deallocate(itk)
c      CALL MEMCIP('MKF04',.false.)
      allocate(lgen(ncf))
      CALL MKF04(zdet,lgen,X(locm(2)))
c      CALL MEMCIP('RELSE',.false.)
      deallocate(lgen)
      if (zenb) WRITE (6,1010) 'ENB ','ENB '
      if (.not.zenb) WRITE (6,1010) 'ENVP','ENVP'
      DO 300 I=1,METAT
      EMPI=EMP(I)*FINV
      EI=EEN(I)*FINV
  300 write (6,1040) i,teste(i),empi,ei,empi+delta,ei+delta
      CALL CPUT(1)
      call flush(6)
      IF (.NOT.ZPERTU) then
      STOP
      endif
      call flush(6)   !tolgo flush
      CALL REIJKL(factor,norb,nsym,itsym,its,nos,jd1)
      finv=1.d0/factor
      CALL MEMCIP('IJKLP',.true.)
      if (ztoul) then
      CALL ITIJKLt(factor,norb,nsym,nos,nijkl,
     * ijkl4,ijkl2,rijkl,ijkl3)
      else
      CALL ITIJKL(factor,norb,nsym,nos,nijkl,
     * ijkl4,ijkl2,rijkl,ijkl3)
      endif
      CALL IJKF(2,TOCOM)
c--renzo      
c      CALL IJKF(1,TOCOM)
      DELTA=DELTA/FINV
      DO 350 M=1,METAT
      E2MP(M)=0.D0
      psimp(M)=0.D0
      psien(M)=0.D0
      do 340 i=1,14
      psimpe(m,i)=0.d0
  340 psiene(m,i)=0.d0
      E2EN(M)=0.D0
      do 350 i=1,12
      e2mpsh(m,i)=0.d0
      psimpsh(m,i)=0.d0
      psiensh(m,i)=0.d0
  350 e2ensh(m,i)=0.d0
      do i=1,14
      ntotex(i)=0
      enddo
      NTOT=0
      NF=NCF+1
      ND(NF)=ND(NCF)+NE(NCF)+1
C 
C     I.C. PERTURBATIVE
C 
      if (zdh) then
      if (ndhb.eq.0) then
      write (6,*) 'Manca la definizione del buco profondo'
      stop
      endif
      endif
C 
      if (zddh) then
      if (ndhb1.eq.0.and.ndhb2.eq.0) then
      write (6,*) 'Manca la definizione del buco profondo'
      stop
      endif
      endif
C
C     Singole sostituzioni
C
c      if (iconf(1).eq.0) then
c      ic=0
c      else
c      ic=1-iconf(1)
c      endif
c--renzo--modif--begin     
      do m=1,ncf
         if(iconf(m).eq.0)iconf(m)=1
      enddo
c--renzo--modif--end      
      ic=1
      mcapos=0
      do m=1,ncf
         if(zcapos(m))then
            icf=ic+iconf(m)-1
            mcapos=mcapos+1
            if(thresh.eq.0.d0)then
               value=1.d0
            else
               value=givval(cconf,ic,icf)
            endif
            if (zdh) then
               if (iocc(ndhb,mcapos).ne.1) then
                  print*,'Problema: det.',m,'non ha deep hole'
c--renzo modif
c                  stop 'no deep hole'
                  goto 4698
               endif
            endif
            if (zddh) then
               if (iocc(ndhb1,mcapos).ne.1.and.iocc(ndhb2,mcapos).ne.1)
     $              then
                  print*,'Problema: det.',m,'non ha deep hole'
c--renzo modif
                  stop 'no deep hole'
                  goto 4698
               endif
            endif
            if(iconf(m).lt.132.and.value.gt.thresh)then
               call singsub(M,IC,mcapos,nconf,iocc,norb)
               if(zdebug)then
                  print*,'singsub: capostipite ',mcapos
                  call flush(6)   !tolgo flush
               endif
            endif
            ic=icf+1
         endif
 4698 enddo
      print*,'Done with single substitutions'
      print*,'E2EN=',e2en(1)/factor
      CALL CPUT(1)
      call flush(6)   !tolgo flush
C
C     Doppie sostituzioni
C
      ic=1
      mcapos=0
      do m=1,ncf
         if(zcapos(m))then
            icf=ic+iconf(m)-1
            mcapos=mcapos+1
            if(thresh.eq.0.d0)then
               value=1.d0
            else
               value=givval(cconf,ic,icf)
            endif
            if (zdh) then
               if (iocc(ndhb,mcapos).ne.1) then
                  print*,'Problema: det.',m,'non ha deep hole'
c--renzo modif
c                  stop 'no deep hole'
                  goto 4699
               endif
            endif
            if(iconf(m).lt.132.and.value.gt.thresh)then
c--renzo debug
c               call doubsub(M,IC,mcapos,nconf,iocc,norb)
cele-provvi
c           if (mcapos.eq.11) then
               call doubsub(M,IC,mcapos,nconf,iocc,norb,e2mp,e2en)
c            print '(a,2f15.9)','doub: e2mp,e2en ',e2mp(5)*finv,e2en(5)
c     $           *finv
c          endif
c               if(mcapos.eq.12) stop
c--renzo debug end
               if(zdebug)then
                  print*,'doubsub: capostipite ',mcapos
                  call flush(6)   !tolgo flush
               endif
            endif
            ic=icf+1
         endif
 4699 enddo
      print*,'Done with double substitutions'
      CALL CPUT(1)
C
      write (6,*) 'Contributi dei vari termini della NEVPT: EN2'
      do i=1,8
      write (6,401) i,psienec(1,i),e2enc(1,i)*finv
      enddo
      write (6,*) 
      write (6,*) 'Contributi dei vari termini della NEVPT: MP2'
      do i=1,8
      write (6,401) i,psimpec(1,i),e2mpc(1,i)*finv
      enddo
  401 format (1x,i2,2f20.15)
      write (6,*)
      WRITE (6,1225) NTOTdet,(i,i=1,14),(ntotex(i),i=1,14)
      write (6,*) 'Ntotconf= ',ntotconf
      write (6,*) 'Ntotcapf= ',ntotcap
      WRITE (6,1400) (i,i=1,14)
      DO 400 II=1,METAT
      write (6,1401) ii,' MPB',psimp(ii),(psimpe(ii,i),i=1,14)
      if (zenb) then
      write (6,1401) ii,' ENB',psien(ii),(psiene(ii,i),i=1,14)
      else
      write (6,1401) ii,'ENVP',psien(ii),(psiene(ii,i),i=1,14)
      endif
  400 continue
      WRITE (6,'(/)')
      WRITE (6,2000)
      WRITE (6,1060)
      WRITE (6,2000)
      WRITE (6,'(/)')
      WRITE (6,1105) ESCF
      if (zenb) WRITE (6,1070) ' ENB'
      if (.not.zenb) WRITE (6,1070) 'ENVP'
      DO 420 M=1,METAT
      E(M)=E(M)*FINV
      E2MP(M)=E2MP(M)*FINV 
      E2EN(M)=E2EN(M)*FINV 
  420 WRITE (6,1100) M,E(M),E2MP(M),E2EN(M)
      WRITE (6,'(//)')
c
c
c     Cele: dati sul file07 per la selezione
c
      if(zthird)goto 1717
      if (zbin) then


        inec=-1
        write (7)  inec,' ',(1.d0,i=1,2*metat)
        write (7) zenb,zmppun
        if (zmppun) then
           write (7) (psimp(i),i=1,metat)
           write (7) (e(i)+ESCF,i=1,metat)
           write (7) (e2mp(i),i=1,metat)
           write (7) (shift(i)*finv,i=1,nshift)
           do m=1,metat
           write (7) (psimpsh(m,i),i=1,nshift)
           enddo
        else
           write (7) (psien(i),i=1,metat)
           write (7) (e(i)+ESCF,i=1,metat)
           write (7) (e2en(i),i=1,metat)
           write (7) (shift(i)*finv,i=1,nshift)
           do m=1,metat
           write (7) (psiensh(m,i),i=1,nshift)
           enddo
        endif


      else

 
        if (zmppun) then
           write (7,'(a5)') 'S MPB'
        elseif (zenb) then
           write (7,'(a5)') 'S ENB'
           else
           write (7,'(a6)') 'S ENVP'
        endif
        if (zmppun) then
           write (7,'(6f14.8)') (psimp(i),i=1,metat)
           write (7,'(6f14.8)') (e(i)+ESCF,i=1,metat)
           write (7,'(6f14.8)') (e2mp(i),i=1,metat)
           write (7,'(12f14.8)') (shift(i)*finv,i=1,nshift)
           do m=1,metat
           write (7,'(12f14.8)') (psimpsh(m,i),i=1,nshift)
           enddo
        else
           write (7,'(6f14.8)') (psien(i),i=1,metat)
           write (7,'(6f14.8)') (e(i)+ESCF,i=1,metat)
           write (7,'(6f14.8)') (e2en(i),i=1,metat)
           write (7,'(12f14.8)') (shift(i)*finv,i=1,nshift)
           do m=1,metat
           write (7,'(12f14.8)') (psiensh(m,i),i=1,nshift)
           enddo
        endif


      endif

c
 1717 if (zenb) WRITE (6,1120) ' ENB'
      if (.not.zenb) WRITE (6,1120) 'ENVP'
      DO 690 M=1,METAT
      E(M)=E(M)+ESCF
      E2MP(M)=E(M)+E2MP(M)
      E2EN(M)=E(M)+E2EN(M)
  690 WRITE (6,1100) M,E(M),E2MP(M),E2EN(M)
c
C     hamiltoniano QDPT sul file 4
      write (4)
     * ((brdmp(i,j)*finv*0.5d0,j=1,i),e2mp(i),i=1,metat),
     * ((brden(i,j)*finv*0.5d0,j=1,i),e2en(i),i=1,metat),
     * (psien(i),i=1,metat),(psimp(i),i=1,metat)
      if (nshift.gt.0) call lshift
      WRITE (6,'(//)')
      IF (ZBRD) CALL WRTBRD(zfa)
      IF (ZBRD) CALL WRTBRD(ztr)
      WRITE (6,'(/)')
      CALL CPUT(1)
      IF (.NOT.ZPRT) STOP
      MNOR1=NORB2+1
      DO 920 K=1,NCF
      L1=ND(K)+1
      L2=ND(K)+NE(K)
  920 WRITE (6,'(1x,30i4)') NE(K),ND(K),(TROU(L),PART(L),L=L1,L2)
      DO 930 I=1,MNOR1
  930 WRITE (6,*) IORB(I),ISPIN(I),ITSYM(I),ZOCS(I)
      STOP
 1010 FORMAT (//' ETAT    SEUIL DE SELECTION          E(ORDRE ZERO)',
     *        '       +DELTA(QUASI-DEGEN)'/
     *        '       DES DETERMINANTS (EN)        MPB      ',a4,
     *        '           MPB      ',a4)
 1040 FORMAT (i4,f16.5,9x,2f11.6,2x,2f11.6)
 1060 FORMAT (5X,'CONTRIBUTIONS A L''ENERGIE TOTALE')
 1070 FORMAT (' ETAT',3X,'DIAGONALISATION        PERTURBATION MP-B',
     * '        PERTURBATION ',a4)
 1100 FORMAT (I4,F17.8,3F24.8)
 1105 FORMAT (' REFERENCE: ',F16.8,//)
 1120 FORMAT (' ETAT',4X,'TOTAL SCF+DIAG          TOTAL MP-B    ',
     * '          TOTAL ',a4)
 1225 FORMAT (//5X,'ONT ETE GENERES:',I15,' DETERMINANTS EFFICACES'/
     * 5x,'NEC =',14i12/5x,'DET =',14i12//)
 1400 FORMAT (5X,'NORME DE LA CORRECTION D''ORDRE 1 A LA FONCTION ',
     * 'D''ONDE'/5x,'ETAT            TOTAL   NEC =',i4,14i10)
 1401 format (i8,2x,a4,f12.6,5x,14f10.6)
 2000 FORMAT (5X,32('*'))
      END
C***********************************************************************
      function givval(c,i,j)
      implicit real*8 (a-h,o-z)
      dimension c(*)
      v=0.d0
      do k=i,j
         v=v+c(k)**2
      enddo
      givval=sqrt(v)
      return
      end
c-----------------------------------------------------------------------      
      SUBROUTINE OFILES(zbin,ndims,ndimtp,zcost)
      implicit character*64(f)
      logical*1 zbin,zsign,zcddet,zcost,zwconf
      logical zexist
      character*10000 namel
      character*100 card
      character*1 cp
      common /fil/ file12,file13,file26,zwconf
      ndims=0
      ndimtp=0
      file03=' '
      file04='CIPPI.FILE04'
      file07=' '
      file09='CIPPI.FILE09'
      file11=' '
      file12=' '
      file13='DIAGO.RESTART'
      file14=' '
      file20='CIPPI.FILE20'
      file25='IJKL.FILE25'
      file26=' '
      file28=' '
      file29=' '
      file50='IJKL.FILE50'
      file05=' '
      filedet=' '
      zbin=.true.
      zcost=.false.
      zwconf=.false.
      call lecnam('&FILES',namel,5)
      call namela(namel,file03,'FILE03')
      call namela(namel,file04,'FILE04')
      call namela(namel,file07,'FILE07')
      call namela(namel,file09,'FILE09')
      call namela(namel,file11,'FILE11')
      call namela(namel,file12,'FILE12')
      call namela(namel,file13,'FILE13')
      call namell(namel,zcost,'ZCOST',1)
      if(zcost)call namela(namel,file14,'FILE14')
      if(zcost)call namela(namel,file28,'FILE28')
      call namela(namel,file20,'FILE20')
      call namela(namel,file25,'FILE25')
      call namela(namel,file26,'FILE26')
      call namell(namel,zwconf,'ZWCONF',1)
      if(zwconf)call namela(namel,file29,'FILE29')
      call namela(namel,file50,'FILE50')
      call namela(namel,file05,'FILE05')
      call namela(namel,filedet,'FILEDET')
      call namell(namel,zbin,'ZBIN',1)
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

      if (file07.ne.' ') then
       if (zbin) then
        open(7,file=file07,form='UNFORMATTED',status='UNKNOWN')
       else
        open(7,file=file07,form='FORMATTED',status='UNKNOWN')
       endif
       rewind 7
      endif

      if (file05.ne.' ')then
         open(5,file=file05,form='FORMATTED',status='OLD')
         rewind 5
      endif

      if (file11.ne.' ') then
      open(11,file=file11,form='UNFORMATTED',status='OLD')
      rewind 11
      endif
      open(9,file=file09,form='UNFORMATTED',status='UNKNOWN')
      rewind 9
      if(zcost)then
         if(file14.eq.' ') stop 'manca il file .Info di cost'
         if(file28.eq.' ') stop 'manca il file .Mono di cost'
         open(14,file=file14,form='FORMATTED',status='OLD')
         rewind 14
         open(28,file=file28,form='UNFORMATTED',status='OLD')
         rewind 28
      endif
      write (6,*) 'ZWCONF=',ZWCONF
      if(zwconf)then
         if(file29.eq.' ') stop 'manca il nome per il file 29'
         open(29,file=file29,form='UNFORMATTED',status='UNKNOWN')
         rewind 29
      endif
c  apro il file 20 con una subroutine da compilare a parte per la
c  compatibilita` ibm-dec
      open(20,file=file20,form='UNFORMATTED',status='UNKNOWN')
c     open(20,file=file20,form='UNFORMATTED',status='UNKNOWN')
      rewind 20
      open(25,file=file25,form='UNFORMATTED',status='UNKNOWN')
      inquire(file=file25,exist=zexist)
      if(zexist) rewind 25
      open(50,file=file50,form='UNFORMATTED',status='OLD')
      maxlong=200000000
      call filesplit('OPEN',50,maxlong,idum,idum2)
      call filesplit('REWIND',50,0,idum,idum2)
      rewind 50
      WRITE (6,'(/''    INPUT DATA:''/)')
      write (6,*) ' FILE04 = ',file04
      write (6,*) ' FILE07 = ',file07
      write (6,*) ' FILE11 = ',file11
      write (6,*) ' FILE25 = ',file25
      write (6,*) ' FILE50 = ',file50
      
      zcddet=.false.
      icount=0
      if (filedet.ne.' ') then
        write (6,*) 'Determinants read from formatted file ', filedet
        open(42,file=filedet,form='FORMATTED',status='OLD')
        zcddet=.true.
      else
        write (6,*) 'Determinants read from standard input (5)'
      endif
      if (zcddet) then
   49 read (5,'(a100)',end=50) card
      write (3,'(a100)') card
      icount=icount+1
      if (icount.le.15) write (6,'(1x,a)') trim(card)
      i=index(card,'ICINP')
      if (i.ne.0) icount=0
      goto 49
      endif
        
   50 if (zcddet) then
       read(42,'(a100)',end=99) card
       select case(card(1:2))
      case('1 ')
         card(1:2)='M '
      case('2 ')
         card(1:2)='D '
      case('3 ')
         card(1:2)='T '
      case('4 ')
         card(1:2)='Q '
      case('5 ')
         card(1:2)='P '
      case('6 ')
         card(1:2)='H '
      case('7 ')
         card(1:2)='E '
      case('8 ')
         card(1:2)='O '
      case('9 ')
         card(1:2)='N '
      case('10')
         card(1:2)='X '
      case('11')
         card(1:2)='U '
      case('12')
         card(1:2)='Z '
      end select
      else
       read (5,'(a100)',end=99) card
      endif
      write (3,'(a100)') card
      cp=card(1:1)
      ndimtmp=0
      ndtptmp=0
      
      call c2EDeg(cp,ndimtmp,ndtptmp)
      ndims = ndims + ndimtmp
      ndimtp = ndimtp + ndtptmp
      icount=icount+1
      if (icount.le.15) write (6,'(1x,a80)') card 
      i=index(card,'ICINP')
      if (i.ne.0) icount=0
      goto 50 

   99 if (icount.gt.15) WRITE (6,'(8('' ..... etc''))')

      ndims=ndims+1

      if (zcddet) write (3,*) ' ' ! end of file 

      call flush(6)
      RETURN
      END
C***********************************************************************
C     C2EDeg : Code to Excitation Degree
C     code (i)   : letter for the excitation degree
C     ndims (o)  : incrementor if a determinant code is found
C     ndimtp (o) : incrementor for the hole/particle count 
      subroutine c2EDeg(code,ndims,ndimtp)
      implicit none
      character code
      integer ndims, ndimtp, i
C     in case of new excitation degrees,
C     remember to add both the new letter and the space for char *n
      character*13 codes
      codes='FMDTQPHEONXUZ'
      ndims=0
      ndimtp=0
      
      i=index(codes,code)
c      write (6,*) 'code, i = ', code, i 

      if (i.eq.0) return

      ndims=1
      ndimtp=i-1
      return
      end
C***********************************************************************
      SUBROUTINE WRTBRD(zen)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      PARAMETER (IEQVC=ID6*ID6+1)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      DIMENSION HEFF(ID6,ID6),VEFF(ID6,ID6)
      allocatable work(:)
c      EQUIVALENCE (C(1),HEFF(1,1)),(C(ieqvc),VEFF(1,1))
 1300 FORMAT (5X,'HAMILTONIEN EFFECTIF QDPT (MOLLER-PLESSET BARIC.)')
 1301 FORMAT (5X,'HAMILTONIEN EFFECTIF QDPT (EPSTEIN-NESBET V.P.)')
 1350 FORMAT (1X,8(1X,F15.8))
 2100 FORMAT (5X,50('*'))
      k662=id6*id6*2
c      if (k662.gt.id7) then
c      write (6,*) '   *** ERREUR DIMENSIONS DANS WRTBRD ***'
c      stop  !but we don't care
c      endif
      WRITE (6,2100)
      if (zen) then
      WRITE (6,1301)
      else
      WRITE (6,1300)
      endif
      WRITE (6,2100)
      WRITE (6,'(/)')
      JJ=0
      DO 510 I=1,metat
      ED=E2MP(I)
      IF (ZEN) ED=E2EN(I)
C     IF (I.GT.METAT) ED=E(I)+ESCF
      HEFF(I,I)=ED
c      IM=I-1
c      IF (IM.EQ.0) GOTO 506
      DO 505 J=1,metat
c      JJ=JJ+1
      if (zen) then
c      HEFF(I,J)=brden(i,j)*FINV*0.5D0
      HEFF(I,J)=brden(i,j)*FINV
      else
c      HEFF(I,J)=brdmp(i,j)*FINV*0.5D0
      HEFF(I,J)=brdmp(i,j)*FINV
      endif
  505 continue
cele 30/9/03
      HEFF(I,I)=HEFF(I,I)+e(i)
  506 WRITE (6,1350) (HEFF(I,J),J=1,metat)
  510 CONTINUE
c         allocate(w(norb))
      allocate(work(3*metat))
      do i=1,metat
         do j=1,metat
            veff(i,j)=(heff(i,j)+heff(j,i))*0.5d0  !!poor man's hermitization
         enddo
      enddo
      call dsyev('V','U',metat,veff,id6,wk,work,3*metat,info)
      deallocate(work)
c      CALL RIERA(metat,HEFF,VEFF,wk,EMP,ID6)
      WRITE (6,520)
  520 FORMAT (///5X,'VALEURS ET VECTEURS PROPRES'/)
      JMAX=0
  530 JMIN=JMAX+1
      JMAX=JMAX+8
      IF (JMAX.GT.metat) JMAX=metat
      WRITE (6,1350) (wk(J),J=JMIN,JMAX)
c      WRITE (6,1350) (wk(J)+e(j),J=JMIN,JMAX) ! e(j) already added 11/03/04
      WRITE (6,'(/)')
      DO 540 I=1,metat
  540 WRITE (6,1350) (VEFF(I,J),J=JMIN,JMAX)
      WRITE (6,'(//)')
      IF (JMAX.LT.METAT) GOTO 530
      RETURN
C-----------------------------------------------------------------------
      ENTRY LSHIFT
      WRITE (6,'(//)')
      WRITE (6,2100)
      WRITE (6,600)
  600 FORMAT (15X,'PERTURBATION AVEC LEVEL SHIFT')
      WRITE (6,2100)
      if (zenb) WRITE (6,610) 'EN-B'
      if (.not.zenb) WRITE (6,610) 'ENVP'
  610 FORMAT (/' ETAT       LEVEL SHIFT          TOTAL MP-B    ',
     * '            TOTAL ',a4,'                 MP NORM',
     * '                EN NORM')
      do 620 m=1,metat
      WRITE (6,'(/)')
      write (6,630) m,0.d0,e2mp(m),e2en(m),psimp(m),psien(m)
      do 620 i=1,nshift
      e2mpsh(m,i)=e(m)+e2mpsh(m,i)*finv
      e2ensh(m,i)=e(m)+e2ensh(m,i)*finv
  620 write (6,630) m,shift(i)*finv,e2mpsh(m,i),e2ensh(m,i),
     * psimpsh(m,i),psiensh(m,i)
  630 FORMAT (I4,F17.8,4F24.8)
      RETURN
      END
C***********************************************************************
      FUNCTION HNTD(II,JJ,JINT,KINT)
C 
C                  SSP DE RECHERCHE DES ELEMENTS D'INTERACTION
C                  ENTRE DEUX DETERMINANTS, NUMEROTES II ET JJ
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION ITA(99),ITB(99),NV1(2),NV2(2),NATU1(2),NATU2(2)
      EQUIVALENCE(NV11,NV1(1)),(NV12,NV1(2)),(NV21,NV2(1)),(NV22,NV2(2))
      KR=II
      KL=JJ
      JINT=0
      KINT=0
C 
C                  LE DETERMINANT - 1 - EST CHOISI COMME LE
C                  DETERMINANT AYANT LE PLUS D'O.M. EXCITEES
C 
      IF (NE(KR)-NE(KL)) 5,7,7
    5 KKKK=KR
      KR=KL
      KL=KKKK
    7 NE1=NE(KR)
C                  TEST SUR LE FONDAMENTAL
      IF (NE1.EQ.0) GOTO 44
      NE2=NE(KL)
C                  NBDIF : EST LE NOMBRE DE SPIN-ORBITALES
C                          DIFFERENTES
      NBDIF=NE1-NE2
      IF (NBDIF.GT.2) GOTO 24
C                  CONSTRUCTION POUR LE DETERMINANT -1- ET -2-
C                  DE IT ET NS, QUI CONTIENNENT LES NUMEROS
C                  DES ORBITALES ET LEURS SPINS
C 
      ND1=ND(KR)
      ND2=ND(KL)
      DO 15 J=1,NE1
      ITA(J)=TROU(J+ND1)
      ITA(J+NE1)=PART(J+ND1)
   15 CONTINUE
C 
C     ZSIG :  EST LA PARITE DU  NOMBRE DE CROISEMENTS
C             DANS LE DIAGRAMME D INTERACTION
      ZSIG=.FALSE.
C                  SI NE2=0 LE DETERMINANT -2- EST LE FONDAMENTAL
      IF (NE2.EQ.0) GOTO 125
  121 CONTINUE
      DO 25 J= 1,NE2
      ITB(J)=TROU(J+ND2)
   25 ITB(J+NE2)=PART(J+ND2)
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
   12 CONTINUE
      NBDIF=NBDIF+1
      IF (NBDIF.GT.2) GOTO 24
      NV1(NBDIF)=NIT1
      NATU1(NBDIF)=K
      NATU2(NBDIF)=I
   10 CONTINUE
    8 CONTINUE
  125 IF (NBDIF.LE.0) GOTO 44
C 
C                   NAR : EXCITATION SUPLEMENTAIRE DU DETERMINANT -1-
C                        PAR RAPPORT AU DETERMINANT -2-
   26 NAR=NE1-NE2
      IF (NAR.LE.0) GOTO 28
   30 DO 32 I=1,NAR
      NET=NE1+1-I
      NV1(I)=ITA(NET)
      NV2(I)=ITA(NET+NE1)
   32 CONTINUE
   28 NBAR=NBDIF-NAR
      IF (NBAR.LE.0) GOTO 34
   36 NAR1=NAR+1
      DO 38 I=NAR1,NBDIF
      K=NATU1(I)
      NI=NATU2(I)
      NU=NV1(I)
      NV2(I)=ITA(NI+NE1*(K-1))
      IF (K.EQ.2) GOTO 38
   40 NV1(I)=NV2(I)
      NV2(I)=NU
      ZSIG=.NOT.ZSIG
   38 CONTINUE
   34 CONTINUE
C                  CALCUL DE L ELEMENT DE MATRICE
    2 HA=0.D0
      N1=IORB(NV11)
      N2=IORB(NV21)
      NS1=ISPIN(NV11)
      NS2=ISPIN(NV21)
      IF (NBDIF.EQ.1) GOTO 58
C 
C     LES -2- DETERMINANTS DIFFERENT PAR -2- SPIN-ORBITALES
C 
   60 N3=IORB(NV12)
      N4=IORB(NV22)
      NS3=ISPIN(NV12)
      NS4=ISPIN(NV22)
      IF (NS1.NE.NS2) GOTO 200
      HA=AI(N1,N2,N3,N4)
      IF (ZSIG) NNJK=-NNJK
      JINT=NNJK
  200 IF (NS1.NE.NS4) GOTO 75
      HA=HA-AI(N1,N4,N3,N2)
      IF (ZSIG) NNJK=-NNJK
      KINT=NNJK
      GOTO 75
C 
C     LES -2- DETERMINANTS DIFFERENT PAR -1- SPIN-ORBITALE
   58 CONTINUE
      IF (N2.GE.N1) GOTO 57
      N12=NUM(N1)+N2
      GOTO 56
   57 N12=NUM(N2)+N1
   56 HA=F(N12)
      DO 72 I=1,NE1
      NP=IORB(ITA(I))
      NSU=ISPIN(ITA(I))
      HA=HA-AI(N1,N2,NP,NP)
      IF (NS1.EQ.NSU) HA=HA+AI(N1,NP,N2,NP)
      NP=IORB(ITA(I+NE1))
      IF (NP.GT.NORB) GOTO 72
      NSU=ISPIN(ITA(I+NE1))
      HA=HA+AI(N1,N2,NP,NP)
      IF (NS1.EQ.NSU) HA=HA-AI(N1,NP,N2,NP)
   72 CONTINUE
   75 HNTD=HA
      IF (ZSIG) HNTD=-HA
      RETURN
   24 HNTD=0.D0
      RETURN
 44   WRITE (6,1006) II,JJ,NE1,NE2
 1006 FORMAT (//5X,'IDENTITE DANS HNTD'/' II,JJ,NE1,NE2 =',4I4)
      WRITE (6,1007) (ITA(K),K=1,NE1)
      WRITE (6,1008) (ITA(K+NE1),K=1,NE1)
      WRITE (6,1007) (ITB(K),K=1,NE2)
      WRITE (6,1008) (ITB(K+NE2),K=1,NE2)
 1007 FORMAT (' TROU',14I4)
 1008 FORMAT (' PART',14I4)
      STOP
      END
C***********************************************************************
      SUBROUTINE IJKF(IPAS,TOCOM)
C 
C     J    INTEGRALES COULOMBIENNES STOCKEES DANS LE TABLEAU AJ
C     K    INTEGRALES D'ECHANGE     STOCKEES DANS LE TABLEAU AK
C     F    ELEMENTS DE L'OPERATEUR DE FOCK
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
c      DIMENSION FDIAG(ID3),TOCOM(*)
      DIMENSION TOCOM(*)
      allocatable fdiag(:)
      allocatable fock(:,:)
      character*80 label
      save fdiag
      dimension nperm(1000) ! da rendere dinamico=id1 (13/01/03)
      common /toul/ ztoul
      common /fil/ file12,file13,file26,zwconf
      character*64 file12,file13,file26
cele 15-11-04      if(.not.allocated(fdiag))allocate(fdiag(norb*2+1))
      if(.not.allocated(fdiag))allocate(fdiag(norb*2+2))
      NTER=NORB*(NORB+1)/2
      allocate(fock(norb,norb))
      if (ztoul) then
        do j=1,norb
        nperm(j)=j
        enddo
        if (file26.eq.' ') then
        write (6,*) ' Errore: manca il file .mono di molcsd'
        call flush(6)   !tolgo flush
        stop
        endif
        open(26,file=file26,form='UNFORMATTED',status='UNKNOWN')
        rewind 26
c      read(26)
c      do
c         read(26) label
c         if(index(label,'==ONEINT MATRIX (TOTAL)').ne.0) exit
c         read(26)
c         read(26)
c      enddo

c      read(26) ! norbttMono,nsymMono,isymMono(1:nsymMono)
c      read(26) F(1:nter)
        READ (26) (F(L),L=1,NTER)
      else
        call filesplit('READ',50,8*nter,ichange,idum2)
        READ (50) (F(L),L=1,NTER)
      endif
cele
      kl=0
      DO LL=1,NORB
c      write (6,444) ll,(f(ij),ij=kl+1,kl+ll)
      kl=kl+ll
      enddo
      if (ztoul) then
      kl=0
      DO LL=1,NORB
      DO KK=1,LL
      kl=kl+1
      fock(nperm(ll),nperm(kk))=f(kl)*factor
      fock(nperm(kk),nperm(ll))=f(kl)*factor
      enddo
      enddo
c      do ij=1,norb
c      write (6,444) ij,(fock(ij,j)/factor,j=1,norb)
c      enddo
      val=0
449   DO 500 ij=1,nocb
      VAL=VAL+Fock(ij,ij)*2.D0
      TT=AI(ij,ij,ij,ij)
      Fock(ij,ij)=Fock(ij,ij)+TT
      VAL=VAL+TT
      IM=ij-1
      IF(IM.EQ.0) GO TO 500
      DO 460 K=1,IM
      TT=AI(ij,ij,K,K)
      TK=AI(ij,K,ij,K)
      TT=TT+TT-TK
      Fock(ij,ij)=Fock(ij,ij)+TT
      Fock(K,K)=Fock(K,K)+TT
      VAL=VAL+TT+TT
      DO 450 J=1,nocb
450   Fock(ij,K)=Fock(ij,K)+AI(ij,K,J,J)*2.D0-AI(ij,J,K,J)
460   continue
500   CONTINUE
      I1=NOCB+1
      DO 600 ij=I1,NORB
      DO 600 J=1,ij
      TT=0.d0
      DO 550 K=1,NOCB
550   TT=TT+2.D0*AI(ij,J,K,K)-AI(ij,K,J,K)
600   Fock(ij,J)=Fock(ij,J)+TT
c      do ij=1,norb
c      write (6,444) ij,(fock(ij,j)/factor,j=1,norb)
c444   format (i4,3x,10F12.8)
c      enddo
c      write (6,'(''Ene. Eletr.'',f15.8)') val/factor
      rewind 25
      read (25)
      read (25)
      read (25) edum,ENUC
      write (6,*) 'Enuc=',enuc
      backspace 25
      escf=enuc+val/factor
      write (25) escf,ENUC,(nperm(j),j=1,norb)
c      write (6,'(''Ene. Nucl.'',f15.8)') enuc
      write (6,'(/''     REFERENCE ENERGY'',f15.8/)') ESCF
      kl=0
      DO LL=1,NORB
      DO KK=1,LL
      kl=kl+1
      f(kl)=fock(ll,kk)/factor
      enddo
      enddo
      call filesplit('READ',50,8*nter,ichange,idum2)
      write (50) (F(L),L=1,NTER)
      endif ! if (ztoul)

      deallocate (fock)


      L=0
      KL=0
      DO 3200 LL=1,NORB
      DO 3100 KK=1,LL
      KL=KL+1
      F(KL)=F(KL)*factor
      AJ(KL)=AI(LL,LL,KK,KK)
 3100 AK(KL)=AI(LL,KK,LL,KK)
c      write (6,'(a,2i3,3f12.6)') 'F,AJ,AK=',ll,kk,f(kl)*finv,aj(kl)*finv
c     $     ,ak(kl)*finv
c      call flush(6)
      L=L+LL
      FDIAG(LL)=F(L)
 3200 FDIAG(LL+NORB)=F(L)
      IF (.NOT.ZPRT.OR.IPAS.GT.1) GOTO 3400
      WRITE (6,3250)
 3250 FORMAT (/5X,'INTEGRALES J')
      MAX=0
      DO 3260 L=1,NORB
      MIN=MAX+1
      MAX=MAX+L
      WRITE (6,3270) L
 3260 WRITE (6,3280) (AJ(KL),KL=MIN,MAX)
 3270 FORMAT (I8)
 3280 FORMAT (1X,10F13.0)
      WRITE (6,3350)
 3350 FORMAT (/5X,'INTEGRALES K')
      MAX=0
      DO 3360 L=1,NORB
      MIN=MAX+1
      MAX=MAX+L
      WRITE (6,3270) L
 3360 WRITE (6,3280) (AK(KL),KL=MIN,MAX)
 3400 IF (.NOT.ZION) GOTO 3600
      NTER=NTER+1
      NN=NTER+NORB+NORB+2
      DO 3500 L=NTER,NN
      F(L)=0.D0
      AJ(L)=0.D0
 3500 AK(L)=0.D0
      FDIAG(NORB2+1)=0.D0
      FDIAG(NORB2+2)=0.D0
 3600 IF (IPAS.GT.1) RETURN
      rewind 25
      read (25) idum,idum,idum,ndegen
      read (25)
      do l=1,ndegen
      read (25)
      enddo
      if (.not.ztoul) then
      if (tocom(1).lt.-1.d5) then
      READ (25) ESCF,(TOCOM(L),L=1,NORB)
      else
      READ (25) ESCF
      endif
      endif ! if(..not.ztoul)
      DO 4300 L=1,NOCB
 4300 TOCOM(L)=TOCOM(L)-2.D0
      DO 4500 L=1,NORB
      T=FDIAG(L)
      LJ=NUM(L)
      DO 4400 J=1,NORB
      LJ=LJ+1
      IF (J.GT.L) LJ=NUM(J)+L
 4400 T=T+TOCOM(J)*(AJ(LJ)-0.5D0*AK(LJ))
      FMPB(L)=T
 4500 FMPB(L+NORB)=T
      IF (ZION) FMPB(NORB2+1)=0.D0
      IF (ZION) FMPB(NORB2+2)=0.D0
      WRITE (6,4600)
 4600 FORMAT (/'    O.M.  ACT.  SYMETRIE   N. D''OCCUP.      ENERGIE',
     * '      INT. COULOMB')
      ntrpa=nd(ncf)+ne(ncf)
      do 4700 l=1,norb
	zact=.false.
	do j=1,ntrpa
	if (trou(j).eq.l) zact=.true.
	if (part(j).eq.l) zact=.true.
	if (trou(j)-norb.eq.l) zact=.true.
	if (part(j)-norb.eq.l) zact=.true.
	enddo
cele 15-11-04      if (l.eq.1) zwww=.true.
cele 15-11-04      if (l.eq.nocb) zwww=.true.
cele 15-11-04      if (l.eq.nocb+1) zwww=.true.
cele 15-11-04      if (l.eq.initb) zwww=.true.
cele 15-11-04      if (l.eq.norb) zwww=.true.
cele 15-11-04      if (l.ge.2.and.itsym(l).ne.itsym(l-1)) zwww=.true.
cele 15-11-04      if (l.lt.norb.and.itsym(l).ne.itsym(l+1)) zwww=.true.
      zwww=.true.
cele
      ll=(l*l+l)/2
      if ((zwww.or.zact).and.l.le.nocb) write (6,4800)
     * 'occ ',l,zact,itsym(l),tocom(l)+2.d0,fmpb(l)*finv,aj(ll)*finv
      if ((zwww.or.zact).and.l.gt.nocb) write (6,4800)
     * 'virt',l,zact,itsym(l),tocom(l),fmpb(l)*finv,aj(ll)*finv
 4700 continue
 4800 format (1x,a4,i4,l4,i8,f15.6,2f16.8)
      RETURN
C-----------------------------------------------------------------------
      ENTRY HDIG(I,ZREF,REF,AB)
C 
C     CE SSP FONCTION CALCULE L'ENERGIE DU DETERMINANT NUMEROTE  I,
C     PAR RAPPORT A L'ENERGIE DE L'ETAT FONDAMENTAL
C 
      NEI=NE(I)
      AB=0.D0
      IF (NEI.EQ.0) RETURN
      L1=ND(I)
      L2=L1+NEI
      L1=L1+1
      LL=L2
      IF (ZREF) LL=L1+1
      IF (ZREF) AB=REF
      DO 1200 L=L1,LL
      N1P=PART(L)
      N1T=TROU(L)
      AB=AB+FDIAG(N1P)-FDIAG(N1T)
      NS1T=ISPIN(N1T)
      NS1P=ISPIN(N1P)
      N1T=IORB(N1T)
      N1P=IORB(N1P)
      N1TP=NUM(N1P)+N1T
      AB=AB-AJ(N1TP)
      IF (NS1T.EQ.NS1P) AB=AB+AK(N1TP)
      IF (L.EQ.L2) GOTO 1200
      LP=L+1
      DO 1100 J=LP,L2
      N2T=TROU(J)
      N2P=PART(J)
      NS2T=ISPIN(N2T)
      NS2P=ISPIN(N2P)
      N2T=IORB(N2T)
      N2P=IORB(N2P)
      IF (N2T.GT.N1T) GOTO 1010
      N12T=NUM(N1T) +N2T
      GOTO 1020
1010  N12T=NUM(N2T)+N1T
1020  IF (N2P.GT.N1P) GOTO 1030
      N12P=NUM(N1P)+N2P
      GOTO 1040
1030  N12P=NUM( N2P)+N1P
1040  N1TP=NUM(N2P)+N1T
      N2TP=NUM(N1P)+N2T
      CA=AJ(N12T)+AJ(N12P)-AJ(N1TP)-AJ(N2TP)
      IF (NS1T.EQ.NS2T) CA=CA-AK(N12T)
      IF (NS1P.EQ.NS2P) CA=CA-AK(N12P)
      IF (NS1T.EQ.NS2P) CA=CA+AK(N1TP)
      IF (NS1P.EQ.NS2T) CA=CA+AK(N2TP)
1100  AB=AB+CA
1200  CONTINUE
      RETURN
C-----------------------------------------------------------------------
      ENTRY HMP(I,AC)
C 
C     ENERGIE DU DETERMINANT IIII DANS LA PARTITION MOLLER-PLESSET
C 
      AC=0.D0
      NEI=NE(I)
      IF (NEI.EQ.0) RETURN
      L1=ND(I)
      L2=L1+NEI
      L1=L1+1
      DO 2200 L=L1,L2
      N1P=PART(L)
      N1T=TROU(L)
 2200 AC=AC+FMPB(N1P)-FMPB(N1T)
      RETURN
      END
C***********************************************************************
      FUNCTION WRTLOG(TF)
      CHARACTER*1 WRTLOG
      LOGICAL*1 TF
      WRTLOG='.'
      IF (TF) WRTLOG='X'
      RETURN
      END
C***********************************************************************
      SUBROUTINE WRTDET(KK,M,MM,EI,AC, ZQUA,zmono)
C 
C     IMPRESSION ET PERFORATION DE DETERMINANTS
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      common/third/zthird,zmp3,zen3
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      CHARACTER*1 ORD(13),cspin(2),arrow(4)
      dimension ecaen(id6),dcaen(id6),ecamp(id6),bcalc(id6),dcaef(id6)
      character*60 ctrou,cpart,cctrou,ccpart
      DATA ORD/'F','M','D','T','Q','P','H','E','O','N','X','U','Z'/
      data cspin/'-','+'/
      data arrow/' ','=','=','>'/
      ZWRT=.TRUE.
      NEC=NE(KK)
      if (nec.gt.12.) return
      MONO1=ND(KK)
      MONO2=MONO1+NEC
      MONO1=MONO1+1
      ctrou=
     * '                                                            '
      cctrou=
     * '                                                            '
      cpart=
     * '                                                            '
      ccpart=
     * '                                                            '
      WRITE (cctrou,40)
     * (IORB(TROU(I)),cspin(ispin(TROU(I))+1),I=MONO1,MONO2)
      WRITE (ccpart,40)
     * (IORB(PART(I)),cspin(ispin(PART(I))+1),I=MONO1,MONO2)
   40 format(12(1x,i3,a1:))
      zet=.true.
      nblt=0
      nblp=0
      lt=0
      lp=0
      do 50 i=1,60
      if (cctrou(i:i).eq.' ') then
      nblt=nblt+1
      if (nblt.eq.1) then
      lt=lt+1
      ctrou(lt:lt)=' '
      endif
      else
      nblt=0
      lt=lt+1
      ctrou(lt:lt)=cctrou(i:i)
      endif
      if (ccpart(i:i).eq.' ') then
      nblp=nblp+1
      if (nblp.eq.1) then
      lp=lp+1
      cpart(lp:lp)=' '
      endif
      else
      nblp=0
      lp=lp+1
      cpart(lp:lp)=ccpart(i:i)
      endif
   50 continue
      if (ctrou(lt:lt).eq.' ') lt=lt-1
      if (cpart(lp:lp).eq.' ') lp=lp-1
      ltot=lt+lp
      IF (KK.LE.NCF) GOTO 200
      IF (.NOT.ZQUA) GOTO 100
C 
C     QUASI DEGENERESCENCE
C 
      EE4=EI*FINV
      EE5=AC*FINV
      do 80 m=1,metat
      EE1=ECAEN(m)*FINV
      EE2=ECAMP(m)*FINV
      IF (ZET) then
         if (ltot.le.64) then
         WRITE (6,60) (ctrou(i:i),i=1,lt),arrow,(cpart(i:i),i=1,lp),
     *   (' ',i=1,64-ltot),M,DCAEN(m),EE1,DCAEF(m),EE2,EE4,EE5
         else
         WRITE (6,60) (ctrou(i:i),i=1,lt),arrow,
     *   ('  ',i=1,64-lt),M,DCAEN(m),EE1,DCAEF(m),EE2,EE4,EE5
         WRITE (6,60) (' ',i=1,18),(cpart(i:i),i=1,lp)
         endif
      else
      WRITE (6,70) M,DCAEN(m),EE1,DCAEF(m),EE2,EE4,EE5
      endif
   60 FORMAT (1X,68A1,1X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,
     *        'QUASI DEG',F6.3,1X,F6.3)
   70 FORMAT (70X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,
     *        'QUASI DEG',F6.3,1X,F6.3)
      zet=.false.
 80   continue
      if(zthird)goto 1717
      if (zpun.and.nec.le.12) then
      if (zmppun) then
        if (nec.lt.7) then
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt),
     *    (cpart(i:i),i=1,lp),(' ',i=1,60-ltot),'Q',
     *    (dabs(dcaef(i)),i=1,metat),(finv*ecamp(i),i=1,metat),
     *    ((dsqrt(shmp(m,i)),m=1,metat),i=1,nshift)
        else
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt)
          write (7,310) '>',(cpart(i:i),i=1,lp),(' ',i=1,60-lp),'Q',
     *   (dabs(dcaef(i)),i=1,metat),(finv*ecamp(i),i=1,metat),
     *   ((dsqrt(shmp(m,i)),m=1,metat),i=1,nshift)
        endif
      else
        if (nec.lt.7) then
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt),
     *    (cpart(i:i),i=1,lp),(' ',i=1,60-ltot),'Q',
     *   (dabs(dcaen(i)),i=1,metat),(finv*ecaen(i),i=1,metat),
     *   ((dsqrt(shen(m,i)),m=1,metat),i=1,nshift)
        else
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt)
          write (7,310) '>',(cpart(i:i),i=1,lp),(' ',i=1,60-lp),'Q',
     *   (dabs(dcaen(i)),i=1,metat),(finv*ecaen(i),i=1,metat),
     *   ((dsqrt(shen(m,i)),m=1,metat),i=1,nshift)
       endif
      endif
      endif
 1717 return
C 
C     CONTRIBUTION PERTURBATIVE IMPORTANTE
C 
  100 EE4=EI*FINV
      EE5=AC*FINV
      do 130 m=1,metat
      ZWRT=dabs(dcaen(m)).GE.TESTE(M).OR.NEC.GT.12
      EE1=ECAEN(m)*FINV
      EE2=ECAMP(m)*FINV
      EE3=BCALC(m)*FINV
      IF (ZWRT) then
      IF (ZET) then
         if (ltot.le.64) then
         WRITE (6,110) (ctrou(i:i),i=1,lt),arrow,(cpart(i:i),i=1,lp),
     *   (' ',i=1,64-ltot),M,DCAEN(m),EE1,DCAEF(m),EE2,EE3,EE4,EE5
         else
         WRITE (6,110) (ctrou(i:i),i=1,lt),arrow,
     *   (' ',i=1,64-lt),M,DCAEN(m),EE1,DCAEF(m),EE2,EE3,EE4,EE5
         WRITE (6,110) (' ',i=1,18),(cpart(i:i),i=1,lp)
         endif
      else
      WRITE (6,120) M,DCAEN(m),EE1,DCAEF(m),EE2,EE3,EE4,EE5
      endif
  110 FORMAT (1X,68A1,
     * 1X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,F8.5,1X,F6.3,1X,F6.3)
  120 FORMAT (70X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,F8.5,
     *        1X,F6.3,1X,F6.3)
      zet=.false.
      endif
  130 continue
c
c     Se zauto=T, in seguito cipselx eseguira` un calcolo con 
c     zauto=T, quindi calcoliamo subito il contributo di tutti i
c     determinanti permutati di spin. Un solo determinante verra` scritto
c     e il suo coefficiente sara` uguale alla radice quadrata della somma
c     dei quadrati di tutti i determinanti con la stessa occupazione
c     orbitalica (cosi quando cipselx fa il quadrato ottiene la somma dei
c     quadrati). Attenzione al cipsel tradizionale ....
c
      if (.not.(zpun.and.nec.le.12)) return
      zstamp=.true.
cele      if (zauto) then
cele      call config(kk,dcaen,ecaen,dcaef,ecamp,zstamp,zmono,zdav,zbin)
cele      endif
c
      if(zthird)goto 1818
      if (zstamp) then
      if (zmppun) then
        if (nec.lt.7) then
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt),
     *    (cpart(i:i),i=1,lp),(' ',i=1,60-ltot),' ',
     *    (dabs(dcaef(i)),i=1,metat),(finv*ecamp(i),i=1,metat),
     *    ((dsqrt(shmp(m,i)),m=1,metat),i=1,nshift)
        else
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt)
          write (7,310) '>',(cpart(i:i),i=1,lp),(' ',i=1,60-lp),' ',
     *    (dabs(dcaef(i)),i=1,metat),(finv*ecamp(i),i=1,metat),
     *    ((dsqrt(shmp(m,i)),m=1,metat),i=1,nshift)
        endif
      else
        if (nec.lt.7) then
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt),
     *     (cpart(i:i),i=1,lp),(' ',i=1,60-ltot),' ',
     *   (dabs(dcaen(i)),i=1,metat),(finv*ecaen(i),i=1,metat),
     *   ((dsqrt(shen(m,i)),m=1,metat),i=1,nshift)
        else
          write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt)
          write (7,310) '>',(cpart(i:i),i=1,lp),(' ',i=1,60-lp),' ',
     *   (dabs(dcaen(i)),i=1,metat),(finv*ecaen(i),i=1,metat),
     *   ((dsqrt(shen(m,i)),m=1,metat),i=1,nshift)
        endif
      endif
      endif
 1818 return
C
C     FONCTION D'ONDE VARIATIONNELLE
C
  200 if (zqua) then
      if (nec.eq.0) then
      write (6,210) kk,(wk(i),i=m,mm)
      else
         if (ltot.le.54) then
         WRITE (6,220) (ctrou(i:i),i=1,lt),arrow,(cpart(i:i),i=1,lp),
     *   (' ',i=1,54-ltot),kk,(wk(i),i=m,mm)
         else
	 if (lt.gt.54) lt=54
	 if (lp.gt.50) lp=54
         WRITE (6,220) (ctrou(i:i),i=1,lt),arrow
         WRITE (6,220) (' ',i=1,2),(cpart(i:i),i=1,lp),
     *   (' ',i=1,56-lp),kk,(wk(i),i=m,mm)
         endif
      endif
  210 FORMAT (1X,'Fondamental',47x,i6,1x,6f11.8)
  220 FORMAT (1X,58A1,i6,1x,6f11.8)
      endif
      if(zthird)return
      if (.not.zpun.or..not.zmono) return
      if (nec.lt.7) then
      write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt),
     * (cpart(i:i),i=1,lp),(' ',i=1,60-ltot),'S',
     * (abs(c(ncf*(k-1)+kk)),k=1,metat),(0.d0,k=1,metat),
     * ((0.d0,i=1,metat),j=1,nshift)
      else
      write (7,310) ord(nec+1),(ctrou(i:i),i=1,lt)
      write (7,310) '>',(cpart(i:i),i=1,lp),(' ',i=1,60-lp),
     * 'S',(abs(c(ncf*(k-1)+kk)),k=1,metat),(0.d0,k=1,metat),
     * ((0.d0,i=1,metat),j=1,nshift)
      endif
  240 continue
      return
  310 FORMAT (61a1,a1,99f15.10)
      END
C***********************************************************************
      SUBROUTINE WRTBIN(KK,M,MM,EI,AC, ZQUA,zmono)
C 
C     IMPRESSION ET PERFORATION DE DETERMINANTS
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      common/third/zthird,zmp3,zen3
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      CHARACTER*1 ORD(13),cspin(2),arrow(4)
      dimension ecaen(id6),dcaen(id6),ecamp(id6),bcalc(id6),dcaef(id6)
      character*60 ctrou,cpart,cctrou,ccpart
      DATA ORD/'F','M','D','T','Q','P','H','E','O','N','X','U','Z'/
      data cspin/'-','+'/
      data arrow/' ','=','=','>'/
      ZWRT=.TRUE.
      NEC=NE(KK)
      if (nec.gt.12.) return
      MONO1=ND(KK)
      MONO2=MONO1+NEC
      MONO1=MONO1+1
      ctrou=
     * '                                                            '
      cctrou=
     * '                                                            '
      cpart=
     * '                                                            '
      ccpart=
     * '                                                            '
      WRITE (cctrou,40)
     * (IORB(TROU(I)),cspin(ispin(TROU(I))+1),I=MONO1,MONO2)
      WRITE (ccpart,40)
     * (IORB(PART(I)),cspin(ispin(PART(I))+1),I=MONO1,MONO2)
   40 format(12(1x,i3,a1:))
      zet=.true.
      nblt=0
      nblp=0
      lt=0
      lp=0
      do 50 i=1,60
      if (cctrou(i:i).eq.' ') then
      nblt=nblt+1
      if (nblt.eq.1) then
      lt=lt+1
      ctrou(lt:lt)=' '
      endif
      else
      nblt=0
      lt=lt+1
      ctrou(lt:lt)=cctrou(i:i)
      endif
      if (ccpart(i:i).eq.' ') then
      nblp=nblp+1
      if (nblp.eq.1) then
      lp=lp+1
      cpart(lp:lp)=' '
      endif
      else
      nblp=0
      lp=lp+1
      cpart(lp:lp)=ccpart(i:i)
      endif
   50 continue
      if (ctrou(lt:lt).eq.' ') lt=lt-1
      if (cpart(lp:lp).eq.' ') lp=lp-1
      ltot=lt+lp
      IF (KK.LE.NCF) GOTO 200
      IF (.NOT.ZQUA) GOTO 100
C 
C     QUASI DEGENERESCENCE
C 
      EE4=EI*FINV
      EE5=AC*FINV
      do 80 m=1,metat
      EE1=ECAEN(m)*FINV
      EE2=ECAMP(m)*FINV
      IF (ZET) then
         if (ltot.le.64) then
         WRITE (6,60) (ctrou(i:i),i=1,lt),arrow,(cpart(i:i),i=1,lp),
     *   (' ',i=1,64-ltot),M,DCAEN(m),EE1,DCAEF(m),EE2,EE4,EE5
         else
         WRITE (6,60) (ctrou(i:i),i=1,lt),arrow,
     *   ('  ',i=1,64-lt),M,DCAEN(m),EE1,DCAEF(m),EE2,EE4,EE5
         WRITE (6,60) (' ',i=1,18),(cpart(i:i),i=1,lp)
         endif
      else
      WRITE (6,70) M,DCAEN(m),EE1,DCAEF(m),EE2,EE4,EE5
      endif
   60 FORMAT (1X,68A1,1X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,
     *        'QUASI DEG',F6.3,1X,F6.3)
   70 FORMAT (70X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,
     *        'QUASI DEG',F6.3,1X,F6.3)
      zet=.false.
 80   continue
      if(zthird)goto 1717
      if (zpun.and.nec.le.12) then
      if (zmppun) then
       write (7) nec,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
     *      'Q',(dabs(dcaef(i)),i=1,metat),(finv*ecamp(i),i=1,metat),
     *    ((dsqrt(shmp(m,i)),m=1,metat),i=1,nshift)
      else
       write (7) nec,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
     *      'Q',(dabs(dcaen(i)),i=1,metat),(finv*ecaen(i),i=1,metat),
     *   ((dsqrt(shen(m,i)),m=1,metat),i=1,nshift)
      endif
      endif
 1717 return
C 
C     CONTRIBUTION PERTURBATIVE IMPORTANTE
C 
  100 EE4=EI*FINV
      EE5=AC*FINV
      do 130 m=1,metat
      ZWRT=dabs(dcaen(m)).GE.TESTE(M).OR.NEC.GT.12
      EE1=ECAEN(m)*FINV
      EE2=ECAMP(m)*FINV
      EE3=BCALC(m)*FINV
      IF (ZWRT) then
      IF (ZET) then
         if (ltot.le.64) then
         WRITE (6,110) (ctrou(i:i),i=1,lt),arrow,(cpart(i:i),i=1,lp),
     *   (' ',i=1,64-ltot),M,DCAEN(m),EE1,DCAEF(m),EE2,EE3,EE4,EE5
         else
         WRITE (6,110) (ctrou(i:i),i=1,lt),arrow,
     *   (' ',i=1,64-lt),M,DCAEN(m),EE1,DCAEF(m),EE2,EE3,EE4,EE5
         WRITE (6,110) (' ',i=1,18),(cpart(i:i),i=1,lp)
         endif
      else
      WRITE (6,120) M,DCAEN(m),EE1,DCAEF(m),EE2,EE3,EE4,EE5
      endif
  110 FORMAT (1X,68A1,
     * 1X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,F8.5,1X,F6.3,1X,F6.3)
  120 FORMAT (70X,I2,2X,F6.3,2X,F9.6,2X,F6.3,1X,F9.6,1X,F8.5,
     *        1X,F6.3,1X,F6.3)
      zet=.false.
      endif
  130 continue
c
c     Se zauto=T, in seguito cipselx eseguira` un calcolo con 
c     zauto=T, quindi calcoliamo subito il contributo di tutti i
c     determinanti permutati di spin. Un solo determinante verra` scritto
c     e il suo coefficiente sara` uguale alla radice quadrata della somma
c     dei quadrati di tutti i determinanti con la stessa occupazione
c     orbitalica (cosi quando cipselx fa il quadrato ottiene la somma dei
c     quadrati). Attenzione al cipsel tradizionale ....
c
      if (.not.(zpun.and.nec.le.12)) return
      zstamp=.true.
cele      if (zauto) then
cele      call config(kk,dcaen,ecaen,dcaef,ecamp,zstamp,zmono,zdav,zbin)
cele      endif
c
      if(zthird)goto 1818
      if (zstamp) then
      if (zmppun) then
       write (7) nec,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
     *      ' ',(dabs(dcaef(i)),i=1,metat),(finv*ecamp(i),i=1,metat),
     *    ((dsqrt(shmp(m,i)),m=1,metat),i=1,nshift)
      else
       write (7) nec,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
     *      ' ',(dabs(dcaen(i)),i=1,metat),(finv*ecaen(i),i=1,metat),
     *   ((dsqrt(shen(m,i)),m=1,metat),i=1,nshift)
      endif
      endif
 1818 return
C
C     FONCTION D'ONDE VARIATIONNELLE
C
  200 if (zqua) then
      if (nec.eq.0) then
      write (6,210) kk,(wk(i),i=m,mm)
      else
         if (ltot.le.54) then
         WRITE (6,220) (ctrou(i:i),i=1,lt),arrow,(cpart(i:i),i=1,lp),
     *   (' ',i=1,54-ltot),kk,(wk(i),i=m,mm)
         else
c
c    Cele: provvisorio
c
         if (lt.lt.58) then
         WRITE (6,220) (ctrou(i:i),i=1,lt),arrow
	 else
         WRITE (6,220) (ctrou(i:i),i=1,56),arrow
	 endif
         if (lp.lt.58) then
         WRITE (6,220) (' ',i=1,2),(cpart(i:i),i=1,lp),
     *   (' ',i=1,56-lp),kk,(wk(i),i=m,mm)
	 else
         WRITE (6,220) (' ',i=1,2),(cpart(i:i),i=1,56),
     *   (' ',i=1,56-lp),kk,(wk(i),i=m,mm)
         endif
         endif
      endif
  210 FORMAT (1X,'Fondamental',47x,i6,1x,6f11.8)
  220 FORMAT (1X,58A1,i6,1x,6f11.8)
      endif
      if(zthird)return
      if (.not.zpun.or..not.zmono) return
       write (7) nec,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
     *      'S',(abs(c(ncf*(k-1)+kk)),k=1,metat),(0.d0,k=1,metat),
     *      ((0.d0,i=1,metat),j=1,nshift)
  240 continue
      return
  310 FORMAT (61a1,a1,99f15.10)
      END
C*********************************************************************
      SUBROUTINE DETER(zdet,TOCOM,INDIS,tgen,zrest,zdav,zdavid2
     $     ,dir,idimorb,iditrpa,ndims)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      character*2000 namel
      character*4 xt
      include 'parameter.cippi'
      PARAMETER (JD2=2*ID5,JD3=3*ID5,JD4=4*ID5,JD5=5*ID5,JD6=6*ID5)
      PARAMETER (JD7=7*ID5,JD8=8*ID5,JD9=9*ID5,JD10=10*ID5)
      PARAMETER (JD11=11*ID5,JD12=12*ID5)
      PARAMETER (JDTOT=156*ID5)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      dimension nt(13),np(13)
c      COMMON /VEC2/ NTM,NPM,NTD,NPD,NTT,NPT,NTQ,NPQ,NTP,NPP,NTH,NPH,
c     * nte,npe,nto,npo,ntn,npn,ntx,npx, 
c     * ntu,npu,ntz,npz, 
c     * KF,KM,KD,KT,KQ,KP,KH,ke,ko,kn,kx,
c     * ku,kz,
c     * NT(13),NP(13),MT(13),MP(13)
c      INTEGER*2 NTM(ID5),NPM(ID5),NTD(JD2),NPD(JD2),NTT(JD3),NPT(JD3),
c     * NTQ(JD4),NPQ(JD4),NTP(JD5),NPP(JD5),NTH(JD6),NPH(JD6),
c     $ nte(jd7),npe(jd7),nto(jd8),npo(jd8),ntn(jd9),npn(jd9),
c     $ ntx(jd10),npx(jd10),ntu(jd11),npu(jd11),ntz(jd12),npz(jd12)    
      common/DVD2/tol,nbx,itmax,ilevel,ngue,nblock
c      INTEGER*2 ITAB(JDTOT)
c      EQUIVALENCE (ITAB(1),NTM(1))
c      EQUIVALENCE (c(1),NTM(1))
      common/third/zthird,zmp3,zen3
      common/debug/zdebug
      common/thresho/thresh
c      dimension c(id7)
      DIMENSION KEX(13)
c      EQUIVALENCE (KEX(1),KF)
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      COMMON /ALE/SPINMAX,SPINMIN,ZSPINSEL,ZINCRE
      character*64 projf,dir
      common /hst/ host
      character*8 host
      INTEGER*2 ietats
      INTEGER*2 IPERM(10)
      DIMENSION INDTR(12),INDPA(12)
cele 15-11-04      DIMENSION IND(ID3),INV(2)
      DIMENSION INV(2)
      allocatable ind(:)
cele
      INTEGER*2 INDIS(IDIMORB,20),idum2
      DIMENSION ISYDEG(20,20)
      DIMENSION TOCOM(*)
      allocatable igel(:)
c      dimension igel(id1)
      common /deephole/ ndhb,zdh,ndhb1,ndhb2,zddh
      common /toul/ ztoul
C 
C     NAMELIST /ICINP/ TEST,ISZ,METAT,ZPRT,ZION,
C    * ZPERTU,ZBRD,ZDET,KGEN,TAU,ZPUN,
C    * ZHMAT,ZRAP,ZDIAB,NROT,ZQD,
C    * NROT,ZQD,FACTOR,DELTA,zhandy,zatom,zspin,zdialz,zsrot
C     NAMELIST /DIAGO/ SFO,SEN,SECR,NITER,IOP,MAXM,ZSS,ISEL,CSEL,ZVEC,
C    * ZWRT,ZFOL
C 
      allocate(ind(idimorb*2+1))
      allocate(igel(idimorb))
      allocate(zfree(idimorb))
      nvecg=id7*8+id14*8+id6*9*8+id3*8+id10*2*8+id6*8*2*8+9*4
c      if (jdtot*2.gt.nvecg) then
c      write (6,*) ' '
c      write (6,*)
c     * '   Attention: common VEC dans DETER/LECNEW/TRANSF plus grand'
c      write (6,*)
c     * '   que dans le reste du programme:',jdtot*2,nvecg,'   bytes'
c      write (6,*) ' '
c      endif
c      INDTR(1)=0
c      INDTR(2)=ID5*2
c      INDTR(3)=ID5*6
c      INDTR(4)=ID5*12
c      INDTR(5)=ID5*20
c      INDTR(6)=ID5*30
c      indtr(7)=id5*42
c      indtr(8)=id5*56
c      indtr(9)=id5*72
c      indtr(10)=id5*90
c      indtr(11)=id5*110
c      indtr(12)=id5*132
c      INDPA(1)=ID5
c      INDPA(2)=ID5*4
c      INDPA(3)=ID5*9
c      INDPA(4)=ID5*16
c      INDPA(5)=ID5*25
c      INDPA(6)=ID5*36
c      indpa(7)=id5*49
c      indpa(8)=id5*64
c      indpa(9)=id5*81
c      indpa(10)=id5*100
c      indpa(11)=id5*121
c      indpa(12)=id5*144
      rewind 25
      allocate(itsym(idimorb*2+1))
      allocate(zgel(idimorb*2+1))
      allocate(iorb(idimorb*2+1))
      allocate(ispin(idimorb*2+1))
      READ (25) NSYM,NORB,NOC,NDEGEN,(ITSYM(I),I=1,NORB),
     *          nao,ngelo,ngelv
      rewind 3
      call cerca(3,'ICINP',itrovo)
      if (itrovo.eq.0) then
      write (6,*) '   *** NAMELIST &ICINP NOT FOUND ***'
      stop 12
      endif
C     read (3,ICINP)
      call lecnam('&ICINP',namel,3)
  50  if (ndegen.gt.0) then
      write (6,*) 
      write (6,*) '   Attenzione!  NDEGEN =',ndegen
      write (6,*)
     * '   Le parti del programma riguardanti il caso NDEGEN>0'
      write (6,*) '   sono state conservate per compatibilita`,'
      write (6,*) '   ma la loro correttezza non e` garantita.'
      write (6,*) '   Ricompilare senza lo STOP che segue e testare.'
      STOP 12
      endif
c
      INV(1)=NORB
      INV(2)=-NORB
      IF (NORB.EQ.1) GOTO 300
      READ (25) ((ITS(I,J),I=1,NSYM),J=1,NSYM)
      DO 10 I=1,NSYM
      DO 9 J=1,NSYM
    9 IF (ITS(1,J).EQ.I) IPERM(I)=J
   10 IF (NDEGEN.GT.0) IPERM(I)=I
      DO 11 I=1,NORB
   11 ITSYM(I+NORB)=ITSYM(I)
      NORB2=NORB+NORB
      MNOR1=NORB2+1
      DO 4 I=1,NSYM
      K=0
      DO 4 J=1,MNOR1
cele      ZOCS(J)=.FALSE.
      IF (ITSYM(J).NE.I) GOTO 4
      K=K+1
      IND(J)=K
      INDIS(K,I)=J
    4 CONTINUE
      DO 5 I=1,NORB
    5 ZFREE(I)=.TRUE.
      NTRSY=1
      DO 15 I=1,NDEGEN
      READ (25) J,((ISYDEG(K,L),K=1,NSYM),L=1,J)
      DO 14 L=1,J
      DO 14 K=1,NORB2
      INSKL=0
      NSKL=ISYDEG(ITSYM(K),L)
      ZSIN=.FALSE.
      IF (NSKL.GE.0) GOTO 13
      NSKL=-NSKL
      ZSIN=.TRUE.
   13 CONTINUE
      IF (NSKL.NE.0) INSKL=INDIS(IND(K),NSKL)
      IF (ZSIN) INSKL=-INSKL
   14 continue
      NTRSY=NTRSY+J
   15 CONTINUE
      ZDET=.FALSE.
      ZION=.FALSE.
      ZPERTU=.TRUE.
      ZENB=.FALSE.
      ZMPPUN=.FALSE.
      ZPUN=.FALSE.
      ZQD=.TRUE.
      ISZ=0
      TEST=0.05D0
      TAU=1.D20
      TDEGEN=1.d-7
      DELTA=0.10D0
      ZPRT=.FALSE.
      ZBRD=.FALSE.
      ZRAP=.FALSE.
      ZDIAB=.FALSE.
      zhmat=.true.
      zhandy=.true.
      zatom=.false.
      zspin=.false.
      zdialz=.false.
      zsrot=.false.
      zauto=.true.
      zfac=.false.
      zdav=.false.
      ztoul=.false.
      NROT=0
      FACTOR=1.D6
      METAT=1
c      KGEN(1)=0
      do i=1,norb
      tocom(i)=0.d0
      enddo
      do i=1,idimorb
      igel(i)=0
      zgel(i)=.false.
      zgel(i+norb)=.false.
      enddo
      do i=1,12
      shift(i)=0.d0
      enddo
      tocom(1)=-1.d10
      tgen=.05D0
      host='        '
      dir='/home/renzo/diago/'
      do i=19,64
         dir(i:i)=' '
      enddo
      zdh=.false.
      zddh=.false.
      ndhb=0
      ndhb1=0
      ndhb2=0
      zcas=.false.
      zthird=.false.
      zmp3=.false.
      zen3=.false.
      thresh=0.0d0
      zdebug=.false.
      call namell(namel,zprt,'ZPRT',1)
      call namell(namel,zcas,'ZCAS',1)
      call namell(namel,zthird,'ZTHIRD',1)
      call namell(namel,zmp3,'ZMP3',1)
      call namell(namel,zen3,'ZEN3',1)
      call namell(namel,zdebug,'ZDEBUG',1)
      if(zthird)then
         if(zmp3.and.zen3)then
            zmp3=.false.
            zen3=.true.
         elseif((.not.zmp3).and.(.not.zen3))then
            zen3=.true.
            zmp3=.false.
         endif
      endif
      zdavid2=.false.
      call namell(namel,zdavid2,'ZDAVID2',1)
      call namell(namel,zion,'ZION',1)
      call namell(namel,zbrd,'ZBRD',1)
      call namell(namel,zpertu,'ZPERTU',1)
      call namell(namel,zenb,'ZENB',1)
      call namell(namel,zmppun,'ZMPPUN',1)
      call namell(namel,zdet,'ZDET',1)
      call namell(namel,zpun,'ZPUN',1)
      call namell(namel,zfac,'ZFAC',1)
      write (6,*) 'zfac',zfac
      call namell(namel,zhmat,'ZHMAT',1)
      call namell(namel,zatom,'ZATOM',1)
      call namell(namel,zspin,'ZSPIN',1)
c      call namell(namel,zauto,'ZAUTO',1)
      call namell(namel,zdialz,'ZDIALZ',1)
      call namell(namel,zsrot,'ZSROT',1)
      call namell(namel,zhandy,'ZHANDY',1)
      call namell(namel,zdh,'ZDH',1)
      call namell(namel,zddh,'ZDDH',1)
      call nameli(namel,ndhb,'NDHB',4)
      call nameli(namel,ndhb1,'NDHB1',4)
      call nameli(namel,ndhb2,'NDHB2',4)
      call namell(namel,zdav,'ZDAV',1)
      call namell(namel,ztoul,'ZTOUL',1)
      call namela(namel,DIR,'DIR')
      if (.not.zhandy) factor=1.d9
      call namell(namel,zrap,'ZRAP',1)
      call namell(namel,zdiab,'ZDIAB',1)
      call namell(namel,zqd,'ZQD',1)
      call nameli(namel,isz,'ISZ',4)
      call nameli(namel,metat,'METAT',4)
      call nameli(namel,nrot,'NROT',4)
      call nameli(namel,igel,'IGEL',4)
      call nameli(namel,NSHIFT,'NSHIFT',4)
      call nameli(namel,kgen,'KGEN',4)
      call nameld(namel,test,'TEST',8)
      call nameld(namel,tau,'TAU',8)
      call nameld(namel,tdegen,'TDEGEN',8)
      call nameld(namel,tgen,'TGEN',8)
      call nameld(namel,factor,'FACTOR',8)
      call nameld(namel,delta,'DELTA',8)
      call nameld(namel,SHIFT,'SHIFT',8)
      call nameld(namel,TOCOM,'TOCOM',8)
      call nameld(namel,THRESH,'THRESH',8)
      call namela(namel,HOST,'HOST')
      if(zdh)then
         print '(a,L4)',' BUCA PROFONDA                ',zdh
         print '(a,I4)',' ORBITALE PROFONDO            ',ndhb
      endif
      if(zddh)then
         print '(a,L4)',' BUCA PROFONDA DOPPIA         ',zddh
         print '(a,I4)',' ORBITALE PROFONDO 1          ',ndhb1
         print '(a,I4)',' ORBITALE PROFONDO 2          ',ndhb2
      endif
      if (zdh.and.zddh) then
      write (6,*) 'Opzioni incompatibili: ZDDH=T, ZDH=T'
      stop
      endif

c
c     Nel caso che NSHIFT sia diverso da zero e non siano stati definiti degli 
c     shift in energia, vengono definiti per default 
c     (servono per la selezione mirata)
c
      if (nshift.ne.0.and.shift(1).eq.0.d0) then
      shift(1)=.05d0
      do i=2,nshift
      shift(i)=shift(i-1)+.05d0
      enddo
      endif
c
      IF (ZDIAB) ZBRD=.TRUE.
      IF (METAT.EQ.1) ZBRD=.FALSE.
      IF (NROT.LT.METAT) NROT=METAT
      IF (.NOT.ZPUN) TAU=TEST
      IF (TAU.GT.TEST) TAU=TEST
      FINV=1.D0/FACTOR
      SFO(1)=1.D-5
      do i=2,5
      SFO(i)=-1.d0
      enddo
      SEN=0.D0
      SECR=5.D-2
      spinmin=-1.d0
      spinmax=20.d0
      NITER=50
      IOP=0
      MAXM=12
      maxgue=200
      ZSS=.FALSE.
      ZVEC=.FALSE.
      ZWRT=.FALSE.
      ZFOL=.FALSE.
      ZSPINSEL=.FALSE.
      ZINCRE=.FALSE.
      zrest=.false.
      projf=' '
      DO 90 I=1,ID12
      CSEL(I)=0.D0
   90 ISEL(I)=0
      do 91 i=1,id6
   91 ietats(i)=i
      rewind 3
      call cerca(3,'DIAGO',itrovo)
      if (itrovo.ne.0) then
C     read (3,DIAGO)
      call lecnam('&DIAGO',namel,3)
      call namela(namel,PROJF,'PROJF')
      call namell(namel,zwrt,'ZWRT',1)
      call namell(namel,zvec,'ZVEC',1)
      call namell(namel,zfol,'ZFOL',1)
      call namell(namel,zrest,'ZREST',1)
      call namell(namel,zss,'ZSS',1)
      call nameli(namel,niter,'NITER',4)
      call nameli(namel,iop,'IOP',4)
      call nameli(namel,maxm,'MAXM',4)
      call nameli(namel,maxgue,'MAXGUE',4)
      call nameli(namel,isel,'ISEL',4)
      call nameli(namel,ietats,'IETATS',2)
      call nameld(namel,sfo,'SFO',8)
      call nameld(namel,sen,'SEN',8)
      call nameld(namel,secr,'SECR',8)
      call nameld(namel,csel,'CSEL',8)
      call namell(namel,ZSPINSEL,'ZSPINSEL',1)
      call nameld(namel,spinmin,'SPINMIN',8)
      call nameld(namel,spinmax,'SPINMAX',8)
c----parte riguardante david2      
      tol=1.0D-6
      nbx=50
      itmax=1000
      ilevel=2
      ngue=400
      nblock=0
      call nameld(namel,tol,'TOL',8)
      call nameli(namel,nbx,'NBX',4)
      call nameli(namel,itmax,'ITMAX',4)
      call nameli(namel,ilevel,'ILEVEL',4)
      call nameli(namel,ngue,'NGUE',4)
      call nameli(namel,nblock,'NBLOCK',4)
      endif
      if (maxgue.gt.id12) then
      write (6,*) 'Attenzione: maxgue maggiore di id12',maxgue,id12
      write (6,*) 'Per aumentare id12 ricordare che al massimo ',
     * 'puo` valere (id7/2)**.5'
      write (6,*) 'MAXGUE e` posto uguale a id12'
      maxgue=id12
      endif
      if (zvec) zss=.false.
      IF (MAXM.GT.ID14) MAXM=ID14
      WRITE (6,3000)
 3000 FORMAT (//)
      WRITE (6,3010)
      WRITE (6,3020)
 3020 FORMAT (' INFORMATIONS GENERALES ET COMMENTAIRES')
      WRITE (6,3010)
 3010 FORMAT (1X,38('*'))
      WRITE (6,3000)
      ISYM=0
      INITB=1
      INITA=NORB+1
      NOCB=NOC
      NOCA=NOC+NORB
      IF (.NOT.ZION) GOTO 109
      ITSYM(MNOR1)=1
      IORB(MNOR1)=NORB+1
      ISPIN(MNOR1)=1
  109 CONTINUE
      DO 150 I=1,NORB
      IORB(I)=I
      IORB(I+NORB)=I
      ISPIN(I)=0
  150 ISPIN(I+NORB)=1
      WRITE (6,3500) NORB
 3500 FORMAT (' NOMBRE D''ORBITALES(NORB):',I3)
      WRITE (6,3520) NOCB
 3520 FORMAT (' NOMBRE D''ORBITALES OCCUPEES :',I3)
      WRITE (6,3580) NSYM
 3580 FORMAT (' NOMBRE DE SYMETRIES POUR LES ORBITALES(NSYM):',I2)
      WRITE (6,3600) ISZ
 3600 FORMAT (' VALEUR DE SZ:',I2)
      XT=' NON'
      IF (ZION) XT=' OUI'
      WRITE (6,3700) XT
 3700 FORMAT (' EXCITATION VERS L''INFINI:',A4)
      IF (METAT.GT.ID6) GOTO 280
      WRITE (6,3590) METAT
 3590 FORMAT (' NOMBRE D''ETATS TRAITES EN PERTURBATION(METAT):',I3)
      WRITE (6,3630) TEST,TAU
 3630 FORMAT (' TEST POUR L''IMPRESSION DES DETERMINANTS (TEST):',F9.6/
     * ' TEST POUR LA PERFORATION DES DETERMINANTS (TAU):',F9.6)
      XT=' NON'
      IF (ZBRD) XT=' OUI'
      WRITE (6,3710) XT
 3710 FORMAT (' HAMILTONIEN EFFECTIF QDPT(ZBRD):',A4)
      WRITE (6,3720)
 3720 FORMAT (' RESULTATS VARIATIONNELS ET ENERGIES GARDES SUR',
     *' LA FILE 04')
      IF (ZRAP) WRITE (6,3730)
 3730 FORMAT (' ON NEGLIGE L''INTERACTION ENTRE DETERMINANTS ',
     * 'NON GENERATEURS ET DETERMINANTS PERTURBATEURS')
      WRITE (6,3000)
      WRITE (6,*)
      if (zenb) WRITE(6,*)' ATTENZIONE !'
      if (zenb) WRITE(6,*)' SI TRATTA DI UN CALCOLO EN BARICENTRICO'
      if (zdav) write(6,*)' Attenzione !'
      if (zdav) write(6,*)' Si tratta di un calcolo EN con energie ',
     * 'mediate sugli stati di spin (ZDAV=T)'
C 
C     LECTURE DES DETERMINANTS DE S ET GENERATION AUTOMATIQUE DES SPINS
C 
c
c   Lettura da file 5
c
      IF (ZDET) GOTO 250
c      DO 180 K=1,ID5
c      NE(K)=0
c      ND(K)=0
c      NDEG(K)=0
cele      ZSAUT(K)=.FALSE.
  180 continue
      CALL LECNEWER(idimorb,ndims)
      IF (NCF.EQ.0) GOTO 290
      IF (KGEN(1).NE.0) GOTO 210
      NCFG=NCF
c      IF (NCFG.GT.ID5) GOTO 320
      DO 200 KG=1,NCFG
  200 KGEN(KG)=KG
      GOTO 230
  210 NCFG=0
c      DO 220 KG=1,ID5
c  220 IF (KGEN(KG).NE.0) NCFG=NCFG+1
      do while (kgen(ncfg+1).ne.0)
         NCFG=NCFG+1
      enddo
 220  continue
  230 WRITE (6,240) NCF,NCFG
  240 FORMAT (5X,'DIMENSION DU SOUS ESPACE VARIATIONNEL, NCF =',I8/
     * 5X,'NOMBRE DE DETERMINANTS GENERATEURS, NCFG =',I8/)
      IF (NCF.eq.NCFG) then
         ZRAP=.TRUE.
      elseif (zprt) then
         write (6,244)
         WRITE (6,245) (KGEN(KG),KG=1,NCFG)
      endif
  244 FORMAT (/5X,'DETERMINANTS GENERATEURS :')
  245 FORMAT (2X,20I6)
      WRITE (4) NORB
c--renzo La scrittura di ndeg e' importante!      
      WRITE (4) NCFG,(KGEN(KG),KG=1,NCFG),ISYM,(NDEG(K),K=1,NCF)
      goto 260
c
c   Lettura da file 4
c
  250 READ (4) NORB,NOCB,NOCA,IDUM,NCF,II
      rewind 4
c      deallocate(ne)
c      deallocate(nd)
c      deallocate(trou)
c      deallocate(part)
c      deallocate(zsaut)
      allocate(ne(ncf+1000))
      allocate(nd(ncf+1000))
      allocate(zsaut(ncf+1000))
      allocate(trou(ii+20000))
      allocate(part(ii+20000))
      iditrpa=ii
      READ (4) NORB,NOCB,NOCA,IDUM,NCF,II,IDUM,(NE(J),ND(J),J=1,NCF),
     * (TROU(J),PART(J),J=1,II)
c      if (II.gt.id8) then
c      write (6,*) 'Attenzione, II maggiore di id8',ii,id8
c      stop 12
c      endif
c
c   Compatibilita` con versioni in cui KGEN e` INTEGER*2
c
      if (zpertu) then
        READ (4) NCFG
        if(associated(kgen))deallocate(kgen)
        allocate(kgen(ncfg+1000))
        rewind 4
        READ (4) NORB,NOCB,NOCA,IDUM,NCF,II,IDUM,(NE(J),ND(J),J=1,NCF),
     * (TROU(J),PART(J),J=1,II)
        READ (4) NCFG,(KGEN(KG),KG=1,NCFG),ISYM
      else
        READ (4) NCFG,(K,KG=1,NCFG),ISYM
      endif
      WRITE (6,240) NCF,NCFG
      IF (NCF.eq.NCFG) then
      ZRAP=.TRUE.
      elseif (zprt) then
      write (6,244)
      WRITE (6,245) (KGEN(KG),KG=1,NCFG)
      endif
c
c     Cele: costruisco ora zsaut, sia per zdet=T che =F.
c
c--renzo bisogna costruire anche ndeg!!
  260 if (zprt) then
      write (6,*) ' determinants (i, nd(i), ne(i), trou, part)'
      call flush(6)
      do i=1,ncf
      write (6,'(4i8,23i4)') i,nd(i),ne(i),
     * (trou(k),k=nd(i)+1,nd(i)+ne(i)),
     * (part(k),k=nd(i)+1,nd(i)+ne(i))
      enddo
      endif
      do i=1,ncf
      zsaut(i)=.false.
      enddo
      if (isz.ne.0.or.zion) goto 270
      znext=.false.
      do 269 i=1,ncf
      if (znext) then
      znext=.false.
      goto 269
      endif
      nec=ne(i)
      do l=1,nec
      nt(l)=trou(nd(i)+l)
      np(l)=part(nd(i)+l)
      enddo
c      if (int(nec/2)*2.ne.nec) goto 265
      if ((nec/2)*2.ne.nec) goto 265
      do l=1,nec/2
      if (nt(l)+norb.ne.nt(l+nec/2)) goto 265
      if (np(l)+norb.ne.np(l+nec/2)) goto 265
      enddo
      goto 269
  265 zsaut(i+1)=.true.
      znext=.true.
 269  continue
c
c     Controllo che nessun orbitale gelato (igel(i).ne.0) sia
c     attivo nei det dello spazio S
c
  270 continue
      ngel=0
      do i=1,norb
      if (igel(i).ne.0) then
         ngel=ngel+1
         zgel(igel(i))=.true.
         zgel(igel(i)+norb)=.true.
      endif
      enddo
      if (ngel.gt.0) then
      do i=1,nd(ncf)+ne(ncf)
      if (zgel(trou(i))) then
         write (6,*) 'Attenzione l''orbitale ',trou(i),
     *    ' non puo` essere gelato (e` attivo in S)'
	 write (6,*) ' buca n. ',i
         stop 12
      endif
      if (zgel(part(i))) then
         write (6,*) 'Attenzione l''orbitale ',part(i),
     *    ' non puo` essere gelato (e` attivo in S)'
	 write (6,*) ' particella n. ',i
         stop 12
      endif
      enddo
      write (6,271) (igel(i),i=1,ngel)
      endif
  271 format(/'     ORBITALI GELATI:'/10(5x,18i4/))
      RETURN
  280 WRITE (6,281)
  281 FORMAT (/' TROP D''ETATS DEMANDES'//)
      STOP 12
  290 WRITE (6,291)
  291 FORMAT (' ERREUR DANS LES DONNEES: AUCUN DETERMINANT DANS LE '
     *,'SOUS ESPACE VARIATIONNEL'/)
      STOP 12
  300 WRITE (6,301)
  301 FORMAT (' ERREUR DANS LES INFORMATIONS GENERALES')
      STOP 12
  320 WRITE (6,321) NCFG
  321 FORMAT (/' TROP DE DETERMINANTS GENERATEURS, NCFG =',I5//)
      STOP 12
      END
C***********************************************************************
      SUBROUTINE DAVID(ncf,N,NVEC,B,AB,zrest,secr,iditrpa)
C-----------------------------------------------------------------------
C 
C     SEARCH OF EIGENVALUES AND EIGENVECTORS FOR LARGE SYMMETRIC
C     MATRICES
C     E. R. DAVIDSON, J. COMP. PHYS 17, 87 (1975)
C 
C-----------------------------------------------------------------------
C 
C     DESCRIPTION OF THE ARGUMENTS
C     ALAM   = EIGENVALUES (MAXIMUM 12)
C     N      = DIMENSION OF THE VECTOR SPACE (MAXIMUM 2000)
C     NVEC   = NUMBER OF DESIRED EIGENVECTORS (AND EIGENVALUES)
C     NITER  = MAXIMUM NUMBER OF ITERATIONS
C     THR1   = THRESHOLD ON THE WAVEFUNCTION 
C     THR2   = THRESHOLD ON THE ENERGY (THE PROGRAM STOPS WHEN
C              EITHER THRESHOLD IS SATISFIED. USUALLY THR2=0)
C     IOP    = 1  ACCUMULATION OF ROUNDOFF ERRORS IS AVOIDED IN
C                 THE STEP OF MATRIX REDUCTION
C            = 0  NORMAL CALCULATION (FASTER BUT LESS ACCURATE)
C     MAXM   = MAXIMUM DIMENSION OF INTERMEDIATE BAB MATRIX
C              (MAXM.LE.ID14)
C 
C 
C     UNIT 20 ...ibuf,jbuf,buffer
C                 ibuf=-1 = end of data
C                 
C 
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      dimension b(n,*),ab(n,*)
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/C,ALAM(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      pointer ort,gue,ibuf,jbuf,a,bnew,abnew,bab,h,alpha,p,zfound
      dimension ort(:,:),gue(:,:),ibuf(:),jbuf(:),a(:),bnew(:),
     $     abnew(:),bab(:,:),h(:,:),alpha(:,:),p(:),zfound(:)
      COMMON /INT2/ NI4,NI3,NI2,NR4,
     * ORT,GUE,IBUF,JBUF,
     * A,BNEW,ABNEW,BAB,
     * H,ALPHA,P,ZFOUND
      allocatable cdet(:)
c      equivalence (cdet(1),ort(1,1))
      INTEGER*2 idum2
      COMMON /DIA/ THR1(5),THR2,THR3,CSEL(ID12),ISEL(ID12),NITER,
     * IOP,MAXM,maxgue,IETATS(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      common /fil/ file12,file13,file26,zwconf
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      character*64 projf,file12,file13,file26
      INTEGER*2 ietats
      allocatable diag(:),q(:),work(:),uvec(:),olsen(:)
c      dimension diag(id4),q(id4)
c      equivalence (diag(1),c(1)),(q(1),c(id4+1))
    2 FORMAT (//5X,'SEARCH OF EIGENVECTORS BY DAVIDSON''S METHOD'/
     * 5X,'MAXIMUM DIMENSION OF THE DIAGONALIZED SUBSPACE ...',I5/
     * 5X,'OPTION FOR THE REDUCTION OF THE SUBSPACE .........',I5/
     * 5X,'DIMENSION OF THE TRIAL VECTORS ...................',I7/
     * 5X,'REQUIRED CONVERGENCE :',D10.3,' FOR THE EIGENVECTOR'/
     * 5X,'OR',D10.3,' FOR THE EIGENVALUE')
    4 FORMAT (5X,'NUMBER OF INTRUDER STATES OF DIFFERENT SYMMETRY ',
     * 'IN THE TRIAL SUBSPACE ',I5)
    5 FORMAT (//10X,'EIGENVALUE NUMBER',I3/
     * 3X,'ITERATION  ROOT     M    EIGENVALUE     CONV. WF.    ',
     * 'CONV. EN.')
    6 FORMAT (3I8,F14.8,2D13.2)
    7 FORMAT (///5x,'OVERALL CONVERGENCE:'/
     * 5x,'SQUARE NORMS OF  (A-E) B(K)  =')
    8 FORMAT (///5X,'THE DAVIDSON PROCEDURE HAS NOT CONVERGED AFTER',
     * I4,' ITERATIONS'///)
    9 FORMAT (60X,'VECTOR B(',I2,') HAS BEEN ELIMINATED')
   10 FORMAT (//5X,'BAD OVERLAP WITH GUESS VECTORS *** STOP'//)
   11 FORMAT (5x,4d16.4)
   12 FORMAT (5X,'<Q|Q> BEFORE AND AFTER ORTHOG.:',2d10.3,
     * ' :  REPEAT ORTHOGONALISATION')
   13 FORMAT (//5X,'REPEAT PROCEDURE WITH NEW WF THRESHOLD =',d10.3)
c      kd42=id4*2
c      if (kd42.gt.id7) then
c      write (6,*) '   WRONG DIMENSIONS IN DAVID 2*ID4>ID7:',id4,id7
c      stop 12
c      endif
      allocate(diag(ncf+1000))
      allocate(q(ncf+1000))
      nintg=id2*3*8+id2*4+599*4+4000*2+id2*3*2+16
      nintl=id12*id6*2*8+5000*2+5000*8+id14*3*8+id14*id14*3*8+id6
c      if (nintl.gt.nintg) then
c      write (6,*) ' '
c      write (6,*)
c     * '   Attention: common INT dans le ssp. DAVID plus grand'
c      write (6,*)
c     * '   que dans le reste du programme:',nintl,nintg,'   bytes'
c      write (6,*) ' '
c      endif
      allocate(ort(id12,id6))
      allocate(gue(id12,id6))
      allocate(ibuf(100000))
      allocate(jbuf(100000))
      allocate(a(100000))
      allocate (bnew(id14))
      allocate (abnew(id14))
      allocate(bab(id14,id14))
      allocate(h(id14,id14))
      allocate(alpha(id14,id14))
      allocate(p(id14))
      allocate(work(3*id14))
      allocate(zfound(id6))
      allocate(uvec(n))
      allocate(olsen(n))
      icount=1
      NNV=N*NVEC
      N2=N+N
      M=NVEC
C
C     INITIAL VECTORS FROM RESTART                               
C
      if (zrest) then
      write (6,*) 'Restart della procedura di DAVIDSON'
      write (6,*) 'Vettore di guess da un calcolo precedente'
      open(13,file=file13,form='UNFORMATTED',status='UNKNOWN')
      rewind 13
      maxdim=n*nvec
      read (13) ((b(i,j),i=1,n),j=1,nvec)
      close (13)
      write (6,*) 'Guess vectors from file ',file13
      call flush(6)
      NORT=0
      goto 15
      endif
C
C     INITIAL VECTORS FROM DIAGONALISATION IN A SELECTED SUBSPACE
C
      IF (.NOT.ZVEC) then
         CALL GUESS(N,NVEC,NGUE,NORT,ab)
         WRITE (6,2) MAXM,IOP,NGUE,THR1(ICOUNT),THR2
         IF (ZSS) WRITE (6,4) NORT
         DO  L=1,M
            DO  I=1,N
               B(I,L)=0.D0
            enddo
            DO  I=1,NGUE
               B(ISEL(I),L)=GUE(I,L)
            enddo
         enddo
      else
C
C     INITIAL VECTORS FROM PREVIOUS RUN WITH THE SAME SUBSPACE
C
         if (projf.eq.' ') then
c            REWIND 4
c            READ (4) IDUM,IDUM,IDUM,IDUM,IDUM,I1,I2,(idum2,IDUM,J=1,N),
c     *           (IDUM,J=1,I1),(B(J),J=1,I2)
c            WRITE (6,2) MAXM,IOP,N,THR1(ICOUNT),THR2
         write (6,*) 'Zvec=T non allowed for CSF calculation'
         stop
         else
C
C     INITIAL VECTORS FROM PREVIOUS RUN WITH A SMALLER SUBSPACE
C
            call projv(n,nvec,ncfpr,b,iditrpa,ncf)
            WRITE (6,2) MAXM,IOP,ncfpr,THR1(ICOUNT),THR2
         endif
      endif
      write (6,*) 'Guess done'
      call flush(6)
C
C     START DAVIDSON PROCEDURE (UP TO 5 TIMES WITH DECREASING THR1)
C
   15 K=1
      KB=1
      KR=1
      WRITE (6,5) K
      ITER=1
   20 EN=1.D9
      ZCONV=.FALSE.
C 
C     CONSTRUCTION OF THE BAB MATRIX IN THE SUBSPACE
C     OF THE TRIAL VECTORS  B(1) ... B(M)
C 
      REWIND 20
      READ (20) (DIAG(IIII),IIII=1,N)
      IL=0
      DO  L=1,M
         do i=1,n
            ab(i,l)=diag(i)*b(i,l)
         enddo
      enddo
      MAX=100000
   35 read (20) (ibuf(i),i=1,size(ibuf)),
     $          (jbuf(i),i=1,size(jbuf)),
     $          (a(i),i=1,size(a))
      if(ibuf(100000).eq.-1)max=jbuf(100000)
      DO  IC=1,MAX
         ih=ibuf(ic)
         jh=jbuf(ic)
c         call daxpy(m,a(ic),b(jh,1),n,ab(ih,1),n)
c         call daxpy(m,a(ic),b(ih,1),n,ab(jh,1),n)
         DO  L=1,M
            AB(ih,l)=AB(ih,l)+A(IC)*B(jh,l)
            AB(jh,l)=AB(jh,l)+A(IC)*B(ih,l)
         enddo
      enddo
      if(ibuf(100000).ne.-1)goto 35
c      call dgemm('T','N',m,m,n,1.d0,b,n,ab,n,0.d0,bab,id14)
      do L=1,M
         do LL=1,L
            sum=ddot(n,b(1,l),1,ab(1,ll),1)
            bab(l,ll)=sum
            bab(ll,l)=sum
         enddo
      enddo
C 
C     DIAGONALIZATION OF BAB ... EIGENVALUES ALAM, EIGENVECTORS ALPHA
C 
      IF (ZWRT) WRITE (6,'(''  BAB MATRIX'')')
 70   DO  L=1,M
         IF (ZWRT) WRITE (6,'(1X,10F13.8)') (BAB(L,LL),LL=1,L)
         call dcopy(m,bab(1,l),1,h(1,l),1)
         call dcopy(m,bab(1,l),1,alpha(1,l),1)
c         DO  LL=1,M
c            H(L,LL)=BAB(L,LL)
c            alpha(l,ll)=H(l,ll)
c         enddo
      enddo
      IF (ITER.GT.1) EN=ALAM(KR)
c      write (6,*) 'dsyev     '
c      call flush(6)
      CALL dsyev('V','U',m,alpha,id14,alam,work,3*id14,info)
      IF (.NOT.ZFOL) GOTO 74
      AMAX=0.D0
      DO 73 L=1,M
      DAL=DABS(ALPHA(KB,L))
      IF (DAL.LT.AMAX) GOTO 73
      AMAX=DAL
      KR=L
   73 CONTINUE
   74 DENE=DABS(ALAM(KR)-EN)
C 
C     SETS UP THE CORRECTION VECTOR Q = B(M+1)
C 
   75 REWIND 20
      READ(20)(DIAG(IIII),IIII=1,N)
      EN=ALAM(KR)
c      qq=0.d0
      DO  I=1,N
         QI=0.D0
         uvec(i)=0.d0
         do l=1,m
            qi=qi+alpha(l,kr)*(ab(i,l)-en*b(i,l))
            uvec(i)=uvec(i)+alpha(l,kr)*b(i,l)
         enddo
         den=en-diag(i)
         if (dabs(den).gt.1.d-25) then
            q(i)=qi/den
            olsen(i)=uvec(i)/den
         else
            q(i)=0.d0
            olsen(i)=0.d0
         endif
c         qq=q(i)*q(i)
      enddo
      uveols=ddot(n,uvec,1,olsen,1)
      if(uveols.gt.1.d-15) then
      cappa=ddot(n,uvec,1,q,1)/ddot(n,uvec,1,olsen,1)
      else
      cappa=0.d0
      endif

      call daxpy(n,-cappa,olsen,1,q,1)
c      call dscal(n,-1.d0,q,1)
c      qq=ddot(n,q,1,q,1)
C 
C     ... ORTHOGONALIZATION WRT  B(1) ... B(M)
C 
 105  continue
      DO  L=1,M
         sum=ddot(n,b(1,l),1,q,1)
         call daxpy(n,-sum,b(1,l),1,q,1)
         IF (ZWRT) WRITE (6,'('' <Q!B'',I3,''> ='',F15.8)') L,SUM
      enddo
C 
C     ... ORTHOGONALIZATION WRT VECTORS OF DIFFERENT SYMMETRY   (ORT)
C 
      IF (NORT.gt.0) then
      DO 136 IORT=1,NORT
      SUM=0.D0
      DO 132 IG=1,NGUE
      I=ISEL(IG)
  132 SUM=SUM+Q(I)*ORT(IG,IORT)
      IF (ZWRT) WRITE (6,'('' <Q!ORT'',I3,''> ='',F15.8)') IORT,SUM
      DO 134 IG=1,NGUE
      I=ISEL(IG)
  134 Q(I)=Q(I)-SUM*ORT(IG,IORT)
  136 CONTINUE
      endif
C 
C     ... NORMALIZATION AND CONVERGENCE TEST
C     IF M OVERFLOWS, THE B SUBSPACE IS THOROUGHLY REDEFINED
C 
      SUM=ddot(n,q,1,q,1)
      cnorm=DSQRT(SUM/N)
      WRITE (6,6) ITER,KR,M,ALAM(KR),cnorm,DENE
      call flush(6)   !tolgo flush
c
      IF (M.GT.NVEC.AND.DENE.LT.THR2) ZCONV=.TRUE.
      IF (cnorm.LT.THR1(ICOUNT)) ZCONV=.TRUE.
      IF (ZCONV) GOTO 200
      IF (ITER.GE.NITER) GOTO 400
      ITER=ITER+1
      IF (M.GE.MAXM) GOTO 200
c      if (sum.lt.qq*1.d-2) then
c      write (6,12) qq,sum
c      goto 105
c      endif
      MX=M*N
      M=M+1
      fac=1.D0/DSQRT(SUM)
      DO  I=1,N
         B(i,m)=fac*Q(I)
      enddo
C 
C     ADDS ONE ROW AND ONE COLUMN TO THE BAB MATRIX
C 
      REWIND 20
      READ(20)(DIAG(IIII),IIII=1,N)
      DO  I=1,N
         AB(i,m)=DIAG(I)*B(i,m)
      enddo
      MAX=100000
 165  read (20)(ibuf(i),i=1,size(ibuf)),
     $          (jbuf(i),i=1,size(jbuf)),
     $          (a(i),i=1,size(a))
      if(ibuf(100000).eq.-1)max=jbuf(100000)
      DO  IC=1,MAX
         ih=ibuf(ic)
         jh=jbuf(ic)
         AB(ih,m)=AB(ih,m)+A(IC)*B(jh,m)
         AB(jh,m)=AB(jh,m)+A(IC)*B(ih,m)
      enddo
      if(ibuf(100000).ne.-1)goto 165
      DO  L=1,M
         SUM=ddot(n,b(1,l),1,ab(1,m),1)
         BAB(L,M)=SUM
         BAB(M,L)=SUM
      enddo
C 
C     BACK TO THE DIAGONALIZATION OF BAB
C 
      GOTO 70
C 
C     THE EIGENVECTORS OF BAB ARE BACK-TRANSFORMED
C     IN THE ORIGINAL BASIS
C 
  200 ZEND=ZCONV.AND.(K.EQ.NVEC)
      NTRAN=M
      IF (ZEND) NTRAN=NVEC
      do 205 l=1,ntran
  205 p(l)=0.d0
      DO  I=1,N
         LI=I-N
         DO  L=1,NTRAN
            BNEW(L)=ddot(m,b(i,1),n,alpha(1,l),1)
            ABNEW(L)=ddot(m,ab(i,1),n,alpha(1,l),1)
         enddo
         DO  L=1,NTRAN
            abab=abnew(l)-alam(l)*bnew(l)
            p(l)=p(l)+abab*abab
            B(I,l)=BNEW(L)
            AB(i,l)=ABNEW(L)
            IF (ZEND) C((l-1)*n+i)=BNEW(L)
         enddo
      enddo
      KB=KR
c
c     Il vettore viene storato sul file13 per un eventuale restart
c
      maxdim=n*nvec
      open(13,file=file13,form='UNFORMATTED',status='UNKNOWN')
      rewind 13
      write (13) ((b(i,j),i=1,n),j=1,nvec)
      close (13)
C 
C     IF THE ALGORITHM HAS CONVERGED FOR THE LAST VECTOR, 
C     CHECK EIGENVALUE PROBLEM AND RETURN
C 
      IF (ZEND) then
      write (6,7)
      write (6,11) (p(i),i=1,nvec)
      call flush(6)   !tolgo flush
      icount=icount+1
c      if (icount.gt.5) return
c      if (thr1(icount).lt.0.d0) return
      secrconf=secr**2
      if(icount.gt.5.or.thr1(icount).lt.0.d0)then
        allocate(cdet(ncf*nvec+1000))
        call zeroe(cdet,ncf*nvec+1000)
         call mkdet(c,cdet,n,nvec,secrconf)
c         do i=1,id7
         do i=1,ncf*nvec
	    cconf(i)=c(i) !DA SISTEMARE (3 vettori quando ne servono 2!!)
            c(i)=cdet(i)
         enddo
         deallocate(cdet)
         deallocate(ort)
         deallocate(gue)
         deallocate(ibuf)
         deallocate(jbuf)
         deallocate(a)
         deallocate (abnew)
         deallocate(bab)
         deallocate(h)
         deallocate(alpha)
         deallocate(p)
         deallocate(olsen)
         deallocate(uvec)
         return
      endif
      write (6,13) thr1(icount)
      m=nvec
      goto 15
      endif
      IF (.not.zss) GOTO 260
C 
C     CHECK THAT THE B VECTORS HAVE THE REQUIRED SYMMETRY
C     (OVERLAP WITH GUESS)
C 
      DO 220 L=1,NVEC
  220 ZFOUND(L)=.FALSE.
      L=1
  230 SMAX=0.D0
      DO 240 LL=1,NVEC
      SUM=0.D0
      IF (ZFOUND(LL)) GOTO 240
      DO 235 IG=1,NGUE
  235 SUM=SUM+B(ISEL(IG),LL)*GUE(IG,LL)
      SUM=DABS(SUM)
      IF (SUM.GT.SMAX) LMAX=LL
      IF (SUM.GT.SMAX) SMAX=SUM
  240 CONTINUE
      IF (SMAX.LT.1.D-3) GOTO 245
      IF (ZWRT) WRITE (6,242) L,LMAX,SMAX
  242 FORMAT (5X,'<B(',I2,')!GUE(',I2,')> =',F12.8)
      ZFOUND(LMAX)=.TRUE.
      IF (L.EQ.NVEC) GOTO 260
      L=L+1
      GOTO 230
  245 WRITE (6,9) L
      M=M-1
      IF (KB.GT.L) KB=KB-1
      IF (M.LT.NVEC) WRITE (6,10)
      IF (M.LT.NVEC) STOP
      IL=LX
      DO 250 LL=L,M
      DO 248 I=1,N
      IL=IL+1
      ILN=IL+N
      B(I,LL)=B(I,LL+1)
  248 AB(I,LL)=AB(I,LL+1)
  250 ALAM(LL)=ALAM(LL+1)
      GOTO 230
  260 M=NVEC
C 
C     GO BACK TO CALCULATE ANOTHER CORRECTION VECTOR Q FOR THE K-TH
C     EIGENVECTOR, OR, IF CONVERGED, FOR THE (K+1)-TH
C     IF IOP=1, RECALCULATE ALSO THE BAB MATRIX
C 
      DO 270 L=1,M
      DO 265 LL=1,M
      BAB(L,LL)=0.D0
  265 ALPHA(L,LL)=0.D0
      BAB(L,L)=ALAM(L)
  270 ALPHA(L,L)=1.D0
      DENE=1.D9
      IF (.NOT.ZCONV.AND.IOP.EQ.1) GOTO 20
      IF (.NOT.ZCONV) GOTO 75
      ITER=1
      K=K+1
      KR=K
      KB=K
      WRITE (6,5) K
      ZCONV=.FALSE.
      IF (IOP.EQ.1) GOTO 20
      GOTO 75
  400 WRITE (6,8) ITER
      STOP 1
      END
c------------------------------------------------------------------
      SUBROUTINE GUESS(N,NVEC,NGUE,NORT,DIAG)
C 
C     SETS UP TRIAL VECTORS FOR DAVIDSON'S METHOD
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      INTEGER*2 ietats
c      COMMON /VEC/ HG(ID12,ID12),VG(ID12,ID12)
      allocatable hg(:,:),vg(:,:)
      pointer ort,gue,ibuf,jbuf,a
      dimension ort(:,:),gue(:,:),ibuf(:),jbuf(:),a(:)
      allocatable inv(:),d(:),w(:)
      COMMON /INT2/ NI4,NI3,NI2,NR4,
     * ORT,GUE,IBUF,JBUF,
     * A
      DIMENSION DIAG(*)
c      dimension ibuf(100000),jbuf(100000),buff(100000)
c      INTEGER*2 II
      allocate(hg(id12,id12))
      allocate(vg(id12,id12))
      allocate(inv(n))
      allocate(d(id12))
      allocate(w(id12*3))
      nvecg=id7*8+id14*8+id6*9*8+id3*8+id10*2*8+id6*8*2*8+9*4
      nvecl=id12*id12*2*8
      nintg=id2*3*8+id2*4+599*4+4000*2+id2*3*2+16
      nintl=id12*id6*2*8+5000*2+5000*8+id12*2*8+id4*4
      do 5 i=1,id12
      do 5 j=1,id12
    5 hg(i,j)=0.d0
      N2=N+N
      NGUE=0
      IF (.NOT.ZSS) GOTO 20
   10 IF (ISEL(NGUE+1).EQ.0) GOTO 20
      NGUE=NGUE+1
      GOTO 10
   20 IF (NGUE.EQ.maxgue) GOTO 80
      IF (NGUE.GT.maxgue) GOTO 200
      REWIND 20
      READ(20)(DIAG(IIII),IIII=1,N)
      IF (NGUE.EQ.0) GOTO 30
      DO 25 I=1,NGUE
      J=ISEL(I)
      HG(I,I)=DIAG(J)
   25 DIAG(J)=4.D10
   30 EMIN=2.D10
      K=1
   32 KK=K
      IF (DIAG(K).GE.EMIN) GOTO 60
   35 IF (KK.GE.N) GOTO 40
      IF (DABS(DIAG(KK+1)-DIAG(K)).GT.1.D-4) GOTO 40
      KK=KK+1
      GOTO 35
   40 IF (NGUE+KK-K+1.LE.maxgue) GOTO 50
      DO 45 L=K,KK
   45 DIAG(L)=4.D10
      GOTO 60
   50 K1=K
      K2=KK
      EMIN=DIAG(K)
   60 K=KK+1
      IF (K.LE.N) GOTO 32
      IF (EMIN.GT.1.D10) GOTO 80
      DO 70 K=K1,K2
      NGUE=NGUE+1
      HG(NGUE,NGUE)=diag(k)
      ISEL(NGUE)=K
   70 DIAG(K)=4.D10
      IF (NGUE.LT.maxgue) GOTO 30
   80 IF (NGUE.EQ.0) GOTO 200
      IF (ZWRT) WRITE (6,85) NGUE
   85 FORMAT (/5X,'CI MATRIX IN THE TRIAL SUBSPACE, DIMENSION =',I3)
      DO 90 K=1,N
      INV(K)=0
      DO 90 I=1,NGUE
   90 IF (ISEL(I).EQ.K) INV(K)=I
      IF (NGUE.EQ.1) GOTO 115
      DO 95 I=2,NGUE
      IM=I-1
      DO 95 J=1,IM
      HG(I,J)=0.D0
   95 HG(J,I)=0.D0
c      ih=1
c      jh=0
c      MAX=5000
c      ibasta=1
      MAX=100000
c      rewind 20
c      do while (ibasta.eq.1)
cele 08/04/2004 98   read (20) (ibuf(i),i=1,size(ibuf))
   98 read (20) (ibuf(i),i=1,size(ibuf)),(jbuf(i),i=1,size(ibuf)),(a(i)
     $     ,i=1,size(a))
c      read (20) ibuf,jbuf,a
      if (ibuf(100000).eq.-1) max=jbuf(100000)
      DO 100 ic=1,MAX
c      ijh=ii(ic)
c      call ijhmat(ijh,ih,jh)
         ih=ibuf(ic)
         jh=jbuf(ic)
         invi=inv(ih)
         invj=inv(jh)
         IF (INVI.EQ.0.OR.INVJ.EQ.0) GOTO 100
         HG(INVI,INVJ)=a(ic)
         HG(INVJ,INVI)=a(ic)
 100  CONTINUE
      if (ibuf(100000).ne.-1) goto 98
c     enddo
      IF (ZWRT) then
      DO 112 I=1,NGUE
      WRITE (6,121) ISEL(I)
  112 WRITE (6,113) (HG(I,J),J=1,I)
      endif
  113 FORMAT (8F14.8)
c  115 CALL RIERA(NGUE,HG,VG,D,W,ID12)
 115  continue
      do i=1,ngue
         do j=1,ngue
            vg(i,j)=hg(i,j)
         enddo
      enddo
      call dsyev('V','U',ngue,vg,id12,d,w,3*id12,info)
      IF (ZWRT) then
      WRITE (6,118)
  118 FORMAT (/5X,'EIGENVALUES AND EIGENVECTORS')
      DO 120 J=1,NGUE
      WRITE (6,121) J,D(J)
  120 WRITE (6,122) (VG(I,J),I=1,NGUE)
      endif
  121 FORMAT (I7,F15.8)
  122 FORMAT (2X,10F12.8)
      K=0
      NORT=0
      L=0
  125 K=K+1
      IF (K.GT.NGUE) GOTO 220
      IF (.NOT.ZSS) GOTO 150
      SUM=0.D0
      DO 130 I=1,NGUE
  130 SUM=SUM+CSEL(I)*VG(I,K)
      IF (DABS(SUM).GT.1.D-6) GOTO 150
      IF (NORT.GE.ID6) GOTO 125
      NORT=NORT+1
      DO 140 I=1,NGUE
  140 ORT(I,NORT)=VG(I,K)
      GOTO 125
  150 L=L+1
      DO 160 I=1,NGUE
  160 GUE(I,L)=VG(I,K)
      IF (L.LT.NVEC) GOTO 125
      deallocate(hg)
      deallocate(vg)
      deallocate(inv)
      deallocate(d)
      deallocate(w)
      RETURN
  200 WRITE (6,210) NGUE
  210 FORMAT (//5X,'ERROR IN SUBROUTINE GUESS, NGUE =',I4//)
      STOP
  220 WRITE (6,230) L
  230 FORMAT (//5X,'ERROR IN SUBROUTINE GUESS, ONLY',I3,
     * ' TRIAL VECTORS OF THE REQUIRED SYMMETRY HAVE BEEN FOUND'//)
      STOP
      END
C***********************************************************************
      SUBROUTINE DAVID2(ncf,N,NVEC,secr,zrest)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer diag,ibuf,jbuf,a,olsen
      dimension diag(:),olsen(:),ibuf(:),jbuf(:),a(:)
      common/mtvc/diag,a,olsen,ibuf,jbuf
      common/DVD2/tol,nbx,itmax,ilevel,ngue,nb
      pointer c,fmpb
      dimension c(:),fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      integer*2 icomp,iconf,isegno
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
c      dimension iseed(5)
      allocatable wrk(:),res(:),iwrk(:),x(:,:),cdet(:),w(:),work(:)
      allocatable hg(:,:),isel(:),inv(:)
      common /fil/ file12,file13,file26,zwconf
      character*64 file12,file13,file26
      external matvec,correc
c      nb=nvec+10
      if(nb.lt.nvec)nb=nvec
      mb=nb
c      nbx=nb+nb+nb+nb+nb+nb+nb+nb
      if(nbx.lt.mb)nbx=mb
      if(nbx.gt.n)nbx=n
      lwrk=2*n*nbx+nbx*nbx+(nbx+n)*mb+8*nbx
      liwrk=6*nbx
      ldx=n
      allocate(wrk(lwrk))
      allocate(iwrk(liwrk))
      allocate(diag(n))
      allocate(olsen(n))
      allocate(res(mb))
      allocate(x(ldx,mb))
      allocate(ibuf(100000))
      allocate(jbuf(100000))
      allocate(a(100000))
      if(ngue.gt.n)ngue=n
      if(ngue.lt.nb)ngue=nb
      allocate(hg(ngue,ngue))
      allocate(isel(n))
      allocate(inv(n))
      do i=1,n
         isel(i)=0
         inv(i)=0
      enddo
      M=NVEC

C
C     INITIAL VECTORS FROM RESTART                               
C
      if (zrest) then
      write (6,*) 'Restart della procedura di DAVIDSON'
      write (6,*) 'Vettore di guess da un calcolo precedente'
      open(13,file=file13,form='UNFORMATTED',status='UNKNOWN')
      rewind 13
      read (13) ((x(i,j),i=1,n),j=1,nvec)
      close (13)
      write (6,*) 'Guess vectors from file ',file13
      call flush(6)
      goto 36
      endif

      REWIND 20
      READ (20) (DIAG(I),I=1,n)
      do i=1,ngue
         imi=idmin(n,diag)
         isel(i)=imi
         hg(i,i)=diag(imi)
         diag(imi)=1.0D90
      enddo
      do i=1,n
         do j=1,ngue
            if(isel(j).eq.i)then
               inv(i)=j
               goto 1
            endif
         enddo
 1    enddo
      MAX=100000
 35   read (20) (ibuf(i),i=1,size(ibuf)),
     $          (jbuf(i),i=1,size(jbuf)),
     $          (a(i),i=1,size(a))

      if(ibuf(100000).eq.-1)max=jbuf(100000)
      DO  IC=1,MAX
         ih=ibuf(ic)
         jh=jbuf(ic)
         invi=inv(ih)
         invj=inv(jh)
         if(invi.ne.0.and.invj.ne.0)then
            hg(invi,invj)=a(ic)
            hg(invj,invi)=a(ic)
         endif
      enddo
      if(ibuf(100000).ne.-1)goto 35
      allocate(w(ngue))
      lwork=5*ngue
      allocate(work(lwork))
      call dsyev('V','U',ngue,hg,ngue,w,work,lwork,info)
      print*,'diagonalization of guess matrix got info=',info
     $     ,' ngue was',ngue
      call flush(6)
      call zeroe(x,n*mb)
      do j=1,mb
         do i=1,ngue
            ii=isel(i)
            x(ii,j)=hg(i,j)
         enddo
      enddo

  36  nev=-nvec
      itmax=1000
      iunit=6
      print*,'    -----------------------------------------------'
      print*,'              CALLING SUBROUTINE DAVSON'
      print*,'    -----------------------------------------------'
      print*,' '
      print '(a,i7)','DIMENSION OF THE MATRIX        ',n
      print '(a,i7)','NUMBER OF EIGENVALUES          ',nvec
      print '(a,i7)','MAXIMUM DIMENSION OF V SPACE   ',nbx
      print '(a,d7.1)','MAXIMUM ERROR on RESIDUAL NORM ',tol
      print '(a,i7)','MAXIMUM NUMBER OF ITERATIONS   ',itmax
      print '(a,i7)','DIMENSION OF THE GUESS MATRIX  ',ngue
      print '(a,i7)','VERBOSITY LEVEL                ',ilevel
         call davson(n,nbx,nb,nev,itmax,iter,nmult,tol,1.0d0,res,e,x,ldx
     $        ,wrk,lwrk,iwrk,liwrk,matvec,correc,ilevel,iunit,ifail)
         call eigsrt(n,abs(nvec),res,e,x,ldx)
         allocate(c(NCF*NVEC+1000))
            do i=1,n
               do j=1,nvec
                  c((j-1)*n+i)=x(i,j)
               enddo
            enddo
C
C     VECTORS STORED IN FILE13 FOR RESTART                               
C
      open(13,file=file13,form='UNFORMATTED',status='UNKNOWN')
      rewind 13
      write (13) ((x(i,j),i=1,n),j=1,nvec)
      close (13)

      deallocate(wrk)
      deallocate(iwrk)
      deallocate(diag)
      deallocate(olsen)
      deallocate(res)
      deallocate(x)
      deallocate(ibuf)
      deallocate(jbuf)
      deallocate(a)
      deallocate(hg)
      deallocate(inv)
      deallocate(isel)
      allocate(cdet(ncf*nvec+1000))
      call zeroe(cdet,ncf*nvec+1000)
      secrconf=secr**2
      call mkdet(c,cdet,n,nvec,secrconf)
      allocate(cconf(ncf*nvec+1000)) 
      do i=1,ncf*nvec
         cconf(i)=c(i)
         c(i)=cdet(i)
      enddo
      deallocate(cdet)
      return
      end
c-----------------------------------------------------
      integer function idmin(n,a)
      real*8 a(*),am
c      integer isel(*)
      am=1.0D90
      do i=1,n
         if(a(i).lt.am)then
            im=i
            am=a(i)
         endif
      enddo
      idmin=im
      return
      end
C***********************************************************************
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
      call flush(6)   !tolgo flush
      RETURN
      END
C***********************************************************************
      SUBROUTINE LECNEWER(idimorb,ndimdum)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
c
      integer*2 trouloc(12),partloc(12),nec
c      integer*1 noccup(id1,id4)
      integer combf
      allocatable noccup(:,:),nec(:)
      integer*1 noccup
      pointer morb(:),mspin(:),isingle(:)
      integer*2 morb,mspin,morbp
      common /ric/ ihalf,nec1,morb,isingle,mspin
      character*1 inexc(13)
      character*100 cnec
      data inexc/'F','M','D','T','Q','P','H','E','O','N',
     *  'X','U','Z'/
      INTEGER*2 NE,TROU,PART
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
c
c
c
      allocate (noccup(1:norb,1:ndimdum))
      allocate (nec(ndimdum))
      do i=1,ndimdum
      do j=1,nocb
      noccup(j,i)=2
      enddo
      do j=nocb+1,norb
      noccup(j,i)=0
      enddo
      enddo
c
      nconf=0
  100 nconf=nconf+1
      read (3,'(A)',END=10000,ERR=10000) cnec
      do
       i=index(cnec,'+')
       if(i.eq.0)exit
       cnec(i:i)=' '
      enddo
      do
       i=index(cnec,'-')
       if(i.eq.0)exit
       cnec(i:i)=' '
      enddo
      do 20 i=1,13
      if (cnec(1:1).eq.inexc(i)) goto 21
   20 continue
c
c     Ho finito di leggere i determinanti salto a 1000
c
      nconf=nconf-1
      goto 1000
c
 21   nec(nconf)=i-1
      read (cnec(2:100),*,END=10002,ERR=10002) 
     * (trouloc(i),i=1,nec(nconf)),(partloc(i),i=1,nec(nconf)) 
      do i=1,nec(nconf)
      noccup(trouloc(i),nconf)=noccup(trouloc(i),nconf)-1
      if (partloc(i).le.norb) then
      noccup(partloc(i),nconf)=noccup(partloc(i),nconf)+1
      endif
      enddo
c
c     Procedo ora all'eliminazione delle occupazioni orbitaliche 
c     doppie 
c
      do iconf=1,nconf-1
      if (znequ(noccup(1,iconf),noccup(1,nconf),norb)) then
       do i=1,nec(nconf)
       noccup(trouloc(i),nconf)=noccup(trouloc(i),nconf)+1
       if (partloc(i).le.norb)noccup(partloc(i),nconf)=noccup(partloc(i)
     $      ,nconf)-1
      enddo 
       nconf=nconf-1
       goto 100
      endif
      enddo

      if (nconf+1.gt.ndimdum) then
      write (6,*) 'Errore: troppi determinanti in lettura'
      write (6,*) 'nconf,id4',nconf,ndimdum
c      stop
      endif
      goto 100
c
c     Si riordinano le occupazioni orbitaliche secondo nec crescente
c
 1000 do iconf=1,nconf
       do jconf=iconf+1,nconf
        if (nec(jconf).lt.nec(iconf)) then
	 do i=1,norb
	  iocc=noccup(i,jconf)
	  noccup(i,jconf)=noccup(i,iconf)
	  noccup(i,iconf)=iocc
       enddo
	 idum=nec(iconf)
	 nec(iconf)=nec(jconf)
	 nec(jconf)=idum
        endif
      enddo
      enddo

c
c     Si generano ora tutti i determinanti dello spazio
c

c----calcolo delle dimensioni----
      isz2=isz*2
      if(zion)isz2=isz2-1
      ncf=0
      ndimtp=0
      do i=1,nconf
         nsing=0
         do j=1,norb
            if(noccup(j,i).eq.1)nsing=nsing+1
         enddo
         nha=(nsing-isz2)/2
         nde=combf(nsing,nha)
         ncf=ncf+nde
         ndimtp=ndimtp+nde*nec(i)
      enddo
      call givdimen(ncf,ndimtp,idimorb)
c------------------------      
      ncf=0
      do iconf=1,nconf
      nf=ncf+1
      ndnf=0
      if (ncf.gt.0) ndnf=nd(ncf)+ne(ncf)
      call scara(noccup(1,iconf),norb,nocb,ncf1,nd(nf),ne(nf),
     * trou(ndnf+1),part(ndnf+1),isz,zion)
      do i=1,ncf1
      nd(ncf+i)=nd(ncf+i)+ndnf
      enddo
      ncf=ncf+ncf1
      enddo
      deallocate (noccup)
      return
       
10000 WRITE (6,10001)
10001 FORMAT (' ERREUR EN LECTURE DETERMINANT')
      STOP
10002 WRITE (6,10003) nconf
10003 FORMAT (' ERROR READING TP IN THE DETERMINANT',I7)
      WRITE (6,'(A)') cnec
      STOP
      END
C***********************************************************************
      SUBROUTINE TRANSF(NEC,IPERM)
C 
C     CETTE SUBROUTINE TESTE LA VALIDITE DES EXCITATIONS
C     ELLE ELIMINE LES DETERMINANTS NON ORTHODOXES
C     ET NE CONSERVE QUE LES EXCITATIONS CONFORMES
C        LEGALITE, CONFORMISME, ORTHODOXIE ! ! !
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      PARAMETER (JDTOT=156*ID5)
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      COMMON /VEC2/ IDUM(JDTOT),JDUM(13),NT(13),NP(13),MT(13),MP(13)
      INTEGER*2 IDUM
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 IPERM(10)
      EQUIVALENCE (K,NCF)
      MNOR1=NORB2+1
      NCFI=K+1
      IF (NEC-1) 10,30,20
   10 K=K+1
      NE(K)=0
      ND(K)=0
      TROU(1)=0
      PART(1)=0
      IF (ISZ.NE.0) K=K-1
      ISYM=1
      write (6,36) isym
      RETURN
C 
C     LES SPINORBITALES BETA SONT NUMEROTEES DE INITB(=1) A NORB
C     LES SPINORBITALES ALPHA DE INITA(=NORB+1) A NORB2(=2*NORB)
C     ORDONNONS TROUS ET PARTICULES PAR NUMEROS CROISSANTS
C 
   20 NECM=NEC-1
      DO 25 I=1,NECM
      IP=I+1
      DO 25 J=IP,NEC
      IF (NT(J).GT.NT(I)) GOTO 24
      N=NT(J)
      NT(J)=NT(I)
      NT(I)=N
   24 IF (NP(J).GT.NP(I)) GOTO 25
      N=NP(J)
      NP(J)=NP(I)
      NP(I)=N
   25 CONTINUE
   30 IS=1
      LS=0
C 
C  TEST SUR LES NUMEROS DES TROUS ET DES PARTICULES
C 
      IF (NT(1).LT.INITB) GOTO 80
      IF (NP(NEC).GT.NORB2.AND..NOT.ZION)GOTO 80
      IF (NT(NEC).GT.NOCA) GOTO 80
      IF (NP(NEC).LE. NOCB) GOTO 80
      DO 35 I=1,NEC
      IF (NT(I).LE.NOCB) GOTO 31
      IF (NT(I).LT.INITA) GOTO 80
   31 IF (NP(I).GT.NOCA) GOTO 32
      IF (NP(I).GT.NORB) GOTO 79
      IF (NP(I).LE.NOCB) GOTO 80
C 
C     TEST SUR LA SYMETRIE ET LE SPIN DES EXCITATIONS
C 
   32 IS=ITS(IPERM(IS),ITSYM(NT(I)))
c     if (is.gt.nsym) write (6,*) 'sr TRANSF, I,J,ITS=',
c    * 999,itsym(nt(i)),is
      IS=ITS(IPERM(IS),ITSYM(NP(I)))
c     if (is.gt.nsym) write (6,*) 'sr TRANSF, I,J,ITS=',
c    * 999,itsym(np(i)),is
      LS=LS+ISPIN(NP(I))
      LS=LS-ISPIN(NT(I))
   35 CONTINUE
      IF (ISYM.EQ.0) then
      ISYM=IS
      write (6,36) isym
   36 format (/'     SYMETRIE DE LA FONCTION D''ONDE:',i4)
      endif
      IF (IS.NE.ISYM) GOTO 90
      IF (LS.NE.ISZ) GOTO 96
C     ON GARDE LES BONS
      IF (K.EQ.0) GOTO 50
      DO 40 I=1,K
      IF (NE(I).NE.NEC) GOTO 40
      NDK=ND(I)
      DO 39 J=1,NEC
      IF (TROU(NDK+J).NE.NT(J)) GOTO 40
      IF (PART(NDK+J).NE.NP(J)) GOTO 40
   39 CONTINUE
      GOTO 70
   40 CONTINUE
      NDK=ND(K)+NE(K)
c  Tolgo questa istruzione che introduce un posto vuoto in trou e part
c  se c'e` il det F (M. Persico, luglio 97)
c     IF (NE(K).EQ.0) NDK=ND(K)+1
      GOTO 51
   50 NDK=0
   51 K=K+1
c      IF (K.GT.ID5) GOTO 56
      NE(K)=NEC
      ND(K)=NDK
      DO 60 J=1,NEC
cele      NTO=NT(J)
cele      IF (NTO.GT.NORB) NTO=NTO-NORB
cele      ZFREE(NTO)=.FALSE.
cele      NPO=NP(J)
cele      IF (NPO.GT.NORB) NPO=NPO-NORB
cele      ZFREE(NPO)=.FALSE.
cele      ZOCS(NT(J))=.TRUE.
cele      ZOCS(NP(J))=.TRUE.
      TROU(NDK+J)=NT(J)
   60 PART(NDK+J)=NP(J)
   70 CONTINUE
      RETURN
   79 IF (ZION) RETURN
   80 WRITE (6,85)
   85 FORMAT (' 1 SPIN ORBITALE EST HORS DE L''ESPACE DE BASE')
      WRITE (6,*) (NT(I),NP(I),I=1,NEC)
      STOP
   90 WRITE (6,95)
      WRITE (6,*) (NT(I),NP(I),I=1,NEC)
   95 FORMAT (' DETERMINANT DE MAUVAISE SYMETRIE')
      STOP
   96 WRITE (6,97)
   97 FORMAT (' MAUVAIS SPIN')
      WRITE (6,*) (NT(I),NP(I),I=1,NEC)
      STOP
   56 WRITE (6,99) ID5
   99 FORMAT (' TROP DE DETERMINANTS LUS DU FICHIER 5: MAXIMUM',I5/)
      STOP
      END
C***********************************************************************
      SUBROUTINE WRTVAR(NVEC,ZVBIS,zbin)
C 
C     RESULTATS DE L'I.C. VARIATIONNELLE
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      common/third/zthird,zmp3,zen3
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      INTEGER*2 ietats
 1000 FORMAT (54X,I4,3X,6F11.8)
 1006 FORMAT (///'     VALEURS ET VECTEURS PROPRES:'/)
 1008 FORMAT (//)
 1010 FORMAT (5X,5F15.8)
 1016 FORMAT (//'     + ENERGIE SCF: ',4f15.8/10(5x,5f15.8/))
 1200 FORMAT (///5X,'VECTEURS TRANSFORMES ET/OU SYMETRISES'/)
      if (zvbis) then
      WRITE (6,1200)
      else
      WRITE (6,1006)
      WRITE (6,1010) (E(J),J=1,nvec)
      DO J=1,nvec
      wk(J)=E(J)+ESCF
      enddo
      WRITE (6,1016) (wk(J),J=1,nvec)
      WRITE (6,1008)
      endif
      if (zpun.and..not.zthird) rewind 7
      DO 15 M=1,METAT
      CMAX=0.D0
      JM=(m-1)*ncf
      DO 10 K=1,NCF
      MK=JM+K
      DC=DABS(C(MK))
   10 IF (CMAX.LT.DC) CMAX=DC
   15 TESTE(M)=CMAX
      
      if(zthird)goto 1717
      if (zbin) then

      zfa=.false.
      if (zpun) then
      WRITE (7) METAT,tau,nshift,zauto,zfa
c      WRITE (7) (TESTE(M),M=1,METAT)
      WRITE (7) (1.d0,M=1,METAT)
      endif

      else

      if (zpun.and..not.zauto) WRITE (7,'(''C''/i3,2X,f14.8,I3)') 
     * METAT,tau,nshift
      if (zpun.and.zauto) WRITE (7,16) METAT,tau,nshift
   16 format ('C',3X,'Attenzione: file preparato per ZAUTO=T in cipselx'
     * /i3,2X,f14.8,I3) 
      if (zpun) WRITE (7,'(6f12.8)') (TESTE(M),M=1,METAT)

      endif
 1717 continue


      DO 18 M=1,METAT
      TAUE(M)=TESTE(M)*TAU
   18 TESTE(M)=TESTE(M)*TEST
      MAX=0
      zpupun=.true.
   20 MIN=MAX+1
      MAX=MAX+6
      IF (max.gt.nvec) max=nvec
      DO 40 K=1,NCF
      JJ=K+NCF*(MIN-2)
      zet=.false.
      DO 30 J=min,max
      JJ=JJ+NCF
      wk(J)=C(JJ)
   30 IF (DABS(wk(J)).GE.SECR) ZET=.TRUE.
      IK=K
      if (zbin) then
         CALL WRTBIN(IK,min,max,FDUM,gdum,zet, zpupun)
      else
         CALL WRTDET(IK,min,max,FDUM,gdum,zet, zpupun)
      endif
   40 continue
      WRITE (6,1008)
      zpupun=.false.
      IF (MAX.LT.nvec) GOTO 20
      RETURN
      END
C***********************************************************************
      SUBROUTINE ELIM(ZVBIS)
C 
C     ELIMINATION DES ETATS INTRUS OU NON DESIRES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      INTEGER*2 ietats
      if (nrot.eq.ietats(nrot)) return
      write (6,20) (ietats(i),i=1,nrot)
   20 format (//' SEUL LES ETATS SUIVANTS SERONT RETENUS:'/26I5)
      zvbis=.true.
      do 50 i=1,nrot
      j=ietats(i)
      e(i)=e(j)
      ncfi=ncf*i-ncf
      ncfj=ncf*j-ncf
      do 50 k=1,ncf
   50 c(ncfi+k)=c(ncfj+k)
      return
      end
C***********************************************************************
      SUBROUTINE IJHMAT(IJ,I,J)
      j=j+ij
   10 continue
      if (j.gt.i) then
      j=j-i
      i=i+1
      goto 10
      endif
      RETURN
      END
C***********************************************************************
      SUBROUTINE PROJV(N,NVEC,NCFPR,B,iditrpa,ncf)
C 
C     SORTS INITIAL COEFFICIENTS FOR DAVIDSON'S METHOD, FROM VECTORS
C     IN A SMALLER SUBSPACE STORED IN FILE 14
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
c      PARAMETER (ID82=ID8*2)
      dimension b(1)
      allocatable trpapr(:)
c      integer*2 trpapr(id82)
      integer*2 trpapr
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      INTEGER*2 ietats
      allocatable ndpr(:),nepr(:)
      INTEGER*2 NEPR
      pointer c
      dimension c(:)
      COMMON /VEC/ C
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      allocate(ndpr(ncf+2))
      allocate(nepr(ncf+2))
      allocate(trpapr(id82))
      open(14,file=projf,form='UNFORMATTED',status='OLD')
      rewind 14
      READ (14) idum,idum,idum,idum,NCFPR,i1,i2,
     * (NEPR(J),NDPR(J),J=1,NCFPR),
     * (TRPAPR(J),TRPAPR(J+IDtrpa),J=1,I1),(c(J),J=1,I2)
      close(14)
      nnn=n*nvec
      do 20 k=1,nnn
   20 b(k)=0.d0
      do 500 ipr=1,ncfpr
      necpr=nepr(ipr)
      i=0
   50 i=i+1
      if (i.gt.n) then
      write (6,60) ipr
   60 format (//'   SUBROUTINE PROJV:  DET N.',I5,' FROM FILE 14'/
     * '   NOT FOUND IN THE CURRENT S SUBSPACE.'/
     * '   THIS CASE NOT TREATED IN THE PRESENT VERSION'/
     * '   (VECTORS WILL NOT BE ORTHONORMALISED):   STOP'/)
      stop 10
      endif
      if (ne(i).ne.necpr) goto 50
      ndi=nd(i)
      ndipr=ndpr(ipr)
      do 100 j=1,necpr
      if (trou(ndi+j).ne.trpapr(ndipr+j)) goto 50
      if (part(ndi+j).ne.trpapr(ndipr+j+iditrpa)) goto 50
  100 continue
      do 200 k=1,nvec
      kncf=(k-1)*n
      kncfpr=(k-1)*ncfpr
  200 b(kncf+i)=c(kncfpr+ipr)
  500 continue
      deallocate(ndpr)
      deallocate(nepr)
      return
      end
C***********************************************************************
      SUBROUTINE DEGEN(zdet,ZVBIS)
C 
C     SYMETRISATION DE VECTEURS DEGENERES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
c      COMMON /INT/ ZOS(ID4),NESP(ID4)
      allocatable zos(:),nesp(:)
      INTEGER*2 NESP
      dimension zmuld(id6)
      nintg=id2*3*8+id2*4+599*4+4000*2+id2*3*2+16
      nintl=id4+id4*2
      allocate(zos(ncf))
      allocate(nesp(ncf))
      REWIND 4
      READ (4)
      if (zdet.and..not.zpertu) then
      READ (4) ii,(idum,KG=1,ii),ISYM,(NESP(K),K=1,NCF)
      else
      READ (4) NCFG,(KGEN(KG),KG=1,NCFG),ISYM,(NESP(K),K=1,NCF)
      endif
      do 50 i=1,nrot
   50 zmuld(i)=.false.
      DO 100 KG=1,NCFG
      K=KGEN(KG)
  100 NDEG(KG)=NESP(K)
      IF (NROT.EQ.1) GOTO 410
      MMAX=NROT-1
      JJ=0
      DO 400 M=1,MMAX
      IF (M.LE.JJ) GOTO 400
      JJ=M
      MN1=M+1
      DO 435 N=MN1,NROT
  435 IF (DABS(E(N)-E(M)).LE.TDEGEN) JJ=N
      IF (JJ.EQ.M) GOTO 400
      IF (JJ.EQ.MN1) GOTO 420
      WRITE (6,415)
  415 FORMAT (//5X,'ATTENTION DEGENERESCENCE MULTIPLE'/
     * 5X,'SYMETRISATION NON PREVUE DANS CE CAS.  POUR LES ENERGIES',
     * 5X,'ON PEUT TRAITER SIMULTANEMENT TOUS LES ETATS DEGENERES'/
     * 5X,'ET DIVISER LA SOMME DE LEURS ENERGIES PERTURBEES PAR LA'/
     * 5X,'DEGENERESCENCE.  POUR LES VECTEURS MODIFIER DEGEN'//)
      do 418 i=m,jj
  418 zmuld(i)=.true.
      GOTO 400
  420 MNCF1=M*NCF-NCF
      MNCF2=MNCF1+NCF
  416 FORMAT (//5X,'LES ETATS ',I2,' ET ',I2,' SONT DEGENERES :',
     *  5X,'SYMETRISATION PAR ROTATION'//)
C     ROTATION DES ETATS
      DO 450 K=1,NCF
      IF (ZSAUT(K)) GOTO 450
      NDEGK=NDEG(K)
      IF (NDEGK.EQ.0) GOTO 450
      DO 440 LL=1,NDEGK
      C11=C(MNCF1+K)
      C12=C(MNCF1+K+LL)
      C1=DABS(C11)
      C2=DABS(C12)
      IF (DABS(C1-C2).LT.1.D-5) GOTO 440
      C21=C(MNCF2+K)
      C22=C(MNCF2+K+LL)
      CC=C21+C22
      IF (DABS(CC).LT.1.D-8) GOTO 440
      TGPHI=-(C11+C12)/CC
      PHI=DATAN(TGPHI)
      CS=DCOS(PHI)
      SI=DSIN(PHI)
      C1=C11*CS+C21*SI
      C2=C12*CS+C22*SI
      C1=DABS(C1)
      C2=DABS(C2)
      C3=DABS(C1-C2)
      IF (C3.GE.1.D-4) GOTO 440
      WRITE (6,416) M,MN1
      ZVBIS=.TRUE.
      DO 480 I=1,NCF
      A=C(MNCF1+I)*CS+C(MNCF2+I)*SI
      C(MNCF2+I)=-C(MNCF1+I)*SI+C(MNCF2+I)*CS
  480 C(MNCF1+I)=A
      GOTO 400
  440 CONTINUE
  450 CONTINUE
  400 CONTINUE
  410 CONTINUE
      IF (ISZ.eq.0.and..not.zion) then
         do k=1,norb+norb+1
            zoc(k)=.false.
         enddo
      DO 510 K=1,NCF
      NDK=ND(K)
      NEK=NE(K)
      ZOS(K)=.FALSE.
      IF (NEK.EQ.0) GOTO 510
      NFK=NDK+NEK
      NDK=NDK+1
      DO 490 I=NDK,NFK
      ZOC(TROU(I))=.TRUE.
  490 ZOC(PART(I))=.TRUE.
      DO 495 I=NDK,NFK
      ITR=TROU(I)
      IFLIP=ITR+NORB
      IF (ITR.GT.NORB) IFLIP=ITR-NORB
      IF (.NOT.ZOC(IFLIP)) ZOS(K)=.TRUE.
      IPA=PART(I)
      IFLIP=IPA+NORB
      IF (IPA.GT.NORB) IFLIP=IPA-NORB
  495 IF (.NOT.ZOC(IFLIP)) ZOS(K)=.TRUE.
      DO 500 I=NDK,NFK
      ZOC(TROU(I))=.FALSE.
  500 ZOC(PART(I))=.FALSE.
  510 CONTINUE
      endif
C
C     SYMETRISATION DES COEFFICIENTS DE DETERMINANTS NE DIFFERANT
C     QUE PAR SPIN FLIP
C
      PIVOT=1.D0/DSQRT(DFLOAT(NCF))
      MNCF=-NCF
      DO 800 M=1,NROT
      MNCF=MNCF+NCF
      if (zmuld(m).or.isz.ne.0.or.zion) goto 600
      K=0
  520 K=K+1
      IF (DABS(C(MNCF+K)).LT.PIVOT) GOTO 520
      ZTRIP=.FALSE.
      IF (.NOT.ZOS(K)) GOTO 530
      DD=DABS(C(MNCF+K)+C(MNCF+K+1))
      IF (DD.LT.DABS(C(MNCF+K))) ZTRIP=.TRUE.
  530 K=0
  540 K=K+1
      IF (K.GT.NCF) GOTO 600
      IF (ZOS(K)) GOTO 550
      IF (ZTRIP) C(MNCF+K)=0.D0
      GOTO 540
  550 CC=C(MNCF+K)
      CCC=C(MNCF+K+1)
      sig=1.d0
      if (dabs(cc).gt.dabs(ccc)) then
      if (cc.lt.0.d0) sig=-sig
      else
      if (ccc.lt.0.d0) sig=-sig
      if (ztrip) sig=-sig
      endif
      CC=sig*dsqrt(0.5D0*(cc*cc+CCC*CCC))
      C(MNCF+K)=CC
      IF (ZTRIP) CC=-CC
      C(MNCF+K+1)=CC
      K=K+1
      GOTO 540
  600 CONTINUE
c     write (6,*) ' dopo spin flip, vettore ',m
c     write (6,2222) (c(mncf+k),k=1,ncf)
C
C     ORTHOGONALISATION
C
      do 700 mm=1,m-1
      mmncf=(mm-1)*ncf
      sum=0.d0
      do 620 k=1,ncf
  620 sum=sum+c(mncf+k)*c(mmncf+k)
c     write (6,*) ' overlap ',mm,m,sum
      if (dabs(sum).gt.1.d-5) write (6,630) mm,m,sum
  630 format ('      SOUSPROGR. DEGEN, ATTENTION: OVERLAP ENTRE LES',
     * ' VECTEURS',I3,' ET',I3,'  =',f15.8)
      do 640 k=1,ncf
  640 c(mncf+k)=c(mncf+k)-c(mmncf+k)*sum
  700 continue
C
C     NORMALISATION
C
      sum=0.d0
      do 720 k=1,ncf
  720 sum=sum+c(mncf+k)*c(mncf+k)
      sum=dsqrt(sum)
c     write (6,*) ' norma ',m,sum
      if (dabs(sum-1.d0).gt.1.d-5) write (6,730) m,sum
  730 format ('      SOUSPROGR. DEGEN, ATTENTION: NORME DU VECTEUR',
     * I3,'  =',f15.8)
      sum=1.d0/sum
      do 740 k=1,ncf
  740 c(mncf+k)=c(mncf+k)*sum
 800  CONTINUE
      deallocate(zos)
      deallocate(nesp)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SPIN(NVEC,ZVBIS)
C 
C     CALCOLO DI S2
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
c      COMMON /INT/ smat(id6,id6),ts(id6,id6),es(id6),wrk(id6),cpro(id6)
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      allocatable smat(:,:),ts(:,:),es(:),wrk(:),cpro(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      COMMON /ALE/SPINMAX,SPINMIN,ZSPINSEL,ZINCRE
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN,IETATS
      WRITE (6,5)
    5 FORMAT (///5X,'CONSTRUCTION DE LA MATRICE DE S2'/)
      call flush(6)   !tolgo flush
      allocate(cpro(id6))
      allocate(wrk(id6*3))
      allocate(es(id6))
      allocate(ts(id6,id6))
      allocate(smat(nvec,nvec))
      do 9 i=1,nvec
      do 9 j=1,nvec
    9 smat(i,j)=0.d0
      do 500 k=1,ncf
c     if (ne(k).ge.nocb) then
c     ndoc=0
c     else
c     ndoc=nocb-ne(k)
c     endif
      sz=float(isz)
      if (zion) sz=float(isz)-0.5d0
      nsoc=2*ne(k)
      if (zion) nsoc=nsoc-1
      do 70 l=1,ne(k)
      do 70 m=1,l-1
      if (iorb(trou(nd(k)+l)).eq.iorb(trou(nd(k)+m))) nsoc=nsoc-2 
   70 if (iorb(part(nd(k)+l)).eq.iorb(part(nd(k)+m))) nsoc=nsoc-2
c     if (iorb(trou(nd(k)+l)).eq.iorb(trou(nd(k)+m))) ndoc=ndoc+1 
c  70 if (iorb(part(nd(k)+l)).eq.iorb(part(nd(k)+m))) ndoc=ndoc+1
c     if (zion) then
c     dssd=float(nocb-ndoc)-0.5d0+(float(isz)-0.5d0)**2
c     else
c     dssd=nocb-ndoc+isz*isz
c     endif
      dssd=sz*sz+0.5d0*float(nsoc)
      do 150 i=1,nvec
      ci=c(ncf*i-ncf+k)
      do 150 j=1,i
      cj=c(ncf*j-ncf+k)
  150 smat(i,j)=smat(i,j)+ci*dssd*cj  
      do 500 kk=1,k-1
      dskd=0.d0
      if (ne(k).ne.ne(kk)) goto 500
      call detdif(k,kk,nbdif,n1,n2,n3,n4,ns1,ns2,ns3,ns4,zsig)
      if (nbdif.gt.2.or.nbdif.le.1) goto 500
      if (ns1.eq.ns3) goto 500
      if (n1.eq.n3) goto 500
      if (n1.eq.n2.and.n3.eq.n4) dskd=1.d0
      if (n1.eq.n4.and.n3.eq.n2) dskd=-1.d0
      do 190 i=1,nvec
      ci=c(ncf*i-ncf+k)
      cii=c(ncf*i-ncf+kk)
      do 190 j=1,i
      cj=c(ncf*j-ncf+kk)
      cjj=c(ncf*j-ncf+k)
      if (zsig) then
      smat(i,j)=smat(i,j)-dskd*(ci*cj+cii*cjj)
      else
      smat(i,j)=smat(i,j)+dskd*(ci*cj+cii*cjj)
      endif
  190 continue
  500 continue
      write(6,*) ' matrice s2 '
      do 600 i=1,nvec
      do 610 l=1,i
  610 smat(l,i)=smat(i,l)
  600 write(6,66) (smat(i,j),j=1,i) 
      if (zvbis) return
      if (nvec.eq.1) return
c      call riera(nvec,smat,ts,es,wrk,id6)
      do i=1,nvec
         do j=1,nvec
            ts(i,j)=smat(i,j)
         enddo
      enddo
      call dsyev('V','U',nvec,ts,id6,es,wrk,3*id6,info)
      write (6,650) 
  650 format (//5X,'VALEURS ET VECTEURS PROPRES DE S2'/)
      write (6,66) (es(i),i=1,nvec)
      write (6,'(/)')
      do 698 k=1,nvec
  698 write (6,66) (ts(k,kk),kk=1,nvec)
C alessandro: luglio 98
      if (ZSPINSEL) then
       zincre=.FALSE.
       ik=1
       do i=1,nvec
        tsmax=0.d0
        if (es(i).lt.spinmax.and.es(i).gt.spinmin) then
 	 do j=1,nvec
 	  do k=1,ik
 	   if (j.eq.ietats(k)) goto 750
 	  enddo
 	  if (abs(ts(j,i)).gt.tsmax) then
	   tsmax=abs(ts(j,i))
	   jmax=j
 	  endif
  750    enddo
	 ietats(ik)=jmax
	 ik=ik+1
        endif
       enddo
       ik=ik-1
       if (ik.lt.metat) then
        zincre=.TRUE.
        return
       else
	ipun=1
  755	minie=id6
	do j=ipun,ik
	 if (ietats(j).lt.minie) then
	  minie=ietats(j)
	  jtar=j
	 endif
	enddo
	idum=ietats(ipun)
	ietats(ipun)=ietats(jtar)
	ietats(jtar)=idum
	if (ipun.lt.ik) then
	 ipun=ipun+1
	 goto 755
	endif
       endif
      endif
C alessandro: luglio 98
      if (.not.zsrot) return
      write (6,800)
  800 format(/5X,'Rotazione degli stati in modo che S2 sia ',
     * 'diagonale'/)
      do 700 i=1,ncf
      do 790 k=1,nvec
      cpro(k)=0.d0
      do 790 l=1,nvec
  790 cpro(k)=c(ncf*l-ncf+i)*ts(l,k)+cpro(k)
      do 795 k=1,nvec
  795 c(ncf*k-ncf+i)=cpro(k)
  700 continue
      IJ=0
      DO 860 IE=1,nvec
      DO 860 JE=1,IE
      IF (JE.NE.IE) IJ=IJ+1
      IH=IJ
      IF (JE.EQ.IE) IH=ID11+IE
      HB=0.D0
      DO 850 K=1,nvec
  850 HB=HB+ts(K,IE)*E(K)*ts(K,JE)
  860 brdmp(Ie,je)=HB
      IE=1
      E(IE)=brdmp(IE,je)
      EPE=E(IE)+ESCF
      WRITE (6,870) IE,EPE
  870 FORMAT (//5X,'MATRICE HAMILTONIENNE APRES ROTATION'/
     * I10/2X,F15.8)
      IB=0
      DO 880 IE=2,nvec
      IM=IE-1
      E(IE)=brdmp(IE,je)
      EPE=E(IE)+ESCF
      WRITE (6,875) IE
  875 FORMAT (I10)
      WRITE (6,890) (brdmp(ie,JE),JE=1,IM),EPE
  880 IB=IB+IM
  890 FORMAT (2X,8F15.8)
   66 format(10f12.8)
      deallocate(cpro)
      deallocate(wrk)
      deallocate(es)
      deallocate(ts)
      deallocate(smat)
      RETURN
      END
C***********************************************************************
      SUBROUTINE DETDIF(II,JJ,NBDIF,N1,N2,N3,N4,NS1,NS2,NS3,NS4,ZSIG)
C 
C                  SSP DE RECHERCHE DES ELEMENTS D'INTERACTION
C                  ENTRE DEUX DETERMINANTS, NUMEROTES II ET JJ
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      DIMENSION ITA(99),ITB(99),NV1(2),NV2(2),NATU1(2),NATU2(2)
      EQUIVALENCE(NV11,NV1(1)),(NV12,NV1(2)),(NV21,NV2(1)),(NV22,NV2(2))
      KR=II
      KL=JJ
      JINT=0
      KINT=0
C 
C                  LE DETERMINANT - 1 - EST CHOISI COMME LE
C                  DETERMINANT AYANT LE PLUS D'O.M. EXCITEES
C 
      IF (NE(KR)-NE(KL)) 5,7,7
    5 KKKK=KR
      KR=KL
      KL=KKKK
    7 NE1=NE(KR)
C                  TEST SUR LE FONDAMENTAL
      IF (NE1.EQ.0) GOTO 44
      NE2=NE(KL)
C                  NBDIF : EST LE NOMBRE DE SPIN-ORBITALES
C                          DIFFERENTES
      NBDIF=NE1-NE2
      IF (NBDIF.GT.2) RETURN
C                  CONSTRUCTION POUR LE DETERMINANT -1- ET -2-
C                  DE IT ET NS, QUI CONTIENNENT LES NUMEROS
C                  DES ORBITALES ET LEURS SPINS
C 
      ND1=ND(KR)
      ND2=ND(KL)
      DO 15 J=1,NE1
      ITA(J)=TROU(J+ND1)
      ITA(J+NE1)=PART(J+ND1)
   15 CONTINUE
C 
C     ZSIG :  EST LA PARITE DU  NOMBRE DE CROISEMENTS
C             DANS LE DIAGRAMME D INTERACTION
      ZSIG=.FALSE.
C                  SI NE2=0 LE DETERMINANT -2- EST LE FONDAMENTAL
      IF (NE2.EQ.0) GOTO 125
  121 CONTINUE
      DO 25 J= 1,NE2
      ITB(J)=TROU(J+ND2)
   25 ITB(J+NE2)=PART(J+ND2)
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
   12 CONTINUE
      NBDIF=NBDIF+1
      IF (NBDIF.GT.2) RETURN
      NV1(NBDIF)=NIT1
      NATU1(NBDIF)=K
      NATU2(NBDIF)=I
   10 CONTINUE
    8 CONTINUE
  125 IF (NBDIF.LE.0) GOTO 44
C 
C                   NAR : EXCITATION SUPLEMENTAIRE DU DETERMINANT -1-
C                        PAR RAPPORT AU DETERMINANT -2-
   26 NAR=NE1-NE2
      IF (NAR.LE.0) GOTO 28
   30 DO 32 I=1,NAR
      NET=NE1+1-I
      NV1(I)=ITA(NET)
      NV2(I)=ITA(NET+NE1)
   32 CONTINUE
   28 NBAR=NBDIF-NAR
      IF (NBAR.LE.0) GOTO 34
   36 NAR1=NAR+1
      DO 38 I=NAR1,NBDIF
      K=NATU1(I)
      NI=NATU2(I)
      NU=NV1(I)
      NV2(I)=ITA(NI+NE1*(K-1))
      IF (K.EQ.2) GOTO 38
   40 NV1(I)=NV2(I)
      NV2(I)=NU
      ZSIG=.NOT.ZSIG
   38 CONTINUE
   34 CONTINUE
C              ORBITALES ET SPIN
C 
C     LES -2- DETERMINANTS DIFFERENT PAR -1- SPIN-ORBITALE
      N1=IORB(NV11)
      N2=IORB(NV21)
      NS1=ISPIN(NV11)
      NS2=ISPIN(NV21)
      IF (NBDIF.EQ.1) RETURN
C 
C     LES -2- DETERMINANTS DIFFERENT PAR -2- SPIN-ORBITALES
C 
      N3=IORB(NV12)
      N4=IORB(NV22)
      NS3=ISPIN(NV12)
      NS4=ISPIN(NV22)
      RETURN
   44 WRITE (6,1006) II,JJ,NE1,NE2
 1006 FORMAT (//5X,'IDENTITE DANS HNTD'/' II,JJ,NE1,NE2 =',4I4)
      WRITE (6,1007) (ITA(K),K=1,NE1)
      WRITE (6,1008) (ITA(K+NE1),K=1,NE1)
      WRITE (6,1007) (ITB(K),K=1,NE2)
      WRITE (6,1008) (ITB(K+NE2),K=1,NE2)
 1007 FORMAT (' TROU',14I4)
 1008 FORMAT (' PART',14I4)
      STOP
      END
C***********************************************************************
      SUBROUTINE ATOM(NVEC,ZVBIS,OM)
C 
C     CALCOLO DI L2
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,ZRAP,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /INT/ alx(id2),aly(id2),alz(id2),all(id1),alz2(id2),
     * almat(id6,id6),alzmat(id6,id6),alblo(id6,id6),tblo(id6,id6),
     * tblot(id6,id6),tldia(id6,id6),almat1(id6,id6),
     * cpro(id6),epro(id6),ell(id6),elz2(id6),wrk(id6),
     * metdr(id6,id6),metdo(id6,id6),ndim(id6),netsel(id6)
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      dimension om(id13,id1)
      nintg=id2*3*8+id2*4+599*4+4000*2+id2*3*2+16
      nintl=id2*4*8+id1*8+id6*id6*7*8+id6*5*8+id6*id6*2*4+id6*2*4
c      if (nintl.gt.nintg) then
c      write (6,*) ' '
c      write (6,*)
c     * '   Attention: common INT dans le ssp. ATOM plus grand'
c      write (6,*)
c     * '   que dans le reste du programme:',nintl,nintg,'   bytes'
c      write (6,*) ' '
c      endif
      WRITE (6,5)
    5 FORMAT (///5X,'CONSTRUCTION DE LA MATRICE DE L2'/)
      call lxyz(zvbis,om)
      ref=0.d0
      refz=0.d0
      do 31 i=1,norb+1
      do 31 j=1,i
      ij=i*(i-1)/2+j
      alz2(ij)=0.d0
      if (i.gt.norb) goto 31
      do 30 k=1,norb
      zsog=.false.
      if (i.ge.k) then
      ik=i*(i-1)/2+k
      else
      zsog=.not.zsog
      ik=k*(k-1)/2+i
      endif
      if (j.ge.k) then
      jk=j*(j-1)/2+k
      else
      zsog=.not.zsog
      jk=k*(k-1)/2+j
      endif
      if (zsog) then
      alz2(ij)=alz2(ij)-alz(ik)*alz(jk)
      else
      alz2(ij)=alz2(ij)+alz(ik)*alz(jk)
      endif
   30 continue
   31 continue
      if (zprt) then
      write(6,*) ' matrice Lz2 negli orbitali '
      max=0
      do 40 i=1,norb
      min=max+1
      max=i*(i-1)/2+i
   40 write(6,66) (alz2(l),l=min,max)
      endif
      do 20 i=1,nocb
      ii=i*(i-1)/2+i
      ref=ref+all(i)*2.d0
      refz=refz+alz2(ii)*2.d0
      do 10 j=1,i
      if (j.eq.i) goto 20
      ij =i*(i-1)/2+j
      refz=refz-alz(ij)*alz(ij)*4.d0
   10 ref=ref-(alx(ij)*alx(ij)+aly(ij)*aly(ij)+alz(ij)*alz(ij))*4.d0
   20 continue
      write (6,*) ' ref= ',ref
      do 11 i=1,nvec
      do 11 j=1,i
      tldia(i,j)=0.d0
   11 tldia(j,i)=tldia(i,j)
    3 do 9 i=1,nvec
      do 9 j=1,i
      alzmat(i,j)=0.d0
      almat(i,j)=0.d0
      tblot(i,j)=0.d0
      tblot(j,i)=tblot(i,j)
      alzmat(j,i)=alzmat(i,j)
    9 almat(j,i)=almat(i,j)
      do 500 k=1,ncf
      dlld=ref 
      dlzd=refz
      do 50 l=1,ne(k)
      ipar=iorb(part(nd(k)+l))
      itro=iorb(trou(nd(k)+l))
      ll=ipar*(ipar-1)/2+ipar
      mm=itro*(itro-1)/2+itro
      dlzd=dlzd+alz2(ll)-alz2(mm)
   50 dlld=dlld+all(ipar)-all(itro)
      do 80 i=1,nocb
      do 70 l=1,ne(k)
      ipar=iorb(part(nd(k)+l))
      itro=iorb(trou(nd(k)+l))
      li=ipar*(ipar-1)/2+i
      dlld=dlld+(-alx(li)*alx(li)-aly(li)*aly(li)-alz(li)*alz(li))*2.d0
      dlzd=dlzd+(-alz(li)*alz(li))*2.d0
      if (i.le.itro) then
      mi=itro*(itro-1)/2+i
      else
      mi=i*(i-1)/2+itro
      endif
      dlld=dlld+(alx(mi)*alx(mi)+aly(mi)*aly(mi)+alz(mi)*alz(mi))*2.d0
      dlzd=dlzd+(alz(mi)*alz(mi))*2.d0
   70 continue
   80 continue
      do 100 l=1,ne(k)
      ipar=iorb(part(nd(k)+l))
      itro2=iorb(trou(nd(k)+l))
      do 90 m=1,ne(k)
      itro=iorb(trou(nd(k)+m))
      ipar2=iorb(part(nd(k)+m))
      if (ispin(part(nd(k)+l)).eq.ispin(trou(nd(k)+m))) then
      lm=ipar*(ipar-1)/2+itro 
      dlld=dlld+(alx(lm)*alx(lm)+aly(lm)*aly(lm)+alz(lm)*alz(lm))*2.d0
      dlzd=dlzd+(alz(lm)*alz(lm))*2.d0
      endif
      if (m.ge.l) goto 90
      if (ispin(part(nd(k)+l)).eq.ispin(part(nd(k)+m))) then
      if (ipar.ge.ipar2) then
      lm=ipar*(ipar-1)/2+ipar2
      else
      lm=ipar2*(ipar2-1)/2+ipar
      endif
      dlld=dlld-(alx(lm)*alx(lm)+aly(lm)*aly(lm)+alz(lm)*alz(lm))*2.d0
      dlzd=dlzd-(alz(lm)*alz(lm))*2.d0
      endif
      if (ispin(trou(nd(k)+l)).eq.ispin(trou(nd(k)+m))) then
      if (itro2.ge.itro) then
      lm=itro2*(itro2-1)/2+itro
      else
      lm=itro*(itro-1)/2+itro2
      endif
      dlld=dlld-(alx(lm)*alx(lm)+aly(lm)*aly(lm)+alz(lm)*alz(lm))*2.d0
      dlzd=dlzd-(alz(lm)*alz(lm))*2.d0
      endif
   90 continue
  100 continue
      do 150 i=1,nvec
      ci=c(ncf*i-ncf+k)
      do 150 j=1,i
      cj=c(ncf*j-ncf+k)
      alzmat(i,j)=alzmat(i,j)+ci*dlzd*cj
  150 almat(i,j)=almat(i,j)+ci*dlld*cj  
      do 500 kk=1,k-1
      dlkd=0.d0
      dlzk=0.d0
      call detdif(k,kk,nbdif,n1,n2,n3,n4,ns1,ns2,ns3,ns4,zsig)
      if (nbdif.eq.1) then
      if (ns1.eq.ns2) then
      if (n1.ge.n2) then
      lm=n1*(n1-1)/2+n2
      else
      lm=n2*(n2-1)/2+n1
      endif
      dlzk=dlzk+alz2(lm)
      endif
      do 160 i=1,nocb
      zsog=.false.
      if (n2.ge.i) then
      im=n2*(n2-1)/2+i
      else
      zsog=.not.zsog
      im=i*(i-1)/2+n2
      endif
      if(n1.le.i) then
      zsog=.not.zsog
      in=i*(i-1)/2+n1
      else
      in=n1*(n1-1)/2+i
      endif
      ctrb=-(alx(in)*alx(im)+aly(in)*aly(im)+alz(in)*alz(im))*2.d0
      ctrbz=-(alz(in)*alz(im))*2.d0
      if (zsog) then 
      ctrb=-ctrb
      ctrbz=-ctrbz
      endif
      dlkd=dlkd+ctrb
      dlzk=dlzk+ctrbz
  160 continue
      do 170 l=1,ne(k)
      ipar=iorb(part(nd(k)+l))
      itro=iorb(trou(nd(k)+l))
      if (ispin(part(nd(k)+l)).ne.ns1) goto 173
      zsog=.false.
      if (n1.ge.ipar) then
      ln=n1*(n1-1)/2+ipar
      else
      zsog=.not.zsog
      ln=ipar*(ipar-1)/2+n1
      endif
      if (ipar.ge.n2) then
      lm=ipar*(ipar-1)/2+n2
      else
      zsog=.not.zsog
      lm=n2*(n2-1)/2+ipar
      endif
      ctrb=(alx(ln)*alx(lm)+aly(ln)*aly(lm)+alz(ln)*alz(lm))*2.d0
      ctrbz=(alz(ln)*alz(lm))*2.d0
      if (zsog) then 
      ctrb=-ctrb
      ctrbz=-ctrbz
      endif
      dlkd=dlkd+ctrb
      dlzk=dlzk+ctrbz
  173 if (ispin(trou(nd(k)+l)).ne.ns1) goto 170
      zsog=.false.
      if (itro.le.n1) then
      ln=n1*(n1-1)/2+itro
      else
      zsog=.not.zsog
      ln=itro*(itro-1)/2+n1
      endif
      if (itro.le.n2) then
      lm=n2*(n2-1)/2+itro
      else
      zsog=.not.zsog
      lm=itro*(itro-1)/2+n2
      endif
      ctrb=+(alx(ln)*alx(lm)+aly(ln)*aly(lm)+alz(ln)*alz(lm))*2.d0
      ctrbz=+(alz(ln)*alz(lm))*2.d0
      if (zsog) then 
      ctrb=-ctrb
      ctrbz=-ctrbz
      endif
      dlkd=dlkd+ctrb
      dlzk=dlzk+ctrbz
  170 continue
      do 180 i=1,nvec
      ci=c(ncf*i-ncf+k)
      cii=c(ncf*i-ncf+kk)
      do 180 j=1,i
      cj=c(ncf*j-ncf+kk)
      cjj=c(ncf*j-ncf+k)
      if (zsig) then
      almat(i,j)=almat(i,j)-dlkd*(ci*cj+cii*cjj)
      alzmat(i,j)=alzmat(i,j)-dlzk*(ci*cj+cii*cjj)
      else
      almat(i,j)=almat(i,j)+dlkd*(ci*cj+cii*cjj)
      alzmat(i,j)=alzmat(i,j)+dlzk*(ci*cj+cii*cjj)
      endif
  180 continue
      endif
      if (nbdif.eq.2) then
      if (ns1.eq.ns2.and.ns3.eq.ns4) then
      zsog=.false.
      if (n2.ge.n1) then
      ln=n2*(n2-1)/2+n1
      else
      zsog=.not.zsog
      ln=n1*(n1-1)/2+n2
      endif
      if (n4.ge.n3) then
      lm=n4*(n4-1)/2+n3
      else
      zsog=.not.zsog
      lm=n3*(n3-1)/2+n4
      endif
      ctrb=-(alx(ln)*alx(lm)+aly(ln)*aly(lm)+alz(ln)*alz(lm))*2.d0
      ctrbz=-(alz(ln)*alz(lm))*2.d0
      if (zsog) then 
      ctrb=-ctrb
      ctrbz=-ctrbz
      endif
      dlkd=dlkd+ctrb
      dlzk=dlzk+ctrbz
      endif
      if (ns1.eq.ns4.and.ns3.eq.ns2) then
      zsog=.false.
      if (n4.ge.n1) then
      ln=n4*(n4-1)/2+n1
      else
      zsog=.not.zsog
      ln=n1*(n1-1)/2+n4
      endif
      if (n2.ge.n3) then
      lm=n2*(n2-1)/2+n3
      else
      zsog=.not.zsog
      lm=n3*(n3-1)/2+n2
      endif
      ctrb=+(alx(ln)*alx(lm)+aly(ln)*aly(lm)+alz(ln)*alz(lm))*2.d0
      ctrbz=+(alz(ln)*alz(lm))*2.d0
      if (zsog) then 
      ctrb=-ctrb
      ctrbz=-ctrbz
      endif
      dlkd=dlkd+ctrb
      dlzk=dlzk+ctrbz
      endif
      do 190 i=1,nvec
      ci=c(ncf*i-ncf+k)
      cii=c(ncf*i-ncf+kk)
      do 190 j=1,i
      cj=c(ncf*j-ncf+kk)
      cjj=c(ncf*j-ncf+k)
      if (zsig) then
      almat(i,j)=almat(i,j)-dlkd*(ci*cj+cii*cjj)
      alzmat(i,j)=alzmat(i,j)-dlzk*(ci*cj+cii*cjj)
      else
      almat(i,j)=almat(i,j)+dlkd*(ci*cj+cii*cjj)
      alzmat(i,j)=alzmat(i,j)+dlzk*(ci*cj+cii*cjj)
      endif
  190 continue
      endif
  500 continue
      write(6,*) ' matrice L2 '
      do 600 i=1,nvec
  600 write(6,66) (almat(i,j),j=1,i) 
   66 format(10f12.8)
      write(6,*) ' matrice Lz2 '
      do 625 i=1,nvec
  625 write(6,66) (alzmat(i,j),j=1,i) 
      if (nvec.eq.1) return
      if (zvbis) return
      do 694 i=1,nvec
      do 694 j=1,i
      almat1(i,j)=almat(i,j)
  694 almat1(j,i)=almat(i,j)
      do 693 i=1,nvec
      do 693 j=1,i-1
      if (dabs(almat1(i,j)).le.1.d-8) almat1(i,j)=0.d0
  693 almat1(j,i)=almat1(i,j)
c      call riera(nvec,almat1,tldia,ell,wrk,id6)
      do i=1,nvec
         do j=1,nvec
            tldia(i,j)=almat1(i,j)
         enddo
      enddo
      call dsyev('V','U',nvec,tldia,id6,ell,wrk,3*id6,info) !mind dim wrk
      write(6,8) 
    8 format (//, ' DIAGONALIZZAZIONE DI L2 '/)
      write(6,*) ' matrice che diagonalizza L2 '
      do 698 k=1,nvec
  698 write(6,66) (tldia(k,kk),kk=1,nvec)
      write(6,77)
   77 FORMAT (//,'autovalori di L2 ')
      write(6,66) (ell(i),i=1,nvec)
      nrd=0
      tdegen1=tdegen*2.d2
      do 299 i=1,nvec
  299 netsel(i)=0
      do 300 i=1,nvec-1
      if (netsel(i).eq.i) goto 300
      metdr(i,1)=i
      netd=1
      do 310 j=i+1,nvec
      aldeg=dabs(almat(i,i)-almat(j,j))
      edeg=dabs(e(i)-e(j))
      if (edeg.lt.tdegen.and.aldeg.lt.tdegen1) then 
      netd=netd+1
      if (netd.eq.2) netsel(i)=i
      netsel(j)=j
      metdr(i,netd)=j
      endif
  310 continue
      if (netd.eq.1) goto 300
      nrd=nrd+1
c     write(6,66) (e(k),k=1,nvec)
      ndim(nrd)=netd
      do 312 l=1,netd
  312 metdo(nrd,l)=metdr(i,l)
  300 continue
      if (nrd.eq.0) then
      write(6,55) tdegen
   55 format(/' non trovati sottospazi degeneri in E ed L2 '/
     * ' tdegen= ',f12.8/)
      endif
      if (nrd.eq.0) return
      if (zprt) then
      write(6,*) ' numero di sottospazi degeneri in E ed L2 ',nrd
      write(6,*) ' suddivisione in blocchi '
      do 351 i=1,nrd
  351 write(6,*) (metdo(i,j),j=1,ndim(i))
      endif
      write(6,6666)
 6666 format(/,' AUTOVALORI DI Lz2 ')
      do 350 i=1,nrd
      do 345 j=1,nvec
      do 345 l=1,j   
      alblo(j,l)=0.d0
      tblo(j,l)=0.d0
      alblo(l,j)=alblo(j,l)
  345 tblo(l,j)=tblo(j,l)
      netd=ndim(i)
      do 340 j=1,netd
      do 340 l=1,j
      alblo(j,l)=alzmat(metdo(i,j),metdo(i,l))
  340 alblo(l,j)=alblo(j,l)
      if (zprt) then
      write(6,*) ' blocco di Lz2 da diagonalizzare '
      do 332 k=1,netd
  332 write(6,66) (alblo(k,kk),kk=1,netd)
      endif
c      call riera(netd,alblo,tblo,elz2,wrk,id6)
      do icp=1,netd
         do j=1,netd
            tblo(icp,j)=alblo(icp,j)
         enddo
      enddo
      call dsyev('V','U',netd,tblo,id6,elz2,wrk,3*id6,info) !mind dim wrk
      write(6,666) (metdo(i,k),k=1,netd)
  666 format(/,' stati ',/15I3)
      write(6,66) (elz2(k),k=1,netd)
      if (zprt) then
      write(6,*) ' matrice di trasformazione '
      do 331 k=1,netd
  331 write(6,66) (tblo(k,kk),kk=1,netd)
      endif
      do 330 j=1,netd
      do 330 l=1,netd
  330 tblot(metdo(i,j),metdo(i,l))=tblo(j,l)
  350 continue 
      do 355 i=1,nvec
  355 if (netsel(i).ne.i) tblot(i,i)=1.d0
      if (zprt) then 
      write(6,665)
  665 format(/,'matrice di trasformazione per la diagonalizzazione ',
     * ' a blocchi di Lz2 ')
      do 360 i=1,nvec
  360 write(6,66) (tblot(i,j),j=1,nvec)
      endif
      if (zdialz) then
      write (6,850)
  850 format(/5X,'rotazione degli stati in modo che lz2 sia 
     * diagonale a blocchi'/)
      zvbis=.true.
      do 400 i=1,ncf
      do 390 k=1,nvec
      cpro(k)=0.d0
      do 390 l=1,nvec
  390 cpro(k)=c(ncf*l-ncf+i)*tblot(l,k)+cpro(k)
      do 395 k=1,nvec
  395 c(ncf*k-ncf+i)=cpro(k)
  400 continue
      endif
      RETURN
      END
C***********************************************************************
      SUBROUTINE LXYZ(ZVBIS,OM)
C 
C     CALCOLO DI LX, LY, LZ, L2 SU BASE MOLECOLARE
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      character*400 namel
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,ZRAP,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      COMMON /INT/ alx(id2),aly(id2),alz(id2),all(id1),
     * wp(id1,10),wd(id1,10)
      dimension ip(10),id(10),om(id13,id1)
      nintg=id2*3*8+id2*4+599*4+4000*2+id2*3*2+16
      nintl=id2*3*8+id1*8+id1*10*2*8
c      if (nintl.gt.nintg) then
c      write (6,*) ' '
c      write (6,*)
c     * '   Attention: common INT dans le ssp. LXYZ plus grand'
c      write (6,*)
c     * '   que dans le reste du programme:',nintl,nintg,'   bytes'
c      write (6,*) ' '
c      endif
      if (.not.zvbis) then 
      do 10 i=1,10
      ip(i)=0
   10 id(i)=0
      attol=1.d-6
      rewind 3
      call cerca(3,'ATOMO',itrovo)
      if (itrovo.eq.0) then
      write (6,*) '   *** NAMELIST &ATOMO NOT FOUND ***'
      stop 12
      endif
      call lecnam('&ATOMO',namel,3)
      call nameli(namel,ip,'IP',4)
      call nameli(namel,id,'ID',4)
      call nameld(namel,attol,'ATTOL',8)
      endif
      rewind 11
      READ (11) ((DUM,I=1,NAO),J=1,NGELO),
     * ((om(I,J),I=1,NAO),J=1,NORB)
      nump=0
      numd=0
      do 30 i=1,10
      if (ip(i).gt.0) nump=nump+1
   30 if (id(i).gt.0) numd=numd+1
      cost=2.d0/dsqrt(3.d0)
      do 500 i=1,norb+1
      all(i)=0.d0
      ij=(i*i-i)/2
      do 40 j=1,i
      ij=ij+1
      alx(ij)=0.d0
      aly(ij)=0.d0
   40 alz(ij)=0.d0
      if (i.gt.norb) goto 500
      zp=.false.
      amx=0.d0
      do 50 k=1,nump
      c1=om(ip(k),i)
      c2=om(ip(k)+1,i)
      c3=om(ip(k)+2,i)
      wp(i,k)=c1*c1+c2*c2+c3*c3
      if (wp(i,k).gt.amx) then 
      amx=wp(i,k)
      kmaxp=k
      endif
   50 if (wp(i,k).ge.attol) zp=.true.
      if (zp.and.zprt) write (6,*) ' ZP=T, I, AMX =',i,amx
      zd=.false.
      amx=0.d0
      do 60 k=1,numd
      c1=om(id(k),i)
      c2=om(id(k)+1,i)
      c3=om(id(k)+2,i)
      c4=om(id(k)+3,i)
      c5=om(id(k)+4,i)
      c6=om(id(k)+5,i)
      wd(i,k)=c1*c1+c2*c2+c3*c3+c4*c4+c5*c5+c6*c6+
     * 2.d0*(c1*c2+c2*c3+c1*c3)/3.d0  
      if (wd(i,k).gt.amx) then 
      amx=wd(i,k)
      kmaxd=k
      endif
      if (wd(i,k).ge.attol) then
      cdif1=c1-c2
      cdif2=c1-c3
      cdif3=c2-c3
      cdif=dabs(cdif1)+dabs(cdif2)+dabs(cdif3)
      cxyz=c4+c5+c6
      if (abs(cxyz).lt.attol.and.cdif.lt.attol) goto 60
      zd=.true.
      endif
   60 continue
      if (zd.and.zprt) write (6,*) ' ZD=T, I, AMX =',i,amx
      if (zp) all(i)=2.d0
      if (zd) all(i)=6.d0
      do 400 j=1,i  
      if (j.eq.i) goto 400
      if (zp) then
      ziden=.true.
      do 70 k=1,nump
   70 if (dabs(wp(i,k)-wp(j,k)).ge.attol) ziden=.false.
      if (ziden) then
      ij=i*(i-1)/2+j
      alx(ij)=om(ip(kmaxp)+2,i)*om(ip(kmaxp)+1,j)-
     *        om(ip(kmaxp)+1,i)*om(ip(kmaxp)+2,j)
      alx(ij)=alx(ij)/wp(i,kmaxp)
      aly(ij)=om(ip(kmaxp),i)*om(ip(kmaxp)+2,j)-
     *        om(ip(kmaxp)+2,i)*om(ip(kmaxp),j)
      aly(ij)=aly(ij)/wp(i,kmaxp)
      alz(ij)=om(ip(kmaxp)+1,i)*om(ip(kmaxp),j)-
     *        om(ip(kmaxp),i)*om(ip(kmaxp)+1,j)
      alz(ij)=alz(ij)/wp(i,kmaxp)
      endif
      endif
      if (zd) then
      ziden=.true.
      do 80 k=1,numd
   80 if (dabs(wd(i,k)-wd(j,k)).ge.attol) ziden=.false.
      if (ziden) then
      ij=i*(i-1)/2+j
      cix2=om(id(kmaxd),i)
      ciy2=om(id(kmaxd)+1,i)
      ciz2=om(id(kmaxd)+2,i)
      cixy=om(id(kmaxd)+3,i)
      cixz=om(id(kmaxd)+4,i)
      ciyz=om(id(kmaxd)+5,i)
      cjx2=om(id(kmaxd),j)
      cjy2=om(id(kmaxd)+1,j)
      cjz2=om(id(kmaxd)+2,j)
      cjxy=om(id(kmaxd)+3,j)
      cjxz=om(id(kmaxd)+4,j)
      cjyz=om(id(kmaxd)+5,j)
      alx(ij)=cost*ciyz*cjy2-cost*ciyz*cjz2+cixz*cjxy-cixy*cjxz
     * +cost*(ciz2-ciy2)*cjyz
      alx(ij)=alx(ij)/wd(i,kmaxd)
      aly(ij)=-cost*cixz*cjx2+cost*cixz*cjz2-ciyz*cjxy
     * +cost*(cix2-ciz2)*cjxz+cixy*cjyz
      aly(ij)=aly(ij)/wd(i,kmaxd)
      alz(ij)=cost*cixy*cjx2-cost*cixy*cjy2+cost*(ciy2-cix2)*cjxy
     * +ciyz*cjxz-cixz*cjyz
      alz(ij)=alz(ij)/wd(i,kmaxd)
      endif
      endif
  400 continue
  500 continue
      if (.not.zprt) return
      write(6,*) ' vettore l2 '
      write(6,'(1p,4d20.12)') (all(i),i=1,norb)
      do 600 k=1,3
      if (k.eq.1) write(6,*) ' matrice -i*lx '
      if (k.eq.2) write(6,*) ' matrice -i*ly '
      if (k.eq.3) write(6,*) ' matrice -i*lz '
      max=0
      do 600 i=1,norb
      min=max+1
      max=i*(i+1)/2
      if (k.eq.1) write(6,13) (alx(ij),ij=min,max)
      if (k.eq.2) write(6,13) (aly(ij),ij=min,max)
      if (k.eq.3) write(6,13) (alz(ij),ij=min,max)
  600 continue
   13 format(10f12.8)
      RETURN
      END
C***********************************************************************
      SUBROUTINE PREPER(zdet,idimorb3,ZEF,ITM,ITK)
C 
C     PREPARATION DU CALCUL PERTURBATIF
C     CONSTRUCTION DU TABLEAU D'INDICES (MUL) GERANT LES REGLES DE
C     SOMMATION DES CONTRIBUTIONS DEGENEREES DANS LA QDPT
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer hdiag,stok,ikk,zocd,zefsp
      dimension hdiag(:),stok(:),ikk(:),zocd(:,:)
     $     ,zefsp(:)
      COMMON /PER/ HDIAG,STOK,IKK,NTI(24),
     * MAXMUL,ZOCD,ZEFSP
      INTEGER*2 IKK,NTI,NTL(12),NPL(12)
      EQUIVALENCE (NTL(1),NTI(1)),(NPL(1),NTI(13))
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      DIMENSION ZEF(ID5,20)
      DIMENSION ZEF(ncfg,20)
      INTEGER*2 ITM(ID6,20),ITK(ncfg,10)
      integer*2 nt(13),npa(13),mp(13),mt(13),itab(26)
C 
C  EMP(M):ENERGIE MP-B D'ORDRE ZERO POUR L'ETAT M
C  LE SEUIL REEL DE SELECTION DES DETERMINANTS POUR CHAQUE ETAT
C  EST AJUSTE SUR LE PLUS GRAND COEFFICIENT DE LA FONCTION
C  MULTICONF DE CET ETAT
C 
c      allocate(hdiag(ncfg))
c      allocate(stok(ncf))
      allocate(ikk(ncf))
cele 14-12-2007      allocate(zocd(idimorb3,ncfg))
      allocate(zefsp(ncfg))
      WRITE (6,1)
    1 FORMAT (/5X,'PREPARATION DU CALCUL PERTURBATIF'/)
c
c     Cele: si costruiscono a questo punto le matrici ZEF, ZOCD
c     ZOCS e ZFREE
c
      MNOR1=NORB2+1
      DO 2 K=1,ncfg
      DO 3 I=1,20
  3   ZEF(K,I)=.FALSE.
      DO 2 I=1,MNOR1
cele 14-12-2007  2    ZOCD(I,K)=.FALSE.
 2    continue
      do i=1,mnor1
      zocs(i)=.false.
      enddo
      do i=1,norb
      zfree(i)=.true.
      enddo
cele 14-12-2007 
c      DO 29 KG=1,NCFG
c      K=KGEN(KG)
c      NEC=NE(K)
c      NDK=ND(K)
c      IF (NEC.EQ.0) GOTO 29
c      DO 28 J=1,NEC
c      NTKJ=TROU(NDK+J)
c      NPKJ=PART(NDK+J)
c      ZOCD(NTKJ,KG)=.TRUE.
c  28  ZOCD(NPKJ,KG)=.TRUE.
c 29   CONTINUE
cele 14-12-2007 
      do i=1,nd(ncf)+ne(ncf)
      zocs(trou(i))=.true.
      zocs(part(i))=.true.
      NTO=trou(i)
      IF (NTO.GT.NORB) NTO=NTO-NORB
      ZFREE(NTO)=.FALSE.
      NPO=part(i)
      IF (NPO.GT.NORB) NPO=NPO-NORB
      ZFREE(NPO)=.FALSE.
      enddo
c
      znext=.false.
      kg=0
1110  kg=kg+1
      ikg=0
      K=KGEN(KG)
      IF (ISZ.NE.0.or.zion) ZEF(kg,1)=.TRUE.
      IF (ISZ.NE.0.or.zion) GOTO 8
      if (ne(kgen(kg)).ne.ne(kgen(kg+1))) goto 8
      nec=ne(k)
      do l=1,nec
      nt(l)=trou(nd(kgen(kg))+l)
      npa(l)=part(nd(kgen(kg))+l)
      mt(l)=nt(l)
      mp(l)=npa(l)
      enddo
      do l=1,nec
      if (nt(l).gt.norb) then
      itab(l)=nt(l)-norb
      else
      itab(l)=nt(l)+norb
      endif
      enddo
      do l=1,nec
      if (npa(l).gt.norb) then
      itab(l+nec)=npa(l)-norb
      else
      itab(l+nec)=npa(l)+norb
      endif
      enddo
      call ordom(nec,itab)
      do l=1,nec
      if (nt(l).ne.itab(l)) goto 665
      if (npa(l).ne.itab(l+nec)) goto 665
      enddo
      goto 8
  665 do l=1,nec
      if (trou(nd(kgen(kg+1))+l).ne.itab(l)) goto 666
      if (part(nd(kgen(kg+1))+l).ne.itab(l+nec)) goto 666
      enddo
      ZEF(kg,1)=.TRUE.
      ZEF(kg+1,1)=.TRUE.
      ikg=ikg+1
      goto 8
  666 write (6,*) 'Errore in preper: il flippato del determinante',kg
      write (6,*) 'non e` un generatore'
      print '(15i4)',(trou(nd(kgen(kg))+l),l=1,nec)
      print '(15i4)',(trou(nd(kgen(kg+1))+l),l=1,nec)
      print '(15i4)',(part(nd(kgen(kg))+l),l=1,nec)
      print '(15i4)',(part(nd(kgen(kg+1))+l),l=1,nec)
      stop 12
    8 kg=kg+ikg
      ikg=0
cele 27.01.2004      if (kg.lt.ncfg) goto 1110
      if (kg.lt.ncfg-1) goto 1110
c
      MNOR1=NORB2+1
      NCF2=NCF*2
c      REWIND 20
c      READ(20)(stok(IIII),IIII=1,NCF)
      DO 5 KG=1,NCFG
    5 HDIAG(KG)=stok(KGEN(KG))
      do 6 ish=1,nshift
    6 shift(ish)=shift(ish)*factor
c   energie di ordine zero, EEN e EMP
      if (zenb) then
         call bari(metat,ncf,stok,c,een)
      else
         do m=1,metat
         een(m)=e(m)
         enddo
      endif
      DO 15 K=1,NCF
      CALL HMP(K,AC)
   15 STOK(K)=AC
      call bari(metat,ncf,stok,c,emp)
c     
      REWIND 4
      READ (4)
      if (zdet.and..not.zpertu) then
      READ (4) ii,(idum,KG=1,ii),ISYM
      else
      READ (4) NCFG,(KGEN(KG),KG=1,NCFG),ISYM
      endif
      DO 42 KG=1,NCFG
      K=KGEN(KG)
   42 ZEFSP(KG)=ZEF(KG,1)
      RETURN
      END
C***********************************************************************
      SUBROUTINE BARI(metat,ncf,stok,c,eb)
C 
C     Calcolo delle energie baricentriche
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      dimension EB(*),C(*),STOK(*)
      DO 30 M=1,METAT
      EBM=0.D0
      JM=(m-1)*ncf
      DO 20 K=1,NCF
      MK=JM+K
      CMK=C(MK)
      CC=CMK*CMK
   20 EBM=EBM+CC*STOK(K)
   30 EB(M)=EBM
      return
      end
C***********************************************************************
      SUBROUTINE MKF04(zdet,lgen,lgen2)
C 
C     Preparazione di alcune variabili da scrivere sul file04 
C     e scrittura del file04.
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      COMMON /LSH/ SHIFT(12),E2MPSH(ID6,12),E2ENSH(ID6,12),
     * SHMP(ID6,24),SHEN(ID6,24),PSIMPSH(ID6,12),PSIENSH(ID6,12),NSHIFT
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      INTEGER*2 IKK,NTI,NTL(12),NPL(12)
c      COMMON /INT/ ZOS(ID4),NESP(ID4)
      allocatable nesp(:)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
c      DIMENSION ZEF(ID5,20)
      INTEGER*2 NESP
      integer*2 nt(13),npa(13),mp(13),mt(13),itab(26)
c      dimension lgen(id4)
      dimension lgen(*)
C 
C     RESULTATS VARIATIONNELS ET PREPARATOIRES DE LA PERTURBATION
C     SUR LE FICHIER 4
C 
      allocate(nesp(ncf))
  800 I=NCF*NROT
c     II=ND(NCF)+NE(NCF)+1
      II=ND(NCF)+NE(NCF)
      III=NORB+NORB
      NNROT=NROT*(NROT+1)/2
      MMETAT=METAT*(METAT-1)/2
      if (.not.zpertu.and.zdet) then
       rewind 4
       read (4)
       READ (4) NCFG,(lGEN(KG),KG=1,NCFG),ISYM,(NESP(K),K=1,NCF)
      endif
      REWIND 4
      WRITE (4) NORB,NOCB,NOCA,METAT,NCF,II,I,(NE(J),ND(J),J=1,NCF),
     * (TROU(J),PART(J),J=1,II),(C(J),J=1,I),(E(J),J=1,METAT),
     * (EMP(J),J=1,METAT),(FMPB(J),J=1,III),
     *     NNROT,(brdum,j=1,nnrot),MMETAT,(ZHEFF(J),J=1
     $     ,MMETAT),zhandy,zthree,FACTOR,escf,zenb,(een(j),j=1,metat)
      if (.not.zpertu.and.zdet) then
      WRITE (4) NCFG,(lgen(KG),KG=1,NCFG),ISYM,(NESP(K),K=1,NCF)
      else
      WRITE (4) NCFG,(KGEN(KG),KG=1,NCFG),ISYM,(NESP(K),K=1,NCF)
      endif
      deallocate(nesp)
      RETURN
      END
c-------------------------------------------------------------------------------
      SUBROUTINE HMAT(nconf,dir,iocc,idimocc)
C 
C     CONSTRUCTION DE LA MATRICE HAMILTONIENNE POUR L'IC VARIATIONNELLE
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer hdiag(:),stok(:)
      COMMON /PER/ HDIAG,STOK
      character*64 dir,filein
c      DIMENSION HOFF(5000),HIC(ID4),ZERO(ID4),ZOS(ID4)
c      dimension hic(id4)
      allocatable hic(:)
      COMMON /CPU/ TCPU,PCPU
#ifdef SINGLET      
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hma(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hma(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hma(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
#endif
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      dimension ibuf(100000),jbuf(100000),buffer(100000)
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      integer*2 icomp,iconf,isegno
      integer*1 iocc(idimocc,*)
crenzo
      integer*8 ilong,ij2long,ijlong,ij1long
crenzo     
      allocate(hic(ncf))
      allocate(hdiag(ncfg))
      allocate(stok(ncf))
c--creazione matrici T
      call cput(0)
      print*,'Calling kotani...'
      T0(1,1)=1.D0
#ifdef SINGLET
      call kotani(0)
#endif
#ifdef DOUBLET
      call kotani(1)
#endif
#ifdef TRIPLET
      call kotani(2)
#endif
      call cput(0)
      print '(a,f5.2,a)','Kotani has taken ',pcpu,' seconds'
c--matrices have been created
c
c      call esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,nocb,norb)
      if(.not.zhmat)then
         inext=1
         iconfi=1
         do i=1,ncf
            if(i.eq.inext)then
               inext=i+icomp(i)
               if(iconf(i).eq.0)then
                  iconff=iconfi
               else
                  iconff=iconfi+iconf(i)-1
               endif
               iconfi=iconff+1
            endif
         enddo
         nconf=iconff
         print*,'There are ',nconf,' configurations'
         deallocate(hic)
         return
      endif
      ZFA=.FALSE.
      WRITE (6,5)
    5 FORMAT (///5X,'CONSTRUCTION DE LA MATRICE DE I.C.'/)
c--open scratch file 21 to be deleted later on-------
      open(21,form='UNFORMATTED',status='UNKNOWN',file='FILE21')
      rewind 21
      inext=1
      indbuf=0
      iconfi=1
      ivolte=0
      DO 15 IIII=1,NCF
      I=IIII
      if(i.eq.inext)then
         call hdigc(i,icomp(i),iconf(i),finv,isegno)
         zcapos(i)=.true.
         inext=i+icomp(i)
         if(iconf(i).eq.0)then
            iconff=iconfi
         else
            iconff=iconfi+iconf(i)-1
         endif
         call putbuf(iconfi,iconff,cmat,ibuf,jbuf,buffer,indbuf,hic
     $        ,100000,ivolte)
         iconfi=iconff+1
      else
         zcapos(i)=.false.
         goto 15
      endif
c      CALL HDIG(I,ZFA,DUM,AB)
c      HIC(Ico)=AB*FINV
 15   continue
      if(indbuf.ne.0)write (21)ibuf,jbuf,buffer
      nconf=iconff
      print*,'There are ',nconf,' configurations'
      REWIND 20
      WRITE(20)(HIC(IIII),IIII=1,NCONF)
      rewind 21
      print*,'ivolte,indbuf=',ivolte,indbuf
      do iii=1,ivolte
         read (21)ibuf,jbuf,buffer
         write (20)ibuf,jbuf,buffer
      enddo
      if(indbuf.ne.0) read (21)ibuf,jbuf,buffer
      close (21,status='DELETE')
      mnext=1
      iconfi=1
      mcapos=0
      DO 300 m=1,NCF
         if (mod(m,1000).eq.0) then
           print '(a,I9)','indice=',m
 	   call flush(6)
         endif
         if(m.ne.mnext)goto 300
         if(iconf(m).eq.0)then
            iconff=iconfi
         else
            iconff=iconfi+iconf(m)-1
         endif
         mcapos=mcapos+1
         mnext=m+icomp(m)
c         call giveocc(iocc,m,nd,ne,trou,part,nocb,norb)
         nnext=1
         jconfi=1
         ncapos=0
         do 301 n=1,m-1
            if(n.eq.nnext)then
               if(iconf(n).eq.0)then
                  jconff=jconfi
               else
                  jconff=jconfi+iconf(n)-1
               endif
               ncapos=ncapos+1
               nedif=ne(m)-ne(n)
               if(nedif.gt.2)goto 302
c               call giveocc(jocc,n,nd,ne,trou,part,nocb,norb)
c               if(ndiff(iocc,jocc,norb).gt.4)goto 302
               if(ndiff(iocc(1,mcapos),iocc(1,ncapos),norb).gt.4)goto
     $              302
               call hntdc(m,n,icomp(m),icomp(n),iconf(m),iconf(n),finv
     $              ,isegno)
               call putbuf2(iconfi,jconfi,iconff,jconff,cmat,ibuf,jbuf
     $              ,buffer,indbuf,100000)
 302           jconfi=jconff+1
               nnext=n+icomp(n)
            else
               goto 301
            endif
c      HIC(n)=HNTD(m,n,JINT,KINT)*FINV  !da pensare
C     STOCKAGE DE LA MATRICE DE IC DANS LE FICHIER 20
C 
 301     continue
         iconfi=iconff+1
 300  CONTINUE
      indbuf=indbuf+1
      ibuf(100000)=-1
      jbuf(100000)=indbuf-1
      write (20)ibuf,jbuf,buffer
      deallocate(hic)
      RETURN
      END
c--------------------------------------------------------
      subroutine esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,nocc
     $     ,norb)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      integer*2 ne,trou,part,icomp,iconf,isegno
      dimension nd(*),ne(*),trou(*),part(*),icomp(*),iconf(*),isegno(*)
      equivalence (ip1,ib1),(ip2,ib2)
#ifdef SINGLET
      nos0=0
      nos2=0
      nos4=0
      nos6=0
      nos8=0
      nos10=0
      nos12=0
#endif
#ifdef DOUBLET
      nos1=0
      nos3=0
      nos5=0
      nos7=0
      nos9=0
      nos11=0
      nos12=0
#endif
#ifdef TRIPLET
      nos2=0
      nos4=0
      nos6=0
      nos8=0
      nos10=0
      nos12=0
#endif
      do m=1,ncf
         nec=ne(m)
         nstart=nd(m)
         iopenh=nec
         iopenp=nec
         do i=1,nec
            ib1=trou(nstart+i)
            if(ib1.gt.norb)ib1=ib1-norb
            do j=i+1,nec
               ib2=trou(nstart+j)
               if(ib2.gt.norb)ib2=ib2-norb
               if(ib2.eq.ib1)then
                  iopenh=iopenh-2
                  goto 1
               endif
            enddo
 1          continue
         enddo
         do i=1,nec
            ip1=part(nstart+i)
            if(ip1.gt.norb)ip1=ip1-norb
            do j=i+1,nec
               ip2=part(nstart+j)
               if(ip2.gt.norb)ip2=ip2-norb
               if(ip2.eq.ip1)then
                  iopenp=iopenp-2
                  goto 2
               endif
            enddo
 2          continue
         enddo
#ifdef SINGLET
         ntot=iopenh+iopenp
#endif
#ifdef DOUBLET
         ntot=iopenh+iopenp-1
#endif
#ifdef TRIPLET
         ntot=iopenh+iopenp
#endif
#ifdef SINGLET
         if(ntot.eq.0)then
            icomp(m)=1
            iconf(m)=0   !diventa 1 nei calcoli
            nos0=nos0+1
         elseif(ntot.eq.2)then
            icomp(m)=2
            iconf(m)=1
            nos2=nos2+1
         elseif(ntot.eq.4)then
            icomp(m)=6
            iconf(m)=2
            nos4=nos4+1
         elseif(ntot.eq.6)then
            icomp(m)=20
            iconf(m)=5
            nos6=nos6+1
         elseif(ntot.eq.8)then
            icomp(m)=70
            iconf(m)=14
            nos8=nos8+1
         elseif(ntot.eq.10)then
            icomp(m)=252
            iconf(m)=42
            nos10=nos10+1
         elseif(ntot.eq.12)then
            icomp(m)=924
            iconf(m)=132
            nos12=nos12+1
         elseif(ntot.gt.12)then
            print*,'sorry, more than 12 electrons not covered yet',ntot
            stop 'too many open shells'
         else
            print*,'something must be really very wrong',ntot,iopenh
     $           ,iopenp
            print*,'determinant m=',m
            print*,'trou'
            print*,(trou(nstart+i),i=1,nec)
            print*,'part'
            print*,(part(nstart+i),i=1,nec)
            stop 'very strange'
         endif
#endif
#ifdef DOUBLET
         if(ntot.eq.1)then
            icomp(m)=1
            iconf(m)=1
            nos1=nos1+1
         elseif(ntot.eq.3)then
            icomp(m)=3
            iconf(m)=2
            nos3=nos3+1
         elseif(ntot.eq.5)then
            icomp(m)=10
            iconf(m)=5
            nos5=nos5+1
         elseif(ntot.eq.7)then
            icomp(m)=35
            iconf(m)=14
            nos7=nos7+1
         elseif(ntot.eq.9)then
            icomp(m)=126
            iconf(m)=42
            nos9=nos9+1
         elseif(ntot.eq.11)then
            icomp(m)=462
            iconf(m)=132
            nos11=nos11+1
         elseif(ntot.gt.11)then
            print*,'sorry, more than 11 electrons not covered yet',ntot
            stop 'too many open shells'
         else
            print*,'something must be really very wrong',ntot,iopenh
     $           ,iopenp
            print*,'determinant m=',m
            print*,'trou'
            print*,(trou(nstart+i),i=1,nec)
            print*,'part'
            print*,(part(nstart+i),i=1,nec)
            stop 'very strange'
         endif
c         print '(a,i4,a,i2)','det. ',m,' isegno ',isegno(m)
#endif
#ifdef TRIPLET
         if(ntot.eq.2)then
            icomp(m)=1
            iconf(m)=1
            nos2=nos2+1
         elseif(ntot.eq.4)then
            icomp(m)=4
            iconf(m)=3
            nos4=nos4+1
         elseif(ntot.eq.6)then
            icomp(m)=15
            iconf(m)=9
            nos6=nos6+1
         elseif(ntot.eq.8)then
            icomp(m)=56
            iconf(m)=28
            nos8=nos8+1
         elseif(ntot.eq.10)then
            icomp(m)=210
            iconf(m)=90
            nos10=nos10+1
         elseif(ntot.eq.12)then
            icomp(m)=792
            iconf(m)=297
            nos12=nos12+1
         elseif(ntot.gt.12)then
            print*,'sorry, more than 12 electrons not covered yet',ntot
            stop 'too many open shells'
         else
            print*,'something must be really very wrong',ntot,iopenh
     $           ,iopenp
            print*,'determinant m=',m
            print*,'trou'
            print*,(trou(nstart+i),i=1,nec)
            print*,'part'
            print*,(part(nstart+i),i=1,nec)
            stop 'very strange'
         endif
#endif
         call mkorb(m,nd,ne,trou,part,nocc,norb,isegno(m))
      enddo
c--Lies, damned lies and statistics
#ifdef SINGLET
      print*,'Configuration statistics:'
      print '(a,i6,a)','# of  0-os ',nos0,' (1x1      det-conf)'
      print '(a,i6,a)','# of  2-os ',nos2/2,' (2x1      det-conf)'
      print '(a,i6,a)','# of  4-os ',nos4/6,' (6x2      det-conf)'
      print '(a,i6,a)','# of  6-os ',nos6/20,' (20x5     det-conf)'
      print '(a,i6,a)','# of  8-os ',nos8/70,' (70x14    det-conf)'
      print '(a,i6,a)','# of 10-os ',nos10/252,' (253x42   det-conf)'
      print '(a,i6,a)','# of 12-os ',nos12/924,' (924x132   det-conf)'
      call flush(6)   !tolgo flush
#endif
#ifdef DOUBLET
      print*,'Configuration statistics:'
      print '(a,i6,a)','# of  1-os ',nos1,' (1x1      det-conf)'
      print '(a,i6,a)','# of  3-os ',nos3/3,' (3x2      det-conf)'
      print '(a,i6,a)','# of  5-os ',nos5/10,' (10x5      det-conf)'
      print '(a,i6,a)','# of  7-os ',nos7/35,' (35x14     det-conf)'
      print '(a,i6,a)','# of  9-os ',nos9/126,' (126x42    det-conf)'
      print '(a,i6,a)','# of 11-os ',nos11/462,' (462x132   det-conf)'
      call flush(6)   !tolgo flush
#endif
#ifdef TRIPLET
      print*,'Configuration statistics:'
      print '(a,i6,a)','# of  2-os ',nos2,' (1x1      det-conf)'
      print '(a,i6,a)','# of  4-os ',nos4/4,' (4x3      det-conf)'
      print '(a,i6,a)','# of  6-os ',nos6/15,' (15x9     det-conf)'
      print '(a,i6,a)','# of  8-os ',nos8/56,' (56x28    det-conf)'
      print '(a,i6,a)','# of 10-os ',nos10/210,' (210x90   det-conf)'
      print '(a,i6,a)','# of 12-os ',nos12/792,' (792x297   det-conf)'
      call flush(6)   !tolgo flush
#endif
      return
      end
c----------------------------------------------------------
      subroutine giveocc(iocc,m,nd,ne,trou,part,nocc,norb)
c      integer*2 iocc(*),ne(*),trou(*),part(*)
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
c      integer*2 iocc(*),jocc(*)
      integer*1 iocc(*),jocc(*)
      ndiff=0
      do i=1,norb
c         ia=abs(iocc(i)-jocc(i))
         ia=iocc(i)-jocc(i)
         if(ia.ne.0)ndiff=ndiff+abs(ia)
c     Da testare la linea qui sotto per la velocita`
c         if(iocc(i).ne.jocc(i))ndiff=ndiff+abs(iocc(i)-jocc(i)
         if(ndiff.gt.4)return
      enddo
      return
      end
c------------------------------------------------------------------
      logical*1 function znequ(iocc,jocc,norb)
c      integer*2 iocc(*),jocc(*)
      integer*1 iocc(*),jocc(*)
      do i=1,norb
         if(iocc(i).ne.jocc(i))then
	 znequ=.false.
	 return
	 endif
      enddo
      znequ=.true.
      return
      end
c----------------------------------------------------------      
#ifdef SINGLET
      SUBROUTINE mkorb(m,nd,ne,trou,part,nocc,norb,isign)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer orba(500),orbb(500),nd(*)
      integer*2 ne(*),trou(*),part(*),isign
      dimension io1(500),io2(500),ipos(500),jpos(500)
      equivalence (io1(1),ipos(1)),(io2(1),jpos(1))
      logical*1 zcli(500),zclj(500)
c---vacuum state creation-----------
      if (nocc.gt.500) then
      write (6,*) 'Too many electrons',nocc
      stop
      endif
      do i=1,nocc
         orbb(i)=i
         orba(i)=i+norb
         zcli(i)=.false.
         zclj(i)=.false.
      enddo
      nstart=nd(m)
      nec=ne(m)
      if(nec.eq.0)goto 21
      do  k=1,nec
         ib=trou(nstart+k)
         ip=part(nstart+k)
         if(ib.le.norb)then
            orbb(ib)=ip
         else
            orba(ib-norb)=ip
         endif
      enddo
   21 idif=0
      jdif=0
      do 30 i=1,nocc
      if(orbb(i).gt.norb)then
      idif=idif+1
      io1(idif)=i
      endif
      if(orba(i).le.norb)then
      jdif=jdif+1
c      io2(jdif)=i+nocc
      io2(jdif)=i
      endif
   30 continue
      isign=1
      do 40 i=1,idif
      isign=-isign
      itamp=orbb(io1(i))
      orbb(io1(i))=orba(io2(i))
      orba(io2(i))=itamp
 40   continue
      do i=1,nocc
         orba(i)=orba(i)-norb
      enddo
c----Looking for closed shells-------
      nclos=0
      do i=1,nocc
         do j=1,nocc
            if(orba(j).eq.orbb(i))then
               nclos=nclos+1
c               ipos(nclos)=i
c               jpos(nclos)=j
               zcli(i)=.true.
               zclj(j)=.true.
               goto 100
            endif
         enddo
 100     continue
      enddo
      do i=1,nocc
         if(.not.zcli(i))then
            do j=i+1,nocc
               if(zcli(j))then
                  itamp=orbb(i)
                  orbb(i)=orbb(j)
                  orbb(j)=itamp
                  zcli(j)=.false.
                  isign=-isign
                  goto 1000
               endif
            enddo
 1000       continue
         endif
      enddo
      do i=1,nocc
         if(.not.zclj(i))then
            do j=i+1,nocc
               if(zclj(j))then
                  itamp=orba(i)
                  orba(i)=orba(j)
                  orba(j)=itamp
                  zclj(j)=.false.
                  isign=-isign
                  goto 2000
               endif
            enddo
 2000       continue
         endif
      enddo
c-----riordino i closed shell      
      do i=1,nclos
         do j=i+1,nclos
            if(orbb(j).lt.orbb(i))then
               itamp=orbb(j)
               orbb(j)=orbb(i)
               orbb(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      do i=1,nclos
         do j=i+1,nclos
            if(orba(j).lt.orba(i))then
               itamp=orba(j)
               orba(j)=orba(i)
               orba(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
c--riordino da nclos+1 a nocc
      do i=nclos+1,nocc
         do j=i+1,nocc
            if(orbb(j).lt.orbb(i))then
               itamp=orbb(j)
               orbb(j)=orbb(i)
               orbb(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      do i=nclos+1,nocc
         do j=i+1,nocc
            if(orba(j).lt.orba(i))then
               itamp=orba(j)
               orba(j)=orba(i)
               orba(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      return
      end
#endif
#ifdef DOUBLET
      SUBROUTINE mkorb(m,nd,ne,trou,part,nocc,norb,isign)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer orba(500),orbb(500),nd(*)
      integer*2 ne(*),trou(*),part(*),isign
      dimension io1(500),io2(500),ipos(500),jpos(500)
      equivalence (io1(1),ipos(1)),(io2(1),jpos(1))
      logical*1 zcli(500),zclj(500)
c---vacuum state creation-----------
      if (nocc.gt.500) then
      write (6,*) 'Too many electrons',nocc
      stop
      endif
      do i=1,nocc
         orbb(i)=i
         orba(i)=i+norb
         zcli(i)=.false.
         zclj(i)=.false.
      enddo
      nstart=nd(m)
      nec=ne(m)
      if(nec.eq.0)goto 21
      do  k=1,nec
         ib=trou(nstart+k)
         ip=part(nstart+k)
         if(ib.le.norb)then
            orbb(ib)=ip
         else
            orba(ib-norb)=ip
         endif
      enddo
   21 idif=0
      jdif=0
      do 30 i=1,nocc
      if(orbb(i).gt.norb)then
      idif=idif+1
      io1(idif)=i
      endif
      if(orba(i).le.norb)then
      jdif=jdif+1
      io2(jdif)=i
      endif
   30 continue
c--idif should be jdif+1------------
c      if((idif-jdif).ne.1)stop 'not equal 1'
      if((idif-jdif).ne.1)then
       print*,'det.',M,'idif - jdif not equal 1'
       stop 'not equal 1'
       call flush(6)
      endif
      isign=1
      do 40 i=1,jdif
      isign=-isign
      itamp=orbb(io1(i))
      orbb(io1(i))=orba(io2(i))
      orba(io2(i))=itamp
 40   continue
c--let's put orbb(idif) into nocc+1
      orba(nocc+1)=orbb(io1(idif))
c--let's account for sign
      if(mod(io1(idif),2).eq.1)isign=-isign
c--let's squeeze orbb
      do i=io1(idif)+1,nocc
         orbb(i-1)=orbb(i)
      enddo
c--let's look for infinity on the right hand side
      do i=1,nocc+1
         if(orba(i).eq.(norb*2+1))goto 101
      enddo
      print*,'Questo non si deve mai vedere'
      stop 'impossibile'
 101  continue
c--let's get rid of infinity
      if(mod(nocc+1-i,2).eq.1)isign=-isign
      do j=i+1,nocc+1
         orba(j-1)=orba(j)
      enddo
c--      
      do i=1,nocc
         orba(i)=orba(i)-norb
      enddo
c      print*,'prima del riordino orbb,orba'
c      print '(7i3)',(orbb(i),i=1,nocc)
c      print '(7i3)',(orba(i),i=1,nocc)
c----Looking for closed shells-------
      nclos=0
      do i=1,nocc-1
         do j=1,nocc
            if(orba(j).eq.orbb(i))then
               nclos=nclos+1
               zcli(i)=.true.
               zclj(j)=.true.
               goto 100
            endif
         enddo
 100     continue
      enddo
      do i=1,nocc-1
         if(.not.zcli(i))then
            do j=i+1,nocc-1
               if(zcli(j))then
                  itamp=orbb(i)
                  orbb(i)=orbb(j)
                  orbb(j)=itamp
                  zcli(j)=.false.
                  isign=-isign
                  goto 1000
               endif
            enddo
 1000       continue
         endif
      enddo
      do i=1,nocc
         if(.not.zclj(i))then
            do j=i+1,nocc
               if(zclj(j))then
                  itamp=orba(i)
                  orba(i)=orba(j)
                  orba(j)=itamp
                  zclj(j)=.false.
                  isign=-isign
                  goto 2000
               endif
            enddo
 2000       continue
         endif
      enddo
c-----riordino i closed shells     
      do i=1,nclos
         do j=i+1,nclos
            if(orbb(j).lt.orbb(i))then
               itamp=orbb(j)
               orbb(j)=orbb(i)
               orbb(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      do i=1,nclos
         do j=i+1,nclos
            if(orba(j).lt.orba(i))then
               itamp=orba(j)
               orba(j)=orba(i)
               orba(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
c--riordino da nclos+1 a nocc-1 e poi a nocc
      do i=nclos+1,nocc-1
         do j=i+1,nocc-1
            if(orbb(j).lt.orbb(i))then
               itamp=orbb(j)
               orbb(j)=orbb(i)
               orbb(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      do i=nclos+1,nocc
         do j=i+1,nocc
            if(orba(j).lt.orba(i))then
               itamp=orba(j)
               orba(j)=orba(i)
               orba(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
c      print*,'dopo il riordino orbb,orba'
c      print '(7i3)',(orbb(i),i=1,nocc)
c      print '(7i3)',(orba(i),i=1,nocc)
      return
      end
#endif
#ifdef TRIPLET
      SUBROUTINE mkorb(m,nd,ne,trou,part,nocc,norb,isign)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer orba(500),orbb(500),nd(*)
      integer*2 ne(*),trou(*),part(*),isign
      dimension io1(500),io2(500),ipos(500),jpos(500)
      equivalence (io1(1),ipos(1)),(io2(1),jpos(1))
      logical*1 zcli(500),zclj(500)
c---vacuum state creation-----------
      if (nocc.gt.500) then
      write (6,*) 'Too many electrons',nocc
      stop
      endif
      do i=1,nocc
         orbb(i)=i
         orba(i)=i+norb
         zcli(i)=.false.
         zclj(i)=.false.
      enddo
      nstart=nd(m)
      nec=ne(m)
      if(nec.eq.0)goto 21 !This should never happen
      do  k=1,nec
         ib=trou(nstart+k)
         ip=part(nstart+k)
         if(ib.le.norb)then
            orbb(ib)=ip
         else
            orba(ib-norb)=ip
         endif
      enddo
   21 idif=0
      jdif=0
      do 30 i=1,nocc
      if(orbb(i).gt.norb)then
      idif=idif+1
      io1(idif)=i
      endif
      if(orba(i).le.norb)then
      jdif=jdif+1
      io2(jdif)=i
      endif
   30 continue
c--idif should be jdif+1------------
      isign=1
      do 40 i=1,jdif
      isign=-isign
      itamp=orbb(io1(i))
      orbb(io1(i))=orba(io2(i))
      orba(io2(i))=itamp
 40   continue
c--let's put orbb(idif) into nocc+1
      orba(nocc+1)=orbb(io1(idif))
c--let's account for sign
      if(mod(io1(idif),2).eq.1)isign=-isign
c--let's squeeze orbb
      do i=io1(idif)+1,nocc
         orbb(i-1)=orbb(i)
      enddo
c--      
      do i=1,nocc+1
         orba(i)=orba(i)-norb
      enddo
c      print*,'prima del riordino orbb,orba'
c      print '(7i3)',(orbb(i),i=1,nocc)
c      print '(7i3)',(orba(i),i=1,nocc)
c----Looking for closed shells-------
      nclos=0
      do i=1,nocc-1
         do j=1,nocc+1
            if(orba(j).eq.orbb(i))then
               nclos=nclos+1
               zcli(i)=.true.
               zclj(j)=.true.
               goto 100
            endif
         enddo
 100     continue
      enddo
      do i=1,nocc-1
         if(.not.zcli(i))then
            do j=i+1,nocc-1
               if(zcli(j))then
                  itamp=orbb(i)
                  orbb(i)=orbb(j)
                  orbb(j)=itamp
                  zcli(j)=.false.
                  isign=-isign
                  goto 1000
               endif
            enddo
 1000       continue
         endif
      enddo
      do i=1,nocc+1
         if(.not.zclj(i))then
            do j=i+1,nocc+1
               if(zclj(j))then
                  itamp=orba(i)
                  orba(i)=orba(j)
                  orba(j)=itamp
                  zclj(j)=.false.
                  isign=-isign
                  goto 2000
               endif
            enddo
 2000       continue
         endif
      enddo
c-----riordino i closed shells     
      do i=1,nclos
         do j=i+1,nclos
            if(orbb(j).lt.orbb(i))then
               itamp=orbb(j)
               orbb(j)=orbb(i)
               orbb(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      do i=1,nclos
         do j=i+1,nclos
            if(orba(j).lt.orba(i))then
               itamp=orba(j)
               orba(j)=orba(i)
               orba(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
c--riordino da nclos+1 a nocc-1 e poi a nocc+1
      do i=nclos+1,nocc-1
         do j=i+1,nocc-1
            if(orbb(j).lt.orbb(i))then
               itamp=orbb(j)
               orbb(j)=orbb(i)
               orbb(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
      do i=nclos+1,nocc+1
         do j=i+1,nocc+1
            if(orba(j).lt.orba(i))then
               itamp=orba(j)
               orba(j)=orba(i)
               orba(i)=itamp
               isign=-isign
            endif
         enddo
      enddo
c      print*,'dopo il riordino orbb,orba'
c      print '(7i3)',(orbb(i),i=1,nocc)
c      print '(7i3)',(orba(i),i=1,nocc)
c      print '(7i3)',isign                   
      return
      end
#endif
c------------------------------------------------------------------------
      subroutine hdigc(i,icomp1,iconf1,finv,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      integer*2 icomp1,iconf1,isegno(*)
      pointer hdiag,stok
      dimension hdiag(:),stok(:)
      COMMON /PER/ HDIAG,STOK
      icomp=icomp1
      iconf=iconf1
      zfa=.false.
      if(iconf.eq.0)then  !closed shell
         CALL HDIG(I,ZFA,DUM,AB)
         cmat(1,1)=ab*finv
         mm=i
	 stok(mm)=ab
         return
      endif
#ifndef SINGLET
      do m=1,icomp
         mm=m+i-1
         CALL HDIG(mm,ZFA,DUM,AB)
         hmat(m,m)=ab*finv
         stok(mm)=ab
         do n=1,m-1
            nn=n+i-1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*FINV*isegno(mm)*isegno(nn)
            hmat(n,m)=hmat(m,n)
         enddo
      enddo
#else
      do m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         CALL HDIG(mm,ZFA,DUM,AB)
         hmat(m,m)=ab*finv
         stok(mm)=ab
         do n=1,m-1
            nn=n+i-1
            if(.not.zmc)znf=mod(n,2).eq.1
            call givecase(zmc,zmc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*FINV*isegno(mm)*isegno(nn)
            hmat(n,m)=hmat(m,n)
            if (ncaso.eq.1)then
               hmat(m,n+1)=hmat(m,n)
               hmat(n+1,m)=hmat(m,n)
            elseif (ncaso.eq.2)then
               hmat(m+1,n)=hmat(m,n)
               hmat(n,m+1)=hmat(m,n)
            elseif (ncaso.eq.3)then
               hmat(m+1,n+1)=hmat(m,n)
               hmat(n+1,m+1)=hmat(m,n)
            elseif (ncaso.eq.4)then
               hmat(m-1,n+1)=hmat(m,n)
               hmat(n+1,m-1)=hmat(m,n)
            endif
 1       enddo
      enddo
#endif
#ifdef SINGLET
      if(iconf.eq.1)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T1,2,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T1,2,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.2)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T2,6,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T2,6,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.5)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T3,20,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T3,20,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.14)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T4,70,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T4,70,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.42)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T5,252,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T5,252,tamp,924,0.d0
     $        ,cmat,132)
         elseif(iconf.eq.132)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T6,924,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T6,924,tamp,924,0.d0
     $        ,cmat,132)
      endif
#endif
#ifdef DOUBLET
      if(iconf.eq.1)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T0,1,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T0,1,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.2)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T1,3,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T1,3,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.5)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T2,10,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T2,10,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.14)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T3,35,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T3,35,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.42)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T4,126,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T4,126,tamp,462,0.d0
     $        ,cmat,132)
         elseif(iconf.eq.132)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T5,462,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T5,462,tamp,462,0.d0
     $        ,cmat,132)
      endif
#endif
#ifdef TRIPLET
      if(iconf.eq.1)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T0,1,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T0,1,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.3)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T1,4,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T1,4,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.9)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T2,15,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T2,15,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.28)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T3,56,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T3,56,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.90)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T4,210,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T4,210,tamp,792,0.d0
     $        ,cmat,297)
         elseif(iconf.eq.297)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T5,792,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T5,792,tamp,792,0.d0
     $        ,cmat,297)
      endif
#endif
      return
      end
c-------------------------------------------------------------------
      subroutine hntdc(i,j,icomp1,jcomp1,iconf1,jconf1,finv,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1))
      integer*2 icomp1,jcomp1,iconf1,jconf1,isegno(*)
      icomp=icomp1
      jcomp=jcomp1
      iconf=iconf1
      jconf=jconf1
#ifndef SINGLET      
      do m=1,icomp
         mm=m+i-1
         do n=1,jcomp
            nn=n+j-1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*FINV*isegno(mm)*isegno(nn)
         enddo
      enddo
#else
      do m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         do n=1,jcomp
            nn=n+j-1
            znc=jcomp.eq.1
            if(.not.znc)znf=mod(n,2).eq.1
            call givecase(zmc,znc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*FINV*isegno(mm)*isegno(nn)
            if (ncaso.eq.1)then
               hmat(m,n+1)=hmat(m,n)
            elseif (ncaso.eq.2)then
               hmat(m+1,n)=hmat(m,n)
            elseif (ncaso.eq.3)then
               hmat(m+1,n+1)=hmat(m,n)
            elseif (ncaso.eq.4)then
               hmat(m-1,n+1)=hmat(m,n)
            endif
 1       enddo
      enddo
#endif
c      call matout(icomp,jcomp,hmat,924)
      indi=ivpos(icomp1)
      indj=ivpos(jcomp1)
      if(iconf.eq.0)then
         ic=1
      else
         ic=iconf
      endif
      if(jconf.eq.0)then
         jc=1
      else
         jc=jconf
      endif
#ifdef SINGLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,924,buff(indj),jcomp,
     $    0.d0,tamp,924)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,924,0.d0
     $     ,cmat,132)
#endif
#ifdef DOUBLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,462,buff(indj),jcomp,
     $    0.d0,tamp,462)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,462,0.d0
     $     ,cmat,132)
#endif
#ifdef TRIPLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,792,buff(indj),jcomp,
     $    0.d0,tamp,792)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,792,0.d0
     $     ,cmat,297)
#endif
      return
      end
c-------------------------------------------------------
      subroutine putbuf(iconfi,iconff,cmat,ibuf,jbuf,buffer,ind,hic,idim
     $     ,ivolte)
#if defined (SINGLET) || defined (DOUBLET)
      dimension ibuf(idim),jbuf(idim),buffer(idim),cmat(132,*),hic(*)
#endif
#ifdef TRIPLET
      dimension ibuf(idim),jbuf(idim),buffer(idim),cmat(297,*),hic(*)
#endif
      real*8 buffer,cmat,hic,thres,temp
      data thres /1.0d-12/
      ii=0
      do i=iconfi,iconff
         ii=ii+1
         jj=0
         do j=iconfi,i
            jj=jj+1
            if(i.eq.j)then
               hic(i)=cmat(ii,jj)
            else
               temp=cmat(ii,jj)
               if(abs(temp).gt.thres)then
               ind=ind+1
               ibuf(ind)=i
               jbuf(ind)=j
               buffer(ind)=cmat(ii,jj)
	       endif
            endif
            if(ind.eq.idim)then
               write(21)ibuf,jbuf,buffer
               ivolte=ivolte+1
               ind=0
            endif
         enddo
      enddo
      return
      end
c-------------------------------------------------------
      subroutine putbuf2(iconfi,jconfi,iconff,jconff,cmat,ibuf,jbuf
     $     ,buffer,ind,idim)
#if defined (SINGLET) || defined (DOUBLET)
      dimension ibuf(idim),jbuf(idim),buffer(idim),cmat(132,*)
#endif
#ifdef TRIPLET
      dimension ibuf(idim),jbuf(idim),buffer(idim),cmat(297,*)
#endif
      real*8 buffer,cmat,thres,temp
      data thres /1.0d-12/
      ii=0
      do i=iconfi,iconff
         ii=ii+1
         jj=0
         do j=jconfi,jconff
            jj=jj+1
            temp=cmat(ii,jj)
            if(abs(temp).gt.thres)then
            ind=ind+1
            ibuf(ind)=i
            jbuf(ind)=j
            buffer(ind)=cmat(ii,jj)
	    endif
            if(ind.eq.idim)then
               write(20)ibuf,jbuf,buffer
               ind=0
            endif
         enddo
      enddo
      return
      end
c----------------------------------------------------
      subroutine mkdet(c,cdet,nconf,nvec,secrconf)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      integer*2 icomp,iconf,isegno,ico
      dimension csq(200)
      common /fil/ file12,file13,file26,zwconf
      character*64 file12,file13,file26
      COMMON /DIA/ SFO(5),SEN,SECR,CSEL(ID12),ISEL(ID12),NITER,IOP,MAXM,
     * maxgue,ietats(ID6),ndb,ndab,PROJF,ZSS,ZVEC,ZWRT,ZFOL
      character*64 projf
      INTEGER*2 ietats
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hma(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hma(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hma(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1))
      dimension c(*),cdet(*)
      if (nvec.gt.200) then
      print *,'Errore, nvec >200',nvec
      stop
      endif
      minit=1
      mfinal=0
      iconfi=1
      ibasta=1
c
c     Scrivo i coefficienti sulle configurazioni nel file29
c
      if (zwconf) then
      write (6,*) 'CSF''s coefficients written on file29'
      write (6,*) 'nconf,nvec',nconf,nvec
      write (6,*) 'nrot,ietats',nrot,(ietats(i),i=1,nrot)
      write (29) nconf,nrot
      write (29) ((c(i+nconf*(ietats(j)-1)),i=1,nconf),j=1,nrot)
      close (29)
      endif

      do while (ibasta.eq.1)
         ico=icomp(minit)
         indi=ivpos(ico)
         icoc=iconf(minit)
         if(icoc.eq.0)icoc=1
         mfinal=minit+ico-1
         iconff=iconfi+icoc-1
         mm=0
	 do i=1,nvec
	 csq(i)=0.d0
	 enddo
         do m=minit,mfinal
            cdet(m)=0.d0
            mm=mm+1
            ii=0
            do i=iconfi,iconff
               ii=ii+1
               ind=ico*(ii-1)+mm+indi-1
               in=-ncf
               inc=-nconf
               do istate=1,nvec
                  in=in+ncf
                  inc=inc+nconf
                  cdet(in+m)=cdet(in+m)+c(inc+i)*buff(ind)*isegno(m)
               enddo
            enddo
	    in=-ncf
	    do istate=1,nvec
	    in=in+ncf
	    csq(istate)=csq(istate)+cdet(in+m)**2
	    enddo
         enddo
	 zpr=.false.
	 do istate=1,nvec
	 if (csq(istate).gt.secrconf) then
	 zpr=.true.
	 goto 98
	 endif
	 enddo
  98     if (zpr) call wrtcap(minit,csq,nvec)
         minit=mfinal+1
         if(mfinal.eq.ncf)ibasta=0
         iconfi=iconff+1
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine singsub(M,IC,mcapos,nconf,iocc,idimocc)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      common/debug/zdebug
      common /deephole/ ndhb,zdh,ndhb1,ndhb2,zddh
      integer*2 icomp,iconf,isegno,ico
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hma(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hma(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hma(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1))
      allocatable jocc(:)
      save jocc
c      integer*1 iocc(idimocc,*),jocc(id1)
      integer*1 iocc(idimocc,*),jocc
c--M e' il determinante capostipite
c--IC e' la prima configurazione che corrisponde a M      
c--iocc e' il vettore di occupazione (2,1,0) relativo al capostipite M      
c--jocc e' il vettore delle singole sostituzioni su M
c--calcoliamo iocc
c      call giveocc(iocc,m,nd,ne,trou,part,nocb,norb)
      if(.not.allocated(jocc))allocate(jocc(norb))
      znexe=.false.
      do i=1,norb
         jocc(i)=iocc(i,mcapos)
      enddo
      nene=2*ne(m)
      do i=1,norb
c         if(zdebug)then
c            print*,'singsub: doing orb. n. ',i
c            call flush(6)
c         endif
c--togliamo un orbitale da jocc
         if(jocc(i).eq.0)goto 1
         if(i.le.nocb)then
            ine=1
         else
            ine=-1
         endif
         jocc(i)=jocc(i)-1
         isy=itsym(i)
         do j=1,norb
c--aggiungiamo un orbitale a jocc
            if(j.eq.i) goto 2 !must not be the one we removed
            if(jocc(j).eq.2) goto 2
            jsy=itsym(j)
            if(its(isy,jsy).ne.1)goto 2 !symmetry check
            if(zdh.and.j.eq.ndhb.and.jocc(j).eq.1)goto 2
         if(j.le.nocb)then
            jne=-1
         else
            jne=1
         endif
            jocc(j)=jocc(j)+1

                   if(zddh) then
c                   zdb1=jocc(ndhb1).eq.2
c                   zdb2=jocc(ndhb2).eq.2
                   zdb1=jocc(ndhb1).eq.1
                   zdb2=jocc(ndhb2).eq.1
c                   if (zdb1.and.zdb2) goto 22
                   if (.not.(zdb1.or.zdb2)) goto 22
                   endif

c--A questo punto jocc e' un capostipite
#if defined (TRIPLET) || defined (QUADRUPLET)
            icountos=0
            do iiii=1,norb
               if(jocc(iiii).eq.1)icountos=icountos+1
            enddo
#endif
#ifdef TRIPLET
            if (icountos.lt.2) goto 22
#endif
#ifdef QUADRUPLET
            if (icountos.lt.3) goto 22
#endif
            if(zcas)then
               zxe=zact(i).and.zact(j)
               if(zxe)goto 22   !true for CAS
            endif
c--controllo che non sia gia' generato
            ncapos=0
            ne2=(nene+ine+jne)/2
            do n=1,m-1          !loop sui capostipiti
               if(zcapos(n))then
                  ncapos=ncapos+1
c                  call giveocc(kocc,n,nd,ne,trou,part,nocb,norb)
                  if(abs(ne(n)-ne2).le.2)then
                     if(ndiff(iocc(1,ncapos),jocc,norb).le.4)goto 22
c--deja genere                     
                  endif
               endif
            enddo
c--controllo che non sia gia' nello spazio S
            if(zcas)goto 1111
            ncapos=mcapos  !m e' capostipite
            do n=m+1,ncf  !loop sui capostipiti
               if(zcapos(n))then
                  ncapos=ncapos+1
c                  call giveocc(kocc,n,nd,ne,trou,part,nocb,norb)
                  if(ne2.eq.ne(n))then
                     if(znequ(iocc(1,ncapos),jocc,norb))goto 22 !deja en S
                  endif
               endif
            enddo
c--   Qui jocc e' utilizzabile per la perturbazione
 1111       kcaso=0
            if (.not.zact(i).and.zact(j)) kcaso=6
            if (zact(i).and..not.zact(j)) kcaso=7
            if (.not.zact(i).and..not.zact(j)) kcaso=8
	    if (kcaso.eq.0) then
	    write (6,*) 'Errore, kcaso=0'
	    write (6,*) i,j,zact(i),zact(j)
	    stop 
	    endif
            call pertu(m,ic,mcapos,jocc,nconf,ne2,iocc,norb,
     *                 znexe,kcaso)
c                          !jocc interacts only with m,m+1,...
 22         jocc(j)=jocc(j)-1
 2          continue
         enddo
         jocc(i)=jocc(i)+1
 1       continue
      enddo
      return
      end
c---------------------------------------------------------------
      subroutine doubsub(M,IC,mcapos,nconf,iocc,idimocc,e2mp,e2en)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
c--renzo debug
      dimension e2mp(*),e2en(*)   !togliere dopo il debug
c--renzo debug end
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      integer*2 icomp,iconf,isegno,ico
      common/debug/zdebug
      common /deephole/ ndhb,zdh,ndhb1,ndhb2,zddh
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hma(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hma(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hma(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1))
      allocatable jocc(:)
      save jocc
      integer*1 iocc(idimocc,*),jocc
c--iocc e' il vettore di occupazione (2,1,0) relativo al capostipite M      
c--jocc e' il vettore delle doppie sostituzioni su M
c--calcoliamo iocc
c      call giveocc(iocc,m,nd,ne,trou,part,nocb,norb)
      if(.not.allocated(jocc))allocate(jocc(norb))
      do i=1,norb
         jocc(i)=iocc(i,mcapos)
      enddo
      nene=ne(m)*2
      do i=1,norb
c         if(zdebug)then
c            print*,'doubsub: doing orb. n. ',i
c            call flush(6)
c         endif
c--tolgo un orbitale da jocc
         if(jocc(i).eq.0) goto 1
         if(i.le.nocb)then
            ine=1
         else
            ine=-1
         endif
         jocc(i)=jocc(i)-1
         isy=itsym(i)
         do j=i,norb
c--tolgo un secondo orbitale da jocc
            if(jocc(j).eq.0)goto 2
            jocc(j)=jocc(j)-1
            if(j.le.nocb)then
               jne=1
            else
               jne=-1
            endif
            jsy=itsym(j)
            ijsy=its(isy,jsy)
c	    write (6,*) 'Faccio',i,j
	    call flush (6)   !tolgo flush
            do k=1,norb
c--aggiungo un orb. a jocc, purche' non sia dove l'ho tolto prima
               if(k.eq.i.or.k.eq.j)goto 3
               if(jocc(k).eq.2)goto 3
               if(zdh.and.k.eq.ndhb.and.jocc(k).eq.1)goto 3

               if(k.le.nocb)then
                  kne=-1
               else
                  kne=1
               endif
               jocc(k)=jocc(k)+1
               ksy=itsym(k)
               ijksy=its(ijsy,ksy)
               do l=k,norb
c--aggiungo un secondo orb. a jocc
                  if(l.eq.i.or.l.eq.j)goto 4
                  if(jocc(l).eq.2)goto 4
                  lsy=itsym(l)
                  if(its(lsy,ijksy).ne.1)goto 4
                  if(zdh.and.l.eq.ndhb.and.jocc(l).eq.1)goto 4

                  if(l.le.nocb)then
                     lne=-1
                  else
                     lne=1
                  endif
                  jocc(l)=jocc(l)+1
c	       write (6,'(''Doppia eccitazione'',4i3,'' M='',i2)') 
c     *                   i,j,k,l,m

                   if(zddh) then
c                   zdb1=jocc(ndhb1).eq.2
c                   zdb2=jocc(ndhb2).eq.2
                   zdb1=jocc(ndhb1).eq.1
                   zdb2=jocc(ndhb2).eq.1
c                   if (zdb1.and.zdb2) goto 44
                   if (.not.(zdb1.or.zdb2)) goto 44
                   endif

#if defined (TRIPLET) || defined (QUADRUPLET)
            icountos=0
            do iiii=1,norb
               if(jocc(iiii).eq.1)icountos=icountos+1
            enddo
#endif
#ifdef TRIPLET
            if (icountos.lt.2) goto 44
#endif
#ifdef QUADRUPLET
            if (icountos.lt.3) goto 44
#endif
c--controllo che jocc non sia gia' generato
            if(zcas)then
               zxe=zact(i).and.zact(j).and.zact(k).and.zact(l)
               if(zxe)goto 44   !true for CAS
	       znexe=.not.(zact(i).or.zact(j).or.zact(k).or.zact(l))
               if(znexe)then
	       if (i.ne.j.and.k.ne.l) icaso=1
	       if (i.eq.j.and.k.ne.l) icaso=2
	       if (i.ne.j.and.k.eq.l) icaso=3
	       if (i.eq.j.and.k.eq.l) icaso=4
c	       write (6,'(''Doppia inattiva'',4i3,'' icaso='',i2)') 
c     *                   i,j,k,l,icaso
	       goto 1111
	       endif
            endif
                  ncapos=0
                  ne2=(nene+ine+jne+kne+lne)/2
                  do n=1,m-1    !loop sui capostipiti
                     if(zcapos(n))then
                        ncapos=ncapos+1
c                        call giveocc(kocc,n,nd,ne,trou,part,nocb,norb)
                        if(abs(ne2-ne(n)).le.2)then
                           if(ndiff(iocc(1,ncapos),jocc,norb).le.4)goto
     $                          44 !deja vu
                        endif
                     endif
                  enddo
c--controllo che non sia gia' nello spazio S
                  if(zcas)goto 1111
                  ncapos=mcapos
            do n=m+1,ncf  !loop sui capostipiti
               if(zcapos(n))then
                  ncapos=ncapos+1
c                  call giveocc(kocc,n,nd,ne,trou,part,nocb,norb)
                  if(ne2.eq.ne(n))then
                     if(znequ(iocc(1,ncapos),jocc,norb))goto 44 !deja en S
                  endif
               endif
            enddo
c--Qui jocc e' utilizzabile per la perturbazione
 1111       kcaso=0
            if (.not.(zact(i).or.zact(j).or.zact(k).or.zact(l))) kcaso=1
            if (.not.(zact(i).or.zact(j)).and.(zact(k).neqv.zact(l)))
     *          kcaso=2
            if ((zact(i).neqv.zact(j)).and..not.(zact(k).or.zact(l)))
     *          kcaso=3
            if (.not.(zact(i).or.zact(j)).and.(zact(k).and.zact(l))) 
     *          kcaso=4
            if ((zact(i).and.zact(j)).and..not.(zact(k).or.zact(l)))
     *          kcaso=5
            if ((zact(i).neqv.zact(j)).and.(zact(k).and.zact(l)))
     *          kcaso=6
            if ((zact(i).and.zact(j)).and.(zact(k).neqv.zact(l))) 
     *          kcaso=7
            if ((zact(i).neqv.zact(j)).and.(zact(k).neqv.zact(l)))
     *          kcaso=8
c            if (kcaso.eq.8) write (6,*) 'indici',i,j,k,l

            if (kcaso.eq.0) then 
	    write (6,*) 'Errore, kcaso=0'
	    write (6,*) i,j,k,l
	    stop 
	    endif
            call pertu(m,ic,mcapos,jocc,nconf,ne2,iocc,norb,znexe,kcaso)
c--renzo debug
c            print '(a,2f15.9)','doub: e2mp,e2en ',e2mp(5)*finv,e2en(5)
c     $           *finv
c--renzo debug end
 44               jocc(l)=jocc(l)-1
 4             enddo
               jocc(k)=jocc(k)-1
 3          enddo
            jocc(j)=jocc(j)+1
 2       enddo
         jocc(i)=jocc(i)+1
 1    enddo
      return
      end
c-----------------------------------------------------------
      subroutine scara(noccup,norb,nocb,ncf,nd,ne,trou,part,isz,zion)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
c      dimension noccup(*),morb(id1),mspin(id1),isingle(id1)
      dimension noccup(*)
      pointer morb(:),mspin(:),isingle(:),iel(:),morbp(:)
c      save morb,mspin,isingle
c      dimension trou(5000),part(5000) 
c      dimension ne(924),nd(924)
      dimension trou(*),part(*) 
      dimension ne(*),nd(*)
      integer*2 ne,trou,part
c      integer*2 noccup,morb,mspin
      integer*2 morb,mspin
      integer*1 noccup
      common /ric/ ihalf,nec1,morb,isingle,mspin,iel,morbp

C     Da passare: noccup,zion,norb,nocb,isz
      ncf=0
      nec1=0
c      do i=1,id1
      do i=1,norb
      morb(i)=0
      mspin(i)=-1
      enddo
      nsingle=0
C  
C     Ciclo sui doppi occupati nel riferimento.
C     nsingle=numero di orb. singoli occupati.
C     Si attribuisce morb che contiene il numero dell'orbitale occupato
C     dalla buca nec-esima.
C
      do 100 i=1,nocb
      if (noccup(i).eq.2) goto 100
      if (noccup(i).eq.1) then
      nec1=nec1+1
      nsingle=nsingle+1
      morb(nec1)=i
      isingle(nsingle)=nec1
      endif
      if (noccup(i).eq.0) then
      nec1=nec1+2
      morb(nec1-1)=i
      mspin(nec1-1)=1
      morb(nec1)=i
      mspin(nec1)=0
      endif
 100  enddo
C  
C     Ciclo sui virtuali nel riferimento.
C     Si attribuisce morb che contiene il numero dell'orbitale occupato
C     dalla particella (nec1)-esima.
C
      do 110 i=nocb+1,norb
      if (noccup(i).eq.0) goto 110
      if (noccup(i).eq.1) then
      nec1=nec1+1
      nsingle=nsingle+1
      morb(nec1)=i
      isingle(nsingle)=nec1
      endif
      if (noccup(i).eq.2) then
      nec1=nec1+2
      morb(nec1-1)=i
      mspin(nec1-1)=1
      morb(nec1)=i
      mspin(nec1)=0
      endif
  110 enddo
C
C     Attribuzione dell'obitale "infinito" se zion=.true.
C
      if (zion) then
      nec1=nec1+1
      morb(nec1)=norb+1
      mspin(nec1)=1
      endif
C
C     nalpha = # di elettroni alpha tra gli nsingle elettroni spaiati
C
      if (zion)  then
      nalpha=isz*2+(nsingle+1-isz*2)/2-1
      else
      nalpha=isz*2+(nsingle-isz*2)/2
      endif
C
C      Scrittura del determinante se non ci sono orbitali singoli occupati
C
      if(nsingle.gt.12)then
         print*,'Warging, troppe open shells',nsingle
cele         stop 'troppe open shells'
         ncf=0
         return
      endif
      if (nsingle.eq.0) then
      nec=nec1/2
      call stora(ncf,nec,norb,isz,nd,ne,trou,part,zion,nsingle)
      else
C   
C      Generazione automatica di tutti i determinanti con stesso
C      sz e stessa occupazione orbitalica.
C
       do i=1,nec1
       if (mspin(i).eq.-1.and.i.le.nec1/2) mspin(i)=1
       if (mspin(i).eq.-1.and.i.gt.nec1/2) mspin(i)=0
       enddo

c       write (6,'(''nec1 ='',I3)') nec1
c       write (6,'(''morb ='',10I3)') (morb(i),i=1,nec1)
c       write (6,'(''mspin='',10I3)') (mspin(i),i=1,nec1)
c       write (6,'(''nsingle ='',I3)') nsingle
c       write (6,'(''isingle='',10I3)') (isingle(i),i=1,nsingle)
C
C      Se SZ=0, allora in ricors verranno calcolati solo meta`
C      dei determinanti, perche` l'altra meta` sono i "flippati".
C      Si calcola quindi ihalf=meta` dei determinanti da generare:
C      ihalf=n!/[(n/2)!*(n/2)!]/2
C
       if (isz.eq.0.and..not.zion) then
        if (nalpha.ne.nsingle/2) then
        write (6,*) 'Con sz=0, deve valere nsingle=2*nalpha'
        write (6,*) 'nsingle=',nsingle,'nalpha=',nalpha
        stop 12
        endif
       ihalf=1
       do ip=nalpha+1,nsingle
       ihalf=ihalf*ip
       enddo
       do ip=1,nalpha
       ihalf=ihalf/ip
       enddo
       ihalf=ihalf/2
       endif
       call genera(ncf,nsingle,nalpha,norb,isz,nd,ne,trou,part,zion)
      endif
      return
      end
C##############################################################
      subroutine genera(ncf,n,k,norb,isz,nd,ne,trou,part,zion)
C
C     Calcola ricorsivamente come disporre kg elettroni alpha 
C     orbitali
C     nec e` due volte il numero di eccitazione, cioe` la somma di
C         buche e particelle.
C
      implicit none
      integer id1,id2,id3,id4,id5,id6,id7,id8,id9,id10,id11,id12,id13
     $     ,id14,id15
      include 'parameter.cippi'
      common /ric/ ihalf,nec1,morb,isingle,mspin,iel,morbp
      integer*2 morb,mspin,morbp
      integer*4 i,n,nj,isingle,j,k,nec,nel,iel,
     *  ihalf,ih,nec1,nd,ncf,isz,kb,metat,ini,norb
      logical*1 zion
      dimension trou(*),part(*) 
      dimension ne(*),nd(*)
      integer*2 ne,trou,part
c      dimension mspin(id1),isingle(id1),morb(id1),iel(id1)
      pointer mspin(:),isingle(:),morb(:),iel(:),morbp(:)
c      save mspin,isingle,morb,iel
c      if(.not.allocated(mspin))then
c         allocate(mspin(norb))
c         allocate(isingle(norb))
c         allocate(morb(norb))
c         allocate(iel(norb))
c      endif
      nec=nec1/2
      ih=0
      nel=1
      if (isingle(1).le.nec) then
      mspin(isingle(1))=0
      else
      mspin(isingle(1))=1
      endif
      iel(nel)=1
c
c     Se nel<k devo ancora sistemare degli elettroni.
c
   10 if (nel.lt.k) then
c
c     Se sto lavorando con il primo elettrone e dalla sua
c     posizione (iel(1)) alla fine (n) ci sono meno orbitali
c     del numero di elettroni che dovro sistemare (k-1)
c     allora ho finito.
c
       if (nel.eq.1.and.(n-iel(nel)).lt.(k-1)) return
c
c     Se il numero di elettroni ancora da sistemare (k-nel)
c     e` maggiore delle posizioni disponibili vado a 300,
c     altrimenti a 100
c
       if ((k-nel).gt.(n-iel(nel))) then
       goto 300
       else
       goto 100
       endif
      else
c 
c     Nel caso che nel = k, l'ultimo elettrone lo metto
c     negli orbitali che stanno da iel(nel-1) (posizione
c     dove ho messo il precedente elettrone) alla fine (n).
c
       if (k.eq.1) then
       ini=0
       else
       ini=iel(nel-1)
       endif
       do i=ini+1,n
       if (isingle(i).le.nec) mspin(isingle(i))=0
       if (isingle(i).gt.nec) mspin(isingle(i))=1
c
c      Scrivo il determinante
c
C       Se ho gia` fatto meta` dei determinanti nel caso isz=0, ritorno
c
        if (isz.eq.0.d0.and..not.zion) then
         ih=ih+1
         if (ih.gt.ihalf) return
        endif
        call stora(ncf,nec,norb,isz,nd,ne,trou,part,zion,n)
       if (isingle(i).le.nec) mspin(isingle(i))=1
       if (isingle(i).gt.nec) mspin(isingle(i))=0
       enddo
       if (nel.eq.1) return
       goto 200
      endif
c
c     Sistemo l'elettrone (nel+1)-esimo nella posizione 
c     iel(nel)+1 e torno a 10
c
  100 nel=nel+1
      iel(nel)=iel(nel-1)+1
      if (isingle(iel(nel)).le.nec) mspin(isingle(iel(nel)))=0
      if (isingle(iel(nel)).gt.nec) mspin(isingle(iel(nel)))=1
      goto 10
c
c     Se arrivo qui, significa che ho gia` sistemato in tutti
c     i modi possibili gli elettroni da nel a k  per una certa
c     configurazione dei precedenti nel elettroni. Allora 
c     sposto il (nel-1)-esimo elettrone in avanti di una 
c     posizione.
c
  200 nel=nel-1
      if (isingle(iel(nel)).le.nec) mspin(isingle(iel(nel)))=1
      if (isingle(iel(nel)).gt.nec) mspin(isingle(iel(nel)))=0
      iel(nel)=iel(nel)+1
      if (isingle(iel(nel)).le.nec) mspin(isingle(iel(nel)))=0  
      if (isingle(iel(nel)).gt.nec) mspin(isingle(iel(nel)))=1
      goto 10
c
c     Rimetto mspin beta a tutte le posizioni da iel(nel) a n
c     e vado a 200
c
  300 do i=iel(nel),n
      if (isingle(i).le.nec) mspin(isingle(i))=1
      if (isingle(i).gt.nec) mspin(isingle(i))=0  
      enddo
      goto 200
c
      end
C#######################################################################
      subroutine stora(ncf,nec,norb,isz,nd,ne,trou,part,zion,nsingle)
      implicit real*8(a-h,o-y),logical*1(z)
      include 'parameter.cippi'
      common /ric/ ihalf,nec1,morb,isingle,mspin,iel,morbp
c      integer*2 morb,mspin,morbp(id1)
      integer*2 morb,mspin,morbp
      dimension trou(*),part(*) 
      dimension ne(*),nd(*)
      integer*2 ne,trou,part
c      dimension mspin(id1),isingle(id1),morb(id1),iel(id1)
      pointer mspin(:),isingle(:),morb(:),iel(:),morbp(:)
c
c        Costruzione del vettore itab degli spin-orbitali
c
      do i=1,2*nec
      morbp(i)=morb(i)+mspin(i)*norb
      enddo
      call ordom(nec,morbp)

      ncf=ncf+1
      if (ncf.eq.1) then
      nd(ncf)=0
      else 
      nd(ncf)=nd(ncf-1)+ne(ncf-1)
      endif
      ne(ncf)=nec
      do i=1,nec
      trou(nd(ncf)+i)=morbp(i)
      part(nd(ncf)+i)=morbp(i+nec)
      enddo
c
c     Se e` il caso storo anche il flippato
c
      if(nsingle.eq.0)return
      if (isz.eq.0.and..not.zion) then
        do i=1,2*nec
           jdu=1-mspin(i)
           jdu=jdu*norb
           idu=morb(i)
           morbp(i)=idu+jdu
c        morbp(i)=morb(i)+(1-mspin(i))*norb !doesn't work with ifc -xK
        enddo
        call ordom(nec,morbp)
        ncf=ncf+1
        nd(ncf)=nd(ncf-1)+ne(ncf-1)
        ne(ncf)=nec
	do i=1,nec
        trou(nd(ncf)+i)=morbp(i)
        part(nd(ncf)+i)=morbp(i+nec)
        enddo
      endif
      return
      end
C***********************************************************************
      subroutine ordom(nec,itab)
c
c     Ordina le buche e le particelle che definiscono il determinante
c     in modo crescente.
c
      implicit none
      integer*2 itab,n
      integer*4 nec,i,j
      dimension itab(*)
c      write (6,*) 'ORDOM: ordino ->',(itab(j),j=1,nec)
      do 24 i=1,nec-1
      do 24 j=i+1,nec
      if (itab(j).gt.itab(i)) goto 24
      n=itab(j)
      itab(j)=itab(i)
      itab(i)=n
   24 continue
      do 25 i=1+nec,nec+nec-1
      do 25 j=i+1,nec+nec
      if (itab(j).gt.itab(i)) goto 25
      N=itab(J)
      itab(j)=itab(I)
      itab(i)=n
   25 continue
      return
      end
c--------------------------------------------------------------
      subroutine pertu(M,IC,mcapos,JOCC,nconf,ne2,iocc,idimocc,
     * znexe,kcaso)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      integer*1 iocc(idimocc,*)
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      pointer cconf,icomp,iconf,isegno,zcapos,zact
      dimension cconf(:),icomp(:),iconf(:),isegno(:),zcapos(:),
     $     zact(:)
      common/detinf/cconf,icomp,iconf,isegno,
     * zcapos,ntotdet,ntotconf,ntotcap,zact
      common/third/zthird,zmp3,zen3
      common/classi/E2ENc(ID6,8),psienec(id6,8),E2MPc(ID6,8),
     * psimpec(id6,8)
      dimension cont(id6,132)
      dimension kaux(4),iposi(500),iposj(500)
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hma(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension vec(ID6,132) !ID6 e' il massimo metat
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hma(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension vec(ID6,132) !ID6 e' il massimo metat
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hma(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension vec(ID6,297) !ID6 e' il massimo metat
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      dimension psiloc(ID6),enloc(ID6),psimploc(id6),emploc(id6),
     * dintloc(id6)
cele 18-11-04
      dimension psisel(ID6),ensel(ID6)
cele 18-11-04
      logical*1 zprint
      integer*2 icomp,iconf,isegno
      integer*1 jocc(*)
#ifdef SINGLET
      ifi=132
#endif
#ifdef DOUBLET
      ifi=132
#endif
#ifdef TRIPLET
      ifi=297
#endif
      do i=1,ifi
         do j=1,metat
            vec(j,i)=0.d0
         enddo
      enddo
      initc=ic
      nf=ncf+1
      ndnf=nd(ncf)+ne(ncf)
      call scara (jocc,norb,nocb,ncf2,
     $     nd(nf),ne(nf),trou(ndnf+1),part(ndnf+1),isz,zion)
      call escla2(nd(nf),ne(nf),trou(ndnf+1),part(ndnf+1),icomp(nf),
     $     iconf(nf),isegno(nf),nocb,norb)
      if(iconf(nf).eq.-1)return   !more than 12 open shells
      ntotdet=ntotdet+ncf2
      ntotconf=ntotconf+iconf(nf)
      ntotcap=ntotcap+1
      do i=1,ncf2
      nd(ncf+i)=nd(ncf+i)+ndnf
      enddo
      iconf2=iconf(nf)
C
C     Energia MP dell'occupazione orbitalica
C
      call hmp(nf,enemp)
      ncapos=mcapos-1
      call giveopen(jocc,norb,iposj,nopenj)
      do n=m,ncf
         if(.not.zcapos(n))goto 1
         ncapos=ncapos+1
c--costruiamo la matrice cmat         
         iconff=iconf(n)
         ifinc=initc+iconff-1
         if(abs(ne2-ne(n)).gt.2)goto 2
         ndiffk=ndiff(iocc(1,ncapos),jocc,norb)
         if(ndiffk.gt.4)goto 2
#ifdef SINGLET
	 if (zcas.and.znexe) then
         call hntdc21(n,nf,icomp(n),icomp(nf),iconf(n),iconf(nf),
     $        isegno)
         else
	   call giveopen(iocc(1,ncapos),norb,iposi,nopeni)
	   ixe=nopenj-nopeni
           if(ixe.eq.4)then
              if(nopeni.eq.0)then
                 do lll=1,4
                    kaux(lll)=lll
                 enddo
              else
                 call seocc(iocc(1,ncapos),jocc,norb,kaux,nopeni,nopenj
     $                ,iposi,iposj)
              endif
              call hntdc24(n,nf,icomp(n),icomp(nf),iconf(n),iconf(nf),
     $                     kaux,iposj,norb,isegno)
           elseif(ixe.eq.2.and.ndiffk.eq.4.and.
     $            zsubset(iposi,nopeni,iposj,nopenj))then
              if(nopeni.eq.0)then
                 kaux(1)=1
                 kaux(2)=2
              else
                 call seocc(iocc(1,ncapos),jocc,norb,kaux,nopeni,nopenj
     $                ,iposi,iposj)
              endif
              call hntdc22(n,nf,icomp(n),icomp(nf),iconf(n),iconf(nf),
     $             kaux,iposj,norb,isegno)
           else
              call hntdc2(n,nf,icomp(n),icomp(nf),iconf(n),iconf(nf),
     $             isegno)
           endif
         endif
#endif
#ifndef SINGLET
c
c          Faccio cosi per semplicita`: in effetti si dovrebbero poter
c          utilizzare le subroutine ottimizzate che uso per il singoletto
c
              call hntdc2(n,nf,icomp(n),icomp(nf),iconf(n),iconf(nf),
     $             isegno)
#endif
c--cmat ha dimensione iconf(n) x iconf(nf)
c--aggiorniamo il vettore vec
c         do ico=initc,ifinc
c            ic2=ico-initc+1
c            do icp=1,iconf2
c               do istate=1,metat
c                  ista=(istate-1)*nconf+ico  !nconf deve essere noto
c                  ista2=(istate-1)*nconf+icp  !nconf deve essere noto
c                  vec(ista2)=vec(ista2)+cconf(ista)*cmat(ic2,icp)
c               enddo
c            enddo
c         enddo
         do istate=1,metat
            call dcopy(iconff,cconf((istate-1)*nconf+initc),1,tamp(1
     $           ,istate),1)
         enddo
#ifdef SINGLET
         call dgemm('T','N',metat,iconf2,iconff,1.d0,tamp,924,cmat,132,
     $        1.d0,vec,id6)
#endif
#ifdef DOUBLET
         call dgemm('T','N',metat,iconf2,iconff,1.d0,tamp,462,cmat,132,
     $        1.d0,vec,id6)
#endif
#ifdef TRIPLET
         call dgemm('T','N',metat,iconf2,iconff,1.d0,tamp,792,cmat,297,
     $        1.d0,vec,id6)
#endif
 2       initc=ifinc+1
 1    enddo
c--calcoliamo le energie delle configurazioni perturbatrici
c      call hdigc(nf,icomp(nf),iconf(nf),finv,isegno)
cele      call hdigc2(nf,icomp(nf),iconf(nf),isegno)
      call hdigc21(nf,icomp(nf),iconf(nf),isegno)
c--i risultati saranno in cmat
c--ora possiamo calcolare i contributi perturbativi      
c--Epstein-Nesbet
      zprint=.false.
cele 18-11-04
      do istate=1,metat
        psisel(istate)=0.d0
        ensel(istate)=0.d0
      enddo
cele 18-11-04
c--renzo debug
c      zcondiz=ntotdet.eq.470960
      do icp=1,iconf2
         ijcou=0
         do istate=1,metat
            psiloc(istate)=0.d0
            enloc(istate)=0.d0
            psimploc(istate)=0.d0
            emploc(istate)=0.d0
            dintloc(istate)=0.d0
c            if(zcondiz.and.istate.eq.5)print*,'icp=',icp
c            if(zcondiz.and.istate.eq.5)print*,'vec,een,cmat',vec(5,icp)
c     $           ,een(5),cmat(icp,icp)
            co=vec(istate,icp)/(een(istate)-cmat(icp,icp))
            if(zen3) cont(istate,icp)=co !Epstein-Nesbet
            psiloc(istate)=psiloc(istate)+co*co

            if (psiloc(istate).ge.test**2)then
            write (6,*) istate,icp,vec(istate,icp),een(istate),
     $                 cmat(icp,icp)
            do ist=1,metat
	     write (6,*) ist,vec(ist,icp)
	    enddo
	    do ist=1,iconff
	     write (6,*) ist,tamp(ist,metat-1),tamp(ist,metat)
	    enddo

       endif
            enloc(istate)=enloc(istate)+co*vec(istate,icp)
cele 18-11-04
	    psisel(istate)=psisel(istate)+psiloc(istate)
	    ensel(istate)=ensel(istate)+enloc(istate)
cele 18-11-04
            co=vec(istate,icp)/(emp(istate)-enemp)
            if(zmp3) cont(istate,icp)=co !Moller-Plesset
            psimp(istate)=psimp(istate)+co*co
            e2mp(istate)=e2mp(istate)+co*vec(istate,icp)
            psimploc(istate)=psimploc(istate)+co*co
            emploc(istate)=emploc(istate)+co*vec(istate,icp)
	    dintloc(istate)=dintloc(istate)+vec(istate,icp)**2
            if(zbrd)then
               do jstate=1,istate
                  ijcou=ijcou+1
                  if(istate.ne.jstate)then
                     brdmp(istate,jstate)=brdmp(istate,jstate)
     $                    +vec(istate,icp)*vec(jstate,icp)/(emp(jstate)
     $                    -enemp)
                     brden(istate,jstate)=brden(istate,jstate)
     $                    +vec(istate,icp)*vec(jstate,icp)/(een(jstate)
     $                    -cmat(icp,icp))
                     brdmp(jstate,istate)=brdmp(jstate,istate)
     $                    +vec(istate,icp)*vec(jstate,icp)/(emp(istate)
     $                    -enemp)
                     brden(jstate,istate)=brden(jstate,istate)
     $                    +vec(istate,icp)*vec(jstate,icp)/(een(istate)
     $                    -cmat(icp,icp))
                  endif
               enddo
            endif
c	 if (kcaso.eq.8) then 
c	 write (6,*) 'Contributo totale'
c	 write (6,*) enloc(istate)*finv
c	 endif
         psienec(istate,kcaso)=psienec(istate,kcaso)+psiloc(istate)
         e2enc(istate,kcaso)=e2enc(istate,kcaso)+enloc(istate)
         psimpec(istate,kcaso)=psimpec(istate,kcaso)+psimploc(istate)
         e2mpc(istate,kcaso)=e2mpc(istate,kcaso)+emploc(istate)     
c     if (kcaso.eq.2) then 
c           write (6,*)
c           write (6,*) 'Contributo ',emploc(1)*finv
c           write (6,*) 'Energia conf',enemp*finv
c           write (6,*) 'Interazione ',sqrt(dintloc(1))*finv
c	   write(6,*) 'E tot',e2mpc(1,1)*finv
c           endif
         psien(istate)=psien(istate)+psiloc(istate)
         e2en(istate)=e2en(istate)+enloc(istate)
         if(zbrd)then
            brdmp(istate,istate)=e2mp(istate)
            brden(istate,istate)=e2en(istate)
         endif
         if(psiloc(istate).ge.tau**2)zprint=.true.
         if(psiloc(istate).ge.test**2)then
            call wrtpert(nf,istate,sqrt(psiloc(istate)),enloc(istate)
     $           *finv)
         endif
         enddo
      enddo
c            if (mcapos.eq.11.or.mcapos.eq.12)then
c	    write (6,*) 'mcapos',mcapos
c	    write (6,*) 'conf. n',icp
c	    write (6,*) 'enloc',enloc(5)               
c	    write (6,*) 'entot',e2en(5)   
c            write (6,*) 'ntotdet',ntotdet
c	    endif
      if(zthird.and.zprint)then
         mono1=nd(nf)+1
         mono2=mono1+ne2-1
         write (7) ne2,(trou(i),i=mono1,mono2),(part(i),i=mono1
     $        ,mono2),iconf2,((cont(istate,icp),istate=1,metat),icp=1
     $        ,iconf2)
      elseif(zprint.and.zpun.and.ne2.le.12)then
         mono1=nd(nf)+1
         mono2=mono1+ne2-1
cele 18-11-04
c         write (7) ne2,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
c     *        ' ',(sqrt(psiloc(i)),i=1,metat),(finv*enloc(i),i=1,metat)
         write (7) ne2,(trou(i),i=mono1,mono2),(part(i),i=mono1,mono2),
     *        ' ',(sqrt(psisel(i)),i=1,metat),(finv*ensel(i),i=1,metat)
cele 18-11-04
      endif
      return
      end
c--------------------------------------------------------
      subroutine escla2(nd,ne,trou,part,icomp,iconf,isegno,nocc,norb)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      integer*2 ne,trou(*),part(*),icomp,iconf,isegno(*),nec
      dimension nd(*),ne(*)
      equivalence (ip1,ib1),(ip2,ib2)
      nec=ne(1)
      iopenh=nec
      iopenp=nec
      do i=1,nec
         ib1=trou(i)
         if(ib1.gt.norb)ib1=ib1-norb
         do j=i+1,nec
            ib2=trou(j)
            if(ib2.gt.norb)ib2=ib2-norb
            if(ib2.eq.ib1)then
               iopenh=iopenh-2
               goto 1
            endif
         enddo
 1       continue
      enddo
      do i=1,nec
         ip1=part(i)
         if(ip1.gt.norb)ip1=ip1-norb
         do j=i+1,nec
            ip2=part(j)
            if(ip2.gt.norb)ip2=ip2-norb
            if(ip2.eq.ip1)then
               iopenp=iopenp-2
               goto 2
            endif
         enddo
 2       continue
      enddo
#if defined (SINGLET) || defined (TRIPLET)
      ntot=iopenh+iopenp
#endif
#ifdef DOUBLET
         ntot=iopenh+iopenp-1
#endif
#ifdef SINGLET
      if(ntot.eq.0)then
         icomp=1
         iconf=1
      elseif(ntot.eq.2)then
         icomp=2
         iconf=1
      elseif(ntot.eq.4)then
         icomp=6
         iconf=2
      elseif(ntot.eq.6)then
         icomp=20
         iconf=5
      elseif(ntot.eq.8)then
         icomp=70
         iconf=14
      elseif(ntot.eq.10)then
         icomp=252
         iconf=42
      elseif(ntot.eq.12)then
         icomp=924
         iconf=132
      elseif(ntot.gt.12)then
c         print*,'sorry, more than 12 electrons not covered yet',ntot
c         stop 'too many open shells'
         iconf=-1
      else
         print*,'something must be really very wrong',ntot,iopenh
     $        ,iopenp
         print*,'determinant m=',m
         print*,'trou'
         print*,(trou(i),i=1,nec)
         print*,'part'
         print*,(part(i),i=1,nec)
         stop 'very strange'
      endif
#endif
#ifdef DOUBLET
         if(ntot.eq.1)then
            icomp=1
            iconf=1
         elseif(ntot.eq.3)then
            icomp=3
            iconf=2
         elseif(ntot.eq.5)then
            icomp=10
            iconf=5
         elseif(ntot.eq.7)then
            icomp=35
            iconf=14
         elseif(ntot.eq.9)then
            icomp=126
            iconf=42
         elseif(ntot.eq.11)then
            icomp=462
            iconf=132
         elseif(ntot.gt.11)then
            print*,'sorry, more than 11 electrons not covered yet',ntot
            stop 'too many open shells'
         else
            print*,'something must be really very wrong',ntot,iopenh
     $           ,iopenp
            print*,'determinant m=',m
            print*,'trou'
            print*,(trou(nstart+i),i=1,nec)
            print*,'part'
            print*,(part(nstart+i),i=1,nec)
            stop 'very strange'
         endif
#endif
#ifdef TRIPLET
         if(ntot.eq.2)then
            icomp=1
            iconf=1
         elseif(ntot.eq.4)then
            icomp=4
            iconf=3
         elseif(ntot.eq.6)then
            icomp=15
            iconf=9
         elseif(ntot.eq.8)then
            icomp=56
            iconf=28
         elseif(ntot.eq.10)then
            icomp=210
            iconf=90
         elseif(ntot.eq.12)then
            icomp=792
            iconf=297
         elseif(ntot.gt.12)then
            print*,'sorry, more than 12 electrons not covered yet',ntot
            stop 'too many open shells'
         else
            print*,'something must be really very wrong',ntot,iopenh
     $           ,iopenp
            print*,'determinant m=',m
            print*,'trou'
            print*,(trou(nstart+i),i=1,nec)
            print*,'part'
            print*,(part(nstart+i),i=1,nec)
            stop 'very strange'
         endif
#endif
      do k=1,icomp
         call mkorb(k,nd,ne,trou,part,nocc,norb,isegno(k))
      enddo
      return
      end
C***********************************************************************
      subroutine wrtpert(kk,istate,psiloc,enloc)
c 
c     impression et perforation de determinants
c 
      implicit real*8(a-h,o-y),logical*1(z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     *     NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     *     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     *     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     *     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     *     NDEG,IORB,ISPIN,ITSYM,
     *   ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
       INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      character*60 ctrou,cpart,cctrou,ccpart
      character*1 ord(13),cspin(2),arrow(4)
      data ord/'F','M','D','T','Q','P','H','E','O','N','X','U','Z'/
      data cspin/'-','+'/
      data arrow/' ','=','=','>'/
      nec=ne(kk)
      if (nec.gt.12.) return
      do i=1,norb
      ispin(i)=0
      ispin(i+norb)=1
      iorb(i)=i
      iorb(i+norb)=i
      enddo
      ispin(2*norb+1)=1
      iorb(2*norb+1)=norb+1
      mono1=nd(kk)
      mono2=mono1+nec
      mono1=mono1+1
      ctrou=
     * '                                                            '
      cctrou=
     * '                                                            '
      cpart=
     * '                                                            '
      ccpart=
     * '                                                            '
      write (cctrou,40)
     * (iorb(trou(i)),cspin(ispin(trou(i))+1),i=mono1,mono2)
      write (ccpart,40)
     * (IORB(PART(I)),cspin(ispin(PART(I))+1),I=MONO1,MONO2)
   40 format(12(1x,i3,a1:))
      nblt=0
      nblp=0
      lt=0
      lp=0
      do 50 i=1,60
      if (cctrou(i:i).eq.' ') then
      nblt=nblt+1
      if (nblt.eq.1) then
      lt=lt+1
      ctrou(lt:lt)=' '
      endif
      else
      nblt=0
      lt=lt+1
      ctrou(lt:lt)=cctrou(i:i)
      endif
      if (ccpart(i:i).eq.' ') then
      nblp=nblp+1
      if (nblp.eq.1) then
      lp=lp+1
      cpart(lp:lp)=' '
      endif
      else
      nblp=0
      lp=lp+1
      cpart(lp:lp)=ccpart(i:i)
      endif
   50 continue
      if (ctrou(lt:lt).eq.' ') lt=lt-1
      if (cpart(lp:lp).eq.' ') lp=lp-1
      ltot=lt+lp
      if (ltot.le.64) then
c      write (99,60) ord(nec+1),(ctrou(i:i),i=1,lt),(cpart(i:i),i=1,lp),
c     * (' ',i=1,61-lt-lp),1.,1.
      write (6,60) ord(nec+1),(ctrou(i:i),i=1,lt),(cpart(i:i),i=1,lp),
     * (' ',i=1,61-lt-lp),istate,psiloc,enloc
      else
c      write (99,60) (ctrou(i:i),i=1,lt),arrow
c      write (99,60) (' ',i=1,18),(cpart(i:i),i=1,lp)
      endif
c      write (6,60) (ctrou(i:i),i=1,lt),arrow
c      write (6,60) (' ',i=1,18),(cpart(i:i),i=1,lp)
      call flush(6)   !tolgo flush
      return
   60 FORMAT (1X,62A1,1X,I2,2X,2f10.6)
c   60 format(62a1,99f15.10)
      END
C***********************************************************************
      subroutine wrtcap(kk,csq,nvec)
c 
c     impression et perforation de determinants
c 
      implicit real*8(a-h,o-y),logical*1(z)
      include 'parameter.cippi'
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     *     NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     *     INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     *     zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     *     ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     *     NDEG,IORB,ISPIN,ITSYM,
     *     ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE
     $     ,ZHEFF
       INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
       dimension csq(200)
      character*60 ctrou,cpart,cctrou,ccpart
      character*1 ord(13),cspin(2),arrow(4)
      data ord/'F','M','D','T','Q','P','H','E','O','N','X','U','Z'/
      data cspin/'-','+'/
      data arrow/' ','=','=','>'/
      nec=ne(kk)
      if (nec.gt.12.) return
      do i=1,norb
      ispin(i)=0
      ispin(i+norb)=1
      iorb(i)=i
      iorb(i+norb)=i
      enddo
      ispin(2*norb+1)=1
      iorb(2*norb+1)=norb+1
      mono1=nd(kk)
      mono2=mono1+nec
      mono1=mono1+1
      ctrou=
     * '                                                            '
      cctrou=
     * '                                                            '
      cpart=
     * '                                                            '
      ccpart=
     * '                                                            '
      write (cctrou,40)
     * (iorb(trou(i)),cspin(ispin(trou(i))+1),i=mono1,mono2)
      write (ccpart,40)
     * (IORB(PART(I)),cspin(ispin(PART(I))+1),I=MONO1,MONO2)
   40 format(12(1x,i3,a1:))
      nblt=0
      nblp=0
      lt=0
      lp=0
      do 50 i=1,60
      if (cctrou(i:i).eq.' ') then
      nblt=nblt+1
      if (nblt.eq.1) then
      lt=lt+1
      ctrou(lt:lt)=' '
      endif
      else
      nblt=0
      lt=lt+1
      ctrou(lt:lt)=cctrou(i:i)
      endif
      if (ccpart(i:i).eq.' ') then
      nblp=nblp+1
      if (nblp.eq.1) then
      lp=lp+1
      cpart(lp:lp)=' '
      endif
      else
      nblp=0
      lp=lp+1
      cpart(lp:lp)=ccpart(i:i)
      endif
   50 continue
      if (ctrou(lt:lt).eq.' ') lt=lt-1
      if (cpart(lp:lp).eq.' ') lp=lp-1
      ltot=lt+lp
      if (ltot.le.61) then
      write (6,60) ord(nec+1),(ctrou(i:i),i=1,lt),(cpart(i:i),i=1,lp),
     * (' ',i=1,61-lt-lp),(csq(i),i=1,nvec)
      else
      write (6,60) ord(nec+1),(ctrou(i:i),i=1,lt)
      write (6,60) '>',(cpart(i:i),i=1,lp),
     * (' ',i=1,61-lp),(csq(i),i=1,nvec)
      endif
      call flush(6)   !tolgo flush
      return
   60 FORMAT (62A1,1X,65f6.3)
c   60 format(62a1,99f15.10)
      END
c-------------------------------------------------------------------
      subroutine hntdc2(i,j,icomp1,jcomp1,iconf1,jconf1,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1))
      integer*2 icomp1,jcomp1,iconf1,jconf1,isegno(*)
      icomp=icomp1
      jcomp=jcomp1
      iconf=iconf1
      jconf=jconf1
#ifndef SINGLET      
      do m=1,icomp
         mm=m+i-1
         do n=1,jcomp
            nn=n+j-1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
         enddo
      enddo
#else
      do m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         do n=1,jcomp
            nn=n+j-1
            znc=jcomp.eq.1
            if(.not.znc)znf=mod(n,2).eq.1
            call givecase(zmc,znc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
            if (ncaso.eq.1)then
               hmat(m,n+1)=hmat(m,n)
            elseif (ncaso.eq.2)then
               hmat(m+1,n)=hmat(m,n)
            elseif (ncaso.eq.3)then
               hmat(m+1,n+1)=hmat(m,n)
            elseif (ncaso.eq.4)then
               hmat(m-1,n+1)=hmat(m,n)
            endif
 1       enddo
      enddo
#endif
      indi=ivpos(icomp1)
      indj=ivpos(jcomp1)
      if(iconf.eq.0)then
         ic=1
      else
         ic=iconf
      endif
      if(jconf.eq.0)then
         jc=1
      else
         jc=jconf
      endif
#ifdef SINGLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,924,buff(indj),jcomp,
     $    0.d0,tamp,924)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,924,0.d0
     $     ,cmat,132)
#endif
#ifdef DOUBLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,462,buff(indj),jcomp,
     $    0.d0,tamp,462)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,462,0.d0
     $     ,cmat,132)
#endif
#ifdef TRIPLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,792,buff(indj),jcomp,
     $    0.d0,tamp,792)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,792,0.d0
     $     ,cmat,297)
#endif
      return
      end
c------------------------------------------------------------------
      subroutine hntdc21(i,j,icomp1,jcomp1,iconf1,jconf1,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1))
      integer*2 icomp1,jcomp1,iconf1,jconf1,isegno(*)
      icomp=icomp1
      jcomp=jcomp1
      iconf=iconf1
      jconf=jconf1
      zclosecc=.false.
      if (icaso.eq.1) then
      ini=kposp4(icomp)
      istep=6
      endif
      if (icaso.eq.2) then
      ini=kposp2v(icomp)
      istep=2
      endif
      if (icaso.eq.3) then
      ini=kposp2c(icomp)
      istep=2
      endif
      if (icaso.eq.4) zclosecc=.true.
      mind=ivpos(icomp1)
      nind=ivpos(jcomp1)
      do m=1,iconf
         do n=1,jconf
            cmat(m,n)=0.0d0
         enddo
      enddo
#ifndef SINGLET      
      do  m=1,icomp
         mm=m+i-1
	 if (zclosecc) then
	    n=m
            nn=n+j-1
            ab=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)

            mi=mind-icomp+m-1
            do ik=1,iconf
               mi=mi+icomp
               buffmi=buff(mi)
               if(buffmi.eq.0.0d0) goto 21
	       ni=nind-jcomp+n-1
               do jk=1,jconf
                ni=ni+jcomp
                buffni=buff(ni)
                if(buffni.ne.0.0d0) 
     *            cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buffni
               enddo
 21         enddo

         else

	    do in=1,istep
	    n=ipos(ini+(m-1)*istep+in-1)
	    nn=n+j-1
	    ab=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
            nind=ivpos(jcomp1)

            mi=mind-icomp+m-1
            do ik=1,iconf
               mi=mi+icomp
               buffmi=buff(mi)
               if(buffmi.eq.0.0d0) goto 22
	       ni=nind-jcomp+n-1
               do jk=1,jconf
                ni=ni+jcomp
                buffni=buff(ni)
                if(buffni.ne.0.0d0) 
     *            cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buffni
               enddo
 22         enddo
            enddo

        endif

      enddo
#else
      do 2 m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         if (zclosecc) then
	    n=m
            nn=n+j-1
            znc=jcomp.eq.1
            if(.not.znc)znf=mod(n,2).eq.1
            call givecase(zmc,znc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 2
            ab=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
	    nind=ivpos(jcomp1)
            mi=mind-icomp+m-1
            do ik=1,iconf
            mi=mi+icomp
            buffmi=buff(mi)
            if(buffmi.eq.0.0d0) goto 23
	       ni=nind-jcomp+n-1
               do jk=1,jconf
               ni=ni+jcomp
               buffni=buff(ni)

               if(buffni.ne.0.0d0)then
               cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buffni
                if (ncaso.eq.1)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buff(ni+1)
                elseif (ncaso.eq.2)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi+1)*buffni
                elseif (ncaso.eq.3)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi+1)*buff(ni+1)
                elseif (ncaso.eq.4)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi-1)*buff(ni+1)
                endif
               endif

            enddo
 23         enddo

         else
	    do 1 in=1,istep
            n=ipos(ini+(m-1)*istep+in-1)
            nn=n+j-1
            znc=jcomp.eq.1
            if(.not.znc)znf=mod(n,2).eq.1
            call givecase(zmc,znc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 1
            ab=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
            nind=ivpos(jcomp1)

            mi=mind-icomp+m-1
            do ik=1,iconf
            mi=mi+icomp
            buffmi=buff(mi)
            if(buffmi.eq.0.0d0) goto 24
	       ni=nind-jcomp+n-1
               do jk=1,jconf
               ni=ni+jcomp
               buffni=buff(ni)

               if(buffni.ne.0.0d0)then
               cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buffni
                if (ncaso.eq.1)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buff(ni+1)
                elseif (ncaso.eq.2)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi+1)*buffni
                elseif (ncaso.eq.3)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi+1)*buff(ni+1)
                elseif (ncaso.eq.4)then
                cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi-1)*buff(ni+1)
                endif
               endif

            enddo
 24         enddo

 1          enddo
        endif
 2    enddo
#endif
      return
      end
c-------------------------------------------------------------------
      subroutine hntdc24(i,j,icomp1,jcomp1,iconf1,jconf1,iaux,iposj,norb
     $     ,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      dimension iaux(*),iposj(*),ivec6(6)
      common/confinfo/ivpos(924),intpos(924)
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      dimension buff(255840)
#endif
      equivalence (t0(1,1),buff(1))
      integer*2 icomp1,jcomp1,iconf1,jconf1,isegno(*)
      icomp=icomp1
      jcomp=jcomp1
      iconf=iconf1
      jconf=jconf1
      do n=1,jconf
         do m=1,iconf
            cmat(m,n)=0.0d0
         enddo
      enddo
#ifndef SINGLET      
c      write (6,*) 'Errore, caso non contemplato in hntdc24'
c      stop
      do m=1,icomp
         mm=m+i-1
         do n=1,jcomp
            nn=n+j-1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
         enddo
      enddo
#else
	 mind=ivpos(icomp1)
	 nind=ivpos(jcomp1)
      do m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
      if(icomp.eq.1)then
         nopen=0
         ii=0
      elseif(icomp.eq.2)then
         nopen=2
         ii=inum1(m)
      elseif(icomp.eq.6)then
         nopen=4
         ii=inum2(m)
      elseif(icomp.eq.20)then
         nopen=6
         ii=inum3(m)
      elseif(icomp.eq.70)then
         nopen=8
         ii=inum4(m)
      elseif(icomp.eq.252)then
         nopen=10
         ii=inum5(m)
         print*,'Warning nopen=',nopen
      elseif(icomp.eq.924)then
         nopen=12
         print*,'Warning nopen=',nopen
         ii=inum6(m)
      endif
c      print*,"chiamo sixnum con ii=",ii,' nopen=',nopen
c      call sixnum(%val(ii),iaux,%val(nopen),ivec6)
      call sixnum(ii,iaux,nopen,ivec6)
      do idu=1,6
         if(jconf.eq.2)then
            n=inv2(ivec6(idu))
         elseif(jconf.eq.5)then
            n=inv3(ivec6(idu))
         elseif(jconf.eq.14)then
            n=inv4(ivec6(idu))
         elseif(jconf.eq.42)then
            n=inv5(ivec6(idu))
         elseif(jconf.eq.132)then
            n=inv6(ivec6(idu))
         else
            print*,'Something is wrong with jconf=',jconf
            stop 'something wrong'
         endif
c         print*,'Sixnum: n=',n,' jconf=',jconf
         nn=n+j-1
         znc=jcomp.eq.1
         if(.not.znc)znf=mod(n,2).eq.1
         call givecase(zmc,znc,zmf,znf,ncaso)
         if (ncaso.lt.0)goto 1
         ab=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
         mi=mind-icomp+m-1
         do ik=1,iconf
           mi=mi+icomp
            buffmi=buff(mi)
            if(buffmi.eq.0.0d0) goto 2
            ni=nind-jcomp+n-1
            do jk=1,jconf
            ni=ni+jcomp
               buffni=buff(ni)
               if(buffni.ne.0.0d0)then
                  cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buffni
            if (ncaso.eq.1)then
               cmat(ik,jk)=cmat(ik,jk)+ab*buffmi*buff(ni+1)
            elseif (ncaso.eq.2)then
               cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi+1)*buffni
            elseif (ncaso.eq.3)then
               cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi+1)*buff(ni+1)
            elseif (ncaso.eq.4)then
               cmat(ik,jk)=cmat(ik,jk)+ab*buff(mi-1)*buff(ni+1)
            endif
         endif
      enddo
 2    enddo
 1    enddo
      enddo
#endif
      return
c      call matout(icomp,jcomp,hmat,924)
      indi=ivpos(icomp1)
      indj=ivpos(jcomp1)
      if(iconf.eq.0)then
         ic=1
      else
         ic=iconf
      endif
      if(jconf.eq.0)then
         jc=1
      else
         jc=jconf
      endif
#ifdef SINGLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,924,buff(indj),jcomp,
     $    0.d0,tamp,924)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,924,0.d0
     $     ,cmat,132)
#endif
#ifdef DOUBLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,462,buff(indj),jcomp,
     $    0.d0,tamp,462)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,462,0.d0
     $     ,cmat,132)
#endif
#ifdef TRIPLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,792,buff(indj),jcomp,
     $    0.d0,tamp,792)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,792,0.d0
     $     ,cmat,297)
#endif
      return
      end
c-------------------------------------------------------------------
      subroutine hntdc22(i,j,icomp1,jcomp1,iconf1,jconf1,iaux,iposj,norb
     $     ,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      dimension iaux(*),iposj(*),ivec2(2)
      common/confinfo/ivpos(924),intpos(924)
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      dimension buff(255840)
#endif
      equivalence (t0(1,1),buff(1))
      integer*2 icomp1,jcomp1,iconf1,jconf1,isegno(*)
      icomp=icomp1
      jcomp=jcomp1
      iconf=iconf1
      jconf=jconf1
c      print*,'Sub. hntdc22: i,j=',i,j
c      print*,'Sub. hntdc22: i-jconf,i-jcomp=',iconf,jconf,icomp,jcomp
      do n=1,jconf
         do m=1,iconf
            cmat(m,n)=0.0d0
         enddo
      enddo
#ifndef SINGLET      
c      write (6,*) 'Errore, caso non contemplato in hntdc22'
c      stop
      do m=1,icomp
         mm=m+i-1
         do n=1,jcomp
            nn=n+j-1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
         enddo
      enddo
#else
      mind=ivpos(icomp1)
      nind=ivpos(jcomp1)
      do m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         if(icomp.eq.1)then
            nopen=0
            ii=0
      elseif(icomp.eq.2)then
         nopen=2
         ii=inum1(m)
      elseif(icomp.eq.6)then
         nopen=4
         ii=inum2(m)
      elseif(icomp.eq.20)then
         nopen=6
         ii=inum3(m)
      elseif(icomp.eq.70)then
         nopen=8
         ii=inum4(m)
      elseif(icomp.eq.252)then
         nopen=10
         ii=inum5(m)
         print*,'Warning nopen=',nopen
      elseif(icomp.eq.924)then
         nopen=12
         print*,'Warning nopen=',nopen
         ii=inum6(m)
      endif
c      print*,"chiamo twonum con ii=",ii,' nopen=',nopen
c      call twonum(%val(ii),iaux,%val(nopen),ivec2)
      call twonum(ii,iaux,nopen,ivec2)
      do idu=1,2
         if(jconf.eq.1)then
            n=inv1(ivec2(idu))
         elseif(jconf.eq.2)then
            n=inv2(ivec2(idu))
         elseif(jconf.eq.5)then
            n=inv3(ivec2(idu))
         elseif(jconf.eq.14)then
            n=inv4(ivec2(idu))
         elseif(jconf.eq.42)then
            n=inv5(ivec2(idu))
         elseif(jconf.eq.132)then
            n=inv6(ivec2(idu))
         else
            print*,'Something is wrong with jconf=',jconf
            stop 'something wrong'
         endif
c         print*,'Twonum: n=',n,' jconf=',jconf
         nn=n+j-1
         znc=jcomp.eq.1
         if(.not.znc)znf=mod(n,2).eq.1
         call givecase(zmc,znc,zmf,znf,ncaso)
         if (ncaso.lt.0)goto 1
         hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
         mi=mind-icomp+m-1
         do ik=1,iconf
         mi=mi+icomp
            buffmi=buff(mi)
            if(buffmi.eq.0.0d0) goto 2
            ni=nind-jcomp+n-1
            do jk=1,jconf
            ni=ni+jcomp
               buffni=buff(ni)
               if(buffni.ne.0.0d0)then
                  cmat(ik,jk)=cmat(ik,jk)+hmat(m,n)*buffmi*buffni
                  if (ncaso.eq.1)then
                     cmat(ik,jk)=cmat(ik,jk)+hmat(m,n)*buffmi*buff(ni+1)
c     hmat(m,n+1)=hmat(m,n)
            elseif (ncaso.eq.2)then
               cmat(ik,jk)=cmat(ik,jk)+hmat(m,n)*buff(mi+1)*buffni
c               hmat(m+1,n)=hmat(m,n)
            elseif (ncaso.eq.3)then
               cmat(ik,jk)=cmat(ik,jk)+hmat(m,n)*buff(mi+1)*buff(ni+1)
c               hmat(m+1,n+1)=hmat(m,n)
            elseif (ncaso.eq.4)then
               cmat(ik,jk)=cmat(ik,jk)+hmat(m,n)*buff(mi-1)*buff(ni+1)
c               hmat(m-1,n+1)=hmat(m,n)
            endif
         endif
      enddo
 2    enddo
 1    enddo
      enddo
#endif
      return
c      call matout(icomp,jcomp,hmat,924)
      indi=ivpos(icomp1)
      indj=ivpos(jcomp1)
      if(iconf.eq.0)then
         ic=1
      else
         ic=iconf
      endif
      if(jconf.eq.0)then
         jc=1
      else
         jc=jconf
      endif
#ifdef SINGLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,924,buff(indj),jcomp,
     $    0.d0,tamp,924)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,924,0.d0
     $     ,cmat,132)
#endif
#ifdef DOUBLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,462,buff(indj),jcomp,
     $    0.d0,tamp,462)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,462,0.d0
     $     ,cmat,132)
#endif
#ifdef TRIPLET
      call dgemm('N','N',icomp,jc,jcomp,1.d0,hmat,792,buff(indj),jcomp,
     $    0.d0,tamp,792)
      call dgemm('T','N',ic,jc,icomp,1.d0,buff(indi),icomp,tamp,792,0.d0
     $     ,cmat,297)
#endif
      return
      end
c------------------------------------------------------------------
      subroutine hdigc2(i,icomp1,iconf1,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      integer*2 icomp1,iconf1,isegno(*)
      pointer hdiag,stok
      dimension hdiag(:),stok(:)
      COMMON /PER/ HDIAG,STOK
      icomp=icomp1
      iconf=iconf1
      zfa=.false.
#ifdef SINGLET
      if(icomp.eq.1)then  !closed shell
         CALL HDIG(I,ZFA,DUM,AB)
         cmat(1,1)=ab
c	 stok(mm)=ab
         return
      endif
#endif
#ifndef SINGLET
      do m=1,icomp
         mm=m+i-1
         CALL HDIG(mm,ZFA,DUM,AB)
         hmat(m,m)=ab
         do n=1,m-1
            nn=n+i-1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
            hmat(n,m)=hmat(m,n)
         enddo
      enddo
#else
      do m=1,icomp
         mm=m+i-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         CALL HDIG(mm,ZFA,DUM,AB)
         hmat(m,m)=ab
         do n=1,m-1
            nn=n+i-1
            if(.not.zmc)znf=mod(n,2).eq.1
            call givecase(zmc,zmc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 1
            hmat(m,n)=HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
            hmat(n,m)=hmat(m,n)
            if (ncaso.eq.1)then
               hmat(m,n+1)=hmat(m,n)
               hmat(n+1,m)=hmat(m,n)
            elseif (ncaso.eq.2)then
               hmat(m+1,n)=hmat(m,n)
               hmat(n,m+1)=hmat(m,n)
            elseif (ncaso.eq.3)then
               hmat(m+1,n+1)=hmat(m,n)
               hmat(n+1,m+1)=hmat(m,n)
            elseif (ncaso.eq.4)then
               hmat(m-1,n+1)=hmat(m,n)
               hmat(n+1,m-1)=hmat(m,n)
            endif
 1       enddo
      enddo
#endif
#ifdef SINGLET
      if(iconf.eq.1)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T1,2,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T1,2,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.2)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T2,6,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T2,6,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.5)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T3,20,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T3,20,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.14)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T4,70,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T4,70,tamp,924,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.42)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T5,252,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T5,252,tamp,924,0.d0
     $        ,cmat,132)
         elseif(iconf.eq.132)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,924,T6,924,0.d0,tamp
     $        ,924)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T6,924,tamp,924,0.d0
     $        ,cmat,132)
      endif
#endif
#ifdef DOUBLET
      if(iconf.eq.1)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T0,1,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T0,1,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.2)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T1,3,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T1,3,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.5)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T2,10,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T2,10,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.14)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T3,35,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T3,35,tamp,462,0.d0
     $        ,cmat,132)
      elseif(iconf.eq.42)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T4,126,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T4,126,tamp,462,0.d0
     $        ,cmat,132)
         elseif(iconf.eq.132)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,462,T5,462,0.d0,tamp
     $        ,462)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T5,462,tamp,462,0.d0
     $        ,cmat,132)
      endif
#endif
#ifdef TRIPLET
      if(iconf.eq.1)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T0,1,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T0,1,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.3)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T1,4,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T1,4,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.9)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T2,15,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T2,15,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.28)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T3,56,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T3,56,tamp,792,0.d0
     $        ,cmat,297)
      elseif(iconf.eq.90)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T4,210,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T4,210,tamp,792,0.d0
     $        ,cmat,297)
         elseif(iconf.eq.297)then
         call dsymm('L','U',icomp,iconf,1.d0,hmat,792,T5,792,0.d0,tamp
     $        ,792)
         call dgemm('T','N',iconf,iconf,icomp,1.d0,T5,792,tamp,792,0.d0
     $        ,cmat,297)
      endif
#endif
      return
      end
c------------------------------------------------------------------
      subroutine hdigc21(k,icomp1,iconf1,isegno)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
#ifdef SINGLET
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hmat(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
      common/zmat/ zint(460673)
      dimension buff(133647)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hmat(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462)
      common/zmat/ zint(115009)
      dimension buff(66823)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hmat(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792)
      common/zmat/ zint(336842)
      dimension buff(255840)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
      equivalence (t0(1,1),buff(1)) 
      integer*2 icomp1,iconf1,isegno(*)
      pointer hdiag,stok
      dimension hdiag(:),stok(:)
      COMMON /PER/ HDIAG,STOK
      icomp=icomp1
      iconf=iconf1
      do i=1,iconf
      cmat(i,i)=0.d0
      enddo
      zfa=.false.
      ini=intpos(icomp)
      mind=ivpos(icomp1)
#ifdef SINGLET
      if(icomp.eq.1)then  !closed shell
         CALL HDIG(k,ZFA,DUM,AB)
         cmat(1,1)=ab
c	 stok(mm)=ab
         return
      endif
#endif
#ifndef SINGLET
      do m=1,icomp
         mm=m+k-1
         CALL HDIG(mm,ZFA,DUM,AB)
	 mi=mind-icomp+m-1
	 do i=1,iconf
	 mi=mi+icomp
	 cmat(i,i)=cmat(i,i)+ab*buff(mi)**2
	 enddo
         do n=1,m-1
	    mn=ini+((m-2)*(m-1))/2+n-1
	    if (zint(mn)) then
            nn=n+k-1
            AB=2.d0*HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
	    mi=mind-icomp+m-1
	    ni=mind-icomp+n-1
	    do i=1,iconf
            mi=mi+icomp
            ni=ni+icomp
	    cmat(i,i)=cmat(i,i)+ab*buff(mi)*buff(ni)
	    enddo
	    endif
         enddo
      enddo
#else
      do m=1,icomp
         mm=m+k-1
         zmc=icomp.eq.1
         if(.not.zmc)zmf=mod(m,2).eq.1
         CALL HDIG(mm,ZFA,DUM,AB)
	 mi=mind-icomp+m-1
	 do i=1,iconf
	 mi=mi+icomp
	 cmat(i,i)=cmat(i,i)+ab*buff(mi)**2
	 enddo
         do n=1,m-1
	    mn=ini+((m-2)*(m-1))/2+n-1
	    if (.not.zint(mn)) goto 1
            nn=n+k-1
            if(.not.zmc)znf=mod(n,2).eq.1
            call givecase(zmc,zmc,zmf,znf,ncaso)
            if (ncaso.lt.0)goto 1
            ab=2.d0*HNTD(mm,nn,JINT,KINT)*isegno(mm)*isegno(nn)
	    mi=mind-icomp+m-1
	    ni=mind-icomp+n-1
            do 2 i=1,iconf
            mi=mi+icomp
            ni=ni+icomp
	    buffmi=buff(mi)
	    buffni=buff(ni)
	    if (buffmi.eq.0.d0.or.buffni.eq.0.d0) goto 2
            cmat(i,i)=cmat(i,i)+ab*buffmi*buffni
	    if (ncaso.eq.1)then
	     cmat(i,i)=cmat(i,i)+ab*buffmi*buff(ni+1)
            elseif (ncaso.eq.2)then
	     cmat(i,i)=cmat(i,i)+ab*buff(mi+1)*buffni
	    elseif (ncaso.eq.3)then
	     cmat(i,i)=cmat(i,i)+ab*buff(mi+1)*buff(ni+1)
            elseif (ncaso.eq.4)then
             if ((m-n).gt.1)
     *         cmat(i,i)=cmat(i,i)+ab*buff(mi-1)*buff(ni+1)
	    endif
 2          enddo
 1       enddo
      enddo
#endif
      return
      end
c------------------------------------------------------------
      subroutine givecase(zmc,znc,zmf,znf,ncaso)
      logical*1 zmc,znc,zmf,znf
      if (zmc.and.znc)then  !both closed
         ncaso=0
         return
      elseif(zmc.and..not.znc)then  !closed-not closed
         if (znf)then               ! n flippable
            ncaso=1
            return
         else
            ncaso=-1                !already dealt with
            return
         endif
      elseif((.not.zmc).and.znc)then  !not closed-close
         if(zmf)then                  ! m flippable
            ncaso=2
            return
         else
            ncaso=-1                  !already dealt with
            return
         endif
      elseif((.not.zmc).and.(.not.znc))then  !neither closed
         if(zmf.and.znf)then                 ! m and n flippable
            ncaso=3
            return
         elseif((.not.zmf).and.znf)then       ! m flipped-n flippable
            ncaso=4
            return
         else
            ncaso=-1
            return
         endif
      endif
      end
c------------------------------------------------------------------
      subroutine definfo(idimorb)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
#ifdef SINGLET
      common/zmat/ zint(460673)
#endif
#ifdef DOUBLET
      common/zmat/ zint(115009)
#endif
#ifdef TRIPLET
      common/zmat/ zint(336842)
#endif
      common/confinfo/ivpos(924),intpos(924),kposp4(924),
     * kposp2v(924),kposp2c(924),ipos(12753),icaso
c      logical*1 zalpha(id1,924)
      logical*1 zalpha
      allocatable zalpha(:,:)
      dimension iel(12),is(6,4),zalphaloc(12)
      data is/1,0,1,0,1,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0/
c
c     Costruisco ivpos
c
      jdimorb=idimorb
      if(idimorb.le.12)jdimorb=12
      allocate(zalpha(jdimorb,924))
#ifdef SINGLET
      ivpos(1)=1
      ivpos(2)=2
      ivpos(6)=4
      ivpos(20)=16
      ivpos(70)=116
      ivpos(252)=1096
      ivpos(924)=11680
#endif 
#ifdef DOUBLET
      ivpos(1)=1
      ivpos(3)=2
      ivpos(10)=8
      ivpos(35)=58
      ivpos(126)=548
      ivpos(462)=5840
#endif
#ifdef TRIPLET
      ivpos(1)=1
      ivpos(4)=2
      ivpos(15)=14
      ivpos(56)=149
      ivpos(210)=1717
      ivpos(792)=20617
#endif
c
c     Costruisco intpos
c
#ifdef SINGLET
      intpos(2)=1
      intpos(6)=2
      intpos(20)=17
      intpos(70)=207
      intpos(252)=2622
      intpos(924)=34248
#endif 
#ifdef DOUBLET
      intpos(3)=1
      intpos(10)=4
      intpos(35)=49
      intpos(126)=644
      intpos(462)=8519
#endif
#ifdef TRIPLET
      intpos(4)=1
      intpos(15)=7
      intpos(56)=122
      intpos(210)=1662
      intpos(792)=23607
#endif
c
c     Preparazione dei tre vettori kposp4(icomp),kposp2v(icomp),
c     kposp2c(icomp), che indicano la posizione di inizio 
C     dell'open-shell icomp nel vettore iposp
c     che permette di sapere quali determinanti perturbatori N
c     interagiscono con un generatore M nel caso di eccitazioni
c     doppie inattive del tipo
c
c     N = a^+_i a^+_j a_k a_l (i,j core, k,l virtuali)
c
c     ij=kposp4(icomp)+6*(m-1)+1
c     iposp4(ij+r,r=0,5)=>  restituisce i 6 determinanti che 
c                 interagiscono con m quando i.ne.j e k.ne.l
c                 (ci sono 4 elettroni spaiati in piu)
c     
c     ij=kposp2v(icomp)+2*(m-1)+1
c     iposp2v(ij+r,r=0,1) restituisce i 2 determinanti che       
c                 interagiscono con m quando i.eq.j e k.ne.l
c                 (ci sono 2 elettroni virtuali spaiati in piu)
c     
c     ij=kposp2c(icomp)+2*(m-1)+1
c     iposp2c(ij+r,r=0,1) restituisce i 2 determinanti che       
c                 interagiscono con m quando i.ne.j e k.eq.l
c                 (ci sono 2 elettroni core spaiati in piu)
c
c     Il vettore iposp sara` costruito dopo
c
#ifdef SINGLET
      kposp4(1)=1
      kposp4(2)=7
      kposp4(6)=19
      kposp4(20)=55
      kposp4(70)=175
      kposp4(252)=595 
      kposp4(924)=2107 
      kposp2v(1)=1+7651
      kposp2v(2)=3+7651
      kposp2v(6)=7+7651
      kposp2v(20)=19+7651
      kposp2v(70)=59+7651
      kposp2v(252)=199+7651
      kposp2v(924)=703+7651
      kposp2c(1)=1+10202
      kposp2c(2)=3+10202
      kposp2c(6)=7+10202
      kposp2c(20)=19+10202
      kposp2c(70)=59+10202
      kposp2c(252)=199+10202
      kposp2c(924)=703+10202
#endif 
#ifdef DOUBLET
      kposp4(1)=1
      kposp4(3)=7
      kposp4(10)=25
      kposp4(35)=85
      kposp4(126)=295
      kposp4(462)=1051
      kposp2v(1)=1+7651
      kposp2v(3)=3+7651
      kposp2v(10)=9 +7651
      kposp2v(35)=29+7651
      kposp2v(126)=99+7651
      kposp2v(426)=351+7651
      kposp2c(1)=1+10202
      kposp2c(3)=3+10202
      kposp2c(10)=9+10202
      kposp2c(35)=29+10202
      kposp2c(126)=99+10202
      kposp2c(426)=351+10202
#endif
#ifdef TRIPLET
      kposp4(1)=1
      kposp4(4)=7
      kposp4(15)=31
      kposp4(56)=121
      kposp4(210)=457
      kposp4(792)=1717
      kposp2v(1)=1+7651
      kposp2v(4)=3+7651
      kposp2v(15)=11+7651
      kposp2v(56)=41+7651
      kposp2v(210)=153+7651
      kposp2v(792)=573+7651
      kposp2c(1)=1+10202
      kposp2c(4)=3+10202
      kposp2c(15)=11+10202
      kposp2c(56)=41+10202
      kposp2c(210)=153+10202
      kposp2c(792)=573+10202
#endif
c
c     Costruisco le matrici di interazione dei determinanti
c     appartenenti ad una stessa occupazione orbitalica
c
#ifdef SINGLET
      isz=0
      zion=.false.
      do 100 k=1,6
      n=k*2
      if (k.eq.1) ihalf=1
      if (k.eq.2) ihalf=3
      if (k.eq.3) ihalf=10
      if (k.eq.4) ihalf=35
      if (k.eq.5) ihalf=126
      if (k.eq.6) ihalf=462
#endif
#ifdef DOUBLET
      isz=1
      zion=.true.
      do 100 k=1,5
      n=k*2+1
      ihalf=0
#endif
#ifdef TRIPLET
      isz=1
      zion=.false.
      do 100 k=1,5
      n=k*2+2
      ihalf=0
#endif
      ncf=0
      call genera2(ncf,n,k,isz,ihalf,zion,zalpha,jdimorb)
c      do i=1,ncf
c      write (6,'(i3,120L1)') i,(zalpha(j,i),j=1,n)
c      enddo
      ini=intpos(ncf)
      do i=2,ncf
       do j=1,i-1
        iint=0
        do l=1,n
	ij=ini+((i-2)*(i-1))/2+j-1
        if (zalpha(l,i).neqv.zalpha(l,j)) iint=iint+1
        if (iint.gt.4) goto 2
        enddo
   2    if (iint.le.4) then
        zint(ij)=.true.
        else
        zint(ij)=.false.
        endif
       enddo
      enddo
c
c     Costruisco il vettore ipos
c     ATTENZIONE: NON VALIDO PER ZION=TRUE
c
#ifdef SINGLET
      if (k.gt.4) goto 100
#endif
#ifdef DOUBLET
      if (k.gt.3) goto 100
#endif
#ifdef TRIPLET
      if (k.gt.3) goto 100
#endif
      do idet=1,ncf
c	write (6,'(''zalpha'',i3,2x,50L1)') idet,(zalpha(l,idet),l=1,n)
	call flush(6)   !tolgo flush
c
       ini=kposp4(ncf)+(idet-1)*6-1
c       write (6,*) 'ncf,idet,ini',ncf,idet,ini
       call flush(6)   !tolgo flush
       do ij=1,n
       zalphaloc(2+ij)=zalpha(ij,idet)
       enddo
       do j=1,6
        zalphaloc(1)=is(j,1).eq.1
        zalphaloc(2)=is(j,2).eq.1
        zalphaloc(n+3)=is(j,3).eq.1
        zalphaloc(n+4)=is(j,4).eq.1
	nloc=n+4
#ifdef SINGLET
	kloc=k+2
#endif
#ifdef DOUBLET
c        write (6,*) 'Attenzione, caso non contemplato'
c        stop
	kloc=k+3 !??????
#endif
#ifdef TRIPLET
	kloc=k+4
#endif
	kl=0
c        write (6,'(''zalphaloc'',2x,50L1)') (zalphaloc(l),l=1,n+4)
	do ij=1,nloc
	if (zalphaloc(ij)) then
	 kl=kl+1
	 iel(kl)=ij
        endif
	enddo
c	write (6,*) (iel(lj),lj=1,kloc)
	call givepos(iel,nloc,kloc,npos)
c       write (6,*) 'j,npos',j,npos
       call flush(6)   !tolgo flush
	ipos(ini+j)=npos
       enddo
c
        ini=kposp2v(ncf)+(idet-1)*2-1
	nloc=n+2
#ifdef SINGLET
	kloc=k+1
#endif
#ifdef DOUBLET
c        write (6,*) 'Attenzione, caso non contemplato'
c        stop
	kloc=k+2 !??????
#endif
#ifdef TRIPLET
	kloc=k+3
#endif
	kl=0
	do ij=1,n
	if (zalpha(ij,idet)) then
	 kl=kl+1
	 iel(kl)=ij
        endif
	enddo
	iel(kloc)=n+1
c	write (6,*) (iel(lj),lj=1,kloc)
	call givepos(iel,nloc,kloc,npos)
c       write (6,*) 'npos',1,npos
       call flush(6)   !tolgo flush
	ipos(ini+1)=npos
	iel(kloc)=n+2
c	write (6,*) (iel(lj),lj=1,kloc)
	call givepos(iel,nloc,kloc,npos)
c       write (6,*) 'npos',2,npos
       call flush(6)   !tolgo flush
	ipos(ini+2)=npos

        ini=kposp2c(ncf)+(idet-1)*2-1
	kl=1
	do ij=1,n
	if (zalpha(ij,idet)) then
	 kl=kl+1
	 iel(kl)=ij+2
        endif
	enddo
	iel(1)=1
c	write (6,*) (iel(lj),lj=1,kloc)
	call givepos(iel,nloc,kloc,npos)
c       write (6,*) 'npos',1,npos
       call flush(6)   !tolgo flush
	ipos(ini+1)=npos
	iel(1)=2
c	write (6,*) (iel(lj),lj=1,kloc)
	call givepos(iel,nloc,kloc,npos)
c       write (6,*) 'npos',2,npos
       call flush(6)   !tolgo flush
	ipos(ini+2)=npos
c
      enddo
c
  100 enddo
c
c     Preparo ipos anche per i casi non contemplati sopra
c     (k=0)
c 
      ini=kposp4(1)-1
      do i=1,6
      ipos(ini+i)=i
      enddo
      ini=kposp2v(1)-1
      do i=1,2
      ipos(ini+i)=i
      enddo
      ini=kposp2c(1)-1
      do i=1,2
      ipos(ini+i)=i
      enddo
      deallocate(zalpha)
      return
      end
C##############################################################
      subroutine genera2(ncf,n,k,isz,ihalf,zion,zalpha,idimorb)
C
C     Calcola come disporre k oggetti in n caselle            
C     ncf=numero di disposizioni ottenute
C
      implicit none
      integer id1,id2,id3,id4,id5,id6,id7,id8,id9,id10,id11,id12,id13
     $     ,id14,id15,idimorb
      include 'parameter.cippi'
      integer*2 mspin
      integer*4 i,n,nj,j,k,nec,nel,iel,npos,kpos,kposold,
     *  ihalf,ih,nec1,nd,ncf,isz,kb,metat,ini,ijk,ndistr,
     *  ncasel,nogg,l,kl,isto
      logical*1 zalpha(idimorb,924)
      logical*1 zion
c      dimension mspin(id1),iel(id1)
      allocatable mspin(:),iel(:)
      allocate(mspin(idimorb))
      allocate(iel(idimorb))

      ih=0
      nel=1
      
      do i=1,n
      mspin(i)=0
      enddo
      iel(nel)=1
      mspin(1)=1
c
c     Se nel<k devo ancora sistemare degli elettroni.
c
   10 if (nel.lt.k) then
c
c     Se sto lavorando con il primo elettrone e dalla sua
c     posizione (iel(1)) alla fine (n) ci sono meno orbitali
c     del numero di elettroni che dovro sistemare (k-1)
c     allora ho finito.
c
       if (nel.eq.1.and.(n-iel(nel)).lt.(k-1)) return
c
c     Se il numero di elettroni ancora da sistemare (k-nel)
c     e` maggiore delle posizioni disponibili vado a 300,
c     altrimenti a 100
c
       if ((k-nel).gt.(n-iel(nel))) then
       goto 300
       else
       goto 100
       endif
      else
c 
c     Nel caso che nel = k, l'ultimo elettrone lo metto
c     negli orbitali che stanno da iel(nel-1) (posizione
c     dove ho messo il precedente elettrone) alla fine (n).
c
       if (k.eq.1) then
       ini=0
       else
       ini=iel(nel-1)
       endif
       do i=ini+1,n
       mspin(i)=1
       iel(nel)=i
c
c      Scrivo il determinante
c
C       Se ho gia` fatto meta` dei determinanti nel caso isz=0, ritorno
c
        if (isz.eq.0.d0.and..not.zion) then
         ih=ih+1
         if (ih.gt.ihalf) return
        endif
c
c        Costruzione di zalpha(i,j)=true se l'orbitale i nel det j e` alpha
c
        do isto=1,n
        zalpha(isto,ncf+1)=(mspin(isto).eq.1)
        enddo
        if (isz.eq.0.d0.and..not.zion) then
         do isto=1,n
         zalpha(isto,ncf+2)=(mspin(isto).eq.0)
         enddo
	 ncf=ncf+1
        endif
	ncf=ncf+1
       mspin(i)=0
       enddo
       if (nel.eq.1) return
       goto 200
      endif
c
c     Sistemo l'elettrone (nel+1)-esimo nella posizione 
c     iel(nel)+1 e torno a 10
c
  100 nel=nel+1
      iel(nel)=iel(nel-1)+1
      mspin(iel(nel))=1
      goto 10
c
c     Se arrivo qui, significa che ho gia` sistemato in tutti
c     i modi possibili gli elettroni da nel a k  per una certa
c     configurazione dei precedenti nel elettroni. Allora 
c     sposto il (nel-1)-esimo elettrone in avanti di una 
c     posizione.
c
  200 nel=nel-1
      mspin(iel(nel))=0
      iel(nel)=iel(nel)+1
      mspin(iel(nel))=1
      goto 10
c
c     Rimetto mspin beta a tutte le posizioni da iel(nel) a n
c     e vado a 200
c
  300 do i=iel(nel),n
      mspin(i)=0  
      enddo
      goto 200
c
      end
C-------------------------------------
       subroutine givepos(iel,n,k,npos)
c
c      n=numero di orbitali singoli occupati
c      k=numero di elettroni spaiati alpha
c      iel(i),i=1,k = posizione dell'i-esimo elettrone 
c                     nelle n caselle
c
       integer iel
       dimension iel(*)
       logical*1 zsz0

       zsz0=n.eq.k*2

       if (zsz0) then
       ntot=1
       do i=k+1,n
       ntot=ntot*i
       enddo
       do i=1,k
       ntot=ntot/i
       enddo
       endif

       if (k.eq.1) then
       npos=iel(1)
       goto 110
       endif

       npos=1
       kpos=0
       kposold=0
       do 111 ijk=1,k-1
       kpos=iel(ijk)-ijk-kposold
       if (kpos.eq.0) goto 111
       nogg=k-ijk
       do 112 kl=ijk+kposold,iel(ijk)-1
       ncasel=n-kl
       ndistr=1
       do l=nogg+1,ncasel
       ndistr=ndistr*l
       enddo
       do l=1,ncasel-nogg
       ndistr=ndistr/l
       enddo
       npos=npos+ndistr
  112  enddo
       kposold=kposold+kpos
  111  enddo
       npos=npos+iel(k)-iel(k-1)-1


  110  if (zsz0) then
        if (npos*2.le.ntot) then
	 npos=npos*2-1
        else
	 npos=(ntot-npos)*2+2	 
        endif
       endif
 
       return
       end
c----------------------------------------------------------------
      subroutine kotani(ispin2)
      implicit real*8(a-h,o-z)
c-----Questo programma produce le funzioni di spin tramite
c-----la formula ricorrente di Kotani
c---- N e' il numero di elettroni spaiati
c---- ispin2 e' il doppio del numero quantico di spin S
c---- vec e' il vettore che contiene i coefficienti delle funzioni di
c---- spin nella base prodotto diretto (alpha,beta)**n
c---- L'ordinamento della base e' lessicografico: ad esempio per n=2 si
c---- ha: aaaa, aaab, aaba, aabb, abaa, abab, abba, abbb, baaa,
c----     baab, baba, babb, bbaa, bbab, bbba, bbbb
c---- Si noti che se a corrispone a 1 e b a 0 tale ordinamento
c---- corrisponde a 15,14,13,12,...,0
c---- Le dimensioni di vec sono 2**n * g dove g e' dato dalla
c---- formula di Wigner
c---- NTRUE da' la dimensione della sottobase della base totale in
c---- cui la funzione di spin e' espansa
c---- Cosi', se n=2 si ha NTRUE=6, corrispondente a :
c---- aabb,abab,abba,baab,baba,bbaa
c---- Il vettore inum (dim. NTRUE) fornisce la posizione dell'i-esimo
c---- elemento della sottobase all'interno della base totale 2**n
c---- I valori inum sono compresi tra 0 e 2**n -1
c      dimension vec(4096,132)
c      dimension inum(4096)
#ifdef SINGLET      
      common/tmat/ T0(1,1),T1(2,1),T2(6,2),T3(20,5),T4(70,14),T5(252,42)
     $     ,T6(924,132),cmat(132,132),hma(924,924),tamp(924,132),
     $     INUM0(1),INUM1(2),INUM2(6),INUM3(20),INUM4(70),INUM5(252)
     $     ,INUM6(924),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
#endif
#ifdef DOUBLET
      common/tmat/ T0(1,1),T1(3,2),T2(10,5),T3(35,14),T4(126,42),T5(462
     $     ,132),cmat(132,132),hma(462,462),tamp(462,132),
     $     INUM0(1),INUM1(3),INUM2(10),INUM3(35),INUM4(126),INUM5(462),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
#endif
#ifdef TRIPLET
      common/tmat/ T0(1,1),T1(4,3),T2(15,9),T3(56,28),T4(210,90),T5(792
     $     ,297),cmat(297,297),hma(792,792),tamp(792,297),
     $     INUM0(1),INUM1(4),INUM2(15),INUM3(56),INUM4(210),INUM5(792),
     $     INV0(1),INV1(4),INV2(16),INV3(64),INV4(256),INV5(1024)
     $     ,INV6(4096)
#endif
c      dimension vec(1),inum(1)
      allocatable vec(:),inum(:)
      integer wigner
c      ispin2=0
#if defined (SINGLET) || defined (TRIPLET)
      do k=1,6
         n=k*2
#endif
#ifdef DOUBLET
         do k=1,6
            n=k*2 - 1
#endif
c      print*,'Dammi n e ispin2'
c      read*,n,ispin2
      ndim=2**n
      ig=wigner(n,ispin2)
c      print*,'Numero di elettroni ',n !renzo
c      print*,'Spin*2 ',ispin2         !renzo
c      print*,'Dimensione g ',ig       !renzo
      ntot=ndim*ig
c      ioffv=iallor8(vec,ntot)
c      ioffi=ialloi4(inum,ndim)
      allocate(vec(ntot))
      allocate(inum(ndim))
c      call F(%val(n),%val(ispin2),vec(ioffv))
      call F(n,ispin2,vec)
c      print '(5f10.6)',(vec(ioffv+i-1),i=1,ntot)
c      call F(%val(n),%val(ispin2),vec)
      call esame(ig,n,ndim,ispin2,ntrue,inum)
c      call esame(%val(ig),%val(n),%val(ndim),%val(ispin2),ntrue
c     $     ,inum(ioffi))
c      call esame(%val(ig),%val(n),%val(ndim),%val(ispin2),ntrue
c     $     ,inum)
c      print*,'ntrue is ',ntrue
c      print*,'Matrice vec'
c      call stampa(%val(ig),%val(n),%val(ndim),%val(ispin2),vec(ioffv)
c     $     ,inum(ioffi),%val(ntrue))
#ifdef SINGLET
      if(k.eq.1)then
         call transfer(vec,t1,inum,inum1,ndim,ig,ntrue)
c         call transfer(vec,t1,inum,inum1,ndim,ig,ntrue)
         call cngsign(t1,inum1,ig,ntrue,k*2,ndim,inv1)
c         call matout(ntrue,ig,t1,ntrue)
c         do iii=1,ntrue
c            print '(6f15.8)',(t1(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.2)then
         call transfer(vec,t2,inum,inum2,ndim,ig,ntrue)
c         call transfer(vec,t2,inum,inum2,ndim,ig,ntrue)
         call transpose(t2,inum2,ig,ntrue)
         call cngsign(t2,inum2,ig,ntrue,k*2,ndim,inv2)
c         do iii=1,ntrue
c            print '(6f15.8)',(t2(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.3)then
         call transfer(vec,t3,inum,inum3,ndim,ig,ntrue)
c         call transfer(vec,t3,inum,inum3,ndim,ig,ntrue)
         call transpose(t3,inum3,ig,ntrue)
         call cngsign(t3,inum3,ig,ntrue,k*2,ndim,inv3)
c         do iii=1,ntrue
c            print '(6f15.8)',(t3(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.4)then
         call transfer(vec,t4,inum,inum4,ndim,ig,ntrue)
c         call transfer(vec,t4,inum,inum4,ndim,ig,ntrue)
         call transpose(t4,inum4,ig,ntrue)
         call cngsign(t4,inum4,ig,ntrue,k*2,ndim,inv4)
      elseif(k.eq.5)then
         call transfer(vec,t5,inum,inum5,ndim,ig,ntrue)
c         call transfer(vec,t5,inum,inum5,ndim,ig,ntrue)
         call transpose(t5,inum5,ig,ntrue)
         call cngsign(t5,inum5,ig,ntrue,k*2,ndim,inv5)
      elseif(k.eq.6)then
         call transfer(vec,t6,inum,inum6,ndim,ig,ntrue)
c         call transfer(vec,t6,inum,inum6,ndim,ig,ntrue)
         call transpose(t6,inum6,ig,ntrue)
         call cngsign(t6,inum6,ig,ntrue,k*2,ndim,inv6)
      endif
#endif
#ifdef DOUBLET
      if(k.eq.1)then
         call transfer(vec,t0,inum,inum0,ndim,ig,ntrue)
         call cngsign(t0,inum0,ig,ntrue,n,ndim,inv1)
c         call matout(ntrue,ig,t0,ntrue)
c         do iii=1,ntrue
c            print '(6f15.8)',(t0(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.2)then
         call transfer(vec,t1,inum,inum1,ndim,ig,ntrue)
         call cngsign(t1,inum1,ig,ntrue,n,ndim,inv2)
c         call matout(ntrue,ig,t1,ntrue)
c         do iii=1,ntrue
c            print '(6f15.8)',(t1(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.3)then
         call transfer(vec,t2,inum,inum2,ndim,ig,ntrue)
         call cngsign(t2,inum2,ig,ntrue,n,ndim,inv3)
c         call matout(ntrue,ig,t2,ntrue)
c         do iii=1,ntrue
c            print '(6f15.8)',(t2(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.4)then
         call transfer(vec,t3,inum,inum3,ndim,ig,ntrue)
         call cngsign(t3,inum3,ig,ntrue,n,ndim,inv4)
c         call matout(ntrue,ig,t3,ntrue)
      elseif(k.eq.5)then
         call transfer(vec,t4,inum,inum4,ndim,ig,ntrue)
         call cngsign(t4,inum4,ig,ntrue,n,ndim,inv5)
      elseif(k.eq.6)then
         call transfer(vec,t5,inum,inum5,ndim,ig,ntrue)
         call cngsign(t5,inum5,ig,ntrue,n,ndim,inv6)
      endif
#endif
#ifdef TRIPLET
      if(k.eq.1)then
         call transfer(vec,t0,inum,inum0,ndim,ig,ntrue)
         call cngsign(t0,inum0,ig,ntrue,k*2,ndim,inv1)
c         do iii=1,ntrue
c            print '(6f15.8)',(t0(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.2)then
         call transfer(vec,t1,inum,inum1,ndim,ig,ntrue)
         call cngsign(t1,inum1,ig,ntrue,k*2,ndim,inv2)
c         do iii=1,ntrue
c            print '(6f15.8)',(t1(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.3)then
         call transfer(vec,t2,inum,inum2,ndim,ig,ntrue)
         call cngsign(t2,inum2,ig,ntrue,k*2,ndim,inv3)
c         do iii=1,ntrue
c            print '(6f15.8)',(t2(iii,jjj),jjj=1,ig)
c         enddo
      elseif(k.eq.4)then
         call transfer(vec,t3,inum,inum3,ndim,ig,ntrue)
         call cngsign(t3,inum3,ig,ntrue,k*2,ndim,inv4)
      elseif(k.eq.5)then
         call transfer(vec,t4,inum,inum4,ndim,ig,ntrue)
         call cngsign(t4,inum4,ig,ntrue,k*2,ndim,inv5)
      elseif(k.eq.6)then
         call transfer(vec,t5,inum,inum5,ndim,ig,ntrue)
         call cngsign(t5,inum5,ig,ntrue,k*2,ndim,inv6)
      endif
#endif
c      call memrel8(vec(ioffv))
c      call memrel4(inum(ioffi))
      deallocate(vec)
      deallocate(inum)
      enddo
      return
      end
c------------------------------------------------------------------
      subroutine transfer(vec,t,inum,inu,ndim,ig,ntrue)
      implicit real*8 (a-h,o-z)
      dimension vec(ndim,ig),t(ntrue,ig),inum(*),inu(*)
      do i=1,ntrue
         do j=1,ig
            inn=inum(i)+1
            t(i,j)=vec(inn,j)
            if(abs(t(i,j)).le.1.d-9)t(i,j)=0.d0
         enddo
         inu(i)=inum(i)
      enddo
      return
      end
c------------------------------------------------------------------
      subroutine transpose(t,inum,ig,ntrue)
      implicit real*8 (a-h,o-z)
c--mette i flippati a coppie
      dimension t(ntrue,ig),inum(*),tamp1(462),tamp2(462)
      dimension ita1(462),ita2(462)
      do j=1,ig
         do i=1,ntrue/2
            tamp1(i)=t(i,j)
            tamp2(i)=t(ntrue-i+1,j)
         enddo
         do i=1,ntrue/2
            i1=2*i-1
            i2=2*i
            t(i1,j)=tamp1(i)
            t(i2,j)=tamp2(i)
         enddo
      enddo
      do i=1,ntrue/2
         ita1(i)=inum(i)
         ita2(i)=inum(ntrue-i+1)
      enddo
      do i=1,ntrue/2
         i1=2*i-1
         i2=2*i
         inum(i1)=ita1(i)
         inum(i2)=ita2(i)
      enddo
      return
      end
c------------------------------------------------------------------
      subroutine cngsign(t,inum,ig,ntrue,nopen,ndim,inv)
      implicit real*8 (a-h,o-z)
      logical*1 zsign
c--riordina la sequenza ababbb... in modo che b preceda a
      dimension t(ntrue,ig),inum(*),tamp1(462),tamp2(462)
      dimension ita1(462),ita2(462)
      dimension iaux(500),ipos(500)
      dimension inv(*)
c      print*,'ntrue=',ntrue
      do i=1,ntrue
         inum(i)=ndim-inum(i)-1
         ii=inum(i)
         call getiocc(nopen,ii,iaux)
c         call getiocc(%val(nopen),%val(ii),iaux)
c         if(ntrue.eq.4)print*,'inum is ',inum(i),' ii is ',ii
c         if(ntrue.eq.4)print '(a,10i1)','Iaux=',(iaux(iiii),iiii=1,nopen
c     $        )
c         if(ntrue.eq.15)print*,'inum is ',inum(i),' ii is ',ii
c         if(ntrue.eq.15)print '(a,10i1)','Iaux=',(iaux(iiii),iiii=1
c     $        ,nopen)
         jj=0
         do j=1,nopen
            if(iaux(j).eq.1)then
               jj=jj+1
               ipos(jj)=j
            endif
         enddo
         nu=0
#ifdef SINGLET
         if(mod(nopen,4).eq.0)then
            do j=1,jj,2
               nu=nu+ipos(j+1)-ipos(j)-1
            enddo
         else
            do j=1,jj-1,2
               nu=nu+ipos(j+1)-ipos(j)-1
            enddo
            nu=nu+nopen-ipos(jj)
         endif
#endif
#ifdef DOUBLET
         nalfa=(nopen+1)/2
c         print*,'nalfa,jj=',nalfa,jj
         if(mod(nalfa,2).eq.0)then
            do j=1,jj,2
               nu=nu+ipos(j+1)-ipos(j)-1
            enddo
         else
            do j=1,jj-1,2
               nu=nu+ipos(j+1)-ipos(j)-1
            enddo
            nu=nu+nopen-ipos(jj)
         endif
#endif
#ifdef TRIPLET
         nalfa=(nopen+2)/2
c         print*,'nalfa,jj=',nalfa,jj
         if(mod(nalfa,2).eq.0)then
            do j=1,jj,2
               nu=nu+ipos(j+1)-ipos(j)-1
            enddo
         else
            do j=1,jj-1,2
               nu=nu+ipos(j+1)-ipos(j)-1
            enddo
            nu=nu+nopen-ipos(jj)
         endif
#endif
         zsign=mod(nu,2).eq.1
         do j=1,ig
            if(zsign)t(i,j)=-t(i,j)
         enddo
         inv(inum(i))=i
c         print*,'inv,inum',i,inum(i)
c         call printinum(%val(inum(i)),%val(nopen))
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine giveopen(iocc,norb,iposi,nopeni)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      dimension iocc(*),iposi(*)
      integer*1 iocc
      nopeni=0
      do i=1,norb
         if(iocc(i).eq.1)then
            nopeni=nopeni+1
            iposi(nopeni)=i
         endif
      enddo
      return
      end
c----------------------------------------------------------
      subroutine seocc(iocc,jocc,norb,iaux,nopeni,nopenj,iposi,iposj
     $     )
      implicit real*8(A-H,O-Y),logical*1(Z)
      dimension iocc(*),jocc(*),iaux(*),iposi(*),iposj(*)
      integer*1 iocc,jocc
c      ii=0
c      do j=1,nopenj
c         if(iposj(j).ne.iposi(j-ii))then
c            ii=ii+1
c            iaux(ii)=j
c         endif
c      enddo
      jj=0
      ii=1
      do j=1,nopenj
         if(iposj(j).ne.iposi(ii))then
            jj=jj+1
            iaux(jj)=j
         else
            if(ii.lt.nopeni)ii=ii+1
         endif
      enddo
      return
      end
c----------------------------------------------------------
      logical*1 function zsubset(iposi,nopeni,iposj,nopenj)
      dimension iposi(*),iposj(*)
      do i=1,nopeni
         do j=1,nopenj
            if(iposi(i).eq.iposj(j))goto 1
         enddo
         zsubset=.false.
         return
 1       continue
      enddo
      zsubset=.true.
      return
      end
c------------------------------------------------------------
      SUBROUTINE REIJKL(factor,norb,nsym,itsym,its,nos,jd1)
C 
C     PREPARATION POUR LE STOCKAGE DES INTEGRALES BIELECTONIQUES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      DIMENSION NOS(20),ITSYM(*),ITS(20,20)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      common/caratte/cara1,cara2
      character*1 cara1,cara2
      dimension norbsym(8)
      common /toul/ ztoul
cele 07/04/04      cara1=char(-1)
      cara2=char(0)
      print*,'REIJKL compiled with intel ifc compiler'
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
      NIJKL=NDEB(NSYM)+(NOS(NSYM)*(NOS(NSYM)+1))/2
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
c_toul
               if (ztoul) then
               call intmaxt(bigint,norb,nsym,nos)
               else
               call intmax(bigint,norb,nsym,nos)
               endif
c_toul
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
         call flush(6)   !tolgo flush
         if (zfac) then
c_toul
            if (ztoul) then
            call intmaxt(bigint,norb,nsym,nos)
            else
            call intmax(bigint,norb,nsym,nos)
            endif
c_toul
            print'(''     MAXIMUM INTEGRAL IS '',D10.4)',bigint
            call flush(6)   !tolgo flush
            factor=2147483646.d0/bigint
         endif
         print'(''     FACTOR WILL BE '',D10.4)',factor
         call flush(6)   !tolgo flush
      endif
   74 FORMAT (//5X,'DIMENSION DU TABLEAU IJKL:',I10,' *',i1,' BYTES')
      RETURN
      END
C***********************************************************************
      SUBROUTINE intmax(bigint,norb,nsym,nos)
C 
C     trova l'integrale piu' grande
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      PARAMETER(MBUF=3255,NBUF=4096)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      bigint=0.d0
      k22=id2*2
c      if (nbuf.gt.k22) then
c      write (6,*) '   *** ERREUR DIMENSIONS DANS INTMAX ***'
c      stop
c      endif
      WRITE (6,1000)
 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      call filesplit('REWIND',50,0,idum,idum2)
      REWIND 50
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
C***********************************************************************
      SUBROUTINE intmaxt(bigint,norb,nsym,nos)
C     
C     trova l'integrale piu' grande
C     
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      PARAMETER(MBUF=3255,NBUF=4096)
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     *     F,AJ,AK,NUM,NDEB(500),NBO(99),
     *     NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      character*1 chin1,chin2,chin3,chin4,chin5
c      dimension chin1(4096),chin2(4096),chin3(4096),chin4(4096)
c     $     ,chin5(4096)
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
      bigint=0.d0
      k22=id2*2
c     if (nbuf.gt.k22) then
c     write (6,*) '   *** ERREUR DIMENSIONS DANS INTMAX ***'
c     stop
c     endif
      WRITE (6,1000)
 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      call filesplit('REWIND',50,0,idum,idum2)
      REWIND 50
 20   continue
c      if(norb.le.255)then
c         READ(50)chin1,chin2,chin3,chin4,xxm
c      else
c         READ(50)chin1,chin2,chin3,chin4,xxm,chin5
         call filesplit('READ',50,16*nbuf,ichange,idum2)
         read(50)ind,jnd,knd,lnd,xxm
c      endif
      do n=1,nbuf
c       if(norb.le.255)then
c          lab=ichar(chin1(N))
       lab=knd(n)
          if (lab.eq.0) return
c       else
c          i=ichar(chin1(n))
c          j=ichar(chin2(n))
c          k=ichar(chin3(n))
c          l=ichar(chin4(n))
c          if (i.lt.0) i=i+256
c          if (j.lt.0) j=j+256
c          if (k.lt.0) k=k+256
c          if (l.lt.0) l=l+256
c          indcorr=ichar(chin5(n))
c          if(indcorr.lt.0)indcorr=indcorr+256
c          lcorr=indcorr/64
c          kdum=mod(indcorr,64)
c          kcorr=kdum/16
c          kdum=mod(kdum,16)
c          jcorr=kdum/4
c          icorr=mod(kdum,4)
c          i=i+icorr*256
c          j=j+jcorr*256
c          k=k+kcorr*256
c          l=l+lcorr*256
c          if(k.eq.0)return
c       endif
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
      include 'parameter.cippi'
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
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
      idummy=0
      k22=id2*2
c      if (mbuf.gt.k22) then
c      write (6,*) '   *** ERREUR DIMENSIONS DANS ITIJKL ***'
c      stop
c      endif
      WRITE (6,1000)
 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      call getenv('HOSTTYPE',host)
      write (6,'(a,a)') '     HOSTTYPE=',host
      write (6,*) '    INTEGRALI DA HONDO + IJKLPI'
      do 500 k=1,ni2
 500     ijkl2(k)=0
      do 501 k=1,ni4
 501     ijkl4(k)=0
      do k=1,ni3
         ijkl3(k)=cara(2:4)
      enddo
      khan=0
      call filesplit('REWIND',50,0,idum,idum2)
      REWIND 50
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
C***********************************************************************
      SUBROUTINE ITIJKLt(factor,norb,nsym,nos,nijkl,ijkl4,
     * ijkl2,rijkl,ijkl3)
C 
C     LECTURE DES INTEGRALES BIELECTRONIQUES
C 
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
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
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      DIMENSION NOS(20),IPAIR(500),XX(MBUF),IX(MBUF),XXm(NBUF),IXm(NBUF)
c_toul
c      character*1 chin1,chin2,chin3,chin4,chin5
c      dimension  chin1(nbuf),chin2(nbuf),chin3(nbuf),chin4(nbuf)
c     $     ,chin5(nbuf)
      integer*2 ind(4096),jnd(4096),knd(4096),lnd(4096)
c_toul
      idummy=0
      k22=id2*2
c      if (mbuf.gt.k22) then
c      write (6,*) '   *** ERREUR DIMENSIONS DANS ITIJKL ***'
c      stop
c      endif
c      if (norb.gt.255) then
c      write (6,*) 'The program does not work with norb>255 and ZTOUL=T'
c      stop
c      endif
      WRITE (6,1000)
 1000 FORMAT (5X,'LECTURE DES INTEGRALES')
      call getenv('HOSTTYPE',host)
      write (6,'(a,a)') '     HOSTTYPE=',host
      write (6,*) '    INTEGRALI LETTI CON IL FORMALISMO DI TOLOSA'
      do 500 k=1,ni2
 500     ijkl2(k)=0
      do 501 k=1,ni4
 501     ijkl4(k)=0
      do k=1,ni3
         ijkl3(k)=cara(2:4)
      enddo
      khan=0
      call filesplit('REWIND',50,0,idum,idum2)
      REWIND 50
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
         call filesplit('READ',50,16*nbuf,ichange,idum2)
         read(50)ind,jnd,knd,lnd,xxm  !renzo new
c      endif
      do indn=1,nbuf
c      i=ichar(chin1(indn))
c      j=ichar(chin2(indn))
c      k=ichar(chin3(indn))
c      l=ichar(chin4(indn))
       i=ind(indn)
       j=jnd(indn)
       k=knd(indn)
       l=lnd(indn)
       if(k.eq.0)return
c      if (i.lt.0) i=i+256
c      if (j.lt.0) j=j+256
c      if (k.lt.0) k=k+256
c      if (l.lt.0) l=l+256
c      if (k.eq.0.and.norb.le.255) return
c      if(norb.gt.255)then
c         indcorr=ichar(chin5(indn))
c         if(indcorr.lt.0)indcorr=indcorr+256
c         lcorr=indcorr/64
c         kdum=mod(indcorr,64)
c         kcorr=kdum/16
c         kdum=mod(kdum,16)
c         jcorr=kdum/4
c         icorr=mod(kdum,4)
c         i=i+icorr*256
c         j=j+jcorr*256
c         k=k+kcorr*256
c         l=l+lcorr*256
c         if(k.eq.0)return
c      endif
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
      print * ,'i,j,k,l',i,j,k,l
      print * ,'ij,kl',ij,kl
      print * ,'Attenzione errore in ITIJKL'
      print * ,'KLS,IJS',KLS,IJS,XXm(indn)
      endif
      KLP=JNDIC(KL)
      N1=MAX0(IJP,KLP)
      N2=MIN0(IJP,KLP)
      NA=NUM(N1)+N2+NDEB(KLS)
      A8=XXm(indn)*factor
c      if (i.eq.j.and.i.eq.k.and.k.eq.l) then
c      print * ,'i,j,k,l',i,j,k,l,XXm(indn)
c      endif
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
      FUNCTION AI(I,J,K,L)
C
C     REPERAGE DE L'INTEGRALE BIELECTRONIQUE IJKL DANS LE TABLEAU IJKL2
C
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      COMMON /HND/ ZHANDY,ZTHREE,zfac
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
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
      if (zhandy) then
         if(.not.zthree)then
      AI=PICK2(nnjk,ijkl2,rijkl)
      else
         ai=pick3(nnjk,ijkl3)
         endif
      else
      AI=PICK4(nnjk,ijkl4)
      endif
      IF (ZSIG) AI=-AI
      RETURN
  300 AI=0.D0
      RETURN
      END
C**********************************************************************
      subroutine givdimen(ncfd,ndim,idimorb)
      IMPLICIT REAL*8(A-H,O-Y),LOGICAL*1(Z)
      include 'parameter.cippi'
      character*2000 namel
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
     * ZHMAT,zrap,zatom,zspin,zdialz,zsrot,zauto,zenb,zmppun,zcas
      pointer nd,ne,trou,part,kgen,ndeg,iorb,ispin,itsym,
     $     zsaut,zoc,zocs,zfree,zheff
      dimension nd(:),ne(:),trou(:),part(:),kgen(:),ndeg(:),iorb(:),
     $     ispin(:),itsym(:),zsaut(:),zoc(:),zocs(:),zfree(:),zheff(:)
      COMMON /DET/ ND,NE,TROU,PART,KGEN,
     * NDEG,IORB,ISPIN,ITSYM,
     * ITS(20,20),ZSAUT,ZOC,ZOCS,ZFREE,ZHEFF
      INTEGER*2 NE,TROU,PART,NDEG,IORB,ISPIN
      pointer c
      dimension c(:)
      pointer fmpb
      dimension fmpb(:)
      COMMON /VEC/ C,E(ID14),wk(ID6),EEN(ID6),EMP(ID6),E2EN(ID6),
     *     E2MP(ID6),TESTE(ID6),TAUE(ID6),FMPB,brdmp(id6,id6),brden(id6
     $     ,id6),psimp(id6),psien(id6),psimpe(id6,14),psiene(id6,14)
     $     ,NTOT,NTOTEX(14)
      pointer f,aj,ak,num,indic,jndic,lndic
      dimension f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      COMMON /INT/ NI4,NI3,NI2,NR4,
     * F,AJ,AK,NUM,NDEB(500),NBO(99),
     * NAD(2000),KT(2000),INDIC,JNDIC,LNDIC
      INTEGER*4 INDIC,JNDIC,LNDIC
      INTEGER*2 NAD,KT
      idimorb2=idimorb*(idimorb+1)/2
      idimorb3=idimorb*2+1
      ndims=ncfd+1000
      ndimtp=ndim+20000
      ndimtp2=ndims+2
c      allocate(zgel(idimorb3))
c      allocate(fmpb(idimorb3))
c      allocate(c(id7))

      allocate(nd(ndimtp2))
      allocate(ne(ndimtp2))
      allocate(trou(ndimtp))
      allocate(part(ndimtp))
      allocate(kgen(ndims))
      allocate(ndeg(ndims))
      do i=1,ndims
         kgen(i)=0
         ndeg(i)=0
      enddo
      rewind 3
      call cerca(3,'ICINP',itrovo)
      if (itrovo.eq.0) then
      write (6,*) '   *** NAMELIST &ICINP NOT FOUND ***'
      stop 12
      endif
      call lecnam('&ICINP',namel,3)
      call nameli(namel,kgen,'KGEN',4)
c      allocate(iorb(idimorb3))
c      allocate(ispin(idimorb3))
c      allocate(itsym(idimorb3))
      allocate(zsaut(ndims))
c      allocate(zoc(idimorb3))
c      allocate(zocs(idimorb3))
c      allocate(zfree(idimorb))
c      allocate(zheff(id11))
c      allocate(f(idimorb2))
c      allocate(aj(idimorb2))
c      allocate(ak(idimorb2))
c      allocate(num(idimorb2))
c      allocate(indic(idimorb2))
c      allocate(jndic(idimorb2))
c      allocate(lndic(idimorb2))
      return
      end
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
      SUBROUTINE MEMCIP(ACTION,doit)
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
      COMMON /INT/ NI4,NI3,NI2,NR4
      pointer zgel
      dimension zgel(:)
      COMMON /CIP/ ESCF,TEST,TAU,TDEGEN,FACTOR,FINV,DELTA,CONVDB,STEPDB,
     * NORB,NORB2,NOCB,NOCA,nao,ngelo,ngelv,
     * INITB,INITA,NCF,NCFG,NSYM,ISYM,ISZ,NTRSY,METAT,NREF,NROT,
     * zgel,ZPRT,ZION,ZPERTU,ZBRD,ZPUN,ZDIAB,ZCOMB,ZRD,ZQD,
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
      norel=.false.
      deallocate (x)
      if(doit)then
      deallocate (ijkl2)
      deallocate (ijkl3)
      deallocate (ijkl4)
      deallocate (rijkl)
      endif
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
      allocate(x(iwant))
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
c     LOCX = loc(x)
C--------------------------------------
      loffs=1
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
      common /big/ px
      pointer x
      dimension x(:)
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
      allocate(x(iwant))
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
c      LOCX = loc(x)
C--------------------------------------
      loffs=1
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
      common /big/ x
      pointer x
      dimension x(:)
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
      allocate(x(iwant))
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
c     LOCX = loc(x)
C--------------------------------------
      loffs=1
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
*DECK LSAME
      LOGICAL FUNCTION LSAME (CA, CB)
C***BEGIN PROLOGUE  LSAME
C***SUBSIDIARY
C***PURPOSE  Test two characters to determine if they are the same
C            letter, except for case.
C***LIBRARY   SLATEC
C***CATEGORY  R, N3
C***TYPE      LOGICAL (LSAME-L)
C***KEYWORDS  CHARACTER COMPARISON, LEVEL 2 BLAS, LEVEL 3 BLAS
C***AUTHOR  Hanson, R., (SNLA)
C           Du Croz, J., (NAG)
C***DESCRIPTION
C
C  LSAME  tests if CA is the same letter as CB regardless of case.
C  CB is assumed to be an upper case letter. LSAME returns .TRUE. if
C  CA is either the same as CB or the equivalent lower case letter.
C
C  N.B. This version of the code is correct for both ASCII and EBCDIC
C       systems.  Installers must modify the routine for other
C       character-codes.
C
C       For CDC systems using 6-12 bit representations, the system-
C       specific code in comments must be activated.
C
C  Parameters
C  ==========
C
C  CA     - CHARACTER*1
C  CB     - CHARACTER*1
C           On entry, CA and CB specify characters to be compared.
C           Unchanged on exit.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   860720  DATE WRITTEN
C   910606  Modified to meet SLATEC prologue standards.  Only comment
C           lines were modified.  (BKS)
C   910607  Modified to handle ASCII and EBCDIC codes.  (WRB)
C   930201  Tests for equality and equivalence combined.  (RWC and WRB)
C***END PROLOGUE  LSAME
C     .. Scalar Arguments ..
      CHARACTER CA*1, CB*1
C     .. Local Scalars ..
      INTEGER IOFF
      LOGICAL FIRST
C     .. Intrinsic Functions ..
      INTRINSIC ICHAR
C     .. Save statement ..
      SAVE FIRST, IOFF
C     .. Data statements ..
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  LSAME
      IF (FIRST) IOFF = ICHAR('a') - ICHAR('A')
C
      FIRST = .FALSE.
C
C     Test if the characters are equal or equivalent.
C
      LSAME = (CA.EQ.CB) .OR. (ICHAR(CA)-IOFF.EQ.ICHAR(CB))
C
      RETURN
C
C  The following comments contain code for CDC systems using 6-12 bit
C  representations.
C
C     .. Parameters ..
C     INTEGER                ICIRFX
C     PARAMETER            ( ICIRFX=62 )
C     .. Scalar Arguments ..
C     CHARACTER*1            CB
C     .. Array Arguments ..
C     CHARACTER*1            CA(*)
C     .. Local Scalars ..
C     INTEGER                IVAL
C     .. Intrinsic Functions ..
C     INTRINSIC              ICHAR, CHAR
C     .. Executable Statements ..
C     INTRINSIC              ICHAR, CHAR
C
C     See if the first character in string CA equals string CB.
C
C     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
C
C     IF (LSAME) RETURN
C
C     The characters are not identical. Now check them for equivalence.
C     Look for the 'escape' character, circumflex, followed by the
C     letter.
C
C     IVAL = ICHAR(CA(2))
C     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
C        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
C     ENDIF
C
C     RETURN
C
C     End of LSAME.
C
      END
c-----------------------------------------------
      subroutine matvec(m,x,n,y,n2)
      implicit real*8(a-h,o-z)
      dimension x(n,m),y(n,m)
      pointer diag,ibuf,jbuf,a,olsen
      dimension diag(:),ibuf(:),jbuf(:),a(:),olsen(:)
      common/mtvc/diag,a,olsen,ibuf,jbuf
      REWIND 20
      READ (20) (DIAG(I),I=1,n)
      do l=1,m
         do i=1,n
            y(i,l)=diag(i)*x(i,l)
         enddo
      enddo
      MAX=100000
 35   read (20) (ibuf(i),i=1,size(ibuf)),
     $          (jbuf(i),i=1,size(jbuf)),
     $          (a(i),i=1,size(a))
      if(ibuf(100000).eq.-1)max=jbuf(100000)
      DO  IC=1,MAX
         ih=ibuf(ic)
         jh=jbuf(ic)
         DO  L=1,M
            y(ih,l)=y(ih,l)+a(ic)*x(jh,l)
            y(jh,l)=y(jh,l)+a(ic)*x(ih,l)
         enddo
      enddo
      if(ibuf(100000).ne.-1)goto 35
      return
      end
c---------------------------------------------
      subroutine correc( m, ritzv, x,ldx, r,ldr )
      implicit real*8(a-h,o-z)
*
      integer m, ldx, ldr
      double precision x(ldx,m), r(ldr,m), ritzv(m)
*
*---  Diagonal correction ...
*
*---  Matrix data and corrector ...
      pointer diag,ibuf,jbuf,a,olsen
      dimension diag(:),ibuf(:),jbuf(:),a(:),olsen(:)
      common/mtvc/diag,a,olsen,ibuf,jbuf
*
      integer i, j
      n=ldx
      do j = 1,m
         do i = 1,n
            olsen(i)=0.d0
            tmp = diag(i) - ritzv(j)
            if ( tmp.ne.0.0d0 )then
               r(i,j) = r(i,j)/tmp
               olsen(i)=x(i,j)/tmp
            endif
         enddo
         cappa=ddot(n,x(1,j),1,r(1,j),1)/ddot(n,x(1,j),1,olsen,1)
         call daxpy(n,-cappa,olsen,1,r(1,j),1)
      enddo
      end
      SUBROUTINE DGEMVI (  M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
      integer*1 A(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      parameter (TRANS='N')
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
c-----------------------------------------------
      subroutine av(n,x,y)
      implicit real*8(a-h,o-z)
      dimension x(*),y(*)
      pointer diag,ibuf,jbuf,a,olsen
      dimension diag(:),ibuf(:),jbuf(:),a(:),olsen(:)
      common/mtvc/diag,a,olsen,ibuf,jbuf
      REWIND 20
      READ (20) (DIAG(I),I=1,n)
      do i=1,n
         y(i)=diag(i)*x(i)
      enddo
      MAX=100000
 35   read (20) (ibuf(i),i=1,size(ibuf)),
     $          (jbuf(i),i=1,size(jbuf)),
     $          (a(i),i=1,size(a))
      if(ibuf(100000).eq.-1)max=jbuf(100000)
      DO  IC=1,MAX
         ih=ibuf(ic)
         jh=jbuf(ic)
         y(ih)=y(ih)+a(ic)*x(jh)
         y(jh)=y(jh)+a(ic)*x(ih)
      enddo
      if(ibuf(100000).ne.-1)goto 35
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
c---------------------------------------------------------------------
*----------------------------------------------------------------------|
* FILE
*     davpack.f - Variable-Block Davidson algorithm with deflation
*
* CONTENTS
*     subroutine davson(...) - Variable-Block Davidson with deflation
*     subroutine mgs(...) - Mofified Gram-Schmidt orthogonalization
*     subroutine ortho(...) - Selective orthonormal expansion
*
* DESCRIPTION
*
* #####################################################################|
* subroutine davson( n, nbx, nb, nev, itmax, iter, nmult, tol,
*    .               anorm, res, eig, x,ldx, wrk,lwrk, iwrk,liwrk,
*    .               matvec, correc, ilevel, ifile, ifail )
* Purpose
*     This is the main computational routine. It computes a few of the 
*     extreme (ie, rightmost or leftmost) eigenpairs of a large sparse
*     symmetric matrix via the variable-block Davidson method with 
*     deflation.  
* Called Routines
*     internal: mgs(...) - Mofified Gram-Schmidt orthogonalization
*     internal: ortho(...) - Selective orthonormal expansion
*     external: matvec(...) - block matrix-vector multiplication routine
*     external: correc(...) - corrector routine
*     blas: dscal(...) - multiplication of a vector by a scalar
*     blas: daxpy(...) - constant times a vector plus a vector
*     blas: ddot(...) - scalar product of two vectors
*     blas: dswap(...) - swaps the contents of two vectors
*     blas: dnrm2(...) - euclidian norm of a vector
*     blas: dcopy(...) - copies a vector into another
*     blas: idamax(...) - returns index of the maximum entry of a vector
*     blas: dgemm(...) - general dense matrix-matrix multiplication
*     lapack: dsyevx(...) - eigenpairs of a dense symmetric matrix
*
* #####################################################################|
* subroutine mgs( n,m,k, v,ldv, x,ldx,mx )
* Purpose
*     Modified Gram-Schmidt acting on v(1:n,k:m). 
*     v(:,k:m) is orthogonalized against x(:,1:mx) and v(:,1:k-1).
* Called Routines
*     blas: dscal(...) - multiplication of a vector by a scalar
*     blas: daxpy(...) - constant times a vector plus a vector
*     blas: ddot(...) - scalar product of two vectors 
*     blas: dnrm2(...) - euclidian norm of a vector
* Calling Routines
*     internal: davson(...) - Variable-Block Davidson with deflation
*
* #####################################################################|
* subroutine ortho( n,m,nb, t,ldt, v,ldv, x,ldx,mx, kb,kbmax, drop)
* Purpose
*     Selective orthonormal expansion of v(1:n,1:m) with vectors from t.
*     Columns of t are orthogonalized and selectively re-orthogonalized
*     against x and the current v. The resulting vector is added in v
*     if its norm is greater than the drop-tolerance drop. 
* Called Routines
*     blas: dscal(...) - multiplication of a vector by a scalar
*     blas: daxpy(...) - constant times a vector plus a vector
*     blas: ddot(...) - scalar product of two vectors
* Calling Routines
*     internal: davson(...) - Variable-Block Davidson with deflation
*
* END DESCRIPTION
*----------------------------------------------------------------------|
*
* AUTHORS
*
* M. Sadkane - sadkane@univ-brest.fr   R. B. Sidje - rbs@maths.uq.edu.au
* Departement de Mathematiques         Department of Mathematics
* Universite de Bretagne Occidentale   University of Queensland
* 6, Avenue Le Gorgeu, B.P. 809        Brisbane QLD 4072
* 29285 Brest Cedex. France            Australia
*
* REVISION DATE: 11/JAN/1999
*----------------------------------------------------------------------|

      subroutine davson( n, nbx, nb, nev, itmax, iter, nmult, tol,
     .                   anorm, res, eig, x,ldx, wrk,lwrk, iwrk,liwrk,
     .                   matvec, correc, ilevel, ifile, ifail )
*
*---  variable-Block Davidson method with deflation ---
*
      implicit none
      integer n, nbx, nb, nev, itmax, ldx, lwrk, liwrk, nmult
      double precision wrk(lwrk), res(*), eig(*), x(ldx,*), tol, anorm
      integer iwrk(liwrk), ilevel, ifile, ifail
      external matvec, correc
*
*-----Purpose----------------------------------------------------------|
*     this subroutine computes a few extreme (ie, rightmost or leftmost) 
*     eigenvalues and the corresponding eigenvectors of a symmetric 
*     matrix by the variable-block Davidson method with deflation.
*     the matrix is referenced implicitly through an external routine 
*     performing the sparse matrix-vector product. Upon completion,
*     the computed eigenpairs are in no particular order. 
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the symmetric matrix.
*
*     nbx    : (input) maximum size allowable for the basis. 
*              nbx .ge. mb, where mb = MAX( nb, |nev| ).
*
*     nb     : (input) initial block-size.
*
*     nev    : (input/output) number of wanted eigenpairs.
*              A positive value will search for the rightmost
*              eigenpairs whereas a negative value will search for the 
*              leftmost eigenpairs. If all of the wanted eigenpairs are
*              found then, on exit nev is set to zero; otherwise |nev|
*              yields the number of eigenpairs that have not converged.
*
*     itmax  : (input) maximum allowable number of iterations.
*
*     iter   : (output) number of iterations used up to convergence.
*
*     nmult  : (output) number of matrix-vector multiplications
*              used by this routine.
*
*     tol    : (input) required tolerance on the residuals. 
*              If tol.le.0.0d0 the square root of eps is used instead.
*              `eps' is the machine unit roundoff (computed internally).
*
*     anorm  : (input) an estimate of some norm of A. this parameter
*              provides a means to select a particular convergence test:
*              - If anorm>0.0d0, an eigenpair (x,lambda) with ||x|| = 1
*                is accepted if the relative residual 
*                ||A*x - lambda*x||/anorm <= tol.
*                If anorm=1.0d0, the test is therefore based on the
*                absolute residual ||A*x - lambda*x|| <= tol.
*              - If anorm=0.0d0, the code uses the relative residual 
*                ||A*x - lambda*x||/MAX(eps^{2/3},ABS(lambda)) <= tol.
*
*     res(mb): (output) the first |nev| entries have the relative
*              residuals with respect to the setting of anorm above, ie.
*              - If anorm > 0.0d0, then
*                   res = ||A*x - lambda*x||/anorm
*                   Thus, if anorm = 1.0d0 then res = ||A*x - lambda*x||
*              - If anorm = 0.0d0, then
*                   res = ||A*x-lambda*x||/MAX(eps^{2/3},ABS(lambda))
*
*     eig(mb): (output) the first |nev| entries are the eigenvalues
*              that have converged.
*
*   x(ldx,mb): (input/output) on entry, x should contain mb starting
*              guesses (i.e., mb starting eigenvector estimations). 
*              On exit, the first |nev| columns contain the final
*              computed eigenvector approximations.
*
*   wrk(lwrk): (workspace) 
*              lwrk .ge. 2*n*nbx+nbx*nbx+(nbx+n)*mb+8*nbx
*              Easier:  lwrk .ge. (n+nbx)*(mb+nbx) + (n+8)*nbx
*              Its utilization is depicted as follows:
*               2*n*nbx + nbx^2 + nbx*mb+n*mb  +    8*nbx   
*              +--------+-------+--------------+--------------+
*               V & W=AV    H    eigv(H) & Residuals. wsp for LAPACK
*
* iwrk(liwrk): (workspace) liwrk .ge. 6*nbx.
*
*     matvec : external subroutine for block matrix-vector product.
*              synopsis: matvec( m, x,ldx, y,ldy )
*                        double precision x(ldx,m), y(ldy,m)
*              computes: y(1:n,1:m) <- A*x(1:n,1:m)
*                        where A is the principal matrix.
*
*     correc : external subroutine for correction (preconditioning).
*              synopsis: correc( m, x,ldx, r,ldr, ritzval ) 
*                        double precision x(ldx,m), r(ldr,m), ritzval(m)
*              computes: for j = 1:m
*                          r(1:n,j) <- M_j( x(1:n,j), r(1:n,j) )
*                        where M_j is the preconditioner operator.
*
*     ilevel : (input, info-level) specifies the amount of runtime 
*              information to be reported:
*              0 -Silent running        3,4 -Convergence info
*              1 -Iteration info          5 -Higly verbose (debug mode).
*              2 -Sub-iteration info
*
*     ifile  : (input, info-file) the logical unit of a file into which
*              information will be reported:
*              0, 6  -Screen
*              k     -Output to file opened with unit k. If such a file
*                     was not opened, information goes into fort.k
*
*     ifail  : (output) exit completion code:
*              0     -the execution was OK
*              not 0 -A problem was encountered:
*
*              *A negative value means bad input arguments:
*               ifail = -1 if mb.gt.n  where mb = MAX( nb, |nev| )
*               ifail = -2 if nbx.gt.n .or. nbx.lt.mb 
*               ifail = -3 if ldx.lt.n 
*               ifail = -4 if liwrk.lt.6*nbx 
*               ifail = -5 if lwrk.lt.n*|nev|+(n+nbx)*(mb+nbx)+(n+8)*nbx
*
*               Note that if inputs are not initialized properly (e.g.
*               if n = 0), one or another error may be reported.
*
*              *A positive value means a runtime failure occurred.
*               Current exit codes are:
*               ifail = 8 if a failure occurred in LAPACK when solving
*                         the Ritz problem.
*               ifail = 9 if the maximum allowable number of iterations
*                         was reached without convergence of all of the
*                         eigenpairs.
*
*     ADDITIONAL NOTES:
*     -----------------
*     When the information level, ilevel > 0, the routine outputs:
*         iter : current number of iterations.
*        basis : current size of the basis.
*        nmult : current number of matrix-vector multiplications used.
*        block : current block size.
*         nevf : number of eigenpairs found so far.
*     residual : maximum `res' (see above).
*
*----------------------------------------------------------------------|
*
      double precision zero, one
      parameter( zero=0.0d0, one=1.0d0 )
      integer iv,iw,i1,i2,i3,i4, i,j,ii,ij, iter, kb,kbnew, nevf, mb, m
      integer iold,ilst, ilower,iupper, ncomp,nconv, isgn, io
      double precision a1, a2, a3, etol, eps, seps, eps23, residual
      logical  endsub
      character fm*255, fm2*255, separator*80
      integer IDAMAX
      double precision DNRM2
      intrinsic ABS, SQRT, MAX, MIN

*
*---  Get working precision...
*
      a1 = 4.0d0/3.0d0
 1    a2 = a1 - 1.0d0
      a3 = a2 + a2 + a2
      eps = ABS( a3-1.0d0 )
      if ( eps.eq.0.0d0 ) goto 1
      seps = SQRT(eps)
      eps23 = eps**(2.0d0/3.0d0)

*
*---  Preliminary checks and initializations...
*
      ifail = 0
      nmult = 0
      m = ABS( nev )
      mb = MAX( nb,m )
      io = ifile
      if ( io.le.0 ) io = 6

      if ( mb.gt.n )                            ifail = -1
      if ( nbx.gt.n .or. nbx.lt.mb )            ifail = -2
      if ( ldx.lt.n )                           ifail = -3
      if ( liwrk.lt.6*nbx )                     ifail = -4
      if ( lwrk.lt.(n+nbx)*(mb+nbx)+(n+8)*nbx ) ifail = -5
 
      if ( ifail.ne.0 ) then
         write(UNIT=io,FMT=*) 'Bad arguments in davson '//
     .       '(see documentation), ifail:',ifail
         return
      endif

      iv = 1
      iw = iv + n*nbx
      i1 = iw + n*nbx 
      i2 = i1 + nbx*nbx
      i3 = i2 + nbx*mb
      i4 = i3 + n*mb

      if ( tol.le.0.0d0 ) then
         etol = seps
      else
         etol = tol 
      endif
      if ( ilevel.ge.1 ) then      
         write(UNIT=io,FMT=*) 'machine epsilon =',eps
         write(UNIT=io,FMT=*) 'tolerance used  =',etol
      endif
      separator = "--------------------------------------------"
      if ( ilevel.gt.0 .and. ilevel.le.2 ) then
         write(UNIT=io,FMT='(/,A)')
     .         'iter   basis  nmult  block  nevf   residual'
         fm  = '(I4,   3x,I4, 3x,I4, 4x,I2, 5x,I2, 2x,1P,E11.3 )'
         fm2 = '(      7x,I4, 3x,I4, 4x,I2, 5x,I2, 2x,1P,E11.3 )'
      else
         fm = '("iter =",I4, ",  basis =",I3, ",  nmult =",I4' // 
     .        ',",  block =",I2, ",  nevf =",I2' //
     .        ',",  residual =",1P,E11.3 )'
      endif

      iter = 0
      nevf = 0
      kb = mb
      call mgs( n,kb,1, x,ldx, x,ldx,0 )
      do j = 1,kb
         call DCOPY( n, x(1,j),1, wrk(iv+(j-1)*n),1 )
      enddo
*----------------------------------------------------------------------|
*---  B E G I N N I N G    O F    T H E    I T E R A T I O N S      ---|

      if ( nev.gt.0 ) then
	 isgn =  1
      else
         isgn = -1
      endif 
 
 10   continue

      ilst = 0
      iter = iter + 1
      if ( iter.gt.itmax ) then
         ifail = 8   
         write(UNIT=io,FMT='("Wanted=",I2," Converged=",I2)') m,nevf 
         write(UNIT=io,FMT=*) 'Number of iterations is reached '//
     .                        'without an overall convergence.'
         return
      endif

*----------------------------------------------------------------------|
*---  B E G I N N I N G    O F   T H E    S U B I T E R A T I O N S ---|

 100  continue

      iold = ilst
      ilst = ilst + kb
      endsub = ilst.eq.nbx
      if ( ilevel.ge.3 ) write(UNIT=io,FMT=*) 'Basis size =',ilst

*     -----------------------------------------------------------------|
*---  Interaction matrix...
*     -----------------------------------------------------------------|
      nmult = nmult + kb
      call matvec( kb, wrk(iv+iold*n),n, wrk(iw+iold*n),n )      
      call DGEMM('t','n',ilst,kb,n,one,wrk(iv),n,
     .            wrk(iw+iold*n),n,zero,wrk(i1+iold*nbx),nbx)
      if ( iold.gt.0 ) then
*---     Restoration: replicate strict upper-part into lower-part...
         ij = i1
         ii = i1
         do j = 2,ilst
            ij = ij + nbx
            ii = ii + 1
            call DCOPY( j-1, wrk(ij),1, wrk(ii),nbx )
         enddo
      endif
      if ( ilevel.ge.5 ) then
         write(UNIT=io,FMT=*) 'Interaction matrix ='
         do i = 1,ilst
            write(UNIT=io,FMT=*) (wrk(i1+(j-1)*nbx+i-1),j=1,ilst)
         enddo
      endif

*     -----------------------------------------------------------------|
*---  Ritzvalues and Ritzvectors...
*     -----------------------------------------------------------------|
*---  Select the range of the spectrum of interest...
      if ( nev.gt.0 ) then
         ilower = MAX( ilst-(mb-nevf)+1,1 )
         iupper = ilst
      else
         ilower = 1
         iupper = MIN( mb-nevf,ilst )
      endif
*---  Save the diagonal of the interaction matrix...
      call DCOPY( ilst, wrk(i1),nbx+1, wrk(i3+nevf*n),1 )
*---  Solve the reduce-sized Ritz problem -- the lower-part of the
*     interaction matrix is destroyed...
      call DSYEVX( 'v','i','l', ilst, wrk(i1),nbx,
     .              a1,a2, ilower,iupper, -etol, ncomp, 
     .              eig(nevf+1), wrk(i2),nbx,
     .              wrk(i4),8*nbx, iwrk(1), iwrk(5*nbx+1), ifail )
*---  Restore diagonal (strict lower-part will be restored at next iter)
      call DCOPY( ilst, wrk(i3+nevf*n),1, wrk(i1),nbx+1 )
      if ( ifail.ne.0 .or. ncomp.ne.iupper-ilower+1 ) then
         ifail = 9
         write(UNIT=io,FMT=*) 'Failure in LAPACK when solving ' //
     .                        'the Ritz problem.'
         return
      endif
      if ( ilevel.ge.3 ) then
         write(UNIT=io,FMT=*) 'Ritzvalues ='
         write(UNIT=io,FMT=*) (eig(nevf+j),j=1,ncomp)
      endif
      if ( ilevel.ge.5 ) then
         write(UNIT=io,FMT=*) 'Ritzvectors ='
         do i = 1,ilst
            write(UNIT=io,FMT=*) (wrk(i2+(j-1)*nbx+i-1),j=1,ncomp)
         enddo
      endif

*     -----------------------------------------------------------------|
*---  Eigenvectors and Residual vectors...
*     -----------------------------------------------------------------|
      call DGEMM('n','n',n,ncomp,ilst,one,wrk(iv),n,
     .            wrk(i2),nbx,zero,x(1,nevf+1),ldx )
      call DGEMM('n','n',n,ncomp,ilst,one,wrk(iw),n,
     .            wrk(i2),nbx,zero,wrk(i3+nevf*n),n )
      do j = 1,ncomp
         i = i3 + (nevf+j-1)*n
         call DAXPY( n, -eig(nevf+j), x(1,nevf+j),1, wrk(i),1 )
         res(nevf+j) = DNRM2( n,wrk(i),1 )
         if ( anorm.gt.0.0d0 ) then         
            res(nevf+j) = res(nevf+j) / anorm
         else
            res(nevf+j) = res(nevf+j) / MAX( eps23,ABS(eig(nevf+j)) )
         endif
      enddo
      residual = res( IDAMAX(mb,res,1) )
      if ( ilevel.ge.3 ) then
         write(UNIT=io,FMT=*) 'Residuals ='
         write(UNIT=io,FMT=*) (res(nevf+j),j=1,ncomp)
      endif

*     -----------------------------------------------------------------|
*---  Deflation - lock the convergent eigenpairs on the left...
*     -----------------------------------------------------------------|
      nconv = 0
      do j = 1,ncomp
         if ( res(nevf+j).le.etol ) then
            nconv = nconv + 1
            if ( nconv.ne.j ) then
               call DSWAP( n, x(1,nevf+nconv),1, x(1,nevf+j),1 )
               call DSWAP( 1, res(nevf+nconv),1, res(nevf+j),1 )
               call DSWAP( 1, eig(nevf+nconv),1, eig(nevf+j),1 )
               call DSWAP( n, wrk(i3+(nevf+nconv-1)*n),1, 
     .                        wrk(i3+(nevf+j-1)*n),1 )
            endif
         endif
      enddo
      nevf = nevf + nconv
      if ( nev.gt.0 ) then
         nev = nev - nconv
      else
         nev = nev + nconv
      endif

*     -----------------------------------------------------------------| 
*---  Check for completion...
*     -----------------------------------------------------------------|
      if ( nevf.ge.m ) then
         if ( ilevel.gt.0 ) then
            residual = res( IDAMAX(m,res,1) )
            write(UNIT=io,FMT=fm) iter, ilst, nmult, kb, nevf, residual
            if ( ilevel.gt.1) write(UNIT=io,FMT='(A)') separator
            write(UNIT=io,FMT=*) ' '
         endif
         nev = 0
         return
      endif

*     -----------------------------------------------------------------|
*---  Update the starting guesses that have not yet converged...
*     -----------------------------------------------------------------|
      iold = ilst
      if ( endsub .or. nconv.ne.0 ) then
         ilst = mb-nevf
         call mgs( n,mb,nevf-nconv+1, x,ldx, x,ldx,0 )
         do j = 1,ilst
            call DCOPY( n, x(1,nevf+j),1, wrk(iv+(j-1)*n),1 )
         enddo
      endif

*     -----------------------------------------------------------------|
*---  Correction of the residuals... 
*     -----------------------------------------------------------------|
      call correc( ncomp-nconv, eig(nevf+1),
     .             x(1,nevf+1),ldx, wrk(i3+nevf*n),n )

*     -----------------------------------------------------------------|
*     Selective incrementation of the basis with the corrections...
*     The maximum allowable number of selections is MIN( nb, nbx-ilst )
*     (nbx-ilst is the remaining room in the basis). The number of
*     corrected residuals actually selected is returned in kbnew.
*     -----------------------------------------------------------------|
      call ortho( n,ilst, isgn*(ncomp-nconv), wrk(i3+nevf*n),n,
     .            wrk(iv),n, x,ldx,nevf, kbnew,MIN(nb,nbx-ilst), 
     .            eps )

*     -----------------------------------------------------------------|
*---  Restart if basis size is reached, or if an eigenpair converged, 
*     or if none of the corrected residuals has been selected.
*     the block-size kb becomes the number of current guesses augmented
*     with the number of selected residuals...
*     -----------------------------------------------------------------|
      if ( endsub .or. nconv.ne.0 .or. kbnew.eq.0 ) then
         if ( ilevel.gt.0 ) then
            write(UNIT=io,FMT=fm) iter, iold, nmult, kb, nevf, residual
            if ( ilevel.gt.1 ) write(UNIT=io,FMT='(A)') separator
         endif
         kb = mb-nevf + kbnew
         call mgs( n,kb,1, wrk(iv),n, x,ldx,nevf )
         goto 10
      endif

*     -----------------------------------------------------------------| 
*---  Otherwise continue with the incremented search subspace...
*     -----------------------------------------------------------------|
      if(ilevel.eq.2) write(UNIT=io,FMT=fm2)iold,nmult,kb,nevf,residual
      kb = kbnew
      goto 100
      end    
*----------------------------------------------------------------------| 
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::| 
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::|
*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::|
*----------------------------------------------------------------------|
      subroutine mgs( n,m,k, v,ldv, x,ldx,mx )
*
*---  Onward Modified Gram-Schmidt acting on v(:,k:m). 
*     v(:,k:m) is orthogonalized against x(:,1:mx) and v(:,1:k-1).
*     It is assumed that x(:,1:mx) and v(:,1:k-1) are already orthogonal
*     If k = 1 and mx = 0, the usual MGS is recovered.
*
      implicit none
      integer n, m, mx, k, ldv, ldx, i, j
      double precision v(ldv,m), x(ldx,mx), s, DDOT, DNRM2

      do i = k,m
*---     orthogonalize against x(:,1:mx)
         do j = 1,mx
            s = -DDOT( n, x(1,j),1, v(1,i),1 )
            call DAXPY( n, s, x(1,j),1, v(1,i),1 )
         enddo
*---     orthogonalize against v(:,1:k-1)
         do j = 1,i-1
            s = -DDOT( n, v(1,j),1, v(1,i),1 )
            call DAXPY( n, s, v(1,j),1, v(1,i),1 )
         enddo
         s = 1.0d0 / DNRM2( n, v(1,i),1 )
         call DSCAL( n, s, v(1,i),1 )
      enddo
      end
*----------------------------------------------------------------------|
      subroutine ortho( n,m,nb, t,ldt, v,ldv, x,ldx,mx, kb,kbmax, drop)
      implicit none
      integer n, m, nb, mx, ldt, ldv, ldx, kb, kbmax
      double precision t(ldt,*), v(ldv,m+kbmax), x(ldx,mx), drop
*
*---  Selective orthonormal expansion of v(1:n,1:m) with vectors from
*     t(1:n,1:|nb|). 
*
*---  The selection process proceeds as follows:
*     Each column of t is a candidate vector which is orthogonalized
*     and selectively re-orthogonalized against x and the current v.
*     The candidate is retained if its norm is greater than drop, the
*     drop-tolerance. A retained candidate is normalized and included
*     in v, thus it subsequently participates in the orthogonalization 
*     and re-orthogonalization processes against the next candidates.
*
*     The sign of nb determines whether the candidates should be taken 
*     from t from the first to the last column, or from the last column
*     to the first. this way of doing so contributes in ensuring that
*     the most meaningful candidates are retained first.
*
*---
*     nb    number of columns in t, nb can be positive or negative:
*           positive value: the interest is in the rightmost eigenpairs, 
*                           scan t from last column to first (backward).
*           negative value: the interest is in the leftmost eigenpairs,
*                           scan t from first column to last (forward).
*
*     kbmax specifies the maximum allowable number of vectors that can
*           be selected.
*
*     kb    returns the actual number of vectors that have been
*           selected.
*---
*
      integer j, next, istep, iorth
      double precision drop2, eta, s, ss
      parameter( eta=0.1d0 )
      double precision DDOT
      intrinsic SQRT

*---  Specify whether to scan the columns of t backward or forward...
      if ( nb.gt.0 ) then
         istep = -1
         next = nb+1
      else
         istep = 1
         next = 0
         nb = -nb
      endif

*---  Orthogonalization + selective-reorthogonalisation + selection ...
      kb = 0
      drop2 = drop*drop
 10   next = next + istep
      if ( next.le.0 .or. next.gt.nb. or. kb.ge.kbmax ) return
      iorth = 2
 11   iorth = iorth - 1
      ss = 0.0d0
*---  Orthogonalization against x...
      do j = 1,mx
         s = -DDOT( n, x(1,j),1, t(1,next),1 )
         call DAXPY( n, s, x(1,j),1, t(1,next),1 )
         ss = ss + s*s
      enddo
*---  Orthogonalization against v...
      do j = 1,m+kb
         s = -DDOT( n, v(1,j),1, t(1,next),1 )
         call DAXPY( n, s, v(1,j),1, t(1,next),1 )
         ss = ss + s*s
      enddo
*---  Selective-reorthogonalization...
      s = DDOT( n, t(1,next),1, t(1,next),1 )
      if ( iorth.gt.0 .and. s.ge.drop2 .and. s.le.eta*ss ) goto 11
*---  Selection...
      if ( s.ge.drop2 ) then
         kb = kb + 1
         s = SQRT( s ) 
         call DCOPY( n, t(1,next),1, v(1,m+kb),1 )
         call DSCAL( n, 1.0d0/s, v(1,m+kb),1 )
      endif
      goto 10
      end
*----------------------------------------------------------------------|
c----------------------------------------------------------------------
*----------------------------------------------------------------------|
* FILE
*     eigsrt.f  - sorts the eigenvalues or residuals and moves the
*                 eigenvectors accordingly
* CONTENTS
*     subroutine eigsrt( n, nev, res, eig, x,ldx )
*
* DESCRIPTION
*
* #####################################################################|
* subroutine eigsrt( n, nev, res, eig, x,ldx )
* Purpose
*     Given a positive integer nev and the arrays eig(1:nev), res(1:nev)
*     x(1:n,1:nev), such that x(:,j) corresponds to eig(j) and res(j), 
*     this routine sorts eig(1:nev) in increasing order and rearranges
*     res(1:nev) and x(1:n,1:nev) so that the initial correspondence is
*     maintained.
*     If nev is negative, the sort of eig(1:|nev|) is done in decreasing 
*     order and the initial correspondence is also maintained.
*     An interested user can sort the eigenpairs w.r.t residuals by
*     simply inverting the role of eig and res in the call, i.e., with:
*     call eigsrt( n, nev, eig, res, x,ldx )
* Called Routines
*     blas: dswap(...) - swaps the contents of two vectors
* Calling Routines
*     runme.f: -main-
*     alone.f: -main-
*
* END DESCRIPTION
*----------------------------------------------------------------------|
*
* AUTHORS
*
* M. Sadkane - sadkane@univ-brest.fr   R. B. Sidje - rbs@maths.uq.edu.au
* Departement de Mathematiques         Department of Mathematics
* Universite de Bretagne Occidentale   University of Queensland
* 6, Avenue Le Gorgeu, B.P. 809        Brisbane QLD 4072
* 29285 Brest Cedex. France            Australia
*
* REVISION DATE: 11/JAN/1999
*----------------------------------------------------------------------|

      subroutine eigsrt( n, nev, res, eig, x,ldx )
      implicit none
      integer          n, nev, ldx
      double precision res(*), eig(*), x(ldx,*)
*
*-----Purpose----------------------------------------------------------|
*     This subroutine sorts the eigenpairs w.r.t. eigenvalues.
*
*-----Arguments--------------------------------------------------------|
*
*     n        : (input) order of x.
*
*     nev      : (input) number of eigenpairs. 
*                A positive value will sort in increasing algebraic
*                order whereas a negative value will sort in decreasing 
*                algebraic order.
*
*     res(*)   : (input/output)
*                On input res(1:|nev|) are the residuals associated
*                to the eigenpairs. On output these are sorted 
*                appropriately.
*
*     eig(*)   : (input/output)
*                On input eig(1:|nev|) are the eigenvalues.
*                On output these are sorted.
*
*     x(ldx,*) : (input/output)
*                On input, x(1:n,1:|nev|) are the eigenvectors.
*                On output these are sorted w.r.t. the eigenvalues.
*
*---  NOTE: 
*     An interested user can sort the eigenpairs w.r.t residuals with:
*     call eigsrt( n, nev, eig, res, x,ldx )
*
*----------------------------------------------------------------------|
*
      integer i, j, m
      logical swap
      intrinsic ABS

      m = ABS( nev )
      do i = 1,m-1
         do j = i+1,m
            if ( nev.gt.0 ) then
               swap = eig(i).gt.eig(j)
            else
               swap = eig(i).lt.eig(j)
            endif
            if ( swap ) then
               call dswap( 1, res(i),1, res(j),1 )
               call dswap( 1, eig(i),1, eig(j),1 )
               call dswap( n, x(1,i),1, x(1,j),1 )
            endif
         enddo
      enddo
      END
*----------------------------------------------------------------------|
      subroutine printtp(ne,nd,trou,part,ncf,norb)
      integer*2 ne(*),trou(*),part(*)
      integer nd(*)
      do i=1,ncf
         ndi=nd(i)
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
      subroutine crea25
      implicit real*8 (a-h,o-z)
      common /fil/ file12,file13,file26,zwconf
      character(len=64)::file12,file13,file26
      logical*1 zwconf
      integer frozen(8),delete(8)
      dimension its(8,8)
      allocatable itsym(:),hmono(:)
      character(len=10000) namel
      character(len=80)aa
      frozen=0
      delete=0
      rewind 14
      call cerca(14,'cost_AO',itrovo)
      if(itrovo.eq.0)then
         print*,'COST_AO not found'
         call flush(6)
      endif
      call lecnam('&COST_AO',namel,14)
      call nameli(namel,ITS,'its',4)
      call nameli(namel,nao,'norb',4)
      rewind 14
      call cerca(14,'cost_MO',itrovo)
      if(itrovo.eq.0)then
         print*,'COST_MO not found'
         call flush(6)
      endif
      write (6,*) 'COST_MO',cost_mo
      call lecnam('&COST_MO',namel,14)
      call nameli(namel,NSYM,'nsym',4)
      call nameli(namel,NORB,'norb',4)
      call nameli(namel,NOC,'noc',4)
      call nameli(namel,FROZEN,'frozen',4)
      call nameld(namel,enuc,'enuc',8)
      delete=nao-norb-frozen
      allocate(itsym(norb))
      call nameli(namel,ITSYM,'itsym',4)
      write(25) NSYM,NORB,NOC,0,(ITSYM(I),I=1,NORB)
     *,nao,sum(frozen),sum(delete)
      write(25) ((ITS(I,J),I=1,NSYM),J=1,NSYM)
      write(25)enuc,enuc
      rewind 25
      nter=NORB*(NORB+1)/2
      allocate(hmono(nter))
      read(28) aa
 18   continue
      read(28) aa
      read(28)
      read(28)
c        essai de lecture en triangle
      if(index(aa,'ONEINT MATRIX (TOTAL)').eq.0) goto 18
      backspace (28)
      backspace (28)
      read(28) 
      read(28) hmono
      if(file26.eq.' ')stop 'manca il file 26!'
      open(26,file=file26,form='UNFORMATTED',status='UNKNOWN')
      rewind 26
      write(26) hmono
      rewind 26
      close(unit=28)
      close(unit=26)
      close(unit=14)
      print*,'End of reading from Cost format'
      return
      end
c----------------------------------------------------------
      recursive subroutine F(n,s2,ivec)
      implicit none
      integer n,s2,nn,n2,i,dime,g,g1,jdummy,n1,wigner,j
      real*8 ivec(0:*)
      real*8,dimension(0:1):: al=(/1.0d0,0.0d0/),be=(/0.0d0,1.0d0/)
      integer*1,allocatable:: SM(:)
      real*8,allocatable::jvec(:),jaux(:),jaux2(:)
      real*8 factor
      nn=dime(n)
      n2=dime(n-1)
      allocate(SM(0:n2**2-1))
      call getSM(n-1,SM)
      if(n.eq.1)then
         ivec(0)=al(0)
         ivec(1)=al(1)
      else
         n1=dime(n-2)
         g=0
         if(s2.gt.0)then
            g=wigner(n-1,s2-1)
            if(g.eq.0)print*,'g is 0 !'
            allocate(jvec(0:n2*g-1))
            call F(n-1,s2-1,jvec)
            do i=0,g-1
               call dirprod(n2,jvec(n2*i),2,al,ivec(nn*i))
            enddo
            deallocate(jvec)
         endif
         if((s2+1).le.(n-1))then
            g1=wigner(n-1,s2+1)
            factor=sqrt(dble((s2+1)*(s2+2)))
            allocate(jvec(0:n2*g1-1))
            allocate(jaux(0:n2-1))
            allocate(jaux2(0:nn-1))
            call F(n-1,s2+1,jvec)
            do i=0,g1-1
               call matveci(n2,SM,jvec(n2*i),jaux)
               call dirprod(n2,jaux,2,al,ivec(nn*(i+g)))
               call dirprod(n2,jvec(n2*i),2,be,jaux2)
               do j=0,nn-1
                  jdummy=nn*(i+g)+j
                  ivec(jdummy)=(ivec(jdummy)+jaux2(j)*dble(s2+1))/factor
               enddo
            enddo
            deallocate(jvec)
         endif
      endif
      deallocate(SM)
      return
      end
c----------------------------------------------------
      recursive integer function dime(n) result(res)
      implicit none
      integer n
      if(n.eq.0)then
         res=1
         return
      endif
      res=2*dime(n-1)
      return
      end
c-----------------------------------------------------------------
      recursive subroutine getSM(n,SM)
      implicit none
      integer n,m,k,i,dime
      integer*1 SM(0:*)
      integer*1, dimension(0:3):: sms=(/0,1,0,0/),one=(/1,0,0,1/)
      integer*1,allocatable,dimension(:)::SM1
      m=dime(n)
      if(n.eq.0)then
         SM(0)=1
      elseif(n.eq.1)then
         call copymat(m,SM,sms)
      else
         k=dime(n-1)
         allocate(SM1(0:k*k-1))
         call getSM(n-1,SM1)
         call matprod(k,SM1,2,one,SM)
         call matprodacc(k,2,sms,SM)
         deallocate(SM1)
      endif
      end
c------------------------------------------------------------
      subroutine matprod(n,a,m,b,c)
      implicit none
      integer n,m,i,j,ij,k,l,kl,nm,nk
      integer*1 a(0:*),b(0:*),c(0:*)
      nm=n*m
      ij=-1
      do i=0,n-1
         do j=0,m-1
            ij=ij+1
            kl=-1
            do k=0,n-1
               nk=n*k
               do l=0,m-1
                  kl=kl+1
                  c(nm*kl+ij)=a(nk+i)*b(m*l+j)
               enddo
            enddo
         enddo
      enddo
      end
c--------------------------------------------------------------
      subroutine matprodacc(n,m,b,c)
      implicit none
      integer n,m,i,j,ij,k,l,kl,nm
      integer*1 b(0:*),c(0:*)
      nm=n*m
      ij=-1
      do i=0,n-1
       do j=0,m-1
        ij=ij+1
        kl=-1
        do k=0,n-1
         do l=0,m-1
          kl=kl+1
          if(i.eq.k)c(nm*kl+ij)=c(nm*kl+ij)+b(m*l+j)
         enddo
        enddo
       enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine dirprod(n,v,m,w,vw)
      implicit none
      integer n,m,i,j,ij
      real*8 v(0:*),w(0:*),vw(0:*)
      ij=-1
      do i=0,n-1
       do j=0,m-1
        ij=ij+1
        vw(ij)=v(i)*w(j)
       enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine matveci(n,s,a,b)
      implicit none
      integer:: n,i,j,inc=1
      integer*1 s(0:*),svalue
      real*8 a(0:*),b(0:*)
      real*8 :: one=-1.0d0,beta=0.d0
      call dgemvi(n,n,one,s,n,a,inc,beta,b,inc)
c      do i=0,n-1
c       b(i)=0.d0
c       do j=0,n-1
c        svalue=s(n*j+i)
c        if(svalue.gt.0)b(i)=b(i)-dble(svalue)*a(j)
c       enddo
c      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copymat(n,a,b)
      implicit none
      integer n,i,j,nj
      integer*1 a(0:*),b(0:*)
      do i=0,n-1
       do j=0,n-1
        nj=n*j+i
        a(nj)=b(nj)
       enddo
      enddo
      return
      end
c------------------------------------------------------------------------
      function wigner(n,s2)
      implicit none
      integer wigner,n,s2,comb
      wigner=comb(n,(n-s2)/2)-comb(n,(n-s2-2)/2)
      return
      end
c------------------------------------------------------------------------
      recursive integer function comb(n,m) result(c)
      implicit none
      integer n,m
      if(m.lt.0)then
         c=0
         return
      endif
      if(m.eq.0)then
         c=1
         return
      endif
      c=(n-m+1)*comb(n,m-1)/m
      return
      end
c-------------------------------------------------------------------
      subroutine esame(g,n,ndim,s2,ntrue,inum)
      implicit none
      integer g,n,ndim,s2,ntrue,inum(0:*),i,j,n1,n0,howmany1
      ntrue=0
      do i=0,ndim-1
       n1=howmany1(ndim-i-1)
       n0=n-n1
       if((n1-n0).eq.s2)then
          inum(ntrue)=i
          ntrue=ntrue+1
       endif
      enddo
      return
      end
c----------------------------------------------------------------
      recursive function howmany1(n) result(h)
      implicit none
      integer howmany,h,n
      if(n.eq.0)then
         h=0
         return
      elseif(n.eq.1)then
         h=1
         return
      else
         h=mod(n,2)+howmany1(n/2)
         return
      endif
      end
c-------------------------------------------------------
      function combf(n,m)
      implicit none
      integer combf,n,m,comb
      combf=comb(n,m)
      return
      end
c--------------------------------------------------------
      subroutine getiocc(n,m,word)
      implicit none
      integer n,m,word(0:*),k
      do k=0,n-1
       if(iand(ishft(1,k),m).ne.0)then
          word(k)=1
       else
          word(k)=0
       endif
      enddo
      call reversi(n,word)
      return
      end
c----------------------------------------------------
      subroutine reversi(n,word)
      implicit none
      integer n,word(0:*),i,c
      i=0
      do while(i.lt.n/2)
       c=word(i)
       word(i)=word(n-i-1)
       word(n-i-1)=c
       i=i+1
      enddo
      return
      end
c-----------------------------------------------
      subroutine sixnum(inum,iaux,nopen,ivec6)
      implicit none
      integer inum,iaux(0:*),nopen,ivec6(0:*)
      integer :: i,j,ii,il,iu=-1,newnum=0
      integer, dimension(0:3) :: iauxp
      integer possibility(0:5,0:3)
      data possibility /1,0,1,0,1,0, 1,0,0,1,0,1, 0,1,1,0,0,1,
     $     0,1,0,1,1,0 /
      iauxp=(/ nopen+5-iaux(3),nopen+5-iaux(2),nopen+5-iaux(1),nopen+5
     $     -iaux(0) /)
      do j=0,3
       il=iu+1
       iu=iauxp(j)-1
       do i=il,iu-1
        if(iand(ishft(1,i-j),inum).ne.0) newnum=ior(newnum,ishft(1,i))
       enddo
      enddo
      il=iu+1
      iu=nopen+4
      do i=il,iu-1
       if(iand(ishft(1,i-4),inum).ne.0) newnum=ior(newnum,ishft(1,i))
      enddo
      do j=0,5
       ivec6(j)=newnum
       do i=0,3
        ii=nopen+4-iaux(i)
        if(possibility(j,i).eq.1)ivec6(j)=ior(ivec6(j),ishft(1,ii))
       enddo
      enddo
      return
      end
c--------------------------------------------------------------
      subroutine twonum(inum,iaux,nopen,ivec2)
      implicit none
      integer inum,iaux(0:*),nopen,ivec2(0:*)
      integer :: i,j,ii,il,iu=-1,newnum=0
      integer, dimension(0:1) :: iauxp
      integer, dimension(0:1,0:1) :: possibility2
      data possibility2/1,0, 0,1 /
      iauxp=(/nopen+3-iaux(1),nopen+3-iaux(0)/)
      do j=0,1
       il=iu+1
       iu=iauxp(j)-1
       do i=il,iu-1
        if(iand(ishft(1,i-j),inum).ne.0) newnum=ior(newnum,ishft(1,i))
       enddo
      enddo
      il=iu+1
      iu=nopen+2
      do i=il,iu-1
       if(iand(ishft(1,i-2),inum).ne.0) newnum=ior(newnum,ishft(1,i))
      enddo
      do j=0,1
       ivec2(j)=newnum
       do i=0,1
        ii=nopen+2-iaux(i)
        if(possibility2(j,i).eq.1) ivec2(j)=ior(ivec2(j),ishft(1,ii))
       enddo
      enddo
      return
      end


