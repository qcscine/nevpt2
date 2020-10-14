!
!-----------------------------------------------------------
! FILESPLIT  FILESPLIT  FILESPLIT  FILESPLIT  FILESPLIT 
!
! Ecriture d'un fichier en plusieurs morceaux
! Lecture correspondante
!
! Fonctionnement : on ouvre normalement le fichier 'bid'
! on ecrit. quan bid devient trop long, on le ferme, et on
! ouvre le fichier bid_1, etc...  limite bid_9999
! idem pour lire
!
! Il faut ajouter dans le code fortran (77) les lignes :
!
! - APRES open :
!   call filesplit('OPEN',num_fil,maxlong,idum,idum2)
!   num_fil = numero de la file
!   maxlong = longueur max d'un fichier 
!   (en general 2Gb, maxlong=2000000)

! - APRES close(num_fil,status='delete')
!   call filesplit('DELETE',num_fil,0,idum,idum2)

! ACCES SEQUENTIEL :

! - avant ou apres rewind :
!   call filesplit('REWIND',num_file,0,idum,idum2)
!   num_fil = numero de la file

! - AVANT read ou write :
!   call filesplit('READ' ou 'WRITE',num_fil,long,ichange,idum2)
!    num_fil = numero de la file
!    long = longueur de l'enregistrement, en bytes
!    ichange=0, =1 si on a fait une operation close/open

! ACCES DIRECT :

! - AVANT write :
!   call filesplit('WRITE',num_fil,long,ichange,idum2)
!    num_fil = numero de la file
!    long = longueur de l'enregistrement, en bytes
!    ichange=0, =1 si on a fait une operation close/open 
!    (variable en sortie)


! AVANT read :
!   call filesplit('READ',num_fil,long,ichange,irec)
!    num_fil = numero de la file
!    long = longueur de l'enregistrement, en bytes
!    ichange=0, =1 si on a fait une operation close/open
!    irec : numero de l'enregistrement a lire
!    (l'instruction qui suit est : "read(num_fil,rec=irec...)"

!-----------------------------------------------------------
!
!
! Voici un exemple de programme principal appelant filesplit
!      dimension a(1010)
!      logical yw
!      character*80 file
!      yw=.false.
!      file='bid'
!      yw=.false.
!      print*,'*************OPEN****************'
!      open(1,file=file,form='unformatted',status='unknown')
!      call filesplit('OPEN',1,20,ichange,idum)
!      print*,'***********ECRITURE**************'
!      do k=1,98
!         do l=1,1010
!            a(l)=1.d0*k
!         enddo
!         call filesplit('WRITE',1,4040,ichange,idum)
!         print*,'ecriture',k
!         write(1) a
!      enddo
!      print*,'************LECTURE**************'
!      do k=1,1010
!         a(k)=0.d0
!      enddo
!      call filesplit('REWIND',1,0,ichange,idum)
!      rewind 1
!      do k=1,98
!         call filesplit('READ',1,4040,ichange,idum)
!         read(1) a
!         print*,'A',a(100)
!      enddo
!      print*,'************DESTRUCTION**********'
!      close(1,status='delete')
!      call filesplit('DELETE',1,0,idum,idum)
!      end
! Voici un exemple de programme principal appelant filesplit 
! (acces direct)
c      dimension a(1000)
c      logical yw
c      yw=.false.
c      open(1,file='bid',form='unformatted',status='unknown',
c     *access='direct',recl=4000)
c      call filesplit('OPEN',1,20,ichange,idum)
c      irec=0
c      do k=1,100
c         do l=1,1000
c            a(l)=1.d0*k
c         enddo
c         call filesplit('WRITE',1,4000,ichange,idum)
c         irec=irec+1
c         if(ichange.eq.1) irec=1
c         print*,'ecriture k,irec=',k,irec
c         write(1,rec=irec) a
c      enddo
c      do k=1,100,9
c         do l=1,1000
c            a(l)=0.d0
c         enddo
c         irec=k
c          if(yw)print*,'irec',irec
c         call filesplit('READ',1,4000,ichange,irec)
c         print*,'irec apres sortie',irec
c         read(1,rec=irec) a
c         print*,'A',a(100)
c      enddo
c      end
! Voici un exemple de programme principal appelant filesplit
! (2 fichiers independants)
c      dimension a(1000),b(300)
c      logical yw
c      yw=.false.
c      open(1,file='bid',form='unformatted',status='unknown')
c      call filesplit('OPEN',1,20,ichange,idum)
c      open(50,file='bid2',form='unformatted',status='unknown')
c      call filesplit('OPEN',50,16,ichange)
c      do k=1,100
c         do l=1,1000
c            a(l)=1.d0*k
c         enddo
c         call filesplit('WRITE',1,4000,ichange,idum)
c         write(1) a
c      enddo
c      do k=1,1000
c         a(k)=0.d0
c      enddo
c      call filesplit('REWIND',1,0,ichange,idum)
c      rewind 1
c      do k=1,100
c         call filesplit('READ',1,4000,ichange,idum)
c         read(1) a
c         do l=1,300
c            B(l)=-A(l)
c         enddo
c         print*,'A',a(100)
c         call filesplit('WRITE',50,1200,ichange,idum)
c         write(50) b
c      enddo
c      print*,'coucou'
c      call filesplit('REWIND',50,0,ichange,idum)
c      rewind 50
c      do k=1,100
c         call filesplit('READ',50,1200,ichange,idum)
c         read(50) b
c         print*,'B',b(100)
c      enddo
c      end
      subroutine filesplit(action,if,long,ichange,irec)
c
c---------------------------------------------------------
c action : OPEN / REWIND / READ / WRITE / DELETE
c if     : (input)   numero du fichier
c long   : action=='OPEN' : longueur max d'un fichier, en bytes
!        : action=='READ/WRITE' : longueur de l'enregistrement, en bytes
! ichange : access=='DIRECT',action=='READ' : donne en sortie
!           irec, numero de l'enregistrement a lire
!         : sinon : 1 si on a fait un close/open, 0 sinon
!---------------------------------------------------------
!
      parameter (maxacf=10)

      dimension maxlong(maxacf),longtot(maxacf),
     *numfile(maxacf),filename(maxacf),access(maxacf)
     *,iactive_files(99)
      integer*8 maxlong,longtot
      save maxlong,longtot,
     *numfile,filename,access,init
     *,iactive_files,nactive_files
      character*300 filename,fileloc
      character*(*) action
      logical yw,ex
      character*80 seq,fm,access
      integer rcl
      character*4 suffix
      data init=0

      yw=.false.
      if(yw) then
         print*,' '
         print*,'entree filesplit'
         print*,'action ',action
         inquire(unit=if,name=fileloc)
         print*,'file',if,',   '//trim(fileloc)
      endif
      ichange=0
c
c Open ou rewind (initialisation des donnees des enregistrements)
c      
      if(action.ne.'OPEN') locacf=iactive_files(if)
      if(action.eq.'OPEN') then
         if(init.eq.0) then
            init=1
            do k=1,99
               iactive_files(k)=0
            enddo
            nactive_files=0
         endif

         nactive_files=nactive_files+1
         iactive_files(if)=nactive_files
         locacf=iactive_files(if)
         maxlong(locacf)=long
         maxlong(locacf)=maxlong(locacf)*1000
!         print*,'LONG(kb),locacf,laxlong(locacf)(bytes)'
!     *   ,long,locacf,maxlong(locacf)
         longtot(locacf)=0
         numfile(locacf)=0
         inquire(if,NAME=filename(locacf))
         inquire(if,access=access(locacf))
         if(yw)print*,'longueur max du fichier',maxlong(locacf),' bytes'
         if(yw)print*,'locacf,filename',locacf,filename(locacf)
      else if(action.eq.'REWIND') then
         longtot(locacf)=0
         numfile(locacf)=0
         inquire(if,NAME=fileloc)
         if(fileloc.ne.filename(locacf)) then
            if(yw)print*,'filename',filename(locacf)
            if(yw)print*,'fileloc',fileloc
            inquire(if,form=fm)
            inquire(if,RECL=rcl)
            if(yw)print*,'filesplit open0',filename(locacf)
            close(if)
c            write (6,*) 'now opening',filename(locacf)
c            call flush (6)
            open(if,form=fm,status='OLD',
     *      file=filename(locacf),access=access(locacf),
     *      recl=rcl)
         endif
c
c ordre read ou write
c
      else if(action.eq.'READ'.or.action.eq.'WRITE') then
         if(access(locacf).eq.'DIRECT') then
               ! DIRECT
               if(yw) print*,' '
               if(yw)print*,'filesplit read DAF.  irec en entree :',irec
c              on cherche dans quel fichier est l'enregistrement
c              et sous quel irec
               inquire(if,recl=rcl)
               if(yw)print*,'longueur de chaque enregistrement',long
               nb_enr_par_fich=maxlong(locacf)/long
               if(yw)print*,'nb_enr_par_fich',nb_enr_par_fich
               nnumfile=(irec-1)/nb_enr_par_fich
               if(yw)print*,'nnumfile',(irec-1),'/',nb_enr_par_fich,
     *         '=',nnumfile
               if(yw)print*,'le irec entre est sur la file n°',nnumfile
               if(yw)print*,'la file ouverte en ce moment est',
     *         numfile(locacf)
               numfile_entree=numfile(locacf)
               numfile_sortie=nnumfile
                if(yw)print*,'numfile_entree',numfile_entree
                if(yw)print*,'numfile_sortie',numfile_sortie
               if(nnumfile.ne.numfile(locacf)) then
                  numfile(locacf)=nnumfile
                  call filesplit_close_open(action,if,numfile(locacf),
     *            filename(locacf))
                  ichange=1
                  if(yw) then
                     print*,'DIR-READ close/open'
                     print*,'action',action
                     print*,'if',if
                     print*,'numfile(locacf)',numfile(locacf)
                     print*,'filename(locacf)',filename(locacf)(1:60)
                  endif
               endif
               if(yw)print*,'nb_enr_par_fich',(maxlong(locacf)),
     *         long,nb_enr_par_fich
               if(yw)print*,'irec',irec
               irec=irec-(numfile_sortie)*nb_enr_par_fich
               if(yw)print*,'numfile(locacf),irec en sortie',
     *         numfile(locacf),irec
         else
            !SEQ
            if(yw) then
               print*,trim(action)//
     *         ' : longueur de l''enregistrement en cours (bytes)'
     *         ,long
               print*,trim(action)//
     *         ' : longueur des enregistrements precedents (bytes)'
     *         ,longtot(locacf)
               print*,trim(action)//
     *         ' : longueur maximale d''un fichier (bytes)'
     *         ,maxlong(locacf)
               print*,trim(action)//
     *         ' : longueur totale lue ou ecrite  (bytes)'
     *         ,longtot(locacf)+long
            endif
            longtot(locacf)=longtot(locacf)+long
            if(long.gt.maxlong(locacf)) then
               print*,'filesplit :'
               print*,'l''enregistrement est plus grand '//
     *         'que la longueur max du fichier'
               print*,'file                        ',if
               print*,'locacf                      ',locacf
               print*,'longueur de l''enregistrement',long
               print*,'maxlong de la file          ',maxlong(locacf)
               print*,'action                      ',action
               print*,'filename(locacf)         ',filename(locacf)(1:60)
               print*,'numfile(locacf)             ',numfile(locacf)
               stop
            endif
            if(longtot(locacf).gt.maxlong(locacf)) then
               if(yw) then
                  print*,trim(action)//' close/open'
                  print*,'if',if
                  print*,'numfile(locacf)',numfile(locacf)
                  print*,'filename(locacf)',filename(locacf)
               endif
               call filesplit_close_open(action,if,numfile(locacf),
     *         filename(locacf))
               if(yw)print*,'numfile(locacf) apres C/O',numfile(locacf)
               if(action.eq.'WRITE') ichange=1
               longtot(locacf)=long
            endif
         endif
c
c delete
c
      else if(action.eq.'DELETE') then
         if(yw)print*,'filesplit DELETE'
         if(yw)print*,'numfile(locacf)',numfile(locacf)
         do kk=1,numfile(locacf)+1
            k=kk-1
            ! trouver le nom du fichier a detruire
            if(k.eq.0) then
               fileloc=filename(locacf)
               if(yw)print*,'a detruire k==0',locacf,fileloc
               goto 1
            else
               if(yw)print*,'a detruire k>0',k,locacf,numfile(locacf)
               numfile(locacf)=numfile(locacf)+1
               l=len_trimf77(fileloc)
               if(yw)print*,'l',l
               suffix='    '
               if(k.lt.10) then
                  write(suffix(2:2),'(I1)') k
               else if(k.lt.100) then
                  write(suffix(2:3),'(I2)') k
               else if(k.lt.100) then
                  write(suffix(2:4),'(I3)') k
               else
                  print*,'trop de files a detruire'
               endif
               suffix(1:1)='_'
              if(yw)print*,'suffix',suffix
               l=len_trimf77(filename(locacf))
              fileloc=filename(locacf)
              fileloc=fileloc(1:l)//suffix
            endif
            ! existe ?
            inquire(FILE=fileloc,EXIST=ex)
            if(yw)print*,'ex ',ex,fileloc
            ! destruction
            if(ex) then
c            write (6,*) 'now opening',fileloc
c            call flush (6)
               open(if,status='OLD',file=fileloc)
               close(if,status='delete')
            endif
 1          continue
         enddo
      endif
      return
      end
c
c
c
       subroutine filesplit_close_open(action,if,numfile,filename)
      character*(*) action
      character*300 fileloc,filename
      character*80 seq,fm,access
      integer rcl
      character*4 suffix
      character*10 format
      logical yw,yww
      yw=.false.
      yww=.false.
      if(yw)print*,'*************entree filesplit_close_open'
      if(yw)print*,'action,if,numfile,filename',
     *action,if,numfile,filename
         
c
! Fermer
c
      inquire(if,NAME=fileloc)
      if(yww)print*,'nom du fichier ouvert, donne par inquire',fileloc
      inquire(if,access=access)
      inquire(if,form=fm)
      inquire(if,RECL=rcl)
      if(yw)print*,'C/O : fermeture de ',fileloc(1:60)
      close(if)
         l=len_trimf77(fileloc)
         do k=l,l-3,-1
            if(fileloc(k:k).eq.'_') then
               fileloc(k:l)=' '
               goto 1
            endif
         enddo
 1       continue
c
! incrementer le nom
c
      if(yww)print*,'action,access',action,access
      if(access.eq.'SEQUENTIAL') then 
         !SEQ
         numfile=numfile+1
         call faisuffix(numfile,suffix)
         l=len_trimf77(fileloc)
         fileloc=fileloc(1:l)//suffix
      else
         !DIRECT
         if(yw)print*,'nom a incrementer',fileloc
         call faisuffix(numfile,suffix)
         if(yw)print*,'suffix',suffix
         l=len_trimf77(fileloc)
         if(yww)print*,'l',l
         fileloc=fileloc(1:l)//suffix
         if(yww)print*,'nouveau nom',fileloc
      endif
c
! ouvrir
c
      if(action.eq.'READ') then
         if(yw)print*,'fm',fm
         if(yw)print*,'access',access
         if(yw)print*,'rcl',rcl
         if(yw)print*,'C/O : ouverture de ',fileloc(1:60)
c            write (6,*) 'now opening',fileloc
c            call flush (6)
         open(if,form=fm,status='OLD',file=fileloc,access=access,
     *   recl=rcl)
      else
         if(yw)print*,'C/O : ouverture de ',fileloc(1:60)
c            write (6,*) 'now opening',fileloc
c            call flush (6)
         open(if,form=fm,status='unknown',file=fileloc,access=access,
     *   recl=rcl)
      endif
      return
      end
       function len_trimf77(namel)
      character*(*)namel
      l=len(namel)
      len_trimf77=0
      do k=l,1,-1
         if(namel(k:k).ne.' ') then
            len_trimf77=k
            return
         endif
      enddo
      return
      end
       subroutine maj2(a2,l)
      character*1 a2(l)
      do k=1,l
c         write(6,*)'k,ichar',k,ichar(a2(k))
         if(ichar(a2(k)).ge.97.and.ichar(a2(k)).le.122)
     *   a2(k)=char(ichar(a2(k))-32)
      enddo
      return
      end
       subroutine faisuffix(numfile,suffix)
      character*4 suffix
c
! Trouver suffix
c
      suffix='    '
      if(numfile.ne.0) suffix(1:1)='_'
      if(numfile.eq.0) then
      else if(numfile.lt.10) then
         write(suffix(2:2),'(i1)') numfile
      else if(numfile.lt.100) then
         write(suffix(2:3),'(i2)') numfile
      else if(numfile.lt.1000) then
         write(suffix(2:4),'(i3)') numfile
      else
         stop'filesplit nom trop long'
      endif
      !print*,'fin faisuffix numfile,suffix',numfile,suffix
      return
      end
