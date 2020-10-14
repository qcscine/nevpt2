
      SUBROUTINE CPUT(IOP)
C
C     DATA, ORA E TEMPO DI CPU
C
C     IOP=-1 ... STAMPA DATA E ORA E INIZIALIZZA I TEMPI
C     IOP= 0 ... CALCOLA I TEMPI TOTALE E PARZIALE (TCPU E PCPU),
C                NON STAMPA
C     IOP= 1 ... STAMPA DATA, ORA E TEMPI
C     IOP= 2 ... STAMPA  TEMPI SOLAMENTE
C
      character*24 dat
      CHARACTER*8 ORA
      character*4 year
      character*3 month
      character*2 day
      REAL*8 PCPU,TCPU
      COMMON /CPU/ TCPU,PCPU,tcpu0
      dimension tarray(2)
      save cpu
c      i=mclock()
      t=etime(tarray)
c     call system('date>file.di.data')
c     open(99,file='file.di.data',form='formatted',status='old')
c     read(99,'(a)')dat
c     close(unit=99,status='delete')
      call gdate(dat)
      month=dat(5:7)
      day=dat(9:10)
      ora=dat(12:19)
      year=dat(21:24)
      IF (IOP.EQ.-1.OR.IOP.EQ.1) WRITE (6,10) month,day,year,ora
   10 FORMAT (/5X,'*****   DATE:  ',A3,', ',A2,', ',A4,
     * '   TIME: ',A8,'   *****')
      if (iop.lt.0) tcpu0=t    !i*1.d-2
      if (iop.lt.0) pcpu=0.d0
      if (iop.ge.0) then
      tcpu=t    !i*1.d-2
      pcpu=tcpu-tcpu0
      tcpu0=tcpu
      if(iop.eq.1)WRITE (6,20) PCPU,TCPU
      if(iop.eq.2)print '(''Elapsed time= '',f9.2,'' total = '',f9.2)'
     *,pcpu,tcpu
   20 FORMAT (5X,'*****   CPU TIME:  PARTIAL...',F9.2,' SEC.',
     * '  TOTAL...',F9.2,' SEC.   *****'/)
      ENDIF
      RETURN
      END
C***********************************************************************
      Subroutine GDate(Date1)
      Implicit Integer(A-Z)
C
C     This wrapper routine either calls FDate (on bsd systems) or
C     builds the 24-character date in some other way.
C
      Character*(*) Date1
C
      Character*26 LDate
      LDate = ' '
      Junk = GCTime(LDate)
      Date1 = LDate
      Return
      End
