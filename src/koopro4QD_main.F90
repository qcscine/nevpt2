program koopro4QD_main

use koopro4QD
use nevpt_header
#ifdef _OPENMP_
use omp_lib
#endif

implicit none

  logical              :: rel_ham
  logical              :: input_found
  integer,dimension(8) :: values
  real*4               :: t13
#ifdef _OPENMP_
  real*8               :: walltime_start,walltime_end
#endif

!
!     This program calculates the second order correction to the energy
!     using the n-electron valence state perturbation theory (NEVPT)
!     in the partially contracted (PC) ad strongly contracted (SC)
!     variants using the Dyall's Hamiltonian for the definition of the
!     zero order energies. All equations are based on the spin free
!     formulation of NEVPT.
!
!     Relevant papers:
!     1) C. Angeli, R. Cimiraglia, J-P Malrieu,
!        Chem. Phys. Lett.,  317, 472-480, (2000)
!        Definition of the inactive (core+virtual) orbital energies.
!     2) C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J-P Malrieu,
!        J. Chem. Phys., 114, 10252, (2001)
!        Introduction of NEVPT
!     3) C. Angeli, R. Cimiraglia, J-P Malrieu,
!        Chem. Phys. Lett., 350, 297-305, (2001)
!        Introduction of the formalism based on the "Koopmans' operators"
!     4) C. Angeli, R. Cimiraglia, J-P Malrieu,
!        J. Chem. Phys., 117, 10525-10533, (2002)
!        Spin free version of NEVPT

      !> set non-relativistic version by default - make it input option in model later
      rel_ham = .false.

      call read_menu_input('input.nevpt', 88, '*MODEL ', input_found)
      call read_menu_input('input.nevpt', 88, '*KOOPRO', input_found)
      if(.not. input_found) stop 'ERROR: no input for NEVPT2 calculation given...'

      call print_nevpt_header(rel_ham)

      call date_and_time(VALUES=values)

      print '(a )',  ' -----------------------------------------------------------------------'
      print '(a  )', '                     Entering koopro4QD_driver                          '
      print '(a,i2,a,i2,a,i4)','                      date: ',values(3),'/',values(2),'/',values(1)
      print '(a,i2,a,i2,a,i2)','                      time: ',values(5),':',values(6),':',values(7)
      print '(a /)', ' -----------------------------------------------------------------------'

      print '(/a/)', ' content of input file: '
      call system('cat input.nevpt')
      print '(/a)', ' '

#ifdef _OPENMP_
      walltime_start = omp_get_wtime()
#endif
      if (no_4rdm_terms) then
        write (6,*) 'WARNING: .NO4RDM has been specified. A-D and A~-D~ matrices will not be calculated. '
      end if
      !> working routine
      ! the last parameter checks whether we're calling the routine from the MOLCAS interface, which we're
      ! not doing here
      call koopro4QD_driver(rel_ham,t13,.false.)

#ifdef _OPENMP_
      walltime_end = omp_get_wtime()
#endif

      ! using keyword arguments
      call date_and_time(VALUES=values)
      print '(//a)',  ' -----------------------------------------------------------------------'
      print '( a,f16.2,a/)','             Final timing of koopro4QD_driver: ',t13,' sec.             '
#ifdef _OPENMP_
      print '( a,f16.2,a/)','             koopro4QD_driver wall time:     : ',walltime_end - walltime_start,' sec.'
#endif
      print '(a,i2,a,i2,a,i4)','                            date: ',values(3),'/',values(2),'/',values(1)
      print '(a,i2,a,i2,a,i2)',    '                            time: ',values(5),':',values(6),':',values(7)
      print '(a  )','                              normal exit                              '
      print '(a//)',  ' -----------------------------------------------------------------------'

end program koopro4QD_main
