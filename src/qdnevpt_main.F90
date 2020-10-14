program qdnevpt_main

use nevpt_header
use qdnevpt_core
use info_state_energy
use nevpt2_cfg
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
      call read_menu_input('input.nevpt', 88, '*QDNEVP', input_found)
      if(.not. input_found) stop 'ERROR: no input for NEVPT2 calculation given...'

      call print_nevpt_header(rel_ham)

      call date_and_time(VALUES=values)

      print '(a )',  ' ------------------------------------------------------------------------'
      print '(a  )', '                      Entering qdnevpt                                   '
      print '(a,i2,a,i2,a,i4)','                       date: ',values(3),'/',values(2),'/',values(1)
      print '(a,i2,a,i2,a,i2)','                       time: ',values(5),':',values(6),':',values(7)
      print '(a /)', ' ------------------------------------------------------------------------'

      print '(/a/)', ' content of input file: '
      call system('cat input.nevpt')
      print '(/a)', ' '

#ifdef _OPENMP_

      print '(a,i4,a)', ' Using OpenMP with max.', omp_get_max_threads(), ' threads'

      walltime_start = omp_get_wtime()
#endif
      !> init energies
      call init_energies(nr_states)
      !> working routine
      ! Leon: the third parameter is whether we're calling it from within MOLCAS, which isn't the case
      call qdnevpt(rel_ham,t13,.false.)

      !> uninitialize routines
      call deinit_energies()
#ifdef _OPENMP_
      walltime_end = omp_get_wtime()
#endif

      call date_and_time(VALUES=values)
      print '( a/)',  ' ------------------------------------------------------------------------'
      print '( a,f16.2,a/)','              Final timing of qdnevpt         : ',t13,' sec.'
#ifdef _OPENMP_
      print '( a,f16.2,a/)','              qdnevpt wall time:              : ',walltime_end - walltime_start,' sec.'
#endif
      print '(a,i2,a,i2,a,i4)','                             date: ',values(3),'/',values(2),'/',values(1)
      print '(a,i2,a,i2,a,i2)',    '                             time: ',values(5),':',values(6),':',values(7)
      print '(a  )','                               normal exit                              '
      print '(a//)',  ' ------------------------------------------------------------------------'


end program qdnevpt_main
