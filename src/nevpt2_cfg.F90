module nevpt2_cfg

! This module contains input options and some important default parameters
  implicit none

  real*8,  public              :: thr        = 1.0d-10
  real*8,  public              :: thr_lindep = 1.0d-6

  integer, public              :: print_level         =  0
  integer, public              :: nr_active_electrons = -1
  integer, public              :: nr_states           = -1
  integer, public              :: ncore               =  0
  integer, public              :: nact                =  0
  integer, public              :: nspin               = -1
  integer, public              :: nr_frozen_orb       =  0
  integer, allocatable, public :: igelo(:)

  logical, public              :: zorder              = .true.

  logical, public              :: do_dmrg_pt          = .false.
  logical, public              :: do_cholesky         = .false.
  logical, public              :: compute_nooccn      = .false.
  logical, public              :: compute_rho1st      = .false.
  logical, public              :: no_pc               = .false.
  logical, public              :: skip_effective_ham  = .false.
  ! Skip terms in Koopmans matrices which require 4-RDM (usually useful only for debugging)
  logical, public              :: no_4rdm_terms       = .false.
  ! Skip calculating Koopmans' matrices if this option is activated
  ! in the MOLCAS interface. No effect in the standalone version.
  logical, public              :: skip_koopro_molcas  = .false.

#ifdef DMRG_NEVPT
  ! Read 4- and 3-RDMs at the beginning of the calculation from QCMaquis
  ! rather than calculating them
  logical, public              :: rdm_read = .false.
  ! read the RDMs from series of files from the distributed calculation
  logical, public              :: rdm_distributed = .false.
  ! Path (relative to WorkDir or absolute) to distributed RDM calculation results
  ! Uses allocation on assignment
  character(len=:),allocatable :: rdm_path

  ! Current working directory with QCMaquis checkpoints
  character(len=255) :: curr_dir
  character(len=255) :: molcas_project
#endif
  character(len=80), public    :: file04    = ' '
  character(len=80), public    :: file50    = ' '

  type states
    integer,            allocatable :: state(:)
    character(len=256), allocatable :: h5_file_name(:)
  end type
  type(states), public         :: MultGroup




end module nevpt2_cfg
