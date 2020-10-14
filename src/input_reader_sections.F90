  subroutine read_input_sections(word, kw_section)

    use input_reader

    implicit none

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: word
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------

    select case (kw_section)

      case ('*KOOPRO')
        call read_input_koopro(word, kw_section)
      case ('*QDNEVP')
        call read_input_qdnevpt(word, kw_section)
      case ('*MODEL')
        call read_input_model(word, kw_section)
      case default
        print *, 'ERROR: section '//kw_section//' not recognized'
        stop 'unknown keyword section'

    end select

  end subroutine

  subroutine read_input_koopro(word, kw_section)

     use input_reader
     use nevpt2_cfg

     implicit none

!    --------------------------------------------------------------------------
     character(kw_length), intent(in) :: word
     character(kw_length), intent(in) :: kw_section
!    --------------------------------------------------------------------------
     integer                          :: i

     call reset_available_kw_list()

     if (kw_matches(word, '.THRESH')) then
         call kw_read(word, thr)
     end if

     if (kw_matches(word, '.SPIN  ')) then
         call kw_read(word, nspin)
     end if

     if (kw_matches(word, '.STATES')) then
         call kw_read(word, nr_states)
     end if

     if (kw_matches(word, '.NACTEL')) then
         call kw_read(word, nr_active_electrons)
     end if

     !> only needed for non-MOLCAS driven (e.g. old interface-based) version
     if (kw_matches(word, '.NCORE ')) then
         call kw_read(word, ncore)
     end if

     !> only needed for non-MOLCAS driven (e.g. old interface-based) version
     if (kw_matches(word, '.NACTOR')) then
         call kw_read(word, nact)
     end if

     if (kw_matches(word, '.FILE04')) then
         call kw_read(word, file04)
     end if

     if (kw_matches(word, '.FILE50')) then
         call kw_read(word, file50)
     end if

     if (kw_matches(word, '.NO ZOR')) then
         zorder = .false.
     end if

     call check_whether_kw_found(word, kw_section)

  end subroutine

  subroutine read_input_qdnevpt(word, kw_section)

     use input_reader
     use nevpt2_cfg

     implicit none

!    --------------------------------------------------------------------------
     character(kw_length), intent(in) :: word
     character(kw_length), intent(in) :: kw_section
!    --------------------------------------------------------------------------

     call reset_available_kw_list()

     if (kw_matches(word, '.LINDEP')) then
         call kw_read(word, thr_lindep)
     end if

     if (kw_matches(word, '.STATES')) then
         call kw_read(word, nr_states)
     end if

     !> freeze means
     !> - the actual frozen orbitals since also integrals
     !    with indices involving frozen orbs have been transformed
     if (kw_matches(word, '.FROZEN')) then
         call kw_read(word, nr_frozen_orb)
         allocate(igelo(nr_frozen_orb))
         igelo = 0; call kw_read(word, igelo, nr_frozen_orb)
     end if

     if (kw_matches(word, '.FILE04')) then
         call kw_read(word, file04)
     end if

     if (kw_matches(word, '.NATORB')) then
         compute_nooccn = .true.
         compute_rho1st = .true.
     end if

     if (kw_matches(word, '.NO PC ')) then
         no_pc = .true.
     end if

     call check_whether_kw_found(word, kw_section)

  end subroutine

  subroutine read_input_model(word, kw_section)

     use input_reader
     use nevpt2_cfg

     implicit none

!    --------------------------------------------------------------------------
     character(kw_length), intent(in) :: word
     character(kw_length), intent(in) :: kw_section
!    --------------------------------------------------------------------------

     call reset_available_kw_list()

     if (kw_matches(word, '.DMRGPT')) then
         do_dmrg_pt  = .true.
     end if

     if (kw_matches(word, '.CHOLES')) then
         do_cholesky = .true.
     end if
     
     if (kw_matches(word, '.NO4RDM')) then
         no_4rdm_terms = .true.
     end if

     call check_whether_kw_found(word, kw_section)

  end subroutine
