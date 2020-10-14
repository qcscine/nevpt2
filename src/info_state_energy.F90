module info_state_energy

  implicit none

  ! Former /VEC/ common block is now here. This allows to access energies also after we quit from qdnevpt()
  real*8, allocatable :: e(:), e2en(:,:), e2mp(:,:), psimp(:), psien(:)
  real*8, allocatable :: heff_sc(:,:), heff_pc(:,:)
  public init_energies,deinit_energies,get_state_energies

  contains

  subroutine init_energies(nstates)
    integer, intent(in) :: nstates
    allocate(e(nstates))
    allocate(e2en(nstates,nstates))
    allocate(e2mp(nstates,nstates))
    allocate(psimp(nstates))
    allocate(psien(nstates))
    allocate(heff_sc(nstates,nstates))
    allocate(heff_pc(nstates,nstates))
    e = 0;e2en = 0; e2mp = 0; psimp = 0; psien = 0
    heff_sc = 0; heff_pc = 0
  end subroutine

  subroutine deinit_energies()

    use nevpt2_cfg, only : MultGroup
    if (allocated(e)) deallocate(e)
    if (allocated(e2en)) deallocate(e2en)
    if (allocated(e2mp)) deallocate(e2mp)
    if (allocated(psimp)) deallocate(psimp)
    if (allocated(psien)) deallocate(psien)
    if (allocated(heff_sc)) deallocate(heff_sc)
    if (allocated(heff_pc)) deallocate(heff_pc)
    ! Leon: MultGroup is deallocated here, in case the energies are still needed
    ! after qdnevpt is finished. It's not very clean though (probably MultGroup
    ! should be here instead of nevpt2_cfg instead)
    if (allocated(MultGroup%State))   deallocate(MultGroup%State)
  end subroutine

  subroutine get_state_energies(energies,nstates,state_ptr)

  ! This routine reads reference energies from the RasOrb files produced by MOLCAS
  ! When we're called from MOLCAS, it's better to read the energies from the JobIph file instead
    integer, intent(in) :: nstates
    integer, intent(in) :: state_ptr(nstates)
    real*8,  intent(out):: energies(nstates)
    real*8              :: state_energy
    integer             :: i, j, k, l, m
    integer             :: ierr, ierr1, ierr2, ierr3
    integer             :: mylen
    character(len=  8)  :: cflag2
    character(len=100)  :: cflag3
    character(len=100)  :: guess
    character(len=100)  :: file_tag

    do m = 1, nstates

      file_tag = ''
      if(state_ptr(m) < 10)then
        write(file_tag,'(a1,i1)') '.',state_ptr(m)
      else if(state_ptr(m) < 100)then
        write(file_tag,'(a1,i2)') '.',state_ptr(m)
      end if

      call system('ls *.RasOrb'//trim(file_tag)//' > STATE_GSS_NEVPT')

      open(unit=101,file="STATE_GSS_NEVPT",iostat=ierr)
      if(ierr.ne.0)then
        write(*,*)" RASSCF orbital file missing in scratch for state # ",state_ptr(m)
        stop
      end if
      read(101,*)guess
      close(101)
      open(unit=102,file=trim(guess))
      cflag2=""

      state_energy = 0.0d0

      do
        read(102,"(A5)",iostat=ierr1)cflag2
        if(ierr1.ne.0)then
          exit
        end if
        if(cflag2(1:5).eq."#INFO")then
          cflag3=""
          read(102,"(A100)",iostat=ierr2)cflag3
          if(ierr2 /= 0) exit
          mylen = len(trim(cflag3))
          do j=1,len(trim(cflag3))
            if(cflag3(j:j).eq."=")then
              read(cflag3(j+1:mylen),*) state_energy
              exit
            else
              cycle
            end if
          end do
        end if
      end do
      close(102)

      !print *, 'state-specific energy for state # ',m,' = ',state_energy

      energies(m) = state_energy
    end do

  end subroutine get_state_energies
! ----------------------------------------------------------------------

end module info_state_energy
