module info_symmetry

  use info_orbital_space 

  implicit none

  public initialize_infsym
  public finalize_infsym

  integer, allocatable, public :: its(:,:)
  integer, allocatable, public :: itsym(:)

  contains

  subroutine initialize_infsym(nsym,norb)

    integer, intent(in) :: nsym
    integer, intent(in) :: norb

    integer             :: i, j, k, l
    integer             :: muld2h

    !> statement function
    muld2h(i,j)= IEOR(i-1,j-1)+1

    allocate(its(nsym,nsym),itsym(norb+norb+1))
    its = 0; itsym = -1

    !print *, 'norb and nsym are ...',norb, nsym

    !> set its array
    do i = 1, nsym
      do j = 1, nsym
        its(i,j) = muld2h(i,j)
      end do
    end do

    k     = 0
    !> inactive
    inforb%nisht = 0
    do i = 1, nsym
      if(inforb%nish(i) == 0) cycle
      do j = 1, inforb%nish(i)
        k             = k + 1
        itsym(k)      = i
        inforb%nisht  = inforb%nisht + 1
      end do
    end do
    !> active
    inforb%nasht = 0
    do i = 1, nsym
      if(inforb%nash(i) == 0) cycle
      do j = 1, inforb%nash(i)
        k             = k + 1
        itsym(k)      = i
        inforb%nasht  = inforb%nasht + 1
      end do
    end do
    !> secondary (virtual)
    inforb%nssht = 0
    do i = 1, nsym
      if(inforb%nssh(i) == 0) cycle
      do j = 1, inforb%nssh(i)
        k             = k + 1
        itsym(k)      = i
        inforb%nssht  = inforb%nssht + 1
      end do
    end do

    if(k /= norb) print *, 'oooppppsss orbitals are not assigned!', k, norb
    if(k /= norb) stop -1
    
  end subroutine initialize_infsym
! ----------------------------------------------------------------------

  subroutine finalize_infsym()
    deallocate(its,itsym)
  end subroutine finalize_infsym

! ----------------------------------------------------------------------

end module info_symmetry
