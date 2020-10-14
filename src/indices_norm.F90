module indices_norm

#ifdef _OPENMP_
  use omp_lib
#endif
  use swap_util

  implicit none

  interface norm
    module procedure norm_int, norm_array, norm_int32, norm_array32
  end interface

  public initialize_indices_norm
  public finalize_indices_norm
  public lindice

  integer, public              :: nwords

  type type_ind_norm
  integer, public, allocatable :: lindi(:)
  integer, public, allocatable :: lindj(:)
  integer, public, allocatable :: lindk(:)
  integer, public, allocatable :: normind(:,:,:,:)
  end type type_ind_norm

  type (type_ind_norm), save   :: ind_norm

  contains

  subroutine initialize_indices_norm(nasht)

    integer, intent(in) :: nasht

    integer             :: i, j, k, l

    allocate(ind_norm%lindi(nasht),ind_norm%lindj(nasht),ind_norm%lindk(nasht))
    allocate(ind_norm%normind(nasht,nasht,nasht,nasht))

    do i = 1, nasht
      ind_norm%lindi(i) = (i-1)*i*(i+1)*(i+2)/24
      ind_norm%lindj(i) = (i-1)*i*(i+1)/6
      ind_norm%lindk(i) = (i-1)*i/2
      do j = 1, nasht
        do k = 1, nasht
          do l = 1, nasht
            ind_norm%normind(i,j,k,l) = norm(i,j,k,l)
          end do
        end do
      end do
    end do

    nwords = nasht*(nasht+1)*(nasht+2)*(nasht+3)/24

  end subroutine initialize_indices_norm
! ----------------------------------------------------------------------

  subroutine finalize_indices_norm()
    deallocate(ind_norm%lindi,ind_norm%lindj,ind_norm%lindk)
    deallocate(ind_norm%normind)
  end subroutine finalize_indices_norm

! ----------------------------------------------------------------------
  integer function lindice(i,j,k,l)

  integer, intent(in) :: i, j, k, l

    lindice = ind_norm%lindi(i) +      &
              ind_norm%lindj(j) +      &
              ind_norm%lindk(k) +      &
              l
  end function lindice

! ----------------------------------------------------------------------

  integer function norm_array(nv_in)
    integer, dimension(4), intent(in) :: nv_in
    integer, dimension(4) :: nv
    integer               :: i,j,jmax,imax

    norm_array = 0


    nv = nv_in ! copy to avoid swapping indices in the input

    do i = 1,4

    imax = nv(i)
    jmax = i

      do j=i+1,4
        if(nv(j) > imax)then
          imax = nv(j)
          jmax = j
        end if
      end do

      if(jmax /= i) call swap(nv(i),nv(jmax))

    end do

    norm_array=(nv(1)-1)*nv(1)*(nv(1)+1)*(nv(1)+2)/24+ &
            (nv(2)-1)*nv(2)*(nv(2)+1)/6+            &
              nv(3)*(nv(3)-1)/2+nv(4)

  end function norm_array

  integer function norm_int(ia,ic,ie,ig)

  integer, intent(in)   :: ia, ic, ie, ig
  norm_int = norm_array((/ ia, ic, ie, ig /))

  end function norm_int


  integer*4 function norm_array32(nv_in)
    integer*4, dimension(4), intent(in) :: nv_in
    integer*4, dimension(4) :: nv
    integer*4               :: i,j,jmax,imax

    norm_array32 = 0
    nv = nv_in ! copy to avoid swapping indices in the input

    do i = 1,4

    imax = nv(i)
    jmax = i

      do j=i+1,4
        if(nv(j) > imax)then
          imax = nv(j)
          jmax = j
        end if
      end do

      if(jmax /= i) call swap(nv(i),nv(jmax))

    end do

    norm_array32=(nv(1)-1)*nv(1)*(nv(1)+1)*(nv(1)+2)/24+ &
            (nv(2)-1)*nv(2)*(nv(2)+1)/6+            &
              nv(3)*(nv(3)-1)/2+nv(4)

  end function norm_array32

  integer*4 function norm_int32(ia,ic,ie,ig)

    integer*4, intent(in)   :: ia, ic, ie, ig
    norm_int32 = norm_array32((/ ia, ic, ie, ig /))

  end function norm_int32

! ----------------------------------------------------------------------

end module indices_norm
