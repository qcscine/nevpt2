module swap_util

! taken from http://rosettacode.org/wiki/Generic_swap

  implicit none

  interface swap
    module procedure swapint, swapreal, swapstring, swapint32
  end interface

contains

  SUBROUTINE swapint(a, b)
    INTEGER, INTENT(IN OUT) :: a, b
    INTEGER :: temp
    temp = a ; a = b ; b = temp
  END SUBROUTINE swapint

  SUBROUTINE swapint32(a, b)
    INTEGER*4, INTENT(IN OUT) :: a, b
    INTEGER*4 :: temp
    temp = a ; a = b ; b = temp
  END SUBROUTINE swapint32

  SUBROUTINE swapreal(a, b)
    REAL, INTENT(IN OUT) :: a, b
    REAL :: temp
    temp = a ; a = b ; b = temp
  END SUBROUTINE swapreal

  SUBROUTINE swapstring(a, b)
    CHARACTER(*), INTENT(IN OUT) :: a, b
    CHARACTER(len(a)) :: temp
    temp = a ; a = b ; b = temp
  END SUBROUTINE swapstring
end module swap_util
