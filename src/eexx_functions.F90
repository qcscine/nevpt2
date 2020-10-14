module eexx_functions

#ifdef _OPENMP_
  use omp_lib
#endif
  use indices_norm
  use ord_utils, only : ord8

  implicit none

  public eeee
  public eeeet
  public eee2
  public eee2t
  public eee2tl
  public ee2
  public ee2t
  public ee2tr
  public ro2t
  public ro3t

  contains

!------------------------------------------------------------------------------
      real*8 function eeee(a,b,c,d,e,f,g,h,d3,d2,d1,ro4,nact,istate,metat)

      integer, intent(in) :: a,b,c,d,e,f,g,h
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: d3(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in) :: d2(metat,nact,nact,nact,nact) 
      real*8 , intent(in) :: d1(metat,nact,nact)
      real*8 , intent(in) :: ro4(nwords,nact,nact,nact,nact,metat)
      integer             :: i1,i2,i3,i4,i5,i6,i7,i8
      integer             :: lindice1,i,j,k,l,ind
!----four-particle density matrix dealt with here
!--a,c,e,g    b,d,f,h

      lindice1(i,j,k,l)=ind_norm%lindi(i)+ind_norm%lindj(j)+ind_norm%lindk(k)+l  ! statement function in Fortran

      call ord8(a,c,e,g,b,d,f,h,i1,i2,i3,i4,i5,i6,i7,i8)

      ind=lindice1(i1,i2,i3,i4)  ! use of the statement function above!
      eeee=ro4(ind,i5,i6,i7,i8,istate)
      if(b.eq.c)eeee=eeee+d3(istate,a,e,g,d,f,h)
      if(b.eq.e)eeee=eeee+d3(istate,a,c,g,f,d,h)
      if(b.eq.g)eeee=eeee+d3(istate,a,c,e,h,d,f)
      if(d.eq.e)eeee=eeee+d3(istate,a,g,c,b,h,f)
      if(d.eq.g)eeee=eeee+d3(istate,a,c,e,b,h,f)
      if(f.eq.g)eeee=eeee+d3(istate,a,e,c,b,h,d)
      if(b.eq.c)then
         if(d.eq.e)eeee=eeee+d2(istate,a,g,f,h)
         if(d.eq.g)eeee=eeee+d2(istate,a,e,h,f)
         if(f.eq.g)eeee=eeee+d2(istate,a,e,d,h)
      endif
      if(b.eq.e.and.d.eq.g)eeee=eeee+d2(istate,a,c,f,h)
      if(b.eq.e.and.f.eq.g)eeee=eeee+d2(istate,a,c,h,d)
      if(b.eq.g.and.d.eq.e)eeee=eeee+d2(istate,a,c,h,f)
      if(d.eq.e.and.f.eq.g)eeee=eeee+d2(istate,a,c,b,h)
      if(b.eq.c.and.d.eq.e.and.f.eq.g)eeee=eeee+d1(istate,a,h)

      end function eeee
!------------------------------------------------------------------------------

      real*8 function eeeet(a,b,c,d,e,f,g,h,d3,d2,d1,ro4,nact,istate,metat)

      integer, intent(in) :: a,b,c,d,e,f,g,h
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: d3(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in) :: d2(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: d1(metat,nact,nact)
      real*8 , intent(in) :: ro4(nwords,nact,nact,nact,nact,metat)

      eeeet=-eeee(a,b,d,c,e,f,g,h,d3,d2,d1,ro4,nact,istate,metat)
      if(c.eq.d)eeeet=eeeet+2.d0*eee2(a,b,e,f,g,h,d3,d2,d1,nact,istate,metat)

      end function eeeet
!------------------------------------------------------------------------------

      real*8 function eee2(a,b,c,d,e,f,daaa,taa,dal,nact,istate,metat)

      integer, intent(in) :: a,b,c,d,e,f
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: daaa(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in) :: taa(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: dal(metat,nact,nact)
!--the same as eee but with the spinless density matrices as input

      eee2=daaa(istate,a,c,e,b,d,f)
      if(e.eq.d)eee2=eee2+taa(istate,a,c,b,f)
      if(e.eq.b)eee2=eee2+taa(istate,a,c,f,d)
      if(b.eq.c)eee2=eee2+taa(istate,a,e,d,f)
      if(b.eq.c.and.e.eq.d)eee2=eee2+dal(istate,a,f)

      end function eee2
!------------------------------------------------------------------------------

      real*8 function eee2t(a,b,c,d,e,f,daaa,taa,dal,nact,istate,metat)

      integer, intent(in) :: a,b,c,d,e,f
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: daaa(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in) :: taa(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: dal(metat,nact,nact)

      eee2t=-eee2(a,b,d,c,e,f,daaa,taa,dal,nact,istate,metat)
      if(c.eq.d)eee2t=eee2t+2.d0*ee2(a,b,e,f,taa,dal,nact,istate,metat)

      end function eee2t
!------------------------------------------------------------------------------

      real*8 function eee2tl(a,b,c,d,e,f,daaa,taa,dal,nact,istate,metat)

      integer, intent(in) :: a,b,c,d,e,f
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: daaa(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in) :: taa(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: dal(metat,nact,nact)

      eee2tl=-eee2(b,a,c,d,e,f,daaa,taa,dal,nact,istate,metat)
      if(a.eq.b)eee2tl=eee2tl+2.d0*ee2(c,d,e,f,taa,dal,nact,istate,metat)

      end function eee2tl
!------------------------------------------------------------------------------

      real*8 function ee2(a,b,c,d,taa,dal,nact,istate,metat)

      integer, intent(in) :: a,b,c,d
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: taa(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: dal(metat,nact,nact)

      ee2=taa(istate,a,c,b,d)
      if(b.eq.c)ee2=ee2+dal(istate,a,d)

      end function ee2
!------------------------------------------------------------------------------

      real*8 function ee2t(a,b,c,d,taa,dal,nact,istate,metat)

      integer, intent(in) :: a,b,c,d
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: taa(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: dal(metat,nact,nact)

      ee2t=-ee2(b,a,c,d,taa,dal,nact,istate,metat)
      if(a.eq.b)ee2t=ee2t+2.d0*dal(istate,c,d)

      end function ee2t
!------------------------------------------------------------------------------

      real*8 function ee2tr(a,b,c,d,taa,dal,nact,istate,metat)

      integer, intent(in) :: a,b,c,d
      integer, intent(in) :: nact, istate, metat
      real*8 , intent(in) :: taa(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: dal(metat,nact,nact)

      ee2tr=-ee2(a,b,d,c,taa,dal,nact,istate,metat)
      if(c.eq.d)ee2tr=ee2tr+2.d0*dal(istate,a,b)

      end function ee2tr
!------------------------------------------------------------------------------

      real*8 function ro3t(a,b,c,ap,bp,cp,d3,d2,d1,nact,istate,metat)

      integer, intent(in) :: a,b,c,ap,bp,cp
      integer, intent(in) :: nact,istate,metat
      real*8 , intent(in) :: d3(metat,nact,nact,nact,nact,nact,nact)
      real*8 , intent(in) :: d2(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: d1(metat,nact,nact)

      ro3t=-d3(istate,ap,bp,cp,a,b,c)
      if(ap.eq.a)ro3t=ro3t+2.d0*d2(istate,bp,cp,b,c)
      if(bp.eq.b)ro3t=ro3t+2.d0*d2(istate,ap,cp,a,c)
      if(cp.eq.c)ro3t=ro3t+2.d0*d2(istate,ap,bp,a,b)
      if(bp.eq.c)ro3t=ro3t-d2(istate,ap,cp,a,b)
      if(ap.eq.c)ro3t=ro3t-d2(istate,cp,bp,a,b)
      if(cp.eq.b)ro3t=ro3t-d2(istate,bp,ap,c,a)
      if(ap.eq.b)ro3t=ro3t-d2(istate,bp,cp,a,c)
      if(cp.eq.a)ro3t=ro3t-d2(istate,ap,bp,c,b)
      if(bp.eq.a)ro3t=ro3t-d2(istate,ap,cp,b,c)
      if(cp.eq.c.and.bp.eq.b)ro3t=ro3t-4.d0*d1(istate,ap,a)
      if(cp.eq.b.and.bp.eq.c)ro3t=ro3t+2.d0*d1(istate,ap,a)
      if(cp.eq.c.and.ap.eq.a)ro3t=ro3t-4.d0*d1(istate,bp,b)
      if(cp.eq.a.and.ap.eq.c)ro3t=ro3t+2.d0*d1(istate,bp,b)
      if(ap.eq.a.and.bp.eq.b)ro3t=ro3t-4.d0*d1(istate,cp,c)
      if(ap.eq.b.and.bp.eq.a)ro3t=ro3t+2.d0*d1(istate,cp,c)
      if(cp.eq.c.and.ap.eq.b)ro3t=ro3t+2.d0*d1(istate,bp,a)
      if(cp.eq.b.and.ap.eq.c)ro3t=ro3t-d1(istate,bp,a)
      if(cp.eq.c.and.bp.eq.a)ro3t=ro3t+2.d0*d1(istate,ap,b)
      if(cp.eq.a.and.bp.eq.c)ro3t=ro3t-d1(istate,ap,b)
      if(bp.eq.b.and.ap.eq.c)ro3t=ro3t+2.d0*d1(istate,cp,a)
      if(bp.eq.c.and.ap.eq.b)ro3t=ro3t-d1(istate,cp,a)
      if(bp.eq.b.and.cp.eq.a)ro3t=ro3t+2.d0*d1(istate,ap,c)
      if(bp.eq.a.and.cp.eq.b)ro3t=ro3t-d1(istate,ap,c)
      if(ap.eq.a.and.cp.eq.b)ro3t=ro3t+2.d0*d1(istate,bp,c)
      if(ap.eq.b.and.cp.eq.a)ro3t=ro3t-d1(istate,bp,c)
      if(ap.eq.a.and.bp.eq.c)ro3t=ro3t+2.d0*d1(istate,cp,b)
      if(ap.eq.c.and.bp.eq.a)ro3t=ro3t-d1(istate,cp,b)
      if(ap.eq.a.and.bp.eq.b.and.cp.eq.c)ro3t=ro3t+8.d0
      if(ap.eq.b.and.bp.eq.a.and.cp.eq.c)ro3t=ro3t-4.d0
      if(ap.eq.a.and.bp.eq.c.and.cp.eq.b)ro3t=ro3t-4.d0
      if(ap.eq.c.and.bp.eq.b.and.cp.eq.a)ro3t=ro3t-4.d0
      if(ap.eq.b.and.bp.eq.c.and.cp.eq.a)ro3t=ro3t+2.d0
      if(ap.eq.c.and.bp.eq.a.and.cp.eq.b)ro3t=ro3t+2.d0

      end function ro3t
!------------------------------------------------------------------------------

      real*8 function ro2t(ap,bp,a,b,d2,d1,nact,istate,metat)

      real*8 , intent(in) :: d2(metat,nact,nact,nact,nact)
      real*8 , intent(in) :: d1(metat,nact,nact)
      integer, intent(in) :: a,b,ap,bp
      integer, intent(in) :: nact,istate,metat

      ro2t=d2(istate,a,b,ap,bp)
      if(ap.eq.a)ro2t=ro2t-2.d0*d1(istate,b,bp)
      if(bp.eq.b)ro2t=ro2t-2.d0*d1(istate,a,ap)
      if(ap.eq.b)ro2t=ro2t+d1(istate,a,bp)
      if(bp.eq.a)ro2t=ro2t+d1(istate,b,ap)
      if(ap.eq.a.and.bp.eq.b)ro2t=ro2t+4.d0
      if(ap.eq.b.and.bp.eq.a)ro2t=ro2t-2.d0

      end function ro2t

!------------------------------------------------------------------------------

end module eexx_functions
