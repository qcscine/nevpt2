module ord_utils

#ifdef _OPENMP_
  use omp_lib
#endif
  use indices_norm

  implicit none

  public ord9
  public ord8

  contains

!------------------------------------------------------------------------------
      subroutine ord9(ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2,i1,i2,i3,i4,i5,i6,i7,i8)
!---  ordina in maniera decrescente ia1, ic1, ie1, ig1 in i1 i2 i3 i4
!---  ed fa la stessa perm anche su ia2, ic2, ie2, ig2 in i5 i6 i7 i8
      integer, intent(in)  :: ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2
      integer, intent(out) :: i1,i2,i3,i4,i5,i6,i7,i8
      integer              :: nv1(4), nv2(4)
      integer              :: i, j, imax, jmax, n
!      write (6,'(''ord9 in'',8I3)') ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2

      if(ind_norm%normind(ia1,ic1,ie1,ig1) >= ind_norm%normind(ia2,ic2,ie2,ig2))then
        nv1(1)=ia1
        nv1(2)=ic1
        nv1(3)=ie1
        nv1(4)=ig1
        nv2(1)=ia2
        nv2(2)=ic2
        nv2(3)=ie2
        nv2(4)=ig2
      else
        nv1(1)=ia2
        nv1(2)=ic2
        nv1(3)=ie2
        nv1(4)=ig2
        nv2(1)=ia1
        nv2(2)=ic1
        nv2(3)=ie1
        nv2(4)=ig1
      endif

      do i=1,4
         imax=nv1(i)
        do j=i+1,4
            if(nv1(j).gt.imax)then
               n=imax
               imax=nv1(j)
               nv1(j)=n
               n=nv2(i)
               nv2(i)=nv2(j)
               nv2(j)=n
            endif
            nv1(i)=imax
         enddo
      enddo
      i1=nv1(1)
      i2=nv1(2)
      i3=nv1(3)
      i4=nv1(4)
      i5=nv2(1)
      i6=nv2(2)
      i7=nv2(3)
      i8=nv2(4)
!      write (6,'(''ord9 ou'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8         
      end subroutine ord9
!------------------------------------------------------------------------------

      subroutine ord8(ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2,i1,i2,i3,i4,i5,i6,i7,i8)
!---  ordina in maniera decrescente ia1, ic1, ie1, ig1 in i1 i2 i3 i4
!---  e fa la stessa perm anche su ia2, ic2, ie2, ig2 in i5 i6 i7 i8
      integer, intent(in)  :: ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2
      integer, intent(out) :: i1,i2,i3,i4,i5,i6,i7,i8
      integer              :: nv1(4), nv2(4)
      integer              :: i, j, imax, jmax, n
!     write (6,'(''ord8 in'',8I3)') ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2

      if(ind_norm%normind(ia1,ic1,ie1,ig1) >= ind_norm%normind(ia2,ic2,ie2,ig2))then
        nv1(1)=ia1
        nv1(2)=ic1
        nv1(3)=ie1
        nv1(4)=ig1
        nv2(1)=ia2
        nv2(2)=ic2
        nv2(3)=ie2
        nv2(4)=ig2
      else
        nv1(1)=ia2
        nv1(2)=ic2
        nv1(3)=ie2
        nv1(4)=ig2
        nv2(1)=ia1
        nv2(2)=ic1
        nv2(3)=ie1
        nv2(4)=ig1
      endif

      do i=1,4
         imax=nv1(i)
         jmax=i
         do j=i+1,4
            if(nv1(j).gt.imax)then
               imax=nv1(j)
               jmax=j  
            endif
         enddo
         if (jmax.ne.i) then
         n=nv1(i)
         nv1(i)=nv1(jmax)
         nv1(jmax)=n
         n=nv2(i)
         nv2(i)=nv2(jmax)
         nv2(jmax)=n
         endif
!      write (6,'(''ord8 ou'',i6,3x,8I3)') i,nv1(1),nv1(2),nv1(3),nv1(4),
!     *                                      nv2(1),nv2(2),nv2(3),nv2(4)
      enddo
      i1=nv1(1)
      i2=nv1(2)
      i3=nv1(3)
      i4=nv1(4)
      i5=nv2(1)
      i6=nv2(2)
      i7=nv2(3)
      i8=nv2(4)
!      write (6,'(''ord8 ou'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8         
      end subroutine ord8
!------------------------------------------------------------------------------

end module ord_utils
