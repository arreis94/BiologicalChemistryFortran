************************************************************************
      subroutine sdist(natoms,mdist,isd)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      parameter (mat=100)
      real*8 mdist,isd
      dimension mdist(mat,mat),isd(mat)
c isd(): distance-sums vector
      do i=1,natoms
         isd(i)=0.d0
         do j=1,natoms
            isd(i)=isd(i)+mdist(i,j)
         end do
      end do
      return
      end
