************************************************************************
      subroutine vdegree(natoms,mcon,ivd)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      parameter (mat=100)
      real*8 ivd,mcon
      dimension mcon(mat,mat),ivd(mat)
c ivd(): vertex-degree vector
      do i=1,natoms
         ivd(i)=0.d0
         do j=1,natoms
            ivd(i)=ivd(i)+mcon(i,j)
         end do
      end do
      return
      end
      
