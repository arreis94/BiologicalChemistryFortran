************************************************************************
      double precision function wiener(natoms,mdist)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      real*8 mdist
      parameter (mat=100)
      dimension mdist(mat,mat)
c
c W index (Wiener,1947)
c
      wiener=0.d0
      do i=1,natoms
         do j=1,natoms
            wiener=wiener+dble(mdist(i,j))
         end do
      end do
      wiener=0.5d0*wiener
      return
      end

