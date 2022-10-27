************************************************************************
      double precision function harary(natoms,mdist)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      parameter (mat=100)
      real*8 mdist
      dimension mdist(mat,mat)
c
c H index (Harary)
c
      harary=0.d0
      do i=1,natoms
         do j=1,natoms
            if(i.ne.j) harary=harary+1.d0/mdist(i,j)
         end do
      end do
      harary=0.5d0*harary
      return
      end

