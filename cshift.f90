      subroutine cshift(natoms, coord, shift)
      implicit real*8 (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (mat=100)
      dimension coord(3,mat), shift(3)
      do i=1,natoms
          do j=1,3
              coord(j,i) = coord(j,i)-shift(j)
          end do
      end do
      return
      end
