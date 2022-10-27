      subroutine cmass(natoms, nat, wt, coord, wmol, com)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (mat=100)
      dimension wt(90), coord(3,mat), nat(mat), com(3)
      sumwx = 0.0
      sumwy = 0.0
      sumwz = 0.0
      wmol = 0.0
      do i=1,natoms
          nati = nat(i)
          wmol = wmol + wt(nati)
          sumwx = sumwx + wt(nati)*coord(1,i)
          sumwy = sumwy + wt(nati)*coord(2,i)
          sumwz = sumwz + wt(nati)*coord(3,i)
      end do
      com(1) = sumwx/wmol
      com(2) = sumwy/wmol
      com(3) = sumwz/wmol
      return
      end
