      subroutine crossp(a,b,c)
      implicit real*8 (a-h, o-z)
      implicit integer*4 (i-n)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
