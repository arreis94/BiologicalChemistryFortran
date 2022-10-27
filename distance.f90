      double precision function distance(a,b)
      implicit real*8 (a-h, o-z)
      implicit integer*4 (i-n)
      dimension a(3),b(3),c(3)
      do i=1,3
          c(i) = a(i)-b(i)
      end do
      distance = vnorm(c)
      return
      end
