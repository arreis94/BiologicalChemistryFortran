      double precision function vnorm(a)
      implicit real*8 (a-h, o-z)
      implicit integer*4 (i-n)
      dimension a(3)
      vnorm=dsqrt(dot(a,a))
      return
      end


