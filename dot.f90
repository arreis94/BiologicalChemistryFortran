      double precision function dot(a,b)
      implicit real*8 (a-h, o-z)
      implicit integer*4 (i-n)
      dimension a(3),b(3)
      dot=0.d0
      do i=1,3,1
         dot=dot+a(i)*b(i)
      end do
      return
      end
      
