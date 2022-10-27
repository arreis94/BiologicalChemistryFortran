************************************************************************
      double precision function randic(nbonds,ivd,ia,ib)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      real*8 ivd
      dimension ivd(*),ia(*),ib(*)
c
c  1khi connectivity index (Randic,1975)
c
      randic=0.d0
      do i=1,nbonds
         randic=randic+dsqrt(1.d0/(ivd(ia(i))*ivd(ib(i))))
      end do
      return
      end
