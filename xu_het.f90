************************************************************************
      double precision function xu(natoms,ivd,isd)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      real*8 ivd,isd
      dimension ivd(*),isd(*)
c
c Xu index (Ren,1998)
c
      sum1=0.d0
      sum2=0.d0
      do i=1,natoms
         sum1=sum1+ivd(i)*isd(i)
         sum2=sum2+ivd(i)*isd(i)*isd(i)
      end do
      xu=dsqrt(dble(natoms))*dlog10(sum2/sum1)
      return
      end

