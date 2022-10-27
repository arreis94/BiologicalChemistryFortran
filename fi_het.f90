************************************************************************
      double precision function fi(nhyd,nbonds,ivd,ia,ib)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      real*8 ivd
      dimension ivd(*),ia(*),ib(*)
c
c  Fi index ()
c
      fi=0.d0
      do i=1,nbonds
         fi=fi+dsqrt(dabs((ivd(ia(i))-1.d0)*(ivd(ib(i))-1.d0)))
      end do
      fi=fi+dble(nhyd)
      return
      end

