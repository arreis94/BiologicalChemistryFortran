************************************************************************
      double precision function balaban(natoms,nbonds,isd,ia,ib)
************************************************************************
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      real*8 isd
      dimension isd(*)
      dimension ia(*),ib(*)
c
c J index (Balaban,1982)
c
      balaban=0.d0
      do i=1,nbonds
         balaban=balaban+
     $           dsqrt(1.d0/(isd(ia(i))*isd(ib(i))))
      end do
      balaban=balaban*dble(nbonds)/(dble(nbonds-natoms+1)+1.d0)
      return
      end
