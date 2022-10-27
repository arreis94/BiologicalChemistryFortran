      block data
      implicit real*8 (a-h, o-z)
      implicit integer*4 (i-n)
      character symb*2
      common/atoms/ wt(90),vdw(90),covrad(90),symb(90)
c
c
      data wt/ 1.00783d0, 4*0.0d0, 12.00000d0, 14.00307d0, 15.99491d0,
     $ 82*0.0d0 /
      data vdw/ 1.185d0, 4*0.0d0, 1.75d0, 1.525d0, 1.4d0,82*0.0d0 /
      data covrad / 0.435d0, 4*0.d0, 0.655d0, 0.75d0, 0.73d0, 82*0.d0 /
      data symb / ' H',4*'  ',' C',' N',' O',82*'  '/
      end
