************************************************************************
      subroutine getmdist(a,n,maxat)
************************************************************************
c
c Floyd-Warshall algorithm for generating the shortest distance matrix
c
      implicit none
      integer*4 n,maxat,i,j,k
      real*8 a
      dimension a(maxat,maxat)

c a(,): distance matrix      
      do k=1,n
         do i=1,n
            do j=1,n
               if(a(i,j).gt.(a(i,k)+a(k,j))) a(i,j)=a(i,k)+a(k,j)
            end do
         end do
      end do
      return
      end

