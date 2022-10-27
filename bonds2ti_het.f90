      program main
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character title*70,symb*2
      logical hetero
      real*8 mcon,ivd,mdist,isd
      parameter(mat=100)
      dimension mcon(mat,mat),mdist(mat,mat)
      dimension ivd(mat),isd(mat)
      dimension ia(mat),ib(mat),nat(mat),br(mat)
      write(*,'(/a)') '*** Molecular topological indices from bonds ***'
      write(*,'(/a$)') 'Number of heavy atoms: '
      read(*,*) natoms
      hetero=.false.
      nhetero=0
      write(*,'(a$)') 'Hetero atoms <yes=1>: '
      read(*,*) nhetero
      if(nhetero.eq.1) hetero=.true.
      if(hetero) then
         do i=1,natoms
            write(*,'(a,i3,a$)') 'Atomic number for atom #',i,': '
            read(*,*) nat(i)
         end do
      else
         do i=1,natoms
            nat(i)=6
         end do
      endif
      write(*,'(a$)') 'Number of H atoms: '
      read(*,*) nhyd
      write(*,'(a$)') 'Number of bonds: '         
      read(*,*) nbonds
      nspcb=0
      write(*,'(a$)') 'Special bonds <yes=1>: '
      read(*,*) nspcb                                   
      do i=1,nbonds
         write(*,'(a,i3,a$)') 'Bond #',i,' atoms [i,j]: '
         read(*,*) ia(i),ib(i)
         br(i)=1.d0
         if(nspcb.eq.1) then
            write(*,'(a$)') 
     $            'Single=1, Double=2, Triple=3, Aromatic=1.5: '
            read(*,*) br(i)
         endif
      end do
      do i=1,natoms
         do j=1,natoms
            mcon(i,j)=0.d0
            mdist(i,j)=dble(nbonds)
         end do
      end do
      do i=1,natoms
         mcon(i,i)=1.d0-6.d0/dble(nat(i))
         mdist(i,i)=0.d0
      end do
      do i=1,nbonds
         k=ia(i)
         l=ib(i)
         mcon(k,l)=36.d0/dble(nat(k))/dble(nat(l))/br(i)
         mcon(l,k)=mcon(k,l)
         mdist(k,l)=mcon(k,l)
         mdist(l,k)=mcon(k,l)
      end do
      write(*,'(/a)') 'Adjacency matrix'
      do i=1,natoms
         write(*,'(20f7.2)') (mcon(i,j),j=1,natoms)
      end do
      call vdegree(natoms,mcon,ivd)     
      write(*,'(/a)') 'Vertex-degree vector'
      write(*,'(20f7.2)') (ivd(i),i=1,natoms) 
      call getmdist(mdist,natoms,mat)
      do i=1,natoms
         mdist(i,i)=mcon(i,i)
      end do
      write(*,'(/a)') 'Distance matrix'
      do i=1,natoms
         write(*,'(20f7.2)') (mdist(i,j),j=1,natoms)
      end do
      call sdist(natoms,mdist,isd)
      write(*,'(/a)') 'Distance-sums vector'
      write(*,'(20f7.2)') (isd(i),i=1,natoms) 
      write(*,'(/a)') 'Topological indices'
      write(*,'(a,f16.4)') 'Wiener (W)  = ',wiener(natoms,mdist)
      write(*,'(a,f16.4)') 'Harary (H)  = ',harary(natoms,mdist)
      write(*,'(a,f16.4)')  'Randic      = ',randic(nbonds,ivd,ia,ib)
      write(*,'(a,f16.4)')  'Randic^1/2  = ',
     $ dsqrt(randic(nbonds,ivd,ia,ib))
      write(*,'(a,f16.4)')  'Balaban (J) = ',
     $ balaban(natoms,nbonds,isd,ia,ib)
      write(*,'(a,f16.4)')  'Xu          = ',xu(natoms,ivd,isd)
      write(*,'(a,f16.4)')  'Fi          = ',fi(nhyd,nbonds,ivd,ia,ib)
      write(*,'(a,f16.4)')  'Fi^1/2      = ',
     $ dsqrt(fi(nhyd,nbonds,ivd,ia,ib))
      write(*,'(a,f16.4)')  'Fix         = ',0.5d0*
     $(dsqrt(fi(nhyd,nbonds,ivd,ia,ib))+dsqrt(randic(nbonds,ivd,ia,ib)))
      end 

