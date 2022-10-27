program rot

implicit none

!Variables
integer :: i,j,k,ii
integer :: Natoms !Number of atoms
integer :: IERR
character (len=2), allocatable :: atoms(:)
double precision, parameter :: degtorad=dacos(-1.0d0)/180.0d0
double precision, parameter :: bohrA = 0.5291772083d0
double precision, parameter :: autocm=219474.63d0
double precision, allocatable :: X0(:,:)   !Cartesian coord.
double precision, allocatable :: xmass(:)
double precision :: Xcom(3), Tmx(3,3), Evec(3), se(3), masstot
double precision :: Imx(3,3), Av, Bv, Cv


open(222,file='CH2FCl.xyz')

read(222,*) Natoms

!Allocate
allocate(atoms(Natoms),X0(Natoms,3),xmass(Natoms))

read(222,*)

!Read atoms and coordinates
do i=1,Natoms
    read(222,*) atoms(i), (X0(i,j),j=1,3)
end do

!Set atomic masses
do i=1,Natoms
  call massselect(atoms(i),xmass(i))
end do

!Calculate total mass
masstot=0.0d0
do i=1,Natoms
  masstot=masstot+xmass(i)
end do

!Center of mass coordinates
do j=1,3
   Xcom(j)=0.0d0
   do i=1,Natoms
      Xcom(j)=Xcom(j)+xmass(i)*X0(i,j)
   end do
   Xcom(j)=Xcom(j)/masstot
end do

!Translate to COM
do j=1,3
   do i=1,Natoms
     X0(i,j)=X0(i,j)-Xcom(j)
   end do
end do

!Compute moment of inertia tensor
do i=1,3
 do j=1,3
  Imx(i,j)=0.0d0
 end do
end do

do i=1,Natoms
 Imx(1,1)=Imx(1,1)+xmass(i)*(X0(i,2)**2+X0(i,3)**2)
 Imx(2,2)=Imx(2,2)+xmass(i)*(X0(i,1)**2+X0(i,3)**2)
 Imx(3,3)=Imx(3,3)+xmass(i)*(X0(i,1)**2+X0(i,2)**2)
 Imx(1,2)=Imx(1,2)-xmass(i)*X0(i,1)*X0(i,2)
 Imx(1,3)=Imx(1,3)-xmass(i)*X0(i,1)*X0(i,3)
 Imx(2,3)=Imx(2,3)-xmass(i)*X0(i,2)*X0(i,3)
end do
Imx(2,1)=Imx(1,2)
Imx(3,1)=Imx(1,3)
Imx(3,2)=Imx(2,3)

!Diagonalize
call SUBTRED2(3,3,Imx,se,Evec,Tmx)
call TQL2(3,3,se,Evec,Tmx,IERR)

!Rotational constants
Av=1.0d0/(2.0d0*se(1))*(bohrA**2)*autocm
Bv=1.0d0/(2.0d0*se(2))*(bohrA**2)*autocm
Cv=1.0d0/(2.0d0*se(3))*(bohrA**2)*autocm

!Results
write(*,*) "Rotational constants: "
write(*,111) "A: ", Av
write(*,111) "B: ", Bv
write(*,111) "C: ", Cv

111 format(A5,F5.2)

close(222)

deallocate(atoms,X0,xmass)

end program

!Atomic mass subroutine
subroutine massselect(atom,mass)
character (len=2), intent(in) :: atom
double precision, intent(out) :: mass

select case(atom)
  case("H") 
     mass=1837.150d0
  case("C")
     mass=21874.66d0
  case("F")
     mass=34631.96d0
  case("Cl")
     mass=35.45270d0*1822.8880d0
end select
end subroutine
