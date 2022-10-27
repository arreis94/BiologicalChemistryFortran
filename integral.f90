PROGRAM integral
	dimension x(3),y(3)
 	write(*,*)'Integral bounds:'
 	read(*,*)a,b
 	n=6

 	h=(b-a)/n
 	x(1)=a
 	y(1)=f(x(1))

	do i=2,n+1
	 x(i)=x(i-1)+h
	 y(i)=f(x(i))
	end do
	s=y(1)+y(n+1)
	
	write(*,*)'0 = Simpson's rule, 1 = Tapezoidal rule'
	read(*,*)modsz
	if (modsz.eq.0) then
		do i=2,n,2
			s=s+4*y(i)
		end do
		do i=3,n-1,2
			s=s+2*y(i)
		end do
		s=h*s/3
		write(*,*)s
	else
		do i=2,n
			s=s+2*y(i)
		end do
		v=h*s/2
		write(*,*)v
	endif
end

function f(x)
	f=1/(1+x*x)
 	return
end
