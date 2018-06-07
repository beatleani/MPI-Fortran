program divcurl
	implicit none
	integer :: i,j,k,ARRAYSIZE
	parameter ( ARRAYSIZE = 512 )
	real,allocatable ,dimension(:,:,:) :: x,y,z,norm,xv,yv,zv,div,curl,Fx,Fy,Fz,P,Q,R
	real  dx,dy,dz,p1,q1,r1,p2,q2,r2
	real ( kind = 8 ) wtime
  	real ( kind = 8 ) wtime1
  	real ( kind = 8 ) wtime2

	call cpu_time ( wtime1 )
	open(unit=1,file='512data.txt')
	 allocate(x(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), y(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(z(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
     allocate(norm(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), xv(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(yv(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), zv(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	!$ parallel do
	do i = 1,ARRAYSIZE
		do j = 1,ARRAYSIZE
			do k = 1,ARRAYSIZE
	  
	  read(1,*) x(i,j,k),y(i,j,k),z(i,j,k),norm(i,j,k),xv(i,j,k),yv(i,j,k),zv(i,j,k)
			end do
		end do
	end do

! Calculating divergence of the vector fields
	 allocate(Fx(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), Fy(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(Fz(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
     allocate(div(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), P(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(Q(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), R(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(curl(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	do i = 1,ARRAYSIZE
		do j = 1,ARRAYSIZE
			do k = 2,ARRAYSIZE-1

				dx = x(i,j,k) - x(i,j,k-1)
				dy = y(i,j,k) - y(i,j,k-1)
				dz = z(i,j,k) - z(i,j,k-1)
			if ( dx .eq. 0 ) then
			Fx(i,j,k) = 0.0 
			q1 = 0.0
			r1 = 0.0
		else 
		Fx(i,j,k) = (xv(i,j,k+1) - xv(i,j,k-1))/((2.0)*dx)
		q1 = (zv(i,j,k+1) - zv(i,j,k-1))/((2.0)*dx)
		r1 = (yv(i,j,k+1) - yv(i,j,k-1))/((2.0)*dx)
		end if

		if ( dy .eq. 0 ) then
			Fy(i,j,k) = 0.0 
			p1 = 0.0
			r2 = 0.0
		else 
		Fy(i,j,k) = (yv(i,j,k+1) - yv(i,j,k-1))/((2.0)*dy)
		p1 = (zv(i,j,k+1) - zv(i,j,k-1))/((2.0)*dy)
		r2 = (xv(i,j,k+1) - xv(i,j,k-1))/((2.0)*dy)
		
		end if

		if ( dz .eq. 0 ) then
			Fz(i,j,k) = 0.0 
			p2 = 0.0
			q2 = 0.0
		else 
		Fz(i,j,k) = (zv(i,j,k+1) - zv(i,j,k-1))/((2.0)*dz)
		 p2 = (yv(i,j,k+1) - yv(i,j,k-1))/((2.0)*dz)
		 q2 = (xv(i,j,k+1) - xv(i,j,k-1))/((2.0)*dz)
		end if
		
		P(i,j,k) = p1 - p2
		Q(i,j,k) = q1 - q2
		R(i,j,k) = r1 - r2
		div(i,j,k) = Fx(i,j,k)+Fy(i,j,k)+Fz(i,j,k)
		curl(i,j,k) = P(i,j,k) - Q(i,j,k) + R(i,j,k)
			end do
		end do
	end do

!				dx = x(i,j,k) - x(i,j-1,k)
!				dy = y(i,j,k) - y(i,j-1,k)
!				dz = z(i,j,k) - z(i,j-1,k)
!					if ( dx .eq. 0 ) then
!			Fx(i,j,k) = 0.0 
!			q1(i,j,k) = 0.0
!			r1(i,j,k) = 0.0
!		else 
!		Fx(i,j,k) = (xv(i,j+1,k) - xv(i,j-1,k))/((2.0)*dx)
!		q1(i,j,k) = (zv(i,j+1,k) - zv(i,j-1,k))/((2.0)*dx)
!		r1(i,j,k) = (yv(i,j+1,k) - yv(i,j-1,k))/((2.0)*dx)
!		end if

!		if ( dy .eq. 0 ) then
!			Fy(i,j,k) = 0.0 
!			p1(i,j,k) = 0.0
!			r2(i,j,k) = 0.0
!		else 
!		Fy(i,j,k) = (yv(i,j+1,k) - yv(i,j-1,k))/((2.0)*dy)
!		p1(i,j,k) = (zv(i,j+1,k) - zv(i,j-1,k))/((2.0)*dy)
!		r2(i,j,k) = (xv(i,j+1,k) - xv(i,j-1,k))/((2.0)*dy)
!		
!		end if

!		if ( dz .eq. 0 ) then
!			Fz(i,j,k) = 0.0 
!			p2(i,j,k) = 0.0
!			q2(i,j,k) = 0.0
!		else 
!		Fz(i,j,k) = (zv(i,j+1,k) - zv(i,j-1,k))/((2.0)*dz)
!		 p2(i,j,k) = (yv(i,j+1,k) - yv(i,j-1,k))/((2.0)*dz)
!		 q2(i,j,k) = (xv(i,j+1,k) - xv(i,j-1,k))/((2.0)*dz)
!		end if
!		
!		P(i,j,k) = p1(i,j,k) - p2(i,j,k)
!		Q(i,j,k) = q1(i,j,k) - q2(i,j,k)
!		R(i,j,k) = r1(i,j,k) - r2(i,j,k)
!		div(i,j,k) = Fx(i,j,k)+Fy(i,j,k)+Fz(i,j,k)
!		curl(i,j,k) = P(i,j,k) - Q(i,j,k) + R(i,j,k)
!		end do
!				dx = x(i,j,k) - x(i-1,j,k)
!				dy = y(i,j,k) - y(i-1,j,k)
!				dz = z(i,j,k) - z(i-1,j,k)
!					if ( dx .eq. 0 ) then
!			Fx(i,j,k) = 0.0 
!			q1(i,j,k) = 0.0
!			r1(i,j,k) = 0.0
!		else 
!		Fx(i,j,k) = (xv(i+1,j,k) - xv(i-1,j,k))/((2.0)*dx)
!		q1(i,j,k) = (zv(i+1,j,k) - zv(i-1,j,k))/((2.0)*dx)
!		r1(i,j,k) = (yv(i+1,j,k) - yv(i-1,j,k))/((2.0)*dx)
!		end if

!		if ( dy .eq. 0 ) then
!			Fy(i,j,k) = 0.0 
!			p1(i,j,k) = 0.0
!			r2(i,j,k) = 0.0
!		else 
!		Fy(i,j,k) = (yv(i+1,j,k) - yv(i-1,j,k))/((2.0)*dy)
!		p1(i,j,k) = (zv(i+1,j,k) - zv(i-1,j,k))/((2.0)*dy)
!		r2(i,j,k) = (xv(i+1,j,k) - xv(i-1,j,k))/((2.0)*dy)
!		
!		end if

!		if ( dz .eq. 0 ) then
!			Fz(i,j,k) = 0.0 
!			p2(i,j,k) = 0.0
!			q2(i,j,k) = 0.0
!		else 
!		Fz(i,j,k) = (zv(i+1,j,k) - zv(i-1,j,k))/((2.0)*dz)
!		 p2(i,j,k) = (yv(i+1,j,k) - yv(i-1,j,k))/((2.0)*dz)
!		 q2(i,j,k) = (xv(i+1,j,k) - xv(i-1,j,k))/((2.0)*dz)
!		end if
!		
!		P(i,j,k) = p1(i,j,k) - p2(i,j,k)
!		Q(i,j,k) = q1(i,j,k) - q2(i,j,k)
!		R(i,j,k) = r1(i,j,k) - r2(i,j,k)
!		div(i,j,k) = Fx(i,j,k)+Fy(i,j,k)+Fz(i,j,k)
!		curl(i,j,k) = P(i,j,k) - Q(i,j,k) + R(i,j,k)
!	end do
		deallocate(x, y, z, norm, xv, yv, zv, Fx, Fy, Fz, P, Q, R)
		
			
	
	!$ end parallel do

	call cpu_time ( wtime2 )

	wtime = wtime2 - wtime1

	write (*,*) 'Time = ', wtime

	 open(unit= 2, file='Resultssequential.txt',form='formatted')

	do i = 1,ARRAYSIZE
		do j = 1,ARRAYSIZE
			do k = 1,ARRAYSIZE
	
	write(2,2) ' The divergences of the vector fields at ', i,j,k ,'are'  ,  div(i,j,k)
	write(2,2) ' The curls of the vector fields at', i,j,k ,'are'  , curl(i,j,k)
	  
			end do
		end do
	end do
	2 format (a,i2,i2,i2,a,e12.4)
deallocate(div,curl)
	 close(2)
end program divcurl













	  

