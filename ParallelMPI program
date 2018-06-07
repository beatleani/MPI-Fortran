program start
	
	!! MPI formalities
	
	use mpi
	implicit none
      integer ierr, pro_num, no_pro, ppt
      integer  MASTER
	  integer  ARRAYSIZE
		
      parameter (ARRAYSIZE = 32)
      parameter (MASTER = 0)
      integer i, mi,  source, chunksize , start1, end1, stop1,start2, stop2, j , k , l 
      real, allocatable, dimension (:,:,:) :: Fx, Fy, Fz, temp_div, P, Q, R, temp_curl,temp_div2,temp_curl2
	  real, allocatable, dimension (:,:,:) :: div,x,y,z,norm,xv,yv,zv,curl
	  !real,dimension(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE)::   div,x,y,z,norm,xv,yv,zv,curl
      real  dx,dy,dz,p1,p2,q1,q2,r1,r2
	  double precision t1,t2,t
      integer  status(MPI_STATUS_SIZE)

	!call cpu_time(t1)
	
    ! use mpi is used for invoking the MPI subroutines
    ! ierr is an error indicator in MPI processes.
    ! ierr is to be zero. If non-zero, then there is an error.
    ! pro_num is the processor number
    ! no_pro is the number of processors

	call MPI_INIT(ierr)
    ! mpi_init starts the MPI processes
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, no_pro, ierr)
    ! mpi_comm_size finds the number of processes under use 9176784438
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, pro_num, ierr)
	 allocate(x(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), y(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(z(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
     allocate(norm(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), xv(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(yv(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), zv(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))
	allocate(div(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE), curl(1:ARRAYSIZE,1:ARRAYSIZE,1:ARRAYSIZE))

		open( unit= 21,file='32data.txt')
		
        do i = 1,ARRAYSIZE
			do j = 1,ARRAYSIZE
				do k = 1,ARRAYSIZE
	  
	  read(21,*) x(i,j,k),y(i,j,k),z(i,j,k),norm(i,j,k),xv(i,j,k),yv(i,j,k),zv(i,j,k)
				end do
			end do
		end do
		print *,' Finished reading data, MPI initializing started'
  

 		ppt = (ARRAYSIZE + no_pro-2)/ ( no_pro-1) 
		

	!***** Master task only ******
      if (pro_num .eq. MASTER) then
		
		print *,'From master processor, no_pro =',no_pro
		
	! Receive Divergence

	 ! call MPI_RECV(start, count, datatype, source, tag, comm, status, ierr)

	 
	 open(unit= 1, file='Results.txt', form ='formatted')
		!div = 0.0
		!curl = 0.0
	t1 = MPI_WTIME()
	do i = 1, no_pro-1
		!do j = 1, no_pro-1
			!do k = 1, no_pro-1
	 
		start1 = (i - 1)*ppt + 1
	end1 = min(arraysize, i*ppt)
	chunksize = end1 - start1 + 1
	 allocate(temp_div2(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE),temp_curl2(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE)) 
		print *, 'From master processor for slave' ,i, 'start , end = ' , start1 , end1
		call MPI_RECV(temp_div2(1,1,1), chunksize*ARRAYSIZE*ARRAYSIZE, MPI_REAL, i, i, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(temp_curl2(1,1,1), chunksize*ARRAYSIZE*ARRAYSIZE, MPI_REAL, i, i+no_pro, MPI_COMM_WORLD, status, ierr)
		
	print *, 'Master processor received data from slave worker ',i
		div(start1:end1,1:ARRAYSIZE,1:ARRAYSIZE) = temp_div2(start1:end1,1:ARRAYSIZE,1:ARRAYSIZE)
		curl(start1:end1,1:ARRAYSIZE,1:ARRAYSIZE) = temp_curl2(start1:end1,1:ARRAYSIZE,1:ARRAYSIZE)
		
	
	deallocate(temp_div2,temp_curl2) 
	t2 = MPI_WTIME()
	end do
		 do l = 1, ARRAYSIZE
		do j = 1,ARRAYSIZE
			do k = 1, ARRAYSIZE
	write(1,1) ' The divergences of the vector fields at ',l,j,k ,'are,', div(l,j,k)
	write(1,1) ' The curls of the vector fields at ',l,j,k,'are,', curl(l,j,k)
	1 format (a,i2,i2,i2,a,e12.4)
			end do
		end do
	end do
	!deallocate(div,curl) 
	t = t2 - t1
	print *,'Time taken = ', t
	end if
	
	!  Non master process
	
	if (pro_num .gt. MASTER) then
	start1 = (pro_num - 1)*ppt + 1
	end1 = min(arraysize, pro_num*ppt)
	chunksize = end1 - start1 + 1
	print *, ' For slave' ,pro_num, 'start , end , chunksize = ' , start1 , end1 , chunksize
   ! calculate divergence
	
	
	 allocate(Fx(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE), Fy(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(Fz(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE))
     allocate(temp_div(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE), P(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(Q(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE), R(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE))
	 allocate(temp_curl(1:chunksize,1:ARRAYSIZE,1:ARRAYSIZE))
	 
	
!	if(pro_num == 1 ) then
!		temp_div(1,1,1) = 0.0
!		temp_curl(1,1,1) = 0.0
!		start2 = 2
!		stop2 = ARRAYSIZE	
!	else if(pro_num == no_pro - 1) then
!		temp_div(chunksize,ARRAYSIZE,ARRAYSIZE) = 0.0
!		temp_curl(chunksize,ARRAYSIZE,ARRAYSIZE) = 0.0
!		start2 = 2
!		stop2 = ARRAYSIZE - 1
!	else
!		start2 = 1
!		stop2 = ARRAYSIZE
!	end if
	print *, 'Process number: ', pro_num ,'started calculation'

	do i = start1,end1
		do j = 1,ARRAYSIZE
			do k = 2,ARRAYSIZE-1
		!mi = start1 + i
		! print *, 'Process number : ', pro_num , ' Y(mi) = ', y(mi)
		! print *, 'Process number : ', pro_num , ' YV(mi) = ', yv(mi)
		dx = x(i,j,k) - x(i,j,k-1)
		dy = y(i,j,k) - y(i,j,k-1)
		dz = z(i,j,k) - z(i,j,k-1)

		! print *, 'Process number : ', pro_num , ' dy(mi) = ', dy

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
		temp_div(i,j,k) = Fx(i,j,k)+Fy(i,j,k)+Fz(i,j,k)
		temp_curl(i,j,k) = P(i,j,k) - Q(i,j,k) + R(i,j,k)
				end do
			end do
		end do


	! print *, 'Process number:' , pro_num ,' divergence = ' , temp_div

	! Send Temp_Divergence to master process
	! call MPI_SEND(start, count, datatype, dest, tag, comm, ierr)
	call MPI_SEND(temp_div(1,1,1), chunksize*ARRAYSIZE*ARRAYSIZE, MPI_REAL, 0,pro_num, MPI_COMM_WORLD, ierr)
	call MPI_SEND(temp_curl(1,1,1), chunksize*ARRAYSIZE*ARRAYSIZE, MPI_REAL, 0, pro_num +no_pro, MPI_COMM_WORLD, ierr)
 	print *,' Process number:' , pro_num ,'finished sending data'


	deallocate(x,y,z,xv,yv,zv,div,curl,Fx, Fy, Fz, temp_div, P, Q, R, temp_curl)

	end if
 	 !call cpu_time(t2)

		
 call MPI_FINALIZE(ierr)

end program start

