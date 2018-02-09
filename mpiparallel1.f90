program start
	
	!! MPI formalities
	
	use mpi
	implicit none
      integer ierr, pro_num, no_pro, ppt
      integer   ARRAYSIZE, MASTER
      parameter (ARRAYSIZE = 32768)
      parameter (MASTER = 0)

      integer i, mi,  source, chunksize , start1, end1, stop1,start2, stop2
      real*8 , allocatable, dimension (:) :: Fx, Fy, Fz, temp_div, P, Q, R, temp_curl,temp_div2,temp_curl2
      real*8 , dimension(1:ARRAYSIZE) ::   div,x,y,z,norm,xv,yv,zv,curl
      real*8  dx,dy,dz,p1,p2,q1,q2,r1,r2
      integer  status(MPI_STATUS_SIZE)
	
    ! use mpi is used for invoking the MPI subroutines
    ! ierr is an error indicator in MPI processes.
    ! ierr is to be zero. If non-zero, then there is an error.
    ! pro_num is the processor number
    ! no_pro is the number of processors

		open( unit= 21,file='32data.txt')
        do i=1, ARRAYSIZE
          read(21,*) x(i) ,y(i) ,z(i) ,norm(i), xv(i), yv(i), zv(i)
        end do

		print *,' Finished reading data, MPI initializing started'
    call MPI_INIT(ierr)
    ! mpi_init starts the MPI processes
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, no_pro, ierr)
    ! mpi_comm_size finds the number of processes under use
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, pro_num, ierr)

 		ppt = (ARRAYSIZE + no_pro-2)/ (no_pro-1)
     

	!***** Master task only ******
      if (pro_num .eq. MASTER) then
		
		print *,'From master processor, no_pro =',no_pro
		
	! Receive Divergence

	 ! call MPI_RECV(start, count, datatype, source, tag, comm, status, ierr)

	 
	 open(unit= 1, file='Results.txt', form ='formatted')
		div = 0.0
		curl = 0.0
	do i = 1, no_pro-1
	  
		start1 = (i - 1)*ppt + 1
	end1 = min(arraysize, i*ppt)
	chunksize = end1 - start1 + 1
	allocate(temp_div2(1:chunksize),temp_curl2(1:chunksize)) 
		print *, 'From master processor for slave' ,i, 'start , end = ' , start1 , end1
		call MPI_RECV(temp_div2(1), chunksize, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
		call MPI_RECV(temp_curl2(1), chunksize, MPI_DOUBLE_PRECISION, i, i+no_pro, MPI_COMM_WORLD, status, ierr)
	print *, 'Master processor received data from slave worker ', i
		div(start1:end1) = temp_div2
		curl(start1:end1) = temp_curl2
	
	write(1,*) ' The divergences of the vector fields calculated by the',i,'th processor are', div(start1:end1)
	write(1,*) ' The curls of the vector fields calculated by the',i,'th processor are', curl(start1:end1)
	deallocate(temp_div2,temp_curl2) 
	end do
		 
	end if
	
	!  Non master process

	if (pro_num .gt. MASTER) then
	start1 = (pro_num - 1)*ppt + 1
	end1 = min(arraysize, pro_num*ppt)
	chunksize = end1 - start1 + 1
	print *, ' For slave' ,pro_num, 'start , end , chunksize = ' , start1 , end1 , chunksize
   ! calculate divergence
	
	 allocate(Fx(1:chunksize), Fy(1:chunksize), Fz(1:chunksize))
     allocate(temp_div(1:chunksize), P(1:chunksize), Q(1:chunksize), R(1:chunksize))
	 allocate(temp_curl(1:chunksize),temp_div2(1:chunksize),temp_curl2(1:chunksize))   
   
	
	if(pro_num == 1 ) then
		temp_div(1) = 0.0
		temp_curl(1) = 0.0
		start2 = 2
		stop2 = chunksize	
	else if(pro_num == no_pro -  1) then
		temp_div(chunksize) = 0.0
		temp_curl(chunksize) = 0.0
		start2 = 1
		stop2 = chunksize - 1
	else
		start2 = 1
		stop2 = chunksize
	end if
	print *, 'Process number: ', pro_num ,'started calculation'
	do i = start2, stop2
		mi = start1 + i
		! print *, 'Process number : ', pro_num , ' Y(mi) = ', y(mi)
		! print *, 'Process number : ', pro_num , ' YV(mi) = ', yv(mi)
		dx = x(mi) - x(mi-1)
		dy = y(mi) - y(mi-1)
		dz = z(mi) - z(mi-1)

		! print *, 'Process number : ', pro_num , ' dy(mi) = ', dy

		if ( dx .eq. 0 ) then
			Fx(i) = 0.0 
			q1 = 0.0
			r1 = 0.0
		else 
		Fx(i) = (xv(mi+1) - xv(mi-1))/((2.0)*dx)
		q1 = (zv(mi+1) - zv(mi-1))/((2.0)*dx)
		r1 = (yv(mi+1) - yv(mi-1))/((2.0)*dx)
		end if
		if ( dy .eq. 0 ) then
			Fy(i) = 0.0 
			p1 = 0.0
			r2 = 0.0
		else 
		Fy(i) = (yv(mi+1) - yv(mi-1))/((2.0)*dy)
		p1 = (zv(mi+1) - zv(mi-1))/((2.0)*dy)
		r2 = (xv(mi+1) - xv(mi-1))/((2.0)*dy)
		
		end if

		if ( dz .eq. 0 ) then
			Fz(i) = 0.0 
			p2 = 0.0
			q2 = 0.0
		else 
		Fz(i) = (zv(mi+1) - zv(mi-1))/((2.0)*dz)
		 p2 = (yv(mi+1) - yv(mi-1))/((2.0)*dz)
		 q2 = (xv(mi+1) - xv(mi-1))/((2.0)*dz)
		end if
		
		P(i) = p1 - p2
		Q(i) = q1 - q2
		R(i) = r1 - r2
		temp_div(i) = Fx(i)+Fy(i)+Fz(i)
		temp_curl(i) = P(i) - Q(i) + R(i)
	end do
	! print *, 'Process number:' , pro_num ,' divergence = ' , temp_div

	! Send Temp_Divergence to master process
	! call MPI_SEND(start, count, datatype, dest, tag, comm, ierr)
	call MPI_SEND(temp_div(1), chunksize, MPI_DOUBLE_PRECISION, 0,pro_num, MPI_COMM_WORLD, ierr)
	call MPI_SEND(temp_curl(1), chunksize, MPI_DOUBLE_PRECISION, 0, pro_num +no_pro, MPI_COMM_WORLD, ierr)
 	print *,' Process number:' , pro_num ,'finished sending data'
	deallocate(Fx, Fy, Fz, temp_div, P, Q, R, temp_curl,temp_div2,temp_curl2)
	end if
		
 call MPI_FINALIZE(ierr)
end program start

