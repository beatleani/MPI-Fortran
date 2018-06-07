# MPI-Fortran
Given a 512 x 512 x 512 grid of data points, it naturally takes a lot of time to run it sequentially. Here, I've devised a code in Fortran, written in MPI to parallelize the sequential program in order to reduce the time taken for computation.
The computation done here is to calculate the divergence and curl of the vector fields in the text file, using central difference method. The file contains 7 columns whose description is provided below.
Columns(1,2,3) are the x,y and z coordinates 
Column 4 is a scalar(norm of the vector in the net three columns) 
Columns(5,6,7) are the x,y and z components of a three dimensional vector
A plot is also submitted where the time taken is on the Y-axis and the number of processors is on the X-axis. 
