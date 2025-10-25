      
      subroutine check_mesh(g)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g
      real :: area_min, dx_error, dy_error, tol, flg
      integer :: ni, nj, i, j 
      real, allocatable :: xvec(:,:), yvec(:,:), perimeter(:,:) 

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Exact checking of floating point numbers never goes well, define a
!     small tolerance to use for a comparative operation instead
      tol = 1e-4 * g%l_min

!     Check that all of the cell areas are positive, either with the intrinsic
!     "minval" function or with nested do loops. Print the output to the screen
!     and flag negative numbers as an error with an if statement to "stop" the
!     program
!     INSERT

      area_min = minval(g%area) 
      write(6,*) "Minimum cell area =", area_min

      flg = 0

      do j = 1, nj-1
            do i = 1, ni-1
                  if (g%area(i,j) < 0.0d0) then
                  write(6,*) "Error, negative area, cell =", g%area
                  flg = 1 
                  end if
            end do
      end do

      write(6,*) "No negative cell areas found"

      if (flg == 1) then  
            write (6,*) "Stopping program"
            stop
      end if
                        

!     Next check that the sum of the edge vectors around every quadrilateral is 
!     very nearly zero in both the x and y-coordinate directions. You can
!     complete this with some elementwise addition of the arrays and use of the
!     "maxval" and "abs" intrinsic functions.
!     INSERT

      allocate(xvec(ni-1,nj-1))
      allocate(yvec(ni-1,nj-1))
      allocate(perimeter(ni-1,nj-1))   

      do j=1, nj-1
            do i =1, ni-1
                  yvec(i,j) = (g%x(i,j+1) - g%x(i,j)) + (g%x(i+1,j) - g%x(i+1,j+1)) &
                            + (g%y(i,j+1) - g%y(i,j)) + (g%y(i+1,j) - g%y(i+1,j+1))

                  xvec(i,j) = (g%x(i+1,j+1)-g%x(i+1,j)) + (g%x(i,j)-g%x(i+1,j)) &
                            + (g%y(i+1,j+1)-g%y(i+1,j)) + (g%y(i,j)-g%y(i+1,j))

                  perimeter(i,j) = yvec(i,j) + xvec(i,j)

                  ! write(6,*) "perimeter =", perimeter(i,j)

            end do
      end do

      write(6,*) "perimeter =", perimeter(1,1)
      write(6,*) "perimeter =", perimeter(5,5)

!     It may be worthwhile to complete some other checks, the prevous call to
!     the "write_output" subroutine has written a file that you can read and
!     postprocess using the Python script plot_mesh.py. This program also has
!     access to all of the mesh parameters used within the solver that you could
!     inspect graphically.

!     Print a blank line
      write(6,*)

      end subroutine check_mesh
