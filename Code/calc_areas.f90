      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g 
      real(8), allocatable :: li(:,:), lj(:,:)
            
      integer :: ni, nj
      
      


!     Declare integers or any extra variables you need here
!     INSERT
      type :: vector2d
            real :: x
            real :: y
      end type vector2d 

      integer :: i, j
      real :: min_i, min_j

      type(vector2d) :: a, b

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
!     INSERT

      do j = 1, nj-1
            do i = 1,ni-1
                  ! a%x = ((g%x(i+1,j) - g%x(i,j))^2 + (g%y(i+1,j) - g%y(i,j))^2)^0.5
                  ! a%y = ((g%x(i+1,j+1) - g%x(i+1,j))^2 + (g%y(i+1,j+1) - g%y(i+1,j))^2)^0.5
                  ! a%x = g%x(i+1,j) - g%x(i,j)
                  ! a%y = g%y(i,j+1) - g%y(i,j)
                  ! b%x = g%x(i,j) - g%x(i+1,j)
                  ! b%y = g%y(i,j+1) - g%y(i,j)
                  ! g%area(i,j) = 0.5 * abs(a%x * b%y - a%y * b%x)

                  ! a = diagonal (i,j) -> (i+1,j+1)
                  a%x = g%x(i+1,j+1) - g%x(i  ,j  )
                  a%y = g%y(i+1,j+1) - g%y(i  ,j  )

                  ! b = other diagonal (i,j+1) -> (i+1,j)
                  b%x = g%x(i  ,j+1) - g%x(i+1,j)
                  b%y = g%y(i  ,j+1) - g%y(i+1,j)

                  g%area(i,j) = 0.5d0 * abs( a%x*b%y - a%y*b%x )

            end do
      end do

      ! write(6,*) 'Calculated cell area', g%area(1,1)


!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
!     INSERT

      do j = 1, nj-1
            do i = 1, ni
                  g%lx_i(i,j) = g%y(i,j+1) - g%y(i,j)
                  g%ly_i(i,j) = g%x(i,j) - g%x(i,j+1)
                  
            end do
      end do      
      
      ! write(6,*) "lx_i vlas", g%lx_i

!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
!     INSERT
      do i = 1, ni-1
            do j = 1, nj
                  g%lx_j(i,j) = g%y(i+1,j) - g%y(i,j)
                  g%ly_j(i,j) = g%x(i+1,j) - g%x(i,j)
            end do
      end do      


!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.
!     INSERT

      allocate(li(ni-1, nj-1))
      allocate(lj(ni-1, nj-1))

      do j = 1, nj-1
            do i = 1, ni-1
                  li(i,j) = hypot(g%lx_i(i,j),g%ly_i(i,j))
                  lj(i,j) = hypot(g%lx_j(i,j),g%ly_j(i,j))
                  
            end do
      end do

      min_i = minval(li)
      min_j = minval(lj)

      g%l_min = min(min_i, min_j)
!
!     Print the overall minimum length size that has been calculated
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', g%l_min
      write(6,*)

      end subroutine calc_areas
