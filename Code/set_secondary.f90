      
      subroutine set_secondary(av,g)

!     This subroutine calculates the secondary flow variables from the primary
!     ones at every node

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g

!     Define any further variables you may need
!     INSERT
      real(8), allocatable :: v(:,:), t_stat(:,:), hstag(:,:)

      

!     The primary flow variables are "ro", "roe", "rovx" and "rovy", these are 
!     the conserved quantities within the Euler equations. Write the code to
!     calculate the secondary flow variables which are the velocity components
!     "vx" and "vy", the static pressure "p" and the stagnation enthalpy
!     "hstag". These are needed at every timestep, there is no need for any 
!     loops as the operations can be performed elementwise, although you may
!     wish to define some intermediate variables to improve readability.
!     INSERT

      allocate(v(g%ni,g%nj))
      allocate(t_stat(g%ni,g%nj))
      allocate(hstag(g%ni,g%nj))

      g%vx = g%rovx/ g%ro
      g%vy = g%rovy/ g%ro
      v = hypot(g%vx, g%vy)

      ! g%p = bcs%pstag - (0.5 * g%ro * g%v **2)

      g%p =g%ro * av%rgas * t_stat
      t_stat=(g%roe/g%ro - (0.5 * v**2) ) / av%cv
      hstag = (av%cp * t_stat) + (1.5 * v**2)
      
      write(6,*) "p=", g%p(5,5) ! why is it printing 5 things.. 
      !write(6,*) "hstag=", hstag


      end subroutine set_secondary


