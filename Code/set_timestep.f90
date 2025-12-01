      
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
      ! real :: astag, v_max
      real, allocatable :: astag(:,:), v(:,:), astag_avrg(:,:), v_avrg(:,:)
      integer :: ni, nj

!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
!     INSERT
      
      ! astag = sqrt(av%gam * av%rgas * bcs%tstag)

      !!!!!!!!!! calculating local speed of sound and local velocity!!!!!!!!!!!

      ni = g%ni; nj = g%nj;

      allocate(astag(1:ni, 1:nj),v(1:ni, 1:nj), astag_avrg(1:ni-1,1:nj-1), v_avrg(1:ni-1,1:nj-1))

      astag=0
      v=0
      astag_avrg=0
      v_avrg=0


      astag(1:ni, 1:nj) = sqrt(av%gam * g%p(1:ni, 1:nj) / g%ro(1:ni, 1:nj))
      v(1:ni, 1:nj) = sqrt(g%vx(1:ni, 1:nj)**2 + g%vy(1:ni, 1:nj)**2)

      !write(6,*) "g%p", g%p(1:10, 1:10)

      astag_avrg(1:ni-1,1:nj-1) = 0.25 *((astag(1:ni-1,1:nj-1) + astag(2:ni,1:nj-1) + astag(1:ni-1,2:nj) &
      + astag(2:ni,2:nj)))


      v_avrg(1:ni-1,1:nj-1) = 0.25 *((v(1:ni-1,1:nj-1) + v(2:ni,1:nj-1) + v(1:ni-1,2:nj) &
      + v(2:ni,2:nj)))

!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number
!     INSERT

      !v_max = astag


!     Calculate the timestep using the CFL number and store it in "av%dt"
!     INSERT

      ! allocate(av%dt(1:ni-1,1:nj-1), av%dt_total(1:ni-1,1:nj-1))

      
      ! New dt for each cell

      av%dt_total(1:ni-1,1:nj-1) = av%cfl * g%l_min_SVT(1:ni-1,1:nj-1) / (v_avrg(1:ni-1,1:nj-1)+astag_avrg(1:ni-1,1:nj-1))
      av%dt(1:ni-1,1:nj-1) = av%dt_total(1:ni-1,1:nj-1)

!     Print the calculated timestep and some intermediate values
!     INSERT
      !write(6,*) "time step = ", av%dt_total

      end subroutine set_timestep


