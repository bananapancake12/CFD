      
      subroutine flow_guess(av,g,bcs,guesstype)

!     This calculates an initial guess of the primary flowfield variables and
!     stores them at the nodes within the mesh dataytype

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid
      
!     Variables required for the crude guess
      real :: t_out, v_out, ro_out, lx, ly, l

!     Variables required for the improved guess, you will need to add to these
      real :: l_i(g%ni)
!     INSERT
      real, allocatable :: ylen(:), v_guess(:), t_guess(:), p_guess(:), ro_guess(:) ! guesses per i-line
      real, allocatable :: lx_new(:,:), ly_new(:,:), l_new(:,:)
      real :: len_end, massflow, mach_lim, t_lim
      


!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Assuming isentropic flow to the the exit plane calculate the static
!     temperature and the exit velocity
      t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
      v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
      ro_out = bcs%p_out / (av%rgas * t_out)

      ! write(6,*) "vout", v_out
      
      ! write(6,*) "ro_out = ", ro_out
      ! write(6,*) "t_out = ", t_out
      ! write(6,*) "cp = ", av%cp
      ! write(6,*) "tstag = ", bcs%tstag


!     Determine which guess calcation method to use by the value of "guesstype"
      if(guesstype == 1) then

!         Store the exit density and internal energy as if they were uniform 
          g%ro = ro_out 
          g%roe  = g%ro * (av%cv * t_out + 0.5 * v_out**2)
          

!         Calculate the gradient of the mesh lines in the centre of the domain
!         to determine the assumed direction of the flow
          j_mid = nj / 2
          do i = 1,ni-1
              lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
              l = hypot(lx,ly)
              g%rovx(i,:) = g%ro(i,:) * v_out * ly / l !pho*u= momentum eqn velocity comp
              g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
          end do

!         Copy the values to the "i = ni" nodes as an approximation
          g%rovx(ni,:) = g%rovx(ni-1,:)
          g%rovy(ni,:) = g%rovy(ni-1,:)

!         Print the guess that has been calculated
          write(6,*) 'Crude flow guess calculated'
          write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
              g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
          write(6,*)

!            else if(guesstype == 2) then
 
! !         Calculate the length of each "i = const" line between the "j = 1" and
! !         "j = nj" boundaries of the domain and store it in the local variable
! !         "l_i". You could calculate the length along each i-facet from the x
! !         and y projected lengths with "hypot" and then sum them up in the
! !         second dimension with "sum".
! !         INSERT
!           l_i = sum(hypot(g%lx_i,g%ly_i),2)
 
! !         Use the exit temperature, density and velocity calculated for the
! !         crude guess with "l_i" to estimate the mass flow rate at the exit
! !         INSERT
!           massflow = ro_out * l_i(ni) * v_out
 
! !         Set a limit to the maximum allowable mach number in the initial
! !         guess, call this "mach_lim", calculate the corresponding temperature,
! !         called "t_lim"
! !         INSERT
!           mach_lim = 1.0
!           t_lim = bcs%tstag / (1.0 + 0.5 * (av%gam-1.0) * mach_lim**2.0)
 
! !         Now estimate the velocity and density at every "i = const" line, call
! !         the velocity "v_guess(i)" and the density "ro_guess(i)":
! !             1. Assume density is constant at the exit value
! !             2. Use continuity and "l_i(i)" to estimate the velocity
! !             3. Assume stagnation temperature is constant for static temp
! !             4. Limit the static temperature, lookup intrinsic "max"
! !             5. Calculate the density throughout "ro_guess(i)"
! !             6. Update the estimate of the velocity "v_guess(i)"
! !         INSERT

!           allocate(v_guess(ni))
!           allocate(t_guess(ni))
!           allocate(p_guess(ni))
!           allocate(ro_guess(ni))
!           ro_guess = ro_out
!           v_guess = massflow / (ro_guess * l_i)
!           t_guess = bcs%tstag - 0.5 * v_guess**2 / av%cp
!           t_guess = max(t_guess,t_lim)
!           p_guess = bcs%pstag * (t_guess / bcs%tstag)**(1/av%fgam)
!           ro_guess = p_guess / (av%rgas * t_guess)
!           v_guess = massflow / (ro_guess * l_i)
 
! !         Direct the calculated velocity to be parallel to the "j = const"
! !         gridlines for all values of i and j. This can be achieved with a
! !         similar calculation to the "j = nj/2" one that was performed in the
! !         crude guess. Then set all of ro, roe, rovx and rovy, note that roe
! !         includes the kinetic energy component of the internal energy.
! !         INSERT
!           do i = 1,ni-1
!               do j = 1,nj
!                  lx = g%lx_j(i,j); ly = g%ly_j(i,j);
!                  l = hypot(lx,ly);
!                  g%rovx(i,j) = ro_guess(i) * v_guess(i) * ly / l
!                  g%rovy(i,j) = -ro_guess(i) * v_guess(i) * lx / l
!               end do
!               g%ro(i,:) = ro_guess(i)
!               g%roe(i,:) = g%ro(i,:) * (av%cv * t_guess(i) + &
!                   0.5 * v_guess(i)**2)
!           end do
             
! !         Make sure the guess has been copied for the "i = ni" values too
! !         INSERT
!           g%ro(ni,:) = g%ro(ni-1,:)
!           g%roe(ni,:) = g%roe(ni-1,:)
!           g%rovx(ni,:) = g%rovx(ni-1,:)
!           g%rovy(ni,:) = g%rovy(ni-1,:)
 
! !         Print the first elements of the guess like for the crude guess
! !         INSERT
!           write(6,*) '1D flow guess calculated'
!           write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
!               g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
!           write(6,*)
 
!       end if



      else if(guesstype == 2) then 

!         Calculate the length of each "i = const" line between the "j = 1" and 
!         "j = nj" boundaries of the domain and store it in the local variable
!         "l_i". You could calculate the length along each i-facet from the x 
!         and y projected lengths with "hypot" and then sum them up in the
!         second dimension with "sum". 
!         INSERT
            l_i = sum(hypot(g%lx_i, g%ly_i), dim=2)

            ! allocate(ylen(nj-1))      
      
            ! do i = 1,ni-1
            !       ylen = 0.0d0 !set to zero before calculating lengths

            !       do j = 1, nj-1
                        
            !             ylen(j) = hypot(g%x(i,j+1) - g%x(i,j), g%y(i,j+1) - g%y(i,j))
                        
            !       end do

            !       l_i(i) =  sum(ylen)
            !       ! write(6,*) "Length of i =", i, "=", l_i(i)
            ! end do


            !len_end = hypot(g%x(i,nj) - g%x(i,1), g%y(i,nj) - g%y(i,1)) ! check they match
            ! write(6,*) "Start Length", len_end

            


!         Use the exit temperature, density and velocity calculated for the 
!         crude guess with "l_i" to estimate the mass flow rate at the exit
!         INSERT

            massflow = l_i(ni) * ro_out * v_out

            ! write(6,*) "li", l_i
            ! write(6,*) "vout", v_out
            ! write(6,*) "length out", l_i(ni)

!         Set a limit to the maximum allowable mach number in the initial
!         guess, call this "mach_lim", calculate the corresponding temperature,
!         called "t_lim"
!         INSERT
            mach_lim = 1.0
            ! t_lim = tstag / (1+ ((av%gam-1)/2)* mach_lim^2)^-1
            t_lim = bcs%tstag / ( 1.0d0 + 0.5d0*(av%gam - 1.0d0)*mach_lim**2.0 )
            

!         Now estimate the velocity and density at every "i = const" line, call 
!         the velocity "v_guess(i)" and the density "ro_guess(i)":
!             1. Assume density is constant at the exit value
!             2. Use continuity and "l_i(i)" to estimate the velocity
!             3. Assume stagnation temperature is constant for static temp
!             4. Limit the static temperature, lookup intrinsic "max"
!             5. Calculate the density throughout "ro_guess(i)"
!             6. Update the estimate of the velocity "v_guess(i)" 
!         INSERT

            allocate(v_guess(ni))
            allocate(t_guess(ni))
            allocate(p_guess(ni))
            allocate(ro_guess(ni))

            

            ro_guess = ro_out

            do i = 1,ni
                  
                  v_guess(i) = massflow / (ro_out * l_i(i))

                  t_guess(i) = max(bcs%tstag - 0.5 * v_guess(i)**2 / av%cp, t_lim) ! energy eqn for isentropic flow
                  p_guess(i) = bcs%pstag*(t_guess(i)/bcs%tstag)**(1/av%fgam)
                  ro_guess(i) = p_guess(i) / (av%rgas * t_guess(i))

                  v_guess(i) = massflow / (ro_guess(i) * l_i(i))
            end do

            
            ! write(6,*) "massflow", massflow
            ! write(6,*) "tguess", t_guess
            ! write(6,*) "li", l_i


!         Direct the calculated velocity to be parallel to the "j = const"
!         gridlines for all values of i and j. This can be achieved with a 
!         similar calculation to the "j = nj/2" one that was performed in the 
!         crude guess. Then set all of ro, roe, rovx and rovy, note that roe 
!         includes the kinetic energy component of the internal energy.
!         INSERT 

            allocate(lx_new(ni-1,nj), ly_new(ni-1,nj), l_new(ni-1,nj))

            do i = 1, ni-1
                  do j = 1, nj
                        lx_new(i,j) = g%lx_j(i,j)
                        ly_new(i,j) = g%ly_j(i,j) !l_i are lengths along i= const
                        l_new(i,j) = hypot(lx_new(i,j), ly_new(i,j))
                        g%rovx(i,j) = ro_guess(i) * v_guess(i) * ly_new(i,j) / l_new(i,j) !pho*u= momentum eqn velocity comp
                        g%rovy(i,j) = -ro_guess(i) * v_guess(i) * lx_new(i,j) / l_new(i,j)
                        g%ro(i,j) = ro_guess(i)
                        g%roe(i,j)  = g%ro(i,j) * (av%cv * t_guess(i) + 0.5 * v_guess(i)**2)
                  end do
            end do 
            
            ! write(6,*) "l hypot",  l_new
            ! write(6,*) "lx",  lx_new
            ! write(6,*) "ly",  ly_new
            ! write(6,*) "nj",  nj

!         Make sure the guess has been copied for the "i = ni" values too
!         INSERT
            g%rovx(ni,:) = g%rovx(ni-1,:)
            g%rovy(ni,:) = g%rovy(ni-1,:)
            g%ro(ni,:) = g%ro(ni-1,:)
            g%roe(ni,:) = g%roe(ni-1,:)



!         Print the first elements of the guess like for the crude guess
!         INSERT
            write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', g%roe(1,1), &
                'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)      
            
            write(6,*) '  At 2nd point ro =', g%ro(2,1), 'roe =', g%roe(2,1), &
                'rovx =', g%rovx(2,1), 'rovy =', g%rovy(2,1)


      end if

      ! write(6,*) "vguess", v_guess
      

!     The initial guess values derived from the boundary conditions are also
!     useful as a reference to non-dimensionalise the convergence
      av%ro_ref = sum(g%ro(1,:)) / nj
      av%roe_ref = sum(g%roe(1,:)) / nj
      av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj

      end subroutine flow_guess


