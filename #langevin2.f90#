program langevin

   use module1, only: dp, pi, r8_normal_01

   integer, parameter :: nruns = 20000, ntrajs=1
   complex(dp), dimension(nruns) ::a,b, a_store, &
        &b_store
   real(dp) :: omega_a, omega_b, g, kappa_prime, kappa_a, kappa_b, &
        &kappa_h, dt, heating_time, dW_hx, n_a, n_b, n_h, &
        start, finish, omega_mod, T
   integer :: i,j, nr_timesteps_heating, seed
   logical :: heating
   character(len=*), parameter :: carriage_return =  char(13)
 !  real(dp), allocatable :: omega_eff(:), Ua(:)
   complex(dp) :: im = complex(0,1), eta_a, eta_b, c0_a, c0_b,&
        vth_a,vth_b

   call cpu_time(start)
   
   omega_a = 2*pi*10! GHz
   omega_b = 2*pi*0.5! GHz
   g=2*pi*0.5 ! GHz
   kappa_a=2*pi*2
   kappa_b = 2*pi*0.05 ! GHz
   kappa_h = kappa_a
   n_a = 0.01
   n_b = 0.01
   n_h = 0.125

   dt=1e-3
   T= nruns*dt

   omega_mod = omega_b
   
   heating_time = pi/omega_mod
   nr_timesteps_heating = floor(heating_time/dt)
 !  allocate(omega_eff(2*nr_timesteps_heating),Ua(2*nr_timesteps_heating))
   a_store = 0
   b_store = 0
   a = 0
   b = 0
       
   heating = .true. 

    do j=1,ntrajs
      seed = j
      write(*,"(2a,i10,$)") carriage_return,"Calculating trajectory: ", j
      a = 0
      b = 0
            
      do i=1,nruns-1
         ! toggle heating
         
       if (mod(i,nr_timesteps_heating) .eq. 0) then
          if (heating .eqv. .false.) then
             heating = .true.
          else
             heating = .false.
          endif
       endif
       if (heating .eqv. .true.) then
          kappa_prime = (kappa_h + kappa_a)/2._dp
          dW_hx = sqrt(kappa_h*n_h)*r8_normal_01(seed)
       else
          kappa_prime = kappa_a/2._dp
          dW_hx = 0
       endif
       
       ! evolve
       eta_a = -(im*omega_a-im*g*(conjg(b(i)) + b(i) ) + kappa_prime/2._dp) 
       eta_b = -(im*omega_b + kappa_b/2._dp)*b(i) + im*g*conjg(a(i))*a(i)

       vth_a = kappa_a*n_a/2._dp!*eta_a)
       vth_b = kappa_b*n_b/2._dp!*eta_b)

       c0_a = exp(-eta_a*dt)
       c0_b = exp(-eta_b*dt)

       a(i+1) = c0_a*a(i) + vth_a*sqrt(1._dp - c0_a**2)*r8_normal_01(seed)+ dW_hx
       b(i+1) = c0_b*b(i) + vth_b*sqrt(1._dp - c0_b**2)*r8_normal_01(seed)
       

      ! print *,
       if ( isnan(real(b(i+1)))) then
          print *, 'Got NaN', i
          call exit(1)
       endif
    enddo
 end do
   print*,
   
   a = a/ntrajs
   b = b/ntrajs!(ntrajs*sqrt(T))
  

  ! omega_eff = omega_eff/ntrajs
   !Ua = Ua/ntrajs

   open(unit=1, file='langevin.dat', action="write")
 !  open(unit=2, file='phase.dat', action="write")
 !  open(unit=3, file='cycle.dat', action="write")
    do i=1,nruns
       
    !   write(1,'(E22.7,A1,E22.7,A1,E22.7)') omega_a*i*dt, char(9), (Xa_store(i)**2 + Ya_store(i)**2)/2._dp,&
       !       &char(9), (Xb_store(i)**2 + Yb_store(i)**2)/2._dp
       
       write(1,'(E22.7,A1,E22.7,A1,E22.7)') i*dt, char(9), real(a(i)*conjg(a(i))), char(9), real(b(i)*conjg(b(i)))

     !  write(2,'(E22.7,A1,E22.7)') Xb_store(i), char(9), Yb_store(i)

    end do

!    do i=1,2*nr_timesteps_heating
!          omega_eff(i) =   omega_a - g*sqrt(2._dp)*Xb_store(nruns-2*nr_timesteps_heating+i)
!         write(3,'(E22.7,A1,E22.7)') omega_eff(i)/omega_a, char(9), &
!              & omega_eff(i)*(Xa_store(i)**2 + Ya_store(i)**2)/(2._dp*omega_a)
!      end do
      
    close(1)
   ! close(2)
 !   close(3)

    call cpu_time(finish)
    print '("Time = ",f10.3," seconds for  solution.")',finish-start
    
end program langevin
