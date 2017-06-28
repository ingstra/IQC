program langevin

   use module1, only: dp, pi, r8_normal_01

   integer, parameter :: nruns = 20000000, ntrajs=100
   real(dp), dimension(nruns) :: Xa, Ya, Xb, Yb, a_store, &
        &b_store, Xa_store, Xb_store, Ya_store, Yb_store
   real(dp) :: omega_a, omega_b, g, kappa_prime, kappa_a, kappa_b, &
        &kappa_h, dt, heating_time, dW_hx, dW_ax, dW_hy, dW_ay, dW_bx, dW_by, n_a, n_b, n_h, &
        start, finish, omega_mod, T
   integer :: i,j, nr_timesteps_heating, seed
   logical :: heating
   character(len=*), parameter :: carriage_return =  char(13)
   real(dp), allocatable :: omega_eff(:), Ua(:)

   call cpu_time(start)
   
   omega_a = 2*pi*10! GHz
   omega_b = 2*pi*0.5! GHz
   g=2*pi*0.5 ! GHz
   kappa_a=2*pi*2
   kappa_b = 2*pi*0.05 ! GHz
   kappa_h = kappa_a
   n_a = 0.01
   n_b = 0.01
   n_h = 0.325

   dt=1e-6
   T= nruns*dt

   omega_mod = omega_b
   
   heating_time = pi/omega_mod
   nr_timesteps_heating = floor(heating_time/dt)
   allocate(omega_eff(2*nr_timesteps_heating),Ua(2*nr_timesteps_heating))
   a_store = 0
   b_store = 0
   Ya_store = 0
   Yb_store = 0
   Xa_store = 0
   Xb_store = 0
   omega_eff=0
   Ua=0
      
   heating = .true. 

    do j=1,ntrajs
      seed = j
      write(*,"(2a,i10,$)") carriage_return,"Calculating trajectory: ", j
      Xa = 0
      Xb = 0
      Ya = 0
      Yb = 0
      
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
          dW_hx = sqrt(kappa_h*n_h*dt)*r8_normal_01(seed)
          dW_hy = sqrt(kappa_h*n_h*dt)*r8_normal_01(seed)
       else
          kappa_prime = kappa_a/2._dp
          dW_hx = 0
          dW_hy = 0
       endif

       dW_ax = sqrt(kappa_a*n_a*dt)*r8_normal_01(seed)
       dW_ay = sqrt(kappa_a*n_a*dt)*r8_normal_01(seed)
       dW_bx = sqrt(kappa_b*n_b*dt)*r8_normal_01(seed)
       dW_by = sqrt(kappa_b*n_b*dt)*r8_normal_01(seed)
       
        ! evolve         
       Xa(i+1) = Xa(i) + (( omega_a*Ya(i) - g*sqrt(2._dp)*Xb(i)*Ya(i) -&
            &kappa_prime*Xa(i) )*dt + dW_hx  + dW_ax)

       Xa_store(i+1) = Xa_store(i+1) + Xa(i+1)
       !
       Ya(i+1) = Ya(i) - (( omega_a*Xa(i)- g*sqrt(2._dp)*Xa(i)*Xb(i)+ &
            &kappa_prime*Ya(i) )*dt + dW_hy + dW_ay)

       Ya_store(i+1) = Ya_store(i+1) + Ya(i+1)
       !
       Xb(i+1) =  Xb(i) + (( omega_b*Yb(i)-kappa_b*Xb(i)/2._dp )*dt + dW_bx)

       Xb_store(i+1) = Xb_store(i+1) + Xb(i+1)
       !       
       Yb(i+1) = Yb(i) - (( omega_b*Xb(i) + kappa_b*Yb(i)/2._dp - &
            &g*(Xa(i)**2 + Ya(i)**2)/sqrt(2._dp) )*dt + dW_by)

       Yb_store(i+1) = Yb_store(i+1) + Yb(i+1)

       ! limit cycle
       if (i .ge. nruns - 2*nr_timesteps_heating) then
          omega_eff(i-nruns +2*nr_timesteps_heating+1) = &
               &omega_eff(i-nruns +2*nr_timesteps_heating+1)+  omega_a - g*sqrt(2._dp)*Xb(i+1)
          Ua(i-nruns +2*nr_timesteps_heating+1) = &
               &Ua(i-nruns +2*nr_timesteps_heating+1) + omega_eff(i-nruns +2*nr_timesteps_heating+1)*(Xa(i+1)**2 + Ya(i+1)**2)/2._dp
       end if

       !
       a_store(i+1) = a_store(i+1) + (Xa(i+1)**2 + Ya(i+1)**2)/2._dp
       b_store(i+1) = b_store(i+1) + (Xb(i+1)**2 + Yb(i+1)**2)/2._dp
       
       if ( isnan(Yb(i+1))) then
          print *, 'Got NaN', i
          call exit(1)
       endif
    enddo
 end do
   print*,
   
   a_store = a_store/ntrajs
   b_store = b_store/(ntrajs*sqrt(T))
   Xa_store = Xa_store/ntrajs
   Ya_store = Ya_store/ntrajs
   Xb_store = Xb_store/ntrajs
   Yb_store = Yb_store/ntrajs

   omega_eff = omega_eff/ntrajs
   Ua = Ua/ntrajs

   open(unit=1, file='langevin.dat', action="write")
   open(unit=2, file='phase.dat', action="write")
   open(unit=3, file='cycle.dat', action="write")
    do i=1,nruns
       
    !   write(1,'(E22.7,A1,E22.7,A1,E22.7)') omega_a*i*dt, char(9), (Xa_store(i)**2 + Ya_store(i)**2)/2._dp,&
       !       &char(9), (Xb_store(i)**2 + Yb_store(i)**2)/2._dp
       
       write(1,'(E22.7,A1,E22.7,A1,E22.7)') i*dt, char(9), a_store(i), char(9), b_store(i)

       write(2,'(E22.7,A1,E22.7)') Xb_store(i), char(9), Yb_store(i)

    end do

    do i=1,2*nr_timesteps_heating
          omega_eff(i) =   omega_a - g*sqrt(2._dp)*Xb_store(nruns-2*nr_timesteps_heating+i)
         write(3,'(E22.7,A1,E22.7)') omega_eff(i)/omega_a, char(9), &
              & omega_eff(i)*(Xa_store(i)**2 + Ya_store(i)**2)/(2._dp*omega_a)
      end do
      
    close(1)
   ! close(2)
    close(3)

    call cpu_time(finish)
    print '("Time = ",f10.3," seconds for  solution.")',finish-start
    
end program langevin
