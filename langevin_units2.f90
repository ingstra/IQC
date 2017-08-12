program langevin

   use module1, only: dp, pi, r8_normal_01

   integer, parameter :: nruns = 10000, ntrajs=1
   real(dp), dimension(nruns) :: Xa, Ya, Xb, Yb, a_store, &
        &b_store, Xa_store, Xb_store, Ya_store, Yb_store
   real(dp) :: omega_a, omega_b, omega_h,g, kappa_prime, kappa_a, kappa_b, &
        &kappa_h, dt, heating_time, dW_hx, dW_ax, dW_hy, dW_ay, dW_bx,  kb,&
        start, finish,  omega_eff, Ua, Ta, Tb, Th, omega_0, &
        & C_c, C_k, L_k, kappa, Q, Z_0, K, k_b, hbar, hbar_over_kb, &
        noise_a, noise_b,noise_h, x_a, x_b,x_h,c_tmp,gamma,phi_0,W, &
        & temp_conv
   integer :: i,j, nr_timesteps_heating, seed
   logical :: heating
   character(len=*), parameter :: carriage_return =  char(13)
   !real(dp), allocatable :: omega_eff(:), Ua(:)

   call cpu_time(start)
   
   omega_a = 2*pi*10d9! GHz
   omega_b = 2*pi*0.5d9! GHz

   dt=1e-4

   !UNITS GHz, kb=1, hbar=1
   temp_conv=130.920295801576

   phi_0 = 2.067833831d-15 ! Wb (kg*m^2*s^-2*A-1) (RESCALED T=1e9)
   Z_0 = 50 ! Ohm
   C_c = 1d-3 ! fF
   C_k = 1d-3 ! pF
   L_k =1d0 ! pH


  Ta = 105d-3*temp_conv ! K
  Tb = 16d-3*temp_conv ! K
  Th = 100*temp_conv! K (n_h=0.125)
   
   kappa = C_c/(C_c+C_k)
   omega_0 = 1d0/sqrt(L_k*(C_c + C_k))
   Q= (C_c+C_k)/(Z_0*C_c**2*omega_0)
   K = sqrt(1._dp - 1._dp/(4._dp*Q**2) )

   gamma=omega_0/(2d0*Q)
   W = sqrt(omega_0**2-gamma**2)

   print *, 'kappa',kappa,'omega_0',omega_0,'Q',Q,'K',K,'gamma',gamma,'W',W

 
   
   x_a = hbar_over_kb *omega_a/(2*Ta)
   noise_a = 4d0*kb*Ta*Z_0 !4*pi*kappa*omega_0*sqrt( Z_0*hbar*omega_a/2d0*cosh(x_a)/sinh(x_a) )

 
   
   ! x_b = hbar_over_kb *omega_b/(2*Tb)
  !  noise_b = sqrt( Z_0*hbar*omega_b/(2._dp * omega_0**2)*cosh(x_b)/sinh(x_b) )

      x_h = hbar_over_kb *omega_a/(2*Th)
      noise_h = 4.d0*kb*Th*Z_0!4*pi*kappa*omega_0*sqrt( Z_0*hbar*omega_a/2d0*cosh(x_h)/sinh(x_h) )

        print *, 'noise_h',noise_h, 'noise_a',noise_a
   !print *, '4*Q*kappa*omega_0**2',4*Q*kappa*omega_0**2
   !print *, 'Ynoise',2.d0*kappa*omega_0/(K*Z_0)

   
    heating_time = 2d-9
    nr_timesteps_heating = floor(heating_time/dt)
   print *, 'nr',nr_timesteps_heating
   if (mod(nr_timesteps_heating,10) .eq. 0) then
      nr_timesteps_heating = nr_timesteps_heating -1
   end if
  !  allocate(omega_eff(2*nr_timesteps_heating),Ua(2*nr_timesteps_heating))
   a_store = 0
   b_store = 0
   Ya_store = 0
   Yb_store = 0
   Xa_store = 0
   Xb_store = 0
   omega_eff=0
   Ua=0
      
   heating = .true.
   
   !open(unit=4, file='noise.dat', action="write")
   
    do j=1,ntrajs
      seed = j
      write(*,"(2a,i10,$)") carriage_return,"Calculating trajectory: ", j
      Xa = 0
      Xb = 0
      Ya = 0
      Yb = 0
      
      Xa(1) =0! sqrt(0.01)
      Xb(1) = 0!sqrt(0.01)
      Ya(1) = 0!sqrt(0.01)
      Yb(1) = 0!sqrt(0.01)
      a_store(1) = 0
      b_store(1) = 0
      
      do i=1,nruns-1
         ! toggle heating
         
       if (mod(i,nr_timesteps_heating) .eq. 0) then
          if (heating .eqv. .false.) then
             heating = .true.
             print *, 'heating on'
           else
             heating = .false.
           endif
       endif
       if (heating .eqv. .true.) then
          kappa_prime = (kappa_h + kappa_a)/2._dp
          dW_hx = sqrt(noise_h*dt)*r8_normal_01(seed)
          dW_hy = sqrt(noise_h*dt)*r8_normal_01(seed)
       else
          dW_hx = 0
          dW_hy = 0
       endif

       dW_ax = sqrt(noise_a*dt)*r8_normal_01(seed)
       dW_ay = sqrt(noise_a*dt)*r8_normal_01(seed)
       dW_bx = sqrt(noise_b*dt)*r8_normal_01(seed)
       dW_by = sqrt(noise_b*dt)*r8_normal_01(seed)

              
        ! evolve         
       Xa(i+1) = Xa(i) + W*Z_0*Ya(i)*dt - gamma*Xa(i)*dt +(dW_hx  + dW_ax)
!print *, '1',W*Z_0*Ya(i)*dt,'2', gamma*Xa(i)*dt,'3',(dW_hx  + dW_ax)
    
       Xa_store(i+1) = Xa_store(i+1) + Xa(i+1)
       !
       Ya(i+1) = Ya(i) -omega_0*Xa(i)*dt/Z_0- gamma*Xa(i)*dt !+(dW_hy + dW_ay)



 !     print *, 'noise' , noise_a,dW_ax
       
       Ya_store(i+1) = Ya_store(i+1) + Ya(i+1)
       !
   !    Xb(i+1) =  Xb(i) + (( omega_b*Yb(i)-kappa_b*Xb(i)/2._dp )*dt + dW_bx)!/sqrt(T)

    !   Xb_store(i+1) = Xb_store(i+1) + Xb(i+1)
       !       
     !  Yb(i+1) = Yb(i) - (( omega_b*Xb(i) + kappa_b*Yb(i)/2._dp - &
      !      &g*(Xa(i)**2 + Ya(i)**2)/sqrt(2._dp) )*dt + dW_by)!/sqrt(T)

      ! Yb_store(i+1) = Yb_store(i+1) + Yb(i+1)

        !write(4,'(E22.7,A1,E22.7)') i*dt, char(9), dW_ax

       !
       a_store(i+1) = a_store(i+1) + (Xa(i+1)**2 + Ya(i+1)**2)/2._dp
!       b_store(i+1) = b_store(i+1) + (Xb(i+1)**2 + Yb(i+1)**2)/2._dp
       
       if ( isnan(Ya(i+1))) then
          print *, 'Got NaN', i
          call exit(1)
       endif
    enddo
 end do
   print*,

   close(4)
   
   a_store = a_store/ntrajs
  ! b_store = b_store/ntrajs!(ntrajs*sqrt(T))
   Xa_store = Xa_store/ntrajs
   Ya_store = Ya_store/ntrajs
  ! Xb_store = Xb_store/ntrajs
   !Yb_store = Yb_store/ntrajs

  ! omega_eff = omega_eff/ntrajs
  ! Ua = Ua/ntrajs

   open(unit=1, file='langevin.dat', action="write")
   open(unit=2, file='phase.dat', action="write")
   open(unit=3, file='cycle.dat', action="write")
   open(unit=5, file='energy.dat',action='write')
    do i=1,nruns
       
    !   write(1,'(E22.7,A1,E22.7,A1,E22.7)') omega_a*i*dt, char(9), (Xa_store(i)**2 + Ya_store(i)**2)/2._dp,&
       !       &char(9), (Xb_store(i)**2 + Yb_store(i)**2)/2._dp

       if (mod(i,1000) .eq. 0) then
          write(1,'(E22.7,A1,E22.7)') i*dt, char(9), a_store(i)
      end if

   !    write(2,'(E22.7,A1,E22.7)') Xb_store(i), char(9), Yb_store(i)


     !  omega_eff =   omega_a - g*sqrt(2._dp)*Xb_store(i)
      ! Ua = omega_eff*a_store(i)

       !write(5,'(E22.7,A1,E22.7)') i*dt, char(9), Ua

    end do

    close(5)

   ! do i=70,2*nr_timesteps_heating
   !    omega_eff =   omega_a - g*sqrt(2._dp)*Xb_store(nruns-2*nr_timesteps_heating+i)
  !     Ua = omega_eff*a_store(i)
  !     write(3,'(E22.7,A1,E22.7)') omega_eff/omega_a, char(9), Ua/omega_a
  !  end do

    
      
    close(1)
   ! close(2)
    close(3)

    call cpu_time(finish)
    print '("Time = ",f10.3," seconds for  solution.")',finish-start
    
end program langevin
