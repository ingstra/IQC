program langevin

   use module1, only: dp, pi, r8_normal_01

   integer, parameter :: nruns = 200000, ntrajs=1
   real(dp), dimension(nruns) :: Xa, Ya, Xb, Yb, a_store, &
        &b_store, Xa_store, Xb_store, Ya_store, Yb_store
   real(dp) :: omega_a, omega_b, omega_h,g, kappa_prime, kappa_a, kappa_b, &
        &kappa_h, dt, heating_time, dW_hx, dW_ax, dW_hy, dW_ay, dW_bx, dW_by, n_a, n_b, n_h, &
        start, finish, omega_mod, omega_eff, Ua, Ta, Tb, Th, omega_0, &
        & C_c, C_k, L_k, kappa, Q, Z_0, K, k_b, hbar, hbar_over_kb, &
        noise_a, noise_b,noise_h, x_a, x_b,x_h,c_tmp
   integer :: i,j, nr_timesteps_heating, seed
   logical :: heating
   character(len=*), parameter :: carriage_return =  char(13)
   !real(dp), allocatable :: omega_eff(:), Ua(:)

   call cpu_time(start)
   
   omega_a = 2*pi*10d9! GHz
   omega_b = 2*pi*0.5d9! GHz
   g=2*pi*0.5!d9 ! GHz
   kappa_a=2*pi*2!d9
   kappa_b = 2*pi*0.05!d9! GHz
   kappa_h = kappa_a
   omega_h = omega_a
   n_a = 0.01
   n_b = 0.01
   n_h = 0.125

   dt=1e-13

   Z_0 = 50 ! Ohm
   C_c = 1d-15 ! fF
   C_k = 1d-12 ! pF
   L_k =1d-9 ! pH

  hbar = 1.0545718d-34 ! Planck's constant m^2 kg s^-1
  kb= 1.38064852d-23 ! Boltzmann's constant m^2 kg s^-2 K^-1

  hbar_over_kb = 1.0545718d-11 / 1.38064852

  Ta = 105d-3 ! K
  Tb = 16d-3 ! K
  Th = 100d-3 ! K
   
   kappa = C_c/(C_c+C_k)
   omega_0 = 1d0/sqrt(L_k*(C_c + C_k))
   Q= (C_c+C_k)/(Z_0*C_c**2*omega_0)
   K = sqrt(1._dp - 1._dp/(4._dp*Q**2) )

   print *, 'kappa',kappa,'omega_0',omega_0,'Q',Q,'K',K

 
   
   x_a = hbar_over_kb *omega_a/(2*Ta)
   noise_a = sqrt( Z_0*hbar*omega_a/(2._dp * omega_0**2)*cosh(x_a)/sinh(x_a) )

  ! print *, 'noise_a',noise_a,'noise_h',noise_h
   !print *, '4*Q*kappa*omega_0**2',4*Q*kappa*omega_0**2
   !print *, 'Ynoise',2.d0*kappa*omega_0/(K*Z_0)

   
    x_b = hbar_over_kb *omega_b/(2*Tb)
    noise_b = sqrt( Z_0*hbar*omega_b/(2._dp * omega_0**2)*cosh(x_b)/sinh(x_b) )

      x_h = hbar_over_kb *omega_a/(2*Th)
    noise_h = sqrt( Z_0*hbar*omega_a/(2._dp * omega_0**2)*cosh(x_h)/sinh(x_h) )  
   
    heating_time = 2
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
      
      Xa(1) = sqrt(0.01)
      Xb(1) = sqrt(0.01)
      Ya(1) = sqrt(0.01)
      Yb(1) = sqrt(0.01)
      a_store(1) = 0.01*ntrajs
      b_store(1) = 0.01*ntrajs
      
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
       Xa(i+1) = Xa(i) + omega_0*K*Z_0*Ya(i)*dt - omega_0/(2._dp*Q)*Xa(i)*dt !+ 4*Q*kappa*omega_0**2*(dW_hx  + dW_ax)!/sqrt(T)

     
  !   print *,'1', omega_0*K*Z_0*Ya(i)*dt, '2',omega_0/(2._dp*Q)*Xa(i) *dt,'3', 4*Q*kappa*omega_0**2*(dW_hx  + dW_ax)

       Xa_store(i+1) = Xa_store(i+1) + Xa(i+1)
       !
       Ya(i+1) = Ya(i) -omega_0*K*Xa(i)*dt/Z_0- omega_0/(2._dp*Q)*Xa(i)*dt !+ 2.d0*kappa*omega_0/(K*Z_0)*(dW_hy + dW_ay)

  !           print *,'1y', omega_0*K*Xa(i)/Z_0*dt, '2y', omega_0/(2._dp*Q)*Xa(i) *dt ,'3',2.d0*kappa*omega_0/(K*Z_0)*(dW_hy + dW_ay)

!       print *, 'noise' , dW_ax
       
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

      ! if (mod(i,10) .eq. 0) then
          write(1,'(E22.7,A1,E22.7)') i*dt, char(9), a_store(i)
       !end if

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
