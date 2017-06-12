
program main

  use omp_lib
  use module1, only: dp, print_matrix_real, trace, dagger, superH, superG,kronecker, expm,&
& liouvillian, delta_rho, r8_normal_01, delta_rho_homodyne, superD, homodyne_detection,&
&  exact_solution, factorial,print_matrix, identity_matrix,pi



  implicit none

  complex(dp), dimension(2,2) :: a,adagger,b,bdagger,rhozero_a,rhozero_b, H_traj
  complex(dp) :: i = complex(0,1),rhovec_zero(16)

 complex(dp), dimension(4,4) :: rhozero,H
  real(dp), parameter ::  dt = 1e-3
  real(dp) ::  gamma, Omega, w, start, finish, ompstart, ompend, const,omega_a,omega_b,T_a,T_b,g,planck,kb,n_c
  integer ::  channels, nruns, ntrajs,nangles
!logical :: milstein = .true.
  logical :: milstein = .false.

  n_c = 0.01
  omega_a = 2*pi*10!e9
  omega_b = 2*pi*0.5!500e6
  T_a = 104e-3
  T_b = 5e-3
  g=2*pi*0.5


  planck = 6.6261e-11
  kb= 1.3806

  a = reshape ( (/double precision :: 0,0,sqrt(1._dp),0/),(/2,2/) )
  b = a
  adagger = reshape ( (/double precision :: 0,sqrt(1._dp),0,0/),(/2,2/) )
  bdagger = adagger

!  rhozero_a = expm(-omega_a*planck/(kb*T_a),matmul(adagger,a))

  rhozero_a = reshape ( (/double precision :: 1,0,0,n_c/),(/2,2/) )
  rhozero_a = rhozero_a/trace(rhozero_a)

 ! rhozero_b = expm(-omega_b*planck/(kb*T_b),matmul(bdagger,b))
 ! rhozero_b = rhozero_b/trace(rhozero_b)
  rhozero_b =  rhozero_a

  rhozero = kronecker(rhozero_a,rhozero_b)
  !rhovec_zero = (/0,0,0,1/)

  call print_matrix_real(real(matmul(adagger,a)*rhozero_a))

  nruns = 20000
  ntrajs = 500
  nangles= 20

   gamma = 1
   Omega = 0.5
   w = 1


   H = omega_a*kronecker(matmul(adagger,a),identity_matrix(2)) + &
        & omega_b*kronecker(identity_matrix(2),matmul(bdagger,b)) + &
        & g*matmul(kronecker(matmul(adagger,a),identity_matrix(2)) , kronecker(identity_matrix(2), b + bdagger))

!H= omega_b*kronecker(identity_matrix(2),matmul(bdagger,b)) + &
 !       & g*n_c*kronecker(matmul(bdagger,b),identity_matrix(2))

   ! resonant fluorescence
   !H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 -i*sqrt(const)*Omega*(sigma_plus-sigma_minus)
  ! H_traj =  w*sigma_z/2

! H_traj = -i*sqrt(const)*Omega*(sigma_plus-sigma_minus) !Omega*(sigma_plus + sigma_minus)/2

call cpu_time(start)

call homodyne_detection(nruns,ntrajs,dt,rhozero,a,adagger,b,bdagger,gamma,'traj.dat',H,milstein,channels,omega_a)

call cpu_time(finish)

print '("Time = ",f10.3," seconds for trajectory solution.")',finish-start


call cpu_time(start)
!call exact_solution(nruns,sigma_minus,sigma_plus,dt,'exact.dat',t_start,rhovec_zero, H)
call cpu_time(finish)
 print '("Time = ",f6.3," seconds for exact solution.")',finish-start

end program main
