
program main

  use omp_lib
  use module1, only: dp, print_matrix_real, trace, dagger, superH, superG,kronecker, expm,excited_probability,&
& liouvillian, delta_rho, r8_normal_01, superD, power_output,&
&  exact_solution, factorial,print_matrix, identity_matrix,pi,get_a



  implicit none

  real(dp), parameter ::  dt = 5e-3
  real(dp) :: start, finish, const,omega_a,omega_b,T_a,T_b,g,hbar_over_kb,n_c,kappa_b
  integer ::   nruns, ntrajs
  integer, parameter :: dim = 6
 
  complex(dp), dimension(dim,dim) ::a,adagger,b,bdagger,rhozero_a,rhozero_b
  complex(dp), dimension(dim**2,dim**2) :: rhozero,H_sys
!logical :: milstein = .true.
  logical :: milstein = .false.
  

  n_c = 0.01
  omega_a = 2*pi*10! GHz
  omega_b = 2*pi*0.5! GHz
  kappa_b = 2*pi*0.05 ! GHz
  T_a = 13.6736048505798
  T_b =  0.683680242528990
  g=2*pi*0.5 ! GHz

  !hbar = 1.0545718e-34 ! Planck's constant m^2 kg s^-1
  !kb= 1.38064852e-23 ! Boltzmann's constant m^2 kg s^-2 K^-1

  !hbar_over_kb = 1.0545718e-11 / 1.38064852
  
  call get_a(dim,a)
  b = a
  adagger = dagger(a)
  bdagger = adagger

  rhozero_a = expm(-omega_a/T_a,matmul(adagger,a))

  rhozero_a = rhozero_a/trace(rhozero_a)

  rhozero_b =  rhozero_a

  rhozero = kronecker(rhozero_a,rhozero_b)

  call print_matrix_real(real(rhozero_a))

  nruns = 10000
  ntrajs = 500


   H_sys = omega_a*kronecker(matmul(adagger,a),identity_matrix(dim)) + &
        & omega_b*kronecker(identity_matrix(dim),matmul(bdagger,b)) - &
        & g*matmul(kronecker(matmul(adagger,a),identity_matrix(dim)) , kronecker(identity_matrix(dim), b + bdagger))


   
call cpu_time(start)

!call power_output(nruns,ntrajs,dt,rhozero,a,adagger,'traj.dat',H_sys,milstein,omega_a,omega_b,dim,n_c,kappa_b)

call cpu_time(finish)

print '("Time = ",f10.3," seconds for trajectory solution.")',finish-start

call cpu_time(start)

call excited_probability(nruns,ntrajs,dt,rhozero,a,adagger,'traj.dat',H_sys,milstein,omega_a,omega_b,dim,n_c,kappa_b)

call cpu_time(finish)

print '("Time = ",f10.3," seconds for trajectory solution.")',finish-start


call cpu_time(start)
!call exact_solution(nruns,sigma_minus,sigma_plus,dt,'exact.dat',t_start,rhovec_zero, H)
call cpu_time(finish)
 print '("Time = ",f6.3," seconds for exact solution.")',finish-start

end program main
