module module1

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: pi = 3.141592653589793

contains

  pure function factorial(n) result(y)
    integer, intent(in) :: n

    integer :: i
    real(dp) :: y

   y = PRODUCT((/(i, i=1,n)/))

  end function factorial

 subroutine print_matrix_real(matrix,long_flag)
    implicit none

    !input
    !>matrix which you want to print
    real(dp),intent(in) :: matrix(:,:)

    !>true value make subroutine print with longer format
    logical, intent(in), optional :: long_flag

    !locals
    !row index of matrix
    integer :: i

    !format string for output
    character(20) :: format

    if (present(long_flag).and.long_flag) then
       format='(20G20.12)'
    else
       format='(20G12.4)'
    end if

    do i=1,size(matrix,1)
       write(*,format) real(matrix(i,:))
    end do

  end subroutine print_matrix_real

  subroutine print_matrix(matrix,long_flag)
    implicit none

    !input
    !>matrix which you want to print
    complex(dp),intent(in) :: matrix(:,:)

    !>true value make subroutine print with longer format
    logical, intent(in), optional :: long_flag

    !locals
    !row index of matrix
    integer :: i

    !format string for output
    character(20) :: format

    if (present(long_flag).and.long_flag) then
       format='(20G20.12)'
    else
       format='(20G12.4)'
    end if

    do i=1,size(matrix,1)
       write(*,format) matrix(i,:)
    end do

  end subroutine print_matrix

  subroutine direct_detection(nruns,ntrajs,dt,rhozero,c,cdagger,filename,H)
    integer, intent(in) :: nruns, ntrajs
    real(dp), intent(in) :: dt
    complex(dp), intent(in) :: rhozero(:,:),c(:,:),cdagger(:,:),H(:,:)
    character(len=*), intent(in) :: filename

    real(dp), allocatable :: store(:)
    integer :: k, j, dN, seed
    complex(dp) :: rho(2,2)
    allocate(store(nruns))

    seed = 123456789
    call srand(seed)
    !open(unit=3, file='trace.dat', action="write")
    store=0
    do k=1,ntrajs
       rho = rhozero
       do j=1,nruns
          !      write(3,'(E22.7,A1,E22.7)') j*dt, char(9), trace(rho)
          store(j) = store(j)+ trace(matmul(cdagger,c)*rho )
          call random_number(random)
          if (random < dt ) then
             dN = 1
          else
             dN = 0
          end if
          rho = rho + delta_rho(rho,H,c,cdagger,dt,dN)
       end do
    end do
    !close(3)

    open(unit=1, file=filename, action="write")
    do j=1,nruns
       write(1,'(E22.7,A1,E22.7)') j*dt, char(9), store(j)/ntrajs
    end do
    close(1)

  end subroutine direct_detection


!!$  subroutine power_output(nruns,ntrajs,dt,rhozero,a,adagger,filename,H,mflag,omega_a,omega_b,dim,n_c,kappa_b)
!!$
!!$    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
!!$    integer, intent(in) :: nruns, ntrajs,dim
!!$    logical, intent(in) :: mflag ! if true, use milstein scheme
!!$    real(dp), intent(in) :: dt,omega_a,omega_b,n_c,kappa_b
!!$    complex(dp), intent(in) :: rhozero(:,:), a(:,:),adagger(:,:), H(:,:)
!!$    character(len=*), intent(in) :: filename
!!$    real(dp), allocatable :: store(:)
!!$
!!$    integer :: k,j,seed
!!$    real(dp) :: dW
!!$    complex(dp), dimension(dim**2,dim**2) :: rho,tmp
!!$    character(len=*), parameter :: carriage_return =  char(13)
!!$    logical :: heating = .false.
!!$    allocate(store(nruns))
!!$
!!$    store = 0
!!$    open(unit=3, file='purity.dat', action="write")
!!$   ! $OMP PARALLEL DO private(k,j) schedule(static)
!!$    do k=1,ntrajs
!!$       rho = rhozero
!!$       seed = k
!!$       do j=1,nruns
!!$          dW = r8_normal_01(seed)*sqrt(dt)
!!$          write(3,'(E22.7,A1,E22.7)') j*dt, char(9), trace(matmul(rho,rho))
!!$
!!$          tmp=kronecker(identity_matrix(dim),matmul(adagger,a))
!!$          store(j) = store(j)+ omega_b*kappa_b*trace(matmul(tmp,rho))
!!$          if (isnan(store(j))) then ! check if invalid value (NaN)
!!$          print *, 'Got Nan. '
!!$          call exit(1)
!!$       end if ! end NaN
!!$          if (mflag .eqv. .true.) then
!!$            ! rho = rho + delta_rho_milstein(rho,a,adagger,H,dt,dW,0._dp)
!!$          else
!!$             rho = rho + delta_rho_homodyne(rho,H,a,adagger,b,bdagger,dt,dW,n_c,kappa_b,dim,heating)
!!$          end if
!!$       end do
!!$       write(stdout,"(2a,i10,$)") carriage_return,"Calculating trajectory: ", k
!!$    end do
!!$! $OMP END PARALLEL DO
!!$
!!$    write(stdout,*) linefeed
!!$    close(3)
!!$print *, 'purity', trace(matmul(rho,rho))
!!$
!!$    open(unit=1, file=filename, action="write")
!!$    do j=1,nruns
!!$       write(1,'(E22.7,A1,E22.7)') omega_a*j*dt, char(9), store(j)/ntrajs
!!$    end do
!!$    close(1)
!!$
!!$  end subroutine power_output


  subroutine excited_probability(nruns,ntrajs,dt,rhozero,a,adagger,&
       &b,bdagger,filename,H,mflag,omega_a,omega_b,dim_a,dim_b,n_c,kappa_b,g)
   use omp_lib
   use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    integer, intent(in) :: nruns, ntrajs,dim_a,dim_b
    logical, intent(in) :: mflag ! if true, use milstein scheme
    real(dp), intent(in) :: dt,omega_a,omega_b,n_c,kappa_b,g
    complex(dp), intent(in) :: rhozero(:,:), a(:,:), adagger(:,:),b(:,:),bdagger(:,:),H(:,:)
    character(len=*), intent(in) :: filename
    real(dp), allocatable :: store_a(:),store_b(:),store_tmp_a(:),store_tmp_b(:)

    integer :: k,j
    real(dp) :: dWb, dWa
    complex(dp), dimension(dim_a*dim_b,dim_a*dim_b) :: rho
    character(len=*), parameter :: carriage_return =  char(13)
    allocate(store_tmp_a(nruns),store_tmp_b(nruns),store_a(nruns),store_b(nruns))

    store = 0
   ! call OMP_SET_NUM_THREADS(4)
   ! open(unit=3, file='purity.dat', action="write")
    ! $OMP PARALLEL DO private(k,store_tmp_a,store_tmp_b) reduction(+:store_a,store_b) schedule(static)

    do k=1,ntrajs
       store_tmp_a = 0
       store_tmp_b = 0
    call evolve_homodyne(nruns,dt,rhozero,a,adagger,b,bdagger,H,&
         &mflag,omega_a,dim_a,dim_b,n_c,kappa_b,k,store_tmp_a,store_tmp_b,g)
       
       write(stdout,"(2a,i10,$)") carriage_return,"Calculating trajectory: ", k
       store_a=store_a+store_tmp_a
       store_b=store_b+store_tmp_b
    end do
    ! $OMP END PARALLEL DO 

    write(stdout,*) linefeed
  !  close(3)
    open(unit=1, file=filename, action="write")
    do j=1,nruns
       write(1,'(E22.7,A1,E22.7,A1,E22.7)') j*dt, char(9), store_a(j)/ntrajs, char(9), store_b(j)/ntrajs
    end do
    close(1)

  end subroutine excited_probability

subroutine evolve_homodyne(nruns,dt,rhozero,a,adagger,b,bdagger,H,mflag,omega_a,dim_a,dim_b,n_c,kappa_b,seed,store_a,store_b,g)

    use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
    integer, intent(in) :: nruns, dim_a,dim_b,seed
    logical, intent(in) :: mflag ! if true, use milstein scheme
    real(dp), intent(in) :: dt,omega_a,n_c,kappa_b,g
    complex(dp), intent(in) :: rhozero(:,:), a(:,:),adagger(:,:), b(:,:),bdagger(:,:),H(:,:)
    real(dp), intent(out) :: store_a(:),store_b(:)
    !real(dp), allocatable :: store(:)

    real(dp) :: y,dWb,dWa,t,const, heating_time,omega_eff,U,omega_mod
    complex(dp), dimension(dim_a*dim_b,dim_a*dim_b) :: rho
    complex(dp) ::  i
    integer :: trueseed_a, trueseed_b,j,nr_timesteps_heating
    logical :: heating
    trueseed_b=seed
    trueseed_a=seed+98765
    store=0

    ! omega_mod is for one full cycle!
    omega_mod = 2*pi*0.5

    ! half cycle
    heating_time = pi/omega_mod
    nr_timesteps_heating = ceiling(heating_time/dt)
    heating = .true.     
    rho = rhozero
    open(unit=5, file='cycle.dat', action="write")
    do j=1,nruns
       ! toggle heating
       if (mod(j,nr_timesteps_heating) .eq. 0) then
          if (heating .eqv. .false.) then
             heating = .true.
          else
             heating = .false.
          endif
       endif

        write(stdout,"(2a,i10,1a,i10,$)") char(13),"Run number ", j, " out of ", nruns
       
        dWb = r8_normal_01(trueseed_b)*sqrt(dt)
        dWa = r8_normal_01(trueseed_a)*sqrt(dt)
       store_b(j) = store_b(j) +  trace(kronecker(identity_matrix(dim_a),matmul(bdagger,b))*rho) ! b
       store_a(j) = store_a(j) +  trace(kronecker(matmul(adagger,a),identity_matrix(dim_b))*rho) ! a
       

       if (isnan(store_b(j))) then ! check if invalid value (NaN)
          print *, 'Got NaN in evolve_homodyne. '
          call exit(1)
       end if ! end NaN
       
          rho = rho + delta_rho_homodyne(rho,H,a,adagger,b,bdagger,dt,dWb,dWa,n_c,kappa_b,dim_a,dim_b,heating)
       
         
    !! *************************************************
   
       if (j .ge. nruns - 2*nr_timesteps_heating) then
          omega_eff = omega_a - g*trace(matmul(kronecker(identity_matrix(dim_a), b+bdagger),rho))
          U=omega_eff* trace(kronecker(matmul(adagger,a),identity_matrix(dim_b))*rho)
        write(5,'(E22.7,A1,E22.7)') omega_eff/omega_a, char(9), U/omega_a
     endif
  end do ! nruns
  print *,
    close(5)
    
    if (seed .eq. 1) then       
       print *, 'purity', trace(matmul(rho,rho))
    endif
   
  end subroutine evolve_homodyne

function delta_rho_homodyne(rho,H,a,adagger,b,bdagger,dt,dWb,dWa,n_c,kappa_b,dim_a,dim_b,heating) result(y)
    complex(dp), intent(in) :: rho(:,:), H(:,:), a(:,:),adagger(:,:),b(:,:),bdagger(:,:)
    real(dp), intent(in) :: dt, dWb, dWa, n_c,kappa_b
    integer, intent(in) :: dim_a, dim_b
    logical, intent(in) :: heating
    complex(dp), allocatable ::y(:,:)


    complex(dp) :: i,id_a(dim_a,dim_a), id_b(dim_b,dim_b)
    integer :: n
    real(dp) :: const, kappa_a, kappa_h, n_h
    kappa_a=2*pi*2
    n_h=0.125
    kappa_h = kappa_a

    i = complex(0,1)
    n = size(rho,1)
    allocate(y(n,n))

    id_a = identity_matrix(dim_a)
    id_b = identity_matrix(dim_b)

    y = -i*(matmul(H,rho)-matmul(rho,H))*dt + &
         &kappa_a*(n_c+1)*superD(kronecker(a,id_b),rho)*dt + kappa_a*n_c*superD(kronecker(adagger,id_b),rho)*dt + &
         &kappa_b*(n_c+1)*superD(kronecker(id_a,b),rho)*dt + kappa_b*n_c*superD(kronecker(id_a,bdagger),rho)*dt + &
         ! 
         & kappa_h*superD(kronecker(a,id_b),rho)*dt + &
         ! stochastic term
         &sqrt(kappa_b)*superH(kronecker(id_a,b),rho)*dWa + &
         &sqrt(kappa_a)*superH(kronecker(a,id_b),rho)*dWa

    if (heating .eqv. .true.) then
       !  y = y + kappa_h*(n_h+1)*superD(kronecker(a,id_b),rho)*dt + kappa_h*n_h*superD(kronecker(adagger,id_b),rho)*dt

        y = y + kappa_h*n_h*superD(kronecker(a,id_b),rho)*dt + kappa_h*n_h*superD(kronecker(adagger,id_b),rho)*dt
    endif

    if (isnan(trace(y))) then
       print *, 'Got NaN in delta_rho_homodyne.'
       call exit(1)
    end if
  end function delta_rho_homodyne


  pure function integrate(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(dp), intent(in)  :: x(:)         !! Variable x
    real(dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(dp)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end associate
end function

  subroutine exact_solution(nruns, c, cdagger, dt ,filename, rhovec_zero, H)
    integer, intent(in) :: nruns
    real(dp), intent(in) :: dt
    complex(dp), intent(in) :: rhovec_zero(:),c(:,:),cdagger(:,:),H(:,:)
   character(len=*), intent(in) :: filename

    complex(dp) :: rho(2,2), rhovec(4)

    open(unit=2, file=filename, action="write")
    do j=1,nruns
          rhovec = matmul( expm(j*dt,liouvillian(c,H)), rhovec_zero)
          rho(1,1) = rhovec(1)
          rho(2,1) = rhovec(2)
          rho(1,2) = rhovec(3)
          rho(2,2) = rhovec(4)
       write(2,'(E22.7,A1,E22.7)') j*dt, char(9), real(trace(matmul(cdagger,c)*rho ))
    end do
    close(2)

  end subroutine exact_solution


  pure function trace(matrix) result(y)
    complex(dp), intent(in) :: matrix(:,:)
    real(dp) :: y
    integer :: n
    n = size(matrix,1)
    y = 0

    do i=1,n
       y = y + real(matrix(i,i))
    enddo
  end function trace

  function superH(r,rho) result(y)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::y(:,:)

    integer :: n
    n=size(r,1)

    allocate(y(n,n))

    y = matmul(r,rho) + matmul(rho,dagger(r)) - trace( matmul(r,rho) + matmul(rho,dagger(r)) )*rho
  end function superH

  pure function superG(r,rho) result(y)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::y(:,:)
    real(dp) :: tracevalue

    integer :: n
    n=size(r,1)

    allocate(y(n,n))

    tracevalue = real( trace( matmul(matmul(r,rho),dagger(r)) ) )

    if (tracevalue /= 0) then
       y = matmul(matmul(r,rho),dagger(r))/tracevalue - rho
    else
       y = 0!- rho
    end if

  end function superG

  pure function superD(c,rho) result(y)

    complex(dp), intent(in) :: c(:,:), rho(:,:)
    complex(dp), allocatable ::y(:,:)

    complex(dp), allocatable :: cdagger(:,:)
    integer :: n
    n=size(c,1)

    allocate(y(n,n),cdagger(n,n))

    cdagger = conjg(transpose(c))

    y = matmul(matmul(c,rho),cdagger) - 0.5_dp*matmul(matmul(cdagger,c),rho) - 0.5_dp*matmul(matmul(rho,cdagger),c)

  end function superD


  pure function dagger(matrix) result(y)
    complex(dp), intent(in) :: matrix(:,:)
    complex(dp), allocatable :: y(:,:)
    integer :: n
    n = size(matrix,1)
    allocate(y(n,n))

    y = conjg(transpose(matrix))
  end function dagger

  pure function identity_matrix(n) result(y)
    integer, intent(in) :: n
    complex(dp), dimension(n,n) :: y

    y = 0
    do i=1,n
       y(i,i) = 1
    end do
  end function identity_matrix

pure SUBROUTINE linspace(d1,d2,n,grid)

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: d1, d2
DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: grid

INTEGER :: indxi


grid(1) = d1
DO indxi= 0,n-2
   grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
END DO
grid(n) = d2

!MATLAB
!grid = [d1+(0:n-2)*(d2-d1)/(floor(n)-1) d2];


END SUBROUTINE

  function delta_rho(rho,H,c,cdagger,delta_t,dN) result(y)
    complex(dp), intent(in) :: rho(:,:), H(:,:), c(:,:), cdagger(:,:)
    real(dp), intent(in) :: delta_t
    integer, intent(in) :: dN
    complex(dp), allocatable ::y(:,:)

    complex(dp) :: i
    integer :: n
    i = complex(0,1)
    n = size(rho,1)
    allocate(y(n,n))

    y = dN*superG(c,rho)-delta_t*superH(i*H + matmul(cdagger,c)/2,rho)

  end function delta_rho

  

!!$  function delta_rho_milstein(rho,c,cdagger,H,dt,dW) result(y)
!!$    complex(dp), intent(in) :: rho(:,:), c(:,:), cdagger(:,:), H(:,:)
!!$    real(dp), intent(in) :: dt, dW
!!$    complex(dp), allocatable ::y(:,:)
!!$
!!$    real(dp) :: const
!!$    complex(dp) :: i
!!$
!!$    complex(dp), dimension(2,2) :: term1, term2
!!$    integer :: n
!!$    i = complex(0,1)
!!$    n = size(rho,1)
!!$    allocate(y(n,n))
!!$
!!$    term1 = matmul(c,c)*rho + 2._dp*matmul(c,rho)*cdagger + matmul(rho,cdagger)*cdagger
!!$    term2 = matmul(c*exp(-i*theta),rho) + matmul(rho,cdagger*exp(i*theta))
!!$
!!$    y = -i*(matmul(H,rho)-matmul(rho,H))*dt + gamma*superD(c,rho)*dt + &
!!$& sqrt(const)*superH(c*exp(-i*theta),rho)*dW + 0.5_dp*(const*(term1-trace(term1)) + &
!!$& 2._dp*const*trace(term2)*(trace(term2)*rho-term2))*(dW**2._dp-dt)
!!$
!!$  end function delta_rho_milstein

  function liouvillian(c, H_eff) result(y)
    complex(dp), intent(in) :: H_eff(:,:), c(:,:)
    complex(dp) :: i = complex(0,1)
    complex(dp), dimension(2,2) :: identity
    complex(dp), allocatable :: y(:,:)
    integer :: n
    n = size(c,1)
    allocate(y(n**2,n**2))

    identity = reshape ( (/1,0,0,1/),(/2,2/) )

    y = -i*kronecker(identity,H_eff) + i*kronecker(conjg(H_eff),identity) +&
         &kronecker(conjg(c),c)

  end function liouvillian

 pure function kronecker(A,B) result(y)
    complex(dp), dimension (:,:), intent(in)  :: A, B
    complex(dp), dimension (:,:), allocatable :: y
    integer :: i, j, m, n, p, q
    i=0
    j=0
    m=0
    n=0
    p=0
    q=0

   allocate(y(size(A,1)*size(B,1),size(A,2)*size(B,2)))
   y=0
   do i = 1,size(A,1)
      do j = 1,size(A,2)
      y(1+size(B,1)*(i-1):size(B,1)+size(B,1)*(i-1),1+size(B,2)*(j-1):size(B,2)+size(B,2)*(j-1))=A(i,j)*B
      enddo
   enddo
end function kronecker


  function expm(t,H) result(expH)
    real(dp), intent(in) :: t
    complex(dp), dimension(:,:), intent(in) :: H
    complex(dp), dimension(size(H,1),size(H,2)) :: expH

    ! Expokit variables
    external :: DGPADM
    integer, parameter :: ideg = 6
    complex(dp), dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
    integer, dimension(size(H,1))  :: iwsp
    integer :: iexp, ns, iflag, n

    if (size(H,1) /= size(H,2)) then
       stop 'expm: matrix must be square'
    end if

    n = size(H,1)
    call ZGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))
  end function expm

  subroutine get_a(dim,a)
    integer, intent(in) ::dim
    complex(dp), intent(out) :: a(:,:)

    integer :: i
    !allocate(a(dim,dim))
    a=0
    do i=1,dim-1
       a(i,i+1) = sqrt(real(i))
    end do
          
  end subroutine get_a

  subroutine ZGPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

    implicit none
    real(dp) :: t
    integer ::         ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
    complex(dp) ::       H(ldh,m), wsp(lwsp)

    !-----Purpose----------------------------------------------------------|
    !
    !     Computes exp(t*H), the matrix exponential of a general complex
    !     matrix in full, using the irreducible rational Pade approximation
    !     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
    !     combined with scaling-and-squaring.
    !
    !-----Arguments--------------------------------------------------------|
    !
    !     ideg      : (input) the degre of the diagonal Pade to be used.
    !                 a value of 6 is generally satisfactory.
    !
    !     m         : (input) order of H.
    !
    !     H(ldh,m)  : (input) argument matrix.
    !
    !     t         : (input) time-scale (can be < 0).
    !
    !     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
    !
    !     ipiv(m)   : (workspace)
    !
    !>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
    !                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
    !                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !                 NOTE: if the routine was called with wsp(iptr),
    !                       then exp(tH) will start at wsp(iptr+iexph-1).
    !
    !     ns        : (output) number of scaling-squaring used.
    !
    !     iflag     : (output) exit flag.
    !                       0 - no problem
    !                      <0 - problem
    !
    !----------------------------------------------------------------------|
    !     Roger B. Sidje (rbs@maths.uq.edu.au)
    !     EXPOKIT: Software Package for Computing Matrix Exponentials.
    !     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
    !----------------------------------------------------------------------|
    !
    integer :: i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
    real(dp) :: hnorm
    complex(dp) :: cp, cq, scale, scale2, ZERO, ONE

    parameter ( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX

    !---  check restrictions on input parameters ...
    mm = m*m
    iflag = 0
    if ( ldh.lt.m ) iflag = -1
    if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
    if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGPADM)'
    !
    !---  initialise pointers ...
    !
    icoef = 1
    ih2 = icoef + (ideg+1)
    ip  = ih2 + mm
    iq  = ip + mm
    ifree = iq + mm
    !
    !---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
    !     and set scale = t/2^ns ...
    !
    do i = 1,m
       wsp(i) = ZERO
    enddo
    do j = 1,m
       do i = 1,m
          wsp(i) = wsp(i) + ABS( H(i,j) )
       enddo
    enddo
    hnorm = 0.0d0
    do i = 1,m
       hnorm = MAX( hnorm,DBLE(wsp(i)) )
    enddo
    hnorm = ABS( t*hnorm )
    if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of ZGPADM.'
    ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
    scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
    scale2 = scale*scale
    !
    !---  compute Pade coefficients ...
    !
    i = ideg+1
    j = 2*ideg+1
    wsp(icoef) = ONE
    do k = 1,ideg
       wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
    enddo
    !
    !---  H2 = scale2*H*H ...
    !
    call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
    !
    !---  initialise p (numerator) and q (denominator) ...
    !
    cp = wsp(icoef+ideg-1)
    cq = wsp(icoef+ideg)
    do j = 1,m
       do i = 1,m
          wsp(ip + (j-1)*m + i-1) = ZERO
          wsp(iq + (j-1)*m + i-1) = ZERO
       enddo
       wsp(ip + (j-1)*(m+1)) = cp
       wsp(iq + (j-1)*(m+1)) = cq
    enddo
    !
    !---  Apply Horner rule ...
    !
    iodd = 1
    k = ideg - 1
100 continue
    iused = iodd*iq + (1-iodd)*ip
    call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,wsp(ih2),m, ZERO,wsp(ifree),m )
    do j = 1,m
       wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
    enddo
    ip = (1-iodd)*ifree + iodd*ip
    iq = iodd*ifree + (1-iodd)*iq
    ifree = iused
    iodd = 1-iodd
    k = k-1
    if ( k.gt.0 )  goto 100
    !
    !---  Obtain (+/-)(I + 2*(p\q)) ...
    !
    if ( iodd.ne.0 ) then
       call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m, H,ldh, ZERO,wsp(ifree),m )
       iq = ifree
    else
       call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,H,ldh, ZERO,wsp(ifree),m )
       ip = ifree
    endif
    call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
    call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
    if ( iflag.ne.0 ) stop 'Problem in ZGESV (within ZGPADM)'
    call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
    do j = 1,m
       wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
    enddo
    iput = ip
    if ( ns.eq.0 .and. iodd.ne.0 ) then
       call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
       goto 200
    endif
    !
    !--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
    !
    iodd = 1
    do k = 1,ns
       iget = iodd*ip + (1-iodd)*iq
       iput = (1-iodd)*ip + iodd*iq
       call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,ZERO,wsp(iput),m )
       iodd = 1-iodd
    enddo
200 continue
    iexph = iput

  END subroutine ZGPADM

  function r8_uniform_01 ( seed )

    !*****************************************************************************80
    !
    !! R8_UNIFORM_01 returns a unit pseudorandom R8.
    !
    !  Discussion:
    !
    !    This routine implements the recursion
    !
    !      seed = 16807 * seed mod ( 2^31 - 1 )
    !      r8_uniform_01 = seed / ( 2^31 - 1 )
    !
    !    The integer arithmetic never requires more than 32 bits,
    !    including a sign bit.
    !
    !    If the initial seed is 12345, then the first three computations are
    !
    !      Input     Output      R8_UNIFORM_01
    !      SEED      SEED
    !
    !         12345   207482415  0.096616
    !     207482415  1790989824  0.833995
    !    1790989824  2035175616  0.947702
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 May 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Second Edition,
    !    Springer, 1987,
    !    ISBN: 0387964673,
    !    LC: QA76.9.C65.B73.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, December 1986, pages 362-376.
    !
    !    Pierre L'Ecuyer,
    !    Random Number Generation,
    !    in Handbook of Simulation,
    !    edited by Jerry Banks,
    !    Wiley, 1998,
    !    ISBN: 0471134031,
    !    LC: T57.62.H37.
    !
    !    Peter Lewis, Allen Goodman, James Miller,
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, 1969, pages 136-143.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
    !    should NOT be 0.
    !    On output, SEED has been updated.
    !
    !    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
    !    strictly between 0 and 1.
    !
    implicit none

    integer, intent(inout) :: seed
    integer ::  k
    real(dp) :: r8_uniform_01


    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
       seed = seed + 2147483647
    end if
    !
    !  Although SEED can be represented exactly as a 32 bit integer,
    !  it generally cannot be represented exactly as a 32 bit real number!
    !
    r8_uniform_01 = seed * 4.656612875D-10


  end function r8_uniform_01

  function r8_normal_01 ( seed )

    !*****************************************************************************80
    !
    !! R8_NORMAL_01 returns a unit pseudonormal R8.
    !
    !  Discussion:
    !
    !    The standard normal probability distribution function (PDF) has
    !    mean 0 and standard deviation 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 August 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, integer  SEED, a seed for the random
    !    number generator.
    !
    !    Output, real  R8_NORMAL_01, a normally distributed
    !    random value.
    !
    implicit none

    real(dp):: r1, r2, r8_normal_01
    real(dp), parameter :: r8_pi = 3.141592653589793D+00

    integer, intent(inout) :: seed

    r1 = r8_uniform_01(seed)
    r2 = r8_uniform_01(seed)
    r8_normal_01 = sqrt( - 2.0_dp * log ( r1 ) ) * cos ( 2.0_dp * r8_pi * r2 )

    !print *, r8_normal_01 , 'bajs', seed

  end function r8_normal_01

  ! From numerical recipes
  SUBROUTINE gasdev_s(harvest)

    REAL(dp), INTENT(OUT) :: harvest
    ! Returns in harvest a normally distributed deviate with zero mean and unit variance, using the intrinsic random_number function as the source of uniform deviates.
    REAL(dp) :: rsq,v1,v2
    REAL(dp), SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.

    integer:: n, clock
    integer, allocatable :: seed(:)
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    call srand(45)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)



    if (gaus_stored) then
       !  We have an extra deviate handy so return it, and unset the flag.
       harvest=g
       gaus_stored=.false.
    else
       ! We don’t have an extra deviate handy, so pick two uniform numbers in the square extending from -1 to +1 in each direction,
       do
          call random_number(v1)
          call random_number(v2)
          !v1 = rand()
          !v2 = rand()
          ! print *, v1, v2
          v1=2.0_dp*v1-1.0_dp
          v2=2.0_dp*v2-1.0_dp
          rsq=v1**2+v2**2
          !see if they are in the unit circle,
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       ! otherwise try again.
       rsq=sqrt(-2.0_dp*log(rsq)/rsq)
       ! Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
       harvest=v1*rsq
       g=v2*rsq
       gaus_stored=.true.
       !Set flag.
    end if
  END SUBROUTINE gasdev_s


end module module1
