program unit_test

  use variable_transformation_mo 

  implicit none

  real, allocatable :: z0(:), x(:), z(:)
  integer n
  
  !=======================================
  ! Test data
  !
  n = 100

  allocate ( z0(n), x(n), z(n), source = -999.0 )

  z0 = rnorm ( n ) ! Standard normal realizations
  x = exp(z0)      ! Log normalized

  !======================================
  ! Test: Yeo-Johnson Transformation 
  !

  ! Yeo-Johnson transformation with MLE of lambda
  z = yeo_johnson ( x )

  call write_csv ( z0, x, z, "yj.csv" )

  ! Yeo-Johnson transformation with user lambda (same as above)
  z = yeo_johnson_lambda ( x, lambda = -1.00446498 )

  call write_csv ( z0, x, z, "yj_lambda.csv" )

  ! Yeo-Johnson transformation with user lambda (obtained from R)
  z = yeo_johnson_lambda ( x, lambda = -1.006137 )

  call write_csv ( z0, x, z, "yj_lambda_r.csv" )

  !======================================
  ! Test: Yeo-Johnson Transformation 
  !

  ! Box-Cox transformation with user lambda (obtained from R)
  z = box_cox_lambda ( x, lambda = -0.1727399 )

  call write_csv ( z0, x, z, "bc_lambda_r.csv" )

contains

  subroutine write_csv ( z0, x, z, file )
    real,         intent(in) :: z0(:), x(:), z(:)
    character(*), intent(in) :: file
    integer u, i
    open ( newunit = u, file = file )
    write ( u, * ) "std_normal, log_normal, transformed" 
    do i = 1, n
      write ( u, '(*(g0, :, ",") )' )  z0(i), x(i), z(i)
!      print *, "Writing: z0 =", z0(i), "x =", x(i), "z =", z(i)
    end do
    close ( u )
  end subroutine

  subroutine read_csv ( z0, x, z, file )
    real,         intent(out) :: z0(:), x(:), z(:)
    character(*), intent(in) :: file
    integer u, i
    open ( newunit = u, file = file )
    read ( u, '()' )
    do i = 1, n
      read ( u, * ) z0(i), x(i), z(i)
!      print *, "Reading: z0 =", z0(i), "x =", x(i), "z =", z(i)
    end do
    close ( u )
  end subroutine

  ! Box-Muller transformation to generate standard normal variables
  function rnorm ( n ) result ( z )
    integer, intent(in) :: n
    real u1, u2, z(n) 
    integer, allocatable :: seeds(:)
    integer i, size
    call random_seed(size = size)
    allocate( seeds(size), source = 777 )
    call random_seed(put = seeds)
    do i = 1, n
      call random_number(u1)
      call random_number(u2)
      z(i) = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2)
    end do
  end function

end program
