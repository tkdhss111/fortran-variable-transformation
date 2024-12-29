program unit_test

  use variable_transformation_mo 

  implicit none

  real, allocatable :: z(:), x(:), z_yj(:)
  integer n

  !======================================
  ! Test: Yeo-Johnson Transformation 
  !
  n = 100

  allocate ( z(n), x(n), z_yj(n), source = -999.0 )

  ! Test data
  z = rnorm ( n ) ! Standard normal realizations
  x = exp(z)      ! Log normalized

  ! Yeo-Johnson transformation with MLE of lambda
  z_yj = yeo_johnson ( x )

  call write_csv ( z, x, z_yj, "data.csv" ) ! for R package comparison

 ! Yeo-Johnson transformation with user lambda
  z_yj = yeo_johnson_lambda ( x, lambda = -1.00446498 )

  call write_csv ( z, x, z_yj, "data_user_lambda.csv" ) ! for R package comparison

contains

  subroutine write_csv ( z, x, z_yj, file )
    real,         intent(in) :: z(:), x(:), z_yj(:)
    character(*), intent(in) :: file
    integer u, i
    open ( newunit = u, file = file )
    write ( u, * ) "std_normal, log_normal, yj_transformed" 
    do i = 1, n
      write ( u, '(*(g0, :, ",") )' )  z(i), x(i), z_yj(i)
!      print *, "Writing: z =", z(i), "x =", x(i), "z_yj =", z_yj(i)
    end do
    close ( u )
  end subroutine

  subroutine read_csv ( z, x, z_yj, file )
    real,         intent(out) :: z(:), x(:), z_yj(:)
    character(*), intent(in) :: file
    integer u, i
    open ( newunit = u, file = file )
    read ( u, '()' )
    do i = 1, n
      read ( u, * ) z(i), x(i), z_yj(i)
!      print *, "Reading: z =", z(i), "x =", x(i), "z_yj =", z_yj(i)
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
