module variable_transformation_mo

  implicit none

  private
  public :: scale
  public :: rnorm
  public :: yeo_johnson_lambda
  public :: opt_yeo_johnson_lambda
  public :: yeo_johnson
  public :: box_cox_lambda
  public :: opt_box_cox_lambda
  public :: box_cox

contains

  !====================================
  ! Scaling Function
  !
  pure function scale ( x, mean, sd ) result ( z )

    real, intent(in) :: x(:)
    real             :: z(size(x))
    real, optional, intent(in) :: mean, sd
    real                       :: mean_, sd_
    real                       :: xbar, sd_x
    integer n

    if ( present( mean ) ) then
      mean_ = mean
    else
      mean_ = 0.0
    end if

    if ( present( sd ) ) then
      sd_ = sd
    else
      sd_ = 1.0
    end if

    n    = size(x)
    xbar = sum(x) / real(n)
    sd_x = sqrt( sum( (x - xbar)**2 )/(n-1) )

    !print *, 'xbar: ', xbar
    !print *, 'sd_x: ', sd_x

    z = ( x - xbar ) / sd_x * sd_ + mean_

  end function

  !====================================
  ! Yeo-Johnson Transformation
  !
  pure function yeo_johnson ( x ) result ( z )

    real, intent(in) :: x(:)
    real             :: z(size(x))
    real             :: lambda_opt

    lambda_opt = opt_yeo_johnson_lambda ( x )

    z = yeo_johnson_lambda ( x, lambda_opt )

  end function yeo_johnson

  pure elemental function yeo_johnson_lambda ( x, lambda ) result ( z )

    real, intent(in) :: x
    real, intent(in) :: lambda
    real             :: z

    if ( x >= 0.0 ) then ! For positive values
      if ( is_eq( lambda, 0.0 ) ) then
        z = log(x + 1.0)
      else
        z = ( (x + 1.0)**lambda - 1.0 ) / lambda
      end if
    else ! For negative values
      if ( is_eq( lambda, 2.0 ) ) then
        z = -log(1.0 - x)
      else
        z = -( (1.0 - x)**(2.0 - lambda) - 1.0) / (2.0 - lambda)
      end if
    end if

  end function yeo_johnson_lambda

  pure function opt_yeo_johnson_lambda ( x ) result ( lambda_opt )

    real, intent(in) :: x(:)
    real             :: lambda_opt
    real             :: const

    const = sum( sign(1.0, x) * log(abs(x) + 1.0) )

    lambda_opt = golden_section_search ( fun = loglike, &
                                         low = -3.0,    &
                                         upp =  3.0,    &
                                         tol = 0.001,   &
                                         maximize = .true. )

    !print *, "Optimal lambda:", lambda_opt

  contains

    pure real function loglike ( lambda_ )

      real, intent(in) :: lambda_
      real             :: z(size(x))
      real             :: zbar, var_z
      integer n

      n = size(x)
      z = yeo_johnson_lambda( x, lambda_ )
      zbar  = sum(z) / real(n)
      var_z = sqrt( sum( (z - zbar)**2 )/(n-1) )
      loglike = -0.5 * n * log(var_z) + (lambda_ - 1.0) * const

    end function

  end function opt_yeo_johnson_lambda 

  !====================================
  ! Box-Cox Transformation
  !

  pure function box_cox ( x, eps ) result ( z )

    real,           intent(in) :: x(:)
    real, optional, intent(in) :: eps
    real                       :: lambda_opt
    real                       :: z(size(x))

    lambda_opt = opt_box_cox_lambda( x, eps )

    z = box_cox_lambda( x, lambda_opt, eps )

  end function box_cox 

  pure function box_cox_lambda ( x, lambda, eps ) result ( z )

    real, intent(in) :: x(:)
    real, intent(in) :: lambda
    real, optional, intent(in) :: eps
    real             :: x_min
    real             :: xp(size(x)) ! Positive value
    real             :: z(size(x))
    real             :: eps_

    if ( present( eps ) ) then
      eps_ = eps
    else
      eps_ = 0.01
    end if

    x_min = minval( x )

    if ( x_min <= 0 ) then
      xp = x - x_min + 1.0 ! Shift data to the positive region
    else
      xp = x
    end if

    if ( abs(lambda) < eps_ ) then
      z = log(xp)
    else
      z = ( sign(1.0, xp) * abs(xp) ** lambda - 1.0 ) / lambda
    end if

  end function box_cox_lambda

  pure function opt_box_cox_lambda ( x, eps ) result ( lambda_opt )

    real,           intent(in) :: x(:)
    real, optional, intent(in) :: eps
    real                       :: xbar, ln_x(size(x))
    real                       :: lambda_opt
    real                       :: eps_
    integer n

    if ( present( eps ) ) then
      eps_ = eps
    else
      eps_ = 0.01
    end if

    n = size(x)
    ln_x = log(x)
    xbar = exp( sum( ln_x ) / real(n) ) ! log-sum mean

    lambda_opt = golden_section_search ( fun = loglike, &
                                         low = -3.0,    &
                                         upp =  3.0,    &
                                         tol = 0.001,   &
                                         maximize = .true. )
    !print *, "Optimal lambda:", lambda_opt

  contains

    pure real function loglike ( lambda_ )

      real, intent(in) :: lambda_
      real             :: z(size(x))
      real             :: eta, zbar, var_z

      eta = xbar**(lambda_ - 1.0)

      if ( abs(lambda_) < eps_ ) then
        z = ln_x / eta
      else
        z = (x**lambda_ - 1.0) / (lambda_ * eta)
      end if

      zbar  = sum(z) / real(n)
      var_z = sqrt( sum( (z - zbar)**2 )/(n-1) )

      loglike = -0.5 * n * log(var_z)

    end function

  end function opt_box_cox_lambda 

  !====================================
  ! Optimizer
  !

  pure function golden_section_search ( fun, low, upp, tol, maximize ) result ( x_opt )

      interface
        pure real function fun(x)
          real, intent(in) :: x
        end function
      end interface
      real,              intent(in) :: low, upp, tol
      logical, optional, intent(in) :: maximize

      real :: low_, upp_
      real :: phi, resphi
      real :: a, b
      real :: f_a, f_b
      real :: x_opt
      real :: s

      s = 1.0 ! Default: minimization

      if ( present( maximize ) ) then
        if ( maximize ) then
          s = -1.0 
        end if
      end if

      phi = (1.0 + sqrt(5.0)) / 2.0 ! Golden ratio
      resphi = 2.0 - phi            ! 1 / phi

      low_ = low
      upp_ = upp

      a = low_ + resphi*(upp_ - low_)
      b = upp_ - resphi*(upp_ - low_)
      f_a = s*fun(a)
      f_b = s*fun(b)

      do while ( abs(upp_ - low_) > tol )
        if (f_a < f_b) then
          upp_ = b
          b = a
          f_b = f_a
          a = low_ + resphi*(upp_ - low_)
          f_a = s*fun(a)
        else
          low_ = a
          a = b
          f_a = f_b
          b = upp_ - resphi*(upp_ - low_)
          f_b = s*fun(b)
        end if
      end do

      ! Objective function f has a minimum (default) at x_opt
      x_opt = (low_ + upp_) / 2.0

  end function golden_section_search

  !====================================
  ! Utilities
  !

  pure elemental logical function is_eq ( x, ref )
    real, intent(in) :: x
    real, intent(in) :: ref
    is_eq = abs(x - ref) < epsilon(ref)
  end function is_eq

  ! Box-Muller
  function rnorm ( n, mean, sd ) result ( z )
    integer, intent(in) :: n
    real,    intent(in) :: mean, sd
    real    :: z(n), u1((n+1)/2), u2((n+1)/2)
    integer :: m, i
    m = (n + 1) / 2
    call random_number( u1 )
    call random_number( u2 )
    do i = 1, m
      z(2*i-1) = sqrt( -2.0 * log(u1(i)) ) * cos(2.0 * acos(-1.0) * u2(i))
      if (2*i <= n) then
        z(2*i) = sqrt( -2.0 * log(u1(i)) ) * sin(2.0 * acos(-1.0) * u2(i))
      end if
    end do
    z = mean + sd * z
  end function rnorm

end module
