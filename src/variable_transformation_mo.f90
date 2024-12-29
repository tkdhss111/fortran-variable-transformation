module variable_transformation_mo

  implicit none

  private
  public :: yeo_johnson
  public :: yeo_johnson_lambda
  public :: box_cox_lambda

contains

  !====================================
  ! Yeo-Johnson Transformation
  !

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

  function yeo_johnson ( x ) result ( z )

    real, intent(in) :: x(:)
    real             :: z(size(x))
    real             :: lambda_opt
    real             :: const

    const = sum( sign(1.0, x) * log(abs(x) + 1.0) )

    lambda_opt = golden_section_search ( fun = loglike, low = -3.0, upp = 3.0, tol = 0.01, maximize = .true. )
    print *, "Optimal lambda:", lambda_opt

    z = yeo_johnson_lambda ( x, lambda_opt )

  contains

    pure real function loglike ( lambda_ )

      real, intent(in) :: lambda_
      real             :: mu, z(size(x))
      real             :: mse
      integer n

      n = size(x)
      z = yeo_johnson_lambda ( x, lambda_ )
      mu = sum(z) / real(n)
      mse = sum( ( z - mu )**2 ) / real(n)
      loglike = -0.5 * n * log(mse) + (lambda_ - 1.0) * const

    end function

  end function yeo_johnson

  !====================================
  ! Box-Cox Transformation
  !

  pure function box_cox_lambda ( x, lambda ) result ( z )

    real, intent(in) :: x(:)
    real, intent(in) :: lambda
    real, parameter  :: EPS = 0.01
    real             :: x_min
    real             :: xp(size(x)) ! Positive value
    real             :: z(size(x))

    x_min = minval(x) 

    if ( x_min <= 0 ) then
      xp = x - x_min + EPS ! Shift data to the positive region
    else
      xp = x
    end if
    
    if ( abs(lambda) < EPS ) then
      z = log(xp)
    else
      z = ( sign(1.0, xp) * abs(xp) ** lambda - 1.0 ) / lambda
    end if

  end function box_cox_lambda

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

  pure elemental logical function is_eq ( x, ref )
    real, intent(in) :: x
    real, intent(in) :: ref
    is_eq = abs(x - ref) < epsilon(ref)
  end function is_eq

end module
