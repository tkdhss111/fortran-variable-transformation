module variable_transformation_mo

  implicit none

  private
  public :: yeo_johnson
  public :: yeo_johnson_lambda

contains

  pure elemental function yeo_johnson_lambda ( x, lambda ) result ( z )

    real, intent(in) :: x
    real, intent(in) :: lambda
    real             :: z
    
    if ( x >= 0.0 ) then ! For positive values
      if ( is_eq( lambda, 0.0 ) ) then
        z = log(x + 1.0)
      else
        z = ( (x + 1)**lambda - 1.0 ) / lambda
      end if
    else ! For negative values
      if ( is_eq( lambda, 2.0 ) ) then
        z = -log(1.0 - x)
      else
        z = -( (1.0 - x)**(2.0 - lambda) - 1.0) / (2.0 - lambda)
      end if
    end if

  end function yeo_johnson_lambda

  pure function yeo_johnson ( x ) result ( z )

    real, intent(in) :: x(:)
    real             :: z(size(x))
    real             :: lambda_opt
    real             :: const

    const = sum( sign(1.0, x) * log(abs(x) + 1.0) )

    lambda_opt = golden_section_search ( low = -3.0, upp = 3.0, tol = 0.01, fun = mle )
    !print *, "Optimal lambda:", lambda_opt

    z = yeo_johnson_lambda ( x, lambda_opt )

  contains

    pure real function mle ( lambda_ )
      real, intent(in) :: lambda_
      mle = - loglike ( x, lambda_, const ) ! Maximization
    end function

  end function yeo_johnson

  pure real function loglike ( x, lambda, const )

    real, intent(in) :: x(:)
    real, intent(in) :: lambda
    real, intent(in) :: const
    real             :: mu, z(size(x))
    real             :: mse
    integer n

    n = size(x)
    z = yeo_johnson_lambda ( x, lambda )
    mu = sum(z) / real(n)
    mse = sum( ( z - mu )**2 ) / real(n)
    loglike = -0.5 * n * log(mse) + (lambda - 1.0) * const

  end function

  pure function golden_section_search ( low, upp, tol, fun ) result ( x_opt )

      real, intent(in) :: low, upp, tol
      interface
        pure real function fun(x)
          real, intent(in) :: x
        end function
      end interface

      real :: low_, upp_
      real :: phi, resphi
      real :: a, b
      real :: f_a, f_b
      real :: x_opt

      phi = (1.0 + sqrt(5.0)) / 2.0 ! Golden ratio
      resphi = 2.0 - phi            ! 1 / phi

      low_ = low
      upp_ = upp

      a = low_ + resphi*(upp_ - low_)
      b = upp_ - resphi*(upp_ - low_)
      f_a = fun(a)
      f_b = fun(b)

      do while ( abs(upp_ - low_) > tol )
        if (f_a < f_b) then
          upp_ = b
          b = a
          f_b = f_a
          a = low_ + resphi*(upp_ - low_)
          f_a = fun(a)
        else
          low_ = a
          a = b
          f_a = f_b
          b = upp_ - resphi*(upp_ - low_)
          f_b = fun(b)
        end if
      end do

      ! Objective function f has a minimum at x_opt
      x_opt = (low_ + upp_) / 2.0

  end function golden_section_search

  pure elemental logical function is_eq ( x, ref )
    real, intent(in) :: x
    real, intent(in) :: ref
    is_eq = abs(x - ref) < epsilon(ref)
  end function is_eq

end module
