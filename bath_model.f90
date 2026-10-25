module bath_model
  use constants
  implicit none
  private

  public :: BathParams, init_bath, nbeta, jjcxx

  type :: BathParams
    real(dp)                :: beta   = 0.0_dp
    real(dp)                :: wc     = 0.0_dp
    real(dp)                :: lambda = 1.0_dp
    complex(dp)             :: aa(n_bath) = zero
    complex(dp)             :: gg(n_bath) = zero
    complex(dp)             :: xx(n_states, n_states) = zero
  end type BathParams

  ! GNU-fit rational decomposition poles
  real(dp), parameter :: c_fit(n_fit) = [ &
      0.144574_dp, 0.136538_dp, 0.0758868_dp, 0.191552_dp ]

  complex(dp), parameter :: z1_fit(n_fit) = [ &
      cmplx(0.967976_dp, 0.255507_dp, kind=dp), &
      cmplx(0.675019_dp, 0.186201_dp, kind=dp), &
      cmplx(0.835287_dp, 0.139491_dp, kind=dp), &
      cmplx(0.404769_dp, 0.324660_dp, kind=dp) ]

  complex(dp), parameter :: z2_fit(n_fit) = [ &
      cmplx(1.655210_dp, 0.317109_dp, kind=dp), &
      cmplx(1.920460_dp, 0.276022_dp, kind=dp), &
      cmplx(2.227760_dp, 0.268375_dp, kind=dp), &
      cmplx(1.348330_dp, 0.340821_dp, kind=dp) ]

contains

  pure function nbeta(beta, w) result(nb)
    real(dp),    intent(in) :: beta
    complex(dp), intent(in) :: w
    complex(dp)             :: nb

    if (real(w, dp) < 0.0_dp) then
      nb = 1.0_dp / (exp(beta * w) - 1.0_dp)
    else
      nb = exp(-beta * w) / (1.0_dp - exp(-beta * w))
    end if
  end function nbeta

  pure function jjcxx(x) result(res)
    complex(dp), intent(in) :: x
    complex(dp)             :: res
    integer(ip)             :: k
    complex(dp)             :: z1k, z2k, zc1k, zc2k

    res = zero
    do k = 1, n_fit
      z1k  = z1_fit(k)
      zc1k = conjg(z1_fit(k))
      z2k  = z2_fit(k)
      zc2k = conjg(z2_fit(k))

      res = res + c_fit(k) * x**3 / &
            ((x - z1k)  * (x + z1k)  * (x - zc1k) * (x + zc1k) * &
             (x - z2k)  * (x + z2k)  * (x - zc2k) * (x + zc2k))
    end do
  end function jjcxx

  subroutine init_bath(bath, wc, beta, lambda)
    type(BathParams), intent(out) :: bath
    real(dp),         intent(in)  :: wc, beta, lambda
    integer(ip)                   :: i, first_matsu, idx
    real(dp)                      :: nu_k
    complex(dp)                   :: zz1, zc1, zz2, zc2, x_matsu

    bath%wc     = wc
    bath%beta   = beta
    bath%lambda = lambda
    bath%aa     = zero
    bath%gg     = zero

    ! System coupling operator
    bath%xx = zero
    bath%xx(2, 2) = one
    bath%xx(4, 4) = 2.0_dp * one

    ! Bath decomposition via fitting poles
    do i = 1, n_fit
      zz1  = z1_fit(i)
      zc1  = conjg(zz1)
      zz2  = z2_fit(i)
      zc2  = conjg(zz2)
      idx  = (i - 1) * n_fit

      bath%gg(idx + 2) =  wc * zz1
      bath%gg(idx + 3) = -wc * zc1
      bath%gg(idx + 4) =  wc * zz2
      bath%gg(idx + 5) = -wc * zc2

      bath%aa(idx + 2) = c_fit(i) * nbeta(beta,  wc * zz1) * im * pi * wc**4 * zz1**2 / &
                         ((zz1**2 - zz2**2) * (zz1**2 - zc1**2) * (zz1**2 - zc2**2))
      bath%aa(idx + 3) = c_fit(i) * nbeta(beta, -wc * zc1) * im * pi * wc**4 * zc1**2 / &
                         ((zc1**2 - zz1**2) * (zc1**2 - zz2**2) * (zc1**2 - zc2**2))
      bath%aa(idx + 4) = c_fit(i) * nbeta(beta,  wc * zz2) * im * pi * wc**4 * zz2**2 / &
                         ((zz2**2 - zz1**2) * (zz2**2 - zc1**2) * (zz2**2 - zc2**2))
      bath%aa(idx + 5) = c_fit(i) * nbeta(beta, -wc * zc2) * im * pi * wc**4 * zc2**2 / &
                         ((zc2**2 - zz1**2) * (zc2**2 - zz2**2) * (zc2**2 - zc1**2))
    end do

    ! Scale with coupling strength
    bath%aa = lambda * bath%aa

    ! Matsubara frequencies
    first_matsu = 4_ip * n_fit + 2_ip
    do i = first_matsu, n_bath
      nu_k         = twopi / beta * real(i - first_matsu + 1, dp)
      x_matsu      = cmplx(0.0_dp, nu_k / wc, kind=dp)
      bath%gg(i)   = cmplx(0.0_dp, nu_k, kind=dp)
      bath%aa(i)   = cmplx(0.0_dp, twopi, kind=dp) * wc**3 / beta * jjcxx(x_matsu)
    end do
  end subroutine init_bath
end module bath_model
