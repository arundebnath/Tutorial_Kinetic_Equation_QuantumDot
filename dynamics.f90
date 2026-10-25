module dynamics
  use constants
  use bath_model
  implicit none
  private

  public :: PulseParams, get_ham, heom_rhs

  type :: PulseParams
    real(dp) :: ep     = 0.0_dp
    real(dp) :: w0     = 0.0_dp
    real(dp) :: taup   = 0.0_dp
    real(dp) :: alpha  = 0.0_dp
    real(dp) :: delta1 = 0.0_dp
    real(dp) :: delta2 = 0.0_dp
  end type PulseParams

contains

  pure subroutine get_ham(t, p, hh)
    real(dp),           intent(in)  :: t
    type(PulseParams),  intent(in)  :: p
    complex(dp),        intent(out) :: hh(n_states, n_states)
    real(dp)                        :: envelope, detuning1, detuning2

    hh = zero
    envelope  = 0.5_dp * p%ep * exp(-(t / p%taup)**2)
    detuning1 = p%w0 + 2.0_dp * p%alpha * t - p%delta1
    detuning2 = p%w0 + 2.0_dp * p%alpha * t - p%delta2

    ! Rotating-wave approximation (RWA)
    hh(1, 1) = cmplx(detuning1,  0.0_dp, kind=dp)
    hh(4, 4) = cmplx(-detuning2, 0.0_dp, kind=dp)

    hh(1, 2) = cmplx(envelope, 0.0_dp, kind=dp)
    hh(2, 1) = cmplx(envelope, 0.0_dp, kind=dp)
    hh(2, 4) = cmplx(envelope, 0.0_dp, kind=dp)
    hh(4, 2) = cmplx(envelope, 0.0_dp, kind=dp)
  end subroutine get_ham

  ! HEOM Right-Hand Side: computes d/dt zz(n_states, n_states, n_bath)
  subroutine heom_rhs(t, zz, zp, bath, pulse)
    real(dp),           intent(in)  :: t
    complex(dp),        intent(in)  :: zz(n_states, n_states, n_bath)
    complex(dp),        intent(out) :: zp(n_states, n_states, n_bath)
    type(BathParams),   intent(in)  :: bath
    type(PulseParams),  intent(in)  :: pulse

    complex(dp) :: hh(n_states, n_states)
    complex(dp) :: zh(n_states, n_states)
    complex(dp) :: commutator, drive
    integer(ip) :: o, m, k, p

    call get_ham(t, pulse, hh)

    ! Auxiliary density matrices (k >= 2)
    do k = 2, n_bath
      do m = 1, n_states
        do o = 1, n_states
          commutator = zero
          drive      = zero
          do p = 1, n_states
            commutator = commutator + hh(o, p) * zz(p, m, k) - zz(o, p, k) * hh(p, m)
            drive      = drive      + bath%xx(o, p) * zz(p, m, 1)
          end do
          zp(o, m, k) = im * bath%gg(k) * zz(o, m, k) - im * commutator + im * bath%aa(k) * drive
        end do
      end do
    end do

    ! Reduced system density matrix (k = 1)
    zp(:, :, 1) = zero
    do m = 1, n_states
      do o = 1, n_states
        commutator = zero
        do p = 1, n_states
          commutator = commutator + hh(o, p) * zz(p, m, 1) - zz(o, p, 1) * hh(p, m)
        end do
        zp(o, m, 1) = -im * commutator
      end do
    end do

    ! Bath dissipation coupling back to system
    do k = 2, n_bath
      do m = 1, n_states
        do o = 1, n_states
          zh(o, m) = zz(o, m, k) + conjg(zz(m, o, k))
        end do
      end do

      do m = 1, n_states
        do o = 1, n_states
          drive = zero
          do p = 1, n_states
            drive = drive + bath%xx(o, p) * zh(p, m) - zh(o, p) * bath%xx(p, m)
          end do
          zp(o, m, 1) = zp(o, m, 1) + im * drive
        end do
      end do
    end do
  end subroutine heom_rhs
end module dynamics
