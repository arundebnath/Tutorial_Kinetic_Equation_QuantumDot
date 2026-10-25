module ode_solver
  use constants
  use bath_model
  use dynamics
  implicit none
  private

  public :: dopri5_step

contains

  ! Fixed-step Dormand-Prince 5(4) integrator step for HEOM
  subroutine dopri5_step(t, h, zz, bath, pulse)
    real(dp),          intent(in)    :: t, h
    complex(dp),       intent(inout) :: zz(n_states, n_states, n_bath)
    type(BathParams),  intent(in)    :: bath
    type(PulseParams), intent(in)    :: pulse

    complex(dp), dimension(n_states, n_states, n_bath) :: k1, k2, k3, k4, k5, k6, y_tmp

    ! Dormand-Prince Butcher tableau coefficients
    real(dp), parameter :: c2 = 0.2_dp, c3 = 0.3_dp, c4 = 0.8_dp, c5 = 8.0_dp/9.0_dp
    real(dp), parameter :: a21 = 0.2_dp
    real(dp), parameter :: a31 = 3.0_dp/40.0_dp,     a32 = 9.0_dp/40.0_dp
    real(dp), parameter :: a41 = 44.0_dp/45.0_dp,    a42 = -56.0_dp/15.0_dp,   a43 = 32.0_dp/9.0_dp
    real(dp), parameter :: a51 = 19372.0_dp/6561.0_dp, a52 = -25360.0_dp/2187.0_dp, &
                           a53 = 64448.0_dp/6561.0_dp, a54 = -212.0_dp/729.0_dp
    real(dp), parameter :: a61 = 9017.0_dp/3168.0_dp,  a62 = -355.0_dp/33.0_dp,     &
                           a63 = 46732.0_dp/5247.0_dp, a64 = 49.0_dp/176.0_dp,     a65 = -5103.0_dp/18656.0_dp
    real(dp), parameter :: b1  = 35.0_dp/384.0_dp,   b3  = 500.0_dp/1113.0_dp, &
                           b4  = 125.0_dp/192.0_dp,  b5  = -2187.0_dp/6784.0_dp, b6  = 11.0_dp/84.0_dp

    call heom_rhs(t, zz, k1, bath, pulse)

    y_tmp = zz + h * a21 * k1
    call heom_rhs(t + c2 * h, y_tmp, k2, bath, pulse)

    y_tmp = zz + h * (a31 * k1 + a32 * k2)
    call heom_rhs(t + c3 * h, y_tmp, k3, bath, pulse)

    y_tmp = zz + h * (a41 * k1 + a42 * k2 + a43 * k3)
    call heom_rhs(t + c4 * h, y_tmp, k4, bath, pulse)

    y_tmp = zz + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    call heom_rhs(t + c5 * h, y_tmp, k5, bath, pulse)

    y_tmp = zz + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    call heom_rhs(t + h, y_tmp, k6, bath, pulse)

    zz = zz + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
  end subroutine dopri5_step
end module ode_solver
