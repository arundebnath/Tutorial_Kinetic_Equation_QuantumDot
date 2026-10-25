program qdot
  use omp_lib
  use constants
  use bath_model
  use dynamics
  use ode_solver
  implicit none

  character(len=64) :: prefix
  real(dp)          :: w0_in, e0_in, tau0_in, phisec_in
  real(dp)          :: lambda_in, delta1_in, delta2_in, temp_k, wc_ev
  real(dp)          :: wc, beta, taup, alpha, ep
  integer(ip)       :: io_status

  integer(ip), parameter :: n_runs = 48_ip
  real(dp)               :: signal(n_runs, 3)
  real(dp)               :: efield_scan(n_runs)

  real(dp), parameter    :: fs_to_au = 1.0_dp / 0.024188843265857_dp
  real(dp), parameter    :: ev_to_au = 1.0_dp / 27.211386245988_dp
  real(dp), parameter    :: cm_to_au = 1.0_dp / 219474.63_dp
  real(dp), parameter    :: k_boltz  = 8.617333262145e-5_dp * ev_to_au

  ! Read configuration
  read(*, '(a)', iostat=io_status) prefix
  if (io_status /= 0) prefix = "sim"
  read(*, *) w0_in, e0_in, tau0_in, phisec_in
  read(*, *) lambda_in
  read(*, *) delta1_in
  read(*, *) delta2_in
  read(*, *) temp_k
  read(*, *) wc_ev

  ! Log parameters
  open(unit=44, file=trim(prefix) // ".param", status="replace", action="write")
  write(44, '(A,A)')     "Simulation Prefix: ", trim(prefix)
  write(44, '(A,ES16.8)') "Central freq w0 [cm^-1]:       ", w0_in
  write(44, '(A,ES16.8)') "Field amplitude e0:            ", e0_in
  write(44, '(A,ES16.8)') "Pulse length tau0 [fs]:        ", tau0_in
  write(44, '(A,ES16.8)') "Chirp parameter phisec [fs^2]: ", phisec_in
  write(44, '(A,ES16.8)') "Bath coupling lambda:          ", lambda_in
  write(44, '(A,ES16.8)') "Detuning delta1 [cm^-1]:       ", delta1_in
  write(44, '(A,ES16.8)') "Detuning delta2 [cm^-1]:       ", delta2_in
  write(44, '(A,ES16.8)') "Temperature [K]:               ", temp_k
  write(44, '(A,ES16.8)') "Cutoff energy [eV]:            ", wc_ev
  close(44)

  ! Convert units to a.u.
  wc   = wc_ev * ev_to_au
  beta = 1.0_dp / (k_boltz * temp_k)

  ! Pulse transform and chirp broadening
  taup  = tau0_in * sqrt(1.0_dp + (2.0_dp * phisec_in / tau0_in**2)**2)
  alpha = 2.0_dp * phisec_in / (tau0_in**4 + (2.0_dp * phisec_in)**2)
  ep    = e0_in * sqrt(tau0_in / taup)

  taup  = taup * fs_to_au
  alpha = alpha / (fs_to_au**2)

  ! Field parameter scan over parallel threads
  !$omp parallel default(none) &
  !$omp& shared(prefix, wc, beta, lambda_in, w0_in, delta1_in, delta2_in, taup, alpha, ep, signal, efield_scan)
  block
    type(BathParams)  :: local_bath
    type(PulseParams) :: local_pulse
    complex(dp)       :: zz(n_states, n_states, n_bath)
    real(dp)          :: dt, t_start, t_end, t_current
    integer(ip)       :: kk, step, n_steps

    call init_bath(local_bath, wc, beta, lambda_in)

    local_pulse%w0     = w0_in * cm_to_au
    local_pulse%delta1 = delta1_in * cm_to_au
    local_pulse%delta2 = delta2_in * cm_to_au
    local_pulse%taup   = taup
    local_pulse%alpha  = alpha

    dt      = 4.0_dp * fs_to_au
    n_steps = int(6.0_dp * taup / dt, kind=ip) + 1_ip
    t_start = -3.0_dp * taup

    !$omp do schedule(dynamic)
    do kk = 1, n_runs
      efield_scan(kk)    = ep * real(kk, dp) / real(n_runs, dp)
      local_pulse%ep     = efield_scan(kk)

      ! Initialize ground state
      zz = zero
      zz(1, 1, 1) = one

      t_current = t_start
      do step = 1, n_steps
        call dopri5_step(t_current, dt, zz, local_bath, local_pulse)
        t_current = t_current + dt
      end do

      ! Record populations (States 1, 2, and 4)
      signal(kk, 1) = real(zz(1, 1, 1), dp)
      signal(kk, 2) = real(zz(2, 2, 1), dp)
      signal(kk, 3) = real(zz(4, 4, 1), dp)
    end do
    !$omp end do
  end block
  !$omp end parallel

  ! Write results
  open(unit=45, file=trim(prefix) // ".result", status="replace", action="write")
  do io_status = 1, n_runs
    write(45, '(4ES24.14)') efield_scan(io_status) * sqrt(taup), &
                           signal(io_status, 1),                 &
                           signal(io_status, 2),                 &
                           signal(io_status, 3)
  end do
  close(45)

end program qdot
