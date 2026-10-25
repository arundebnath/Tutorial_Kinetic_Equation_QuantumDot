module constants
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  private

  public :: dp, ip, pi, twopi, im, zero, one
  public :: n_states, n_bath, bign, n_fit

  integer, parameter :: ip = int32
  integer, parameter :: dp = real64

  real(dp),    parameter :: pi    = 3.14159265358979323846264338327950288_dp
  real(dp),    parameter :: twopi = 2.0_dp * pi
  complex(dp), parameter :: im    = (0.0_dp, 1.0_dp)
  complex(dp), parameter :: zero  = (0.0_dp, 0.0_dp)
  complex(dp), parameter :: one   = (1.0_dp, 0.0_dp)

  ! Physical model dimensions (N=4 matches the Hamiltonian dimension in RWA)
  integer(ip), parameter :: n_states = 4_ip
  integer(ip), parameter :: n_bath   = 20_ip
  integer(ip), parameter :: bign     = 2_ip * n_states * n_states * n_bath
  integer(ip), parameter :: n_fit    = 4_ip
end module constants
