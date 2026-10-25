# QuantumDotHEOM

Basic tutorial for open quantum system dynamics (quantum dots driven by chirped ultrafast laser pulses). Dissipative exciton-phonon interactions are treated with a rational pole-decomposition of the spectral density.

---

## Directory Structure

```text
.
├── CMakeLists.txt     # CMake build configuration
├── constants.f90      # Physical constants, precision (real64), array bounds
├── bath_model.f90     # Bath parametrization, pole decomposition, Matsubara modes
├── dynamics.f90       # Laser Hamiltonian & HEOM derivative (RHS) evaluation
├── ode_solver.f90     # High-order Dormand-Prince (DOPRI5) explicit integrator
└── main.f90           # Parameter parsing, unit conversions, OpenMP parameter scan

```

---

## Prerequisites

* Modern Fortran compiler supporting Fortran 2008+ (e.g., `gfortran` $\ge$ 9.0, `ifx`, or `flang`)
* OpenMP $\ge$ 4.5
* CMake $\ge$ 3.18 (or GNU Make)

---

## Compilation

### Using CMake

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j$(nproc)

```

### Manual Compilation (`gfortran`)

```bash
gfortran -O3 -fopenmp -std=f2008 \
  constants.f90 bath_model.f90 dynamics.f90 ode_solver.f90 main.f90 \
  -o qdot

```

---

## Usage

The program expects inputs via standard input (`stdin`) in the following order:

```text
<prefix>
<w0 [cm^-1]>  <e0 [a.u.]>  <tau0 [fs]>  <phisec [fs^2]>
<lambda>
<delta1 [cm^-1]>
<delta2 [cm^-1]>
<temperature [K]>
<wc [eV]>

```

### Example

Create an `input.dat` file:

```text
dot_sim
1000.0   0.05    100.0   500.0
0.1
1000.0
980.0
4.2
0.002

```

Run the executable:

```bash
./qdot < input.dat

```

---

## Output

* `<prefix>.param`: Echo of converted parameters and run metadata.
* `<prefix>.result`: Four-column table generated across the parallelized field amplitude scan:
* **Col 1:** Effective pulse area parameter ($\mathcal{E} \sqrt{\tau_p}$)
* **Col 2:** Population of ground state $\rho_{11}$
* **Col 3:** Population of exciton state $\rho_{22}$
* **Col 4:** Population of biexciton state $\rho_{44}$



---

## License

MIT License. See [LICENSE](https://www.google.com/search?q=LICENSE) for details.
