# NSE solver

A Fortran solver for nuclear statistical equilibrium (NSE) in hot astrophysical matter.

The code computes an NSE composition at fixed density, temperature, and electron fraction by solving for the neutron and proton chemical potentials subject to

\[
\sum_i X_i = 1,
\qquad
\sum_i \frac{Z_i}{A_i}X_i = Y_e.
\]

The current implementation stores nuclear-set-dependent quantities in `type(nse_network_t)`, so different NSE species sets can be initialized and solved independently in the same program. Currently supported sets include a large WinVNE-based set, an `aprox21`-like reduced set, and a small four-species test set.

> **Status:** research code under active development. Some auxiliary drivers in `src/` are development tools and may not always track the latest module interface.

## Features

- NSE solution at fixed `rho`, `T`, and `Ye`
- damped Newton iteration in neutron/proton chemical potentials
- logarithmic abundance evaluation for numerical stability
- WinVNE nuclear masses, spins, and partition functions
- optional Rauscher partition functions with WinVNE fallback
- Hempel-type Coulomb correction
- explicit `nse_network_t` objects for multiple nuclear sets
- large-NSE and reduced `aprox21` NSE calculations
- diagnostics for `Abar`, heavy-nucleus moments, light-particle abundances, mass excess, and Coulomb energy

## Repository layout

```text
.
├── data/
│   └── winvn_v2.0.dat
├── src/
│   ├── module_nse.f90
│   ├── module_ptf_reaclib.f90
│   ├── module_ptf_rauscher.f90
│   ├── const_mod.f90
│   ├── nse_single.f90
│   ├── nse_aprox21.f90
│   └── ...
├── utils/
│   ├── pl_nse.py
│   └── pl_conv.py
└── makefile
```

The default `makefile` currently builds `src/nse_single.f90`.

## Build

The supplied makefile uses `gfortran`.

```bash
git clone https://github.com/sfujibayashi/nse.git
cd nse
make
```

The executable is created as

```text
bin/a.out
```

No external numerical library is required for the default build.

## Nuclear data

### WinVNE data

The repository includes

```text
data/winvn_v2.0.dat
```

which provides nuclear masses, ground-state spins, and tabulated partition functions.

### Rauscher partition functions

The code can optionally use a Rauscher partition-function table. This file is not currently bundled in the repository and must be supplied separately.

`module_ptf_rauscher.f90` expects a 23-line header followed by one record per nucleus of the form

```text
Z  A  J0  PF(T9_1) ... PF(T9_72)
```

where the temperature grid is defined in the module.

At present the single-point driver initializes the Rauscher table even if the Rauscher partition functions are disabled, so a valid file path is still required.

## Quick start: single-point NSE

Create a parameter file, for example `single.para`:

```text
# WinVNE data
data/winvn_v2.0.dat
# Use Rauscher partition functions?
T
# Rauscher partition-function table
/path/to/rauscher_partition_functions.dat
# Output file
nse.dat
# rho [g cm^-3], T [K], Ye
1.0e8 9.0e9 0.50
```

Run

```bash
./bin/a.out single.para
```

The output file contains

```text
index  name  A  Z  N  Xi  Yi  gi
```

where

- `Xi` is the mass fraction,
- `Yi = Xi/Ai` is the abundance per baryon,
- `gi` is the nuclear statistical weight including the partition function.

## General Fortran usage

The main API is provided by

```fortran
use module_nse
```

A typical large-NSE calculation is

```fortran
program example_nse
  use module_nse
  use module_ptf_reaclib
  use module_ptf_rauscher
  implicit none

  type(nse_network_t) :: net
  real(8), allocatable :: x(:)

  real(8) :: rho, temp, ye
  integer, parameter :: itrlim = 300
  real(8), parameter :: tol = 1.d-10
  logical :: nsefail

  call init_ptf_reaclib("data/winvn_v2.0.dat")
  call init_ptf_rauscher("/path/to/rauscher_partition_functions.dat")

  call nse_init_reaclib(net, .true.)
  allocate(x(net%n_spec))

  rho  = 1.d8
  temp = 9.d9
  ye   = 0.30d0

  call calc_nse(net, rho, temp, ye, itrlim, tol, &
       x, nsefail, .false.)

  ! Robust fallback for difficult points
  if (nsefail) then
     call calc_nse(net, rho, temp, ye, itrlim, tol, &
          x, nsefail, .true.)
  endif

  if (nsefail) error stop "NSE did not converge"

  call output_nse_full(net, rho, temp, ye, x, "nse.dat")

end program example_nse
```

For robust calculations, the recommended strategy is to first solve with `use_TNAguess=.false.` and retry with `use_TNAguess=.true.` only if `nsefail=.true.`. This was found to be more reliable than always using the two-nuclei approximation as the initial guess.

## Nuclear sets

### Large NSE set

Initialize with

```fortran
call nse_init_reaclib(net, use_rauscher_ptf)
```

Despite the historical `reaclib` naming in the source, nuclear properties are read from the WinVNE-format input table.

With the current large-set construction:

- nuclei matched to the Rauscher table use the Rauscher partition function when enabled;
- unmatched nuclei with `Z < 87` are retained using the WinVNE partition-function data;
- unmatched nuclei with `Z >= 87` are excluded.

The exact number of retained species therefore depends on the input data files.

### aprox21 NSE set

Initialize with

```fortran
call nse_init_aprox21(net, use_rauscher_ptf)
```

The NSE set contains 20 physically distinct nuclei:

```text
n, p, He3, He4,
C12, N14, O16, Ne20, Mg24, Si28, S32, Ar36,
Ca40, Ti44, Cr48, Cr56, Fe52, Fe54, Fe56, Ni56
```

The Microphysics `aprox21` reaction network contains both `H1` and `p` as evolved network variables. They correspond to the same physical nucleus `(A,Z)=(1,1)` and therefore must not be counted twice in an NSE partition sum. The NSE implementation consequently contains a single physical proton species.

The large and reduced sets can coexist:

```fortran
type(nse_network_t) :: net_large, net_aprox21

call nse_init_reaclib(net_large, .true.)
call nse_init_aprox21(net_aprox21, .true.)
```

and can then be solved independently at the same `(rho,T,Ye)` point.

### Four-species test set

A small built-in test set is available through

```fortran
call nse_init_four(net)
```

and contains

```text
n, p, He4, Ni56
```

## Main solver interface

The central routine is

```fortran
call calc_nse(net, rho, temp, ye, itrlim, tol, &
     xnse, nsefail, use_TNAguess)
```

with

```text
rho             density [g cm^-3]
temp            temperature [K]
ye              electron fraction
itrlim          maximum number of Newton iterations
tol             convergence tolerance
xnse(:)         returned mass fractions Xi
nsefail         convergence flag
use_TNAguess    use the two-nuclei approximation for the initial guess
```

Additional optional arguments provide convergence histories, iteration counts, residuals, and final neutron/proton chemical-potential variables.

Always check `nsefail` before using the returned composition.

## Diagnostics

Basic composition moments can be obtained with

```fortran
call statistic(net, x, mexc_ave, z_heavy, a_heavy, &
     y_heavy, ytot, xsum, yesum)
```

For comparisons with CompOSE-like quantities,

```fortran
type(stat_t) :: stat
call statistic_compose(net, rho, x, stat)
```

returns quantities including

```text
stat%yn, stat%yp
stat%yh2, stat%yh3, stat%yhe3, stat%yhe4
stat%a_n, stat%z_n, stat%y_n
stat%abar
stat%mexc
stat%ecoul
```

## Mass-excess convention

Care is required when comparing with external EOS or nuclear-data tables.

The WinVNE input contains atomic mass excesses. Internally they are converted to bare-nuclear mass excesses as

\[
m_{\rm exc,nuc}
=
m_{\rm exc,atomic} - Zm_ec^2.
\]

`statistic_compose()` adds the rest mass of the balancing electrons,

\[
Y_e m_e c^2,
\]

when reporting the average mass excess per baryon.

Mass-excess quantities in `module_nse` are in MeV unless otherwise noted.

## Partition functions

For a WinVNE/Rauscher-based network, the statistical weight is

\[
g_i(T)=(2J_i+1)G_i(T),
\]

where `J_i` is the ground-state spin and `G_i(T)` is the nuclear partition function.

When `net%use_rauscher_ptf=.true.` and a matching Rauscher nucleus is available, the Rauscher spin and partition function are used. Otherwise the WinVNE values are used as a fallback.

## Coulomb correction

For the WinVNE-based NSE sets, `calc_nse()` includes the Coulomb correction implemented by

```fortran
calc_coulomb_HS()
fcoulomb_HS()
```

The default nuclear saturation density is

```text
n0 = 0.1583 fm^-3
```

stored in `net%n0_fm`.

## Numerical method

The abundance of each species is evaluated in logarithmic form,

\[
\ln X_i = \ln G_i + Z_i\eta_p + N_i\eta_n,
\]

with the density, temperature, mass, and Coulomb contributions included in the species-dependent prefactor.

The two nonlinear constraints

\[
\sum_i X_i = 1,
\qquad
\sum_i \frac{Z_i}{A_i}X_i=Y_e
\]

are solved with a damped Newton iteration. Exponentials are evaluated relative to the largest logarithmic abundance to reduce overflow and underflow.

For sparse nuclear sets, the Jacobian can become poorly conditioned when a trial composition is dominated by a single nucleus. The two-nuclei initial-guess fallback is useful in this regime.

## Validation

Useful checks for a converged solution are

\[
\sum_i X_i \simeq 1,
\qquad
\sum_i \frac{Z_i}{A_i}X_i \simeq Y_e.
\]

Additional validation quantities include

```text
Abar
mexc/b
Coulomb energy/b
free n/p/alpha abundances
heavy-nucleus moments
dominant nuclei
```

The refactor from module-global nuclear-set data to `nse_network_t` was regression-tested against the previous large-NSE implementation, reproducing the tested output exactly.

## Development notes

The repository also contains drivers used to compare

- the large NSE set with the reduced `aprox21` NSE set,
- NSE quantities with CompOSE/tabulated-EOS quantities,
- NSE compositions with Helmholtz-EOS thermodynamics.

These are currently research/development utilities rather than a stable command-line interface.
