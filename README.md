# NSE solver

A Fortran code for nuclear statistical equilibrium (NSE) calculations in hot astrophysical matter.

The code solves for the equilibrium nuclear composition at fixed density, temperature, and electron fraction,

$$
(\rho,T,Y_e),
$$

subject to baryon-number and charge conservation,

$$
\sum_i X_i = 1,
\qquad
\sum_i \frac{Z_i}{A_i}X_i = Y_e.
$$

The repository is currently being developed not only as an NSE solver, but also as a tool for constructing thermodynamically consistent NSE + Helmholtz-EOS tables in a CompOSE-compatible HDF5 format.

> **Status:** research code under active development. Interfaces and development drivers may change.

## Main features

* NSE composition at fixed \((\rho,T,Y_e)\)
* large WinVNE-based nuclear species set
* reduced `aprox21`-like NSE set
* selectable nuclear-mass prescriptions
* selectable nuclear statistical-weight prescriptions
* Coulomb corrections
* conventional two-dimensional Newton NSE solver
* nested one-dimensional NSE solver with safeguarded root finding
* continuation of chemical-potential variables across thermodynamic grids
* Helmholtz EOS interface
* generation of CompOSE-style thermodynamic and composition tables
* diagnostics for composition moments, mass excess, and Coulomb energy

## Repository structure

```text
.
├── data/
│   └── winvn_v2.0.dat
├── src/
│   ├── module_nse.f90
│   ├── module_nuclear_data_winvne.f90
│   ├── module_nuclear_data_HS.f90
│   ├── module_nuclear_mass_policy.f90
│   ├── module_stat_weight_policy.f90
│   ├── module_stat_weight_HS.f90
│   ├── module_nse_species_policy.f90
│   ├── module_ptf_rauscher.f90
│   ├── module_eos_helmholtz.f90
│   ├── make_compose_helmholtz.f90
│   ├── make_nse_table.f90
│   ├── compose_extraction.f90
│   └── ...
├── utils/
│   ├── pl_nse.py
│   └── pl_conv.py
└── makefile
```

The central NSE implementation is contained in

```text
src/module_nse.f90
```

while the current main development driver is

```text
src/make_compose_helmholtz.f90
```

which combines an NSE composition with the Helmholtz EOS and constructs a CompOSE-style table.

## Nuclear data

### WinVNE

The repository contains

```text
data/winvn_v2.0.dat
```

with nuclear masses, ground-state spins, and tabulated nuclear partition functions.

The mass excess in the WinVNE file is an atomic mass excess. Internally the electron rest masses are removed,

$$
m_{\mathrm{exc,nuc}}
=
m_{\mathrm{exc,atomic}}-Z m_e c^2,
$$

so that the NSE calculation uses bare nuclear masses.

### Rauscher partition functions

A Rauscher partition-function table can also be read by

```text
module_ptf_rauscher.f90
```

The Rauscher table is not distributed with this repository and must be supplied separately.

### HS nuclear data

The code can also use nuclear masses and statistical weights associated with the Hempel-Schaffner-Bielich-type nuclear treatment.

The corresponding nuclear-data binary file is not included in this repository and must be supplied separately.

## Nuclear-data policies

The current implementation separates three choices that were previously implicit in the NSE initialization.

### Nuclear masses

The nuclear-mass source is selected using

```fortran
type(nuclear_mass_policy_t) :: nuclear_mass_policy
```

Available sources are

```fortran
NUCLEAR_MASS_WINVNE
NUCLEAR_MASS_HS
```

with an optional fallback source.

For example,

```fortran
nuclear_mass_policy%primary  = NUCLEAR_MASS_WINVNE
nuclear_mass_policy%fallback = NUCLEAR_MASS_NONE
```

uses WinVNE masses only.

### Statistical weights

The statistical-weight source is controlled by

```fortran
type(stat_weight_policy_t) :: stat_weight_policy
```

Available prescriptions are

```fortran
STAT_WEIGHT_WINVNE
STAT_WEIGHT_RAUSCHER
STAT_WEIGHT_HS
```

with an optional fallback.

For example,

```fortran
stat_weight_policy%primary  = STAT_WEIGHT_HS
stat_weight_policy%fallback = STAT_WEIGHT_WINVNE
```

uses the HS prescription whenever requested and falls back to the WinVNE statistical weight otherwise.

### NSE species set

The large WinVNE NSE set is selected using

```fortran
type(nse_species_policy_t) :: species_policy
```

Current choices include

```fortran
NSE_SPECIES_LEGACY
NSE_SPECIES_ALL_WINVNE
```

The legacy selection keeps nuclei available in the Rauscher set together with WinVNE nuclei having \(Z<87\).

## Initializing an NSE network

A typical initialization is

```fortran
use module_nse
use module_nuclear_data_winvne
use module_nuclear_data_HS
use module_ptf_rauscher
use module_stat_weight_policy
use module_nuclear_mass_policy
use module_nse_species_policy

type(nse_network_t) :: net
type(stat_weight_policy_t) :: stat_weight_policy
type(nuclear_mass_policy_t) :: nuclear_mass_policy
type(nse_species_policy_t) :: species_policy

call init_winvne("data/winvn_v2.0.dat")
call init_nuclear_data_HS("/path/to/HS-nuclear-data")
call init_ptf_rauscher("/path/to/rauscher-partition-functions")

nuclear_mass_policy%primary  = NUCLEAR_MASS_WINVNE
nuclear_mass_policy%fallback = NUCLEAR_MASS_NONE

stat_weight_policy%primary  = STAT_WEIGHT_HS
stat_weight_policy%fallback = STAT_WEIGHT_WINVNE

species_policy%mode = NSE_SPECIES_LEGACY

call nse_init_winvne( &
     net, stat_weight_policy, species_policy, nuclear_mass_policy)
```

The nuclear properties belonging to a particular NSE set are stored in

```fortran
type(nse_network_t)
```

so that multiple NSE networks can coexist in the same program.

## NSE solvers

Two solver implementations are currently available.

### `calc_nse`

The original solver solves simultaneously for the neutron and proton chemical-potential variables,

```fortran
call calc_nse(net, rho, temp, ye, itrlim, tol, &
     xnse, nsefail, use_TNAguess)
```

where

```text
rho       rest-mass density [g cm^-3]
temp      temperature [K]
ye        electron fraction
xnse      returned mass fractions
nsefail   convergence flag
```

Optional arguments can be used to supply an initial guess and obtain iteration diagnostics.

### `calc_nse_nested_1d`

A newer solver is provided by

```fortran
call calc_nse_nested_1d( &
     net, rho, temp, ye, itrlim, tol, &
     xnse, nsefail, &
     xn_guess=xn_guess, xp_guess=xp_guess, &
     xn_out=xn_out, xp_out=xp_out)
```

It introduces

$$
u=\eta_n,
\qquad
v=\eta_p-\eta_n.
$$

For a fixed \(v\), the inner one-dimensional problem determines \(u\) from baryon normalization. The outer problem then determines \(v\) from the required \(Y_e\).

The outer solve uses bracketing together with a safeguarded Newton iteration.

This formulation is particularly useful when constructing EOS tables because the converged chemical-potential variables at one density can be used as the initial guess at the neighboring density.

## Reduced aprox21 NSE set

A reduced NSE species set corresponding to the physical nuclei represented by the Microphysics `aprox21` network is available through

```fortran
call nse_init_aprox21( &
     net, stat_weight_policy, nuclear_mass_policy)
```

The set contains 20 physically distinct nuclei,

```text
n, p, He3, He4,
C12, N14, O16, Ne20, Mg24, Si28,
S32, Ar36, Ca40, Ti44, Cr48, Cr56,
Fe52, Fe54, Fe56, Ni56
```

The reaction-network variables `H1` and `p` correspond to the same physical proton and are therefore represented by only one NSE species.

## Coulomb correction

The large NSE calculation includes the Coulomb correction implemented through

```fortran
calc_coulomb_HS()
fcoulomb_HS()
```

The saturation-density parameter used by the current network object is stored as

```fortran
net%n0_fm
```

and can be changed if needed.

## Helmholtz EOS

The repository contains a Fortran Helmholtz-EOS implementation in

```text
src/module_eos_helmholtz.f90
```

with the main interfaces

```fortran
init_eos
eos_all
eos_get_misc
eos_get_temp_from_pres
eos_get_temp_from_eps
```

The EOS can be evaluated using NSE-derived quantities such as

```text
Ye
total ion abundance
average nuclear mass excess per baryon
```

in addition to density and temperature.

The Helmholtz table itself is not distributed with this repository and must be supplied separately.

## CompOSE-style NSE + Helmholtz table

The current development driver

```text
src/make_compose_helmholtz.f90
```

constructs a three-dimensional table on

$$
(n_b,T,Y_q).
$$

The density and temperature grids are logarithmic, while the \(Y_q\) grid is linear.

The current implementation uses

```text
nuclear masses:       WinVNE
statistical weights:  HS with WinVNE fallback
species set:           legacy large NSE set
```

at every grid point.

For each point it

1. solves NSE,
2. checks baryon-number and charge conservation,
3. computes composition moments,
4. evaluates the Helmholtz EOS,
5. constructs the CompOSE thermodynamic quantities,
6. writes thermodynamic and composition arrays to HDF5.

The density conversion is

$$
\rho=m_u n_b,
$$

with \(n_b\) expressed in \(\mathrm{fm}^{-3}\).

### Thermodynamic quantities

The table generator currently constructs the CompOSE quantities

```text
Q1 ... Q7
cs2
```

including

$$
Q_1=\frac{P}{n_b},
$$

the entropy per baryon,

$$
Q_2=S,
$$

and the normalized internal/free-energy quantities.

The chemical-potential quantities are reconstructed consistently from the free energy. In particular,

$$
\mu_l =
\left.
\frac{\partial F_b}{\partial Y_q}
\right|_{n_b,T}
$$

is evaluated numerically along the equally spaced \(Y_q\) grid.

The composition output includes

```text
Ye
Yn
Yp
YH2
YH3
YHe3
YHe4
Ynuc
Anuc
Znuc
Abar
```

## Parameter file for `make_compose_helmholtz`

The driver expects a text parameter file containing alternating comment/header and value lines.

A template is

```text
# WinVNE nuclear data
data/winvn_v2.0.dat

# Rauscher partition-function table
/path/to/rauscher.dat

# HS nuclear-data file
/path/to/hs_nuclear_data.bin

# Helmholtz EOS table
/path/to/helm_table.dat

# output HDF5 file
helmholtz_nse.h5

# nnb, nb_min [fm^-3], nb_max [fm^-3]
326 1.0e-12 1.0

# nt, T_min [MeV], T_max [MeV]
80 0.1 100.0

# nyq, Yq_min, Yq_max
60 0.01 0.60
```

The actual ranges and resolutions should be chosen for the intended EOS application.

The \(Y_q\) grid must contain at least three points because the current calculation of the chemical potentials uses a three-point derivative.

Run with

```bash
./bin/a.out parameters.dat
```

## Build

The current development makefile uses the HDF5 Fortran compiler wrapper

```text
h5fc
```

and therefore requires a Fortran compiler and an HDF5 installation with Fortran support.

The intended build command is

```bash
make
```

with the executable written to

```text
bin/a.out
```

### Current development-tree note

The makefile is presently a development makefile rather than a stable build system. The active source driver is selected directly in the `SRC` definitions.

At the current state of the repository, the default target also references

```text
src/module_compose_hdf5.f90
```

which is not yet tracked in the repository. Consequently the current `main` branch is not yet a fully self-contained build of `make_compose_helmholtz`.

Several older development drivers also predate the recent NSE API refactor and may require updating before compilation.

## Numerical checks

For every converged NSE solution, the primary consistency checks are

$$
\sum_i X_i \simeq 1
$$

and

$$
\sum_i \frac{Z_i}{A_i}X_i \simeq Y_e.
$$

The CompOSE-table driver checks these conditions explicitly and aborts if their errors exceed the prescribed tolerance.

Other useful diagnostics include

```text
Abar
total ion abundance
mass excess per baryon
free neutron abundance
free proton abundance
alpha-particle abundance
heavy-nucleus moments
Coulomb energy
```

## Energy convention

The NSE calculation uses bare nuclear masses.

When the average nuclear mass excess is passed to the Helmholtz EOS or compared with charge-neutral tabulated EOS data, the electron rest-mass contribution

$$
Y_e m_e c^2
$$

is added back to the nuclear mass excess.

This convention is important when comparing NSE + Helmholtz quantities with tables such as DD2/CompOSE.

## Development status

This repository is currently being used to develop and test

* robust NSE solution methods over large EOS grids,
* alternative nuclear-mass prescriptions,
* alternative nuclear statistical-weight prescriptions,
* consistency between NSE and tabulated nuclear EOS data,
* Helmholtz-EOS treatment outside the tabulated-NSE regime,
* generation of CompOSE-format EOS tables.

The code should therefore currently be regarded as research/development software rather than a stable public library.
