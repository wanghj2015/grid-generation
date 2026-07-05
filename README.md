# Grid Generation for Curvilinear Field-Line-Following Coordinates

Fortran code for generating curvilinear magnetic **field-line-following coordinate
systems** for ionosphere–plasmasphere modeling.

It builds two complementary programs:

| Program | What it does |
|---------|--------------|
| **`dipole_coordinates`** | Generates coordinates **analytically** for an axis-centered dipole field (closed-form formulas and an exact algebraic coordinate inversion — no field-line tracing). This is the fast, exact reference case. |
| **`general_coordinates`** | Generates coordinates **numerically** by tracing field lines through the IGRF geomagnetic field (Runge–Kutta integration + numerically computed metric). With `use_dipole=1` it reduces the field to a pure centered dipole, reproducing the analytical result numerically. |

The paper compares the analytical dipole against the numerically traced "IGRF
dipole," which is why both programs are kept.

## Paper

This code accompanies:

> **A general curvilinear magnetic field-line-following coordinate system for ionosphere-plasmasphere modeling**
> Houjun Wang - wanghoujun@gmail.com

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5768674.svg)](https://doi.org/10.5281/zenodo.5768674)

## Repository layout

```
src/            Fortran sources for both programs (shared modules + per-program code)
data/           igrf13coeffs.txt — IGRF-13 spherical-harmonic coefficients
scripts/run/    Ready-to-run shell scripts that build a run directory and launch a program
scripts/plot/   Python plotting scripts for the output
CMakeLists.txt
```

## Requirements

- A Fortran compiler (tested with `gfortran`)
- CMake ≥ 3.16
- Python 3 with NumPy + Matplotlib (only for the plotting scripts)

## Building

```bash
cmake -S . -B build          # configure (Release/-O2 by default)
cmake --build build -j       # compile both programs
```

The executables are written to `build/bin/`:

```
build/bin/dipole_coordinates
build/bin/general_coordinates
build/bin/igrf13coeffs.txt   # copied here for convenient direct runs
```

For a debug build with runtime checks:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build -j
```

## Running

Each program reads two files **from the current working directory**:

- `grid_params.namelist` — run parameters (see below)
- `igrf13coeffs.txt` — the IGRF coefficient table (in `data/`)

and writes unformatted `fort.*` output files there (see
[Output files](#output-files) for what each one contains).

### Quick start (recommended)

The run scripts set up a clean run directory, drop in the coefficient file, write
a namelist, and launch the program:

```bash
bash scripts/run/run_dipole.sh       # analytical dipole      -> run/dipole/
bash scripts/run/run_igrf_dipole.sh  # numerical IGRF dipole  -> run/igrf/dipole/
bash scripts/run/run_igrf_full.sh    # full IGRF field        -> run/igrf/full/
```

### Running by hand

```bash
mkdir -p run/mycase && cd run/mycase
cp ../../data/igrf13coeffs.txt .
# create grid_params.namelist (see the table below), then:
../../build/bin/general_coordinates > grid.output 2>&1
```

## Namelist parameters (`grid_params.namelist`)

All parameters live in a single `&grid_params ... /` group. Example:

```fortran
&grid_params
epoch   = 2000.0,   ! IGRF epoch (year)
nptx    = 101,      ! points along each field line (odd)
nlp     = 45,       ! number of field lines / flux tubes
nmp     = 1,        ! number of magnetic-longitude planes
theta1  = 45.0,     ! inner boundary colatitude at r = 1 Re [deg]
theta2  = 82.0,     ! outer boundary colatitude at r = 1 Re [deg]
hmin    = 90.0,     ! minimum altitude above Re [km]
ds_step = 5.0,      ! field-line integration step [km] (general only)
use_psi = 0,        ! 0 = mu coordinate, 1 = modified psi coordinate
aa      = 3.0,      ! psi stretching parameter (used when use_psi = 1)
use_dipole = 1,     ! 1 = pure dipole, 0 = full IGRF (general only)
/
```

| Parameter | Type | Meaning |
|-----------|------|---------|
| `epoch`      | real | IGRF epoch year used to evaluate the geomagnetic coefficients (e.g. `2000.0`). |
| `nptx`       | int  | Number of grid points **along** each field line. Should be **odd** — the midpoint `(nptx+1)/2` sits on the magnetic equator. |
| `nlp`        | int  | Number of field lines (flux tubes) spanning the L-shell range. |
| `nmp`        | int  | Number of magnetic-longitude (φ) planes. Use `1` for a single meridian. |
| `theta1`     | real | Colatitude [deg] at `r = 1 Re` of the **inner** (high-latitude) boundary field line. |
| `theta2`     | real | Colatitude [deg] at `r = 1 Re` of the **outer** (low-latitude) boundary field line. `theta1 < theta2`. |
| `hmin`       | real | Minimum altitude above the reference radius `Re` [km]; field lines are cut off below `Re + hmin` (ionospheric base, e.g. `90.0`). |
| `ds_step`    | real | Step size [km] for integrating the field-line equations. Used by `general_coordinates`; ignored by the analytical `dipole_coordinates`. |
| `use_psi`    | int  | Field-aligned coordinate choice: `0` = standard dipole coordinate `mu`; `1` = modified coordinate `psi` (concentrates resolution near the ionosphere). |
| `aa`         | real | Stretching factor for the `psi` coordinate; larger values push more points toward low altitude. Only used when `use_psi = 1`. |
| `use_dipole` | int  | `1` = reduce the field to a pure centered dipole (zeroes every IGRF coefficient except g₁₀); `0` = use the full IGRF field. Used by `general_coordinates`. For `dipole_coordinates` the field is always a dipole regardless. |

`Re = 6371.2 km` (mean spherical reference radius) is fixed in the code.

## Output files

Both programs write Fortran **unformatted** `fort.*` files into the run directory.
Each is written per magnetic longitude (`nmp`) and per field line (`nlp`).
`general_coordinates` outputs SI units (metres, tesla); the analytical
`dipole_coordinates` outputs normalized quantities.

| File | Written by | Contents |
|------|-----------|----------|
| `fort.20` | both | Grid: the grid dimensions (`nlp, nmp, npts`) and index arrays, then per field line the position (radius `r`, colatitude `theta`, and longitude `phi` for the general program), the field magnitude `|B|`, and the arc length `s`. |
| `fort.21` | both | Metric scale factors `h1, h2, h3` along each field line. |
| `fort.22` | both | Coordinate curves running across the field lines (constant-`mu`/`psi` surfaces). |
| `fort.24` | both | The magnetic coordinates themselves: `mu` (or `psi`), `chi`, `phi`. |
| `fort.51` | dipole only | Reference equatorial arc length of each field line. |
| `fort.23` | general only | Contravariant and covariant metric basis vectors. |
| `fort.25` | general only | Diagnostics: angles between the basis vectors and a field-magnitude ratio. |
| `fort.26` | general only | Additional basis-vector angle diagnostics. |

The plotting scripts in `scripts/plot/` read these files back.

## Plotting

The scripts in `scripts/plot/` read the `fort.*` output and produce figures. Each
takes the run directory (or directories) as an argument and defaults to the
directories created by the run scripts. Pass `-h` to any of them for usage. All
figures are written to **`run/plots/`** (the script prints the saved path).

```bash
python scripts/plot/plot_grid.py          # grid geometry        (default: run/igrf/full)
python scripts/plot/plot_basis_angles.py  # basis-vector angles  (default: run/igrf/full)
python scripts/plot/plot_euler_ratio.py   # Euler-potential ratios (default: run/igrf/full)
python scripts/plot/plot_mertic.py        # metric terms: dipole vs IGRF dipole
                                          #   (defaults: run/dipole  run/igrf/dipole)

# override the run directory, e.g. plot a different case:
python scripts/plot/plot_grid.py run/igrf/dipole
```

## License

See [LICENSE](LICENSE).
