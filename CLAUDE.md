# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CPU port of the CUDAKerr GPU code for simulating null geodesics (light paths) in Kerr spacetime around black holes. The primary scientific application is modeling X-ray reflection and reverberation from accretion discs.

## Build Commands

```bash
# Configure build system (run once, or after CMakeLists.txt changes)
cmake .

# Build a specific executable
cmake --build . --target <target_name>

# Examples
cmake --build . --target emissivity
cmake --build . --target imageplane_disc_image
cmake --build . --target gsl_test
```

Executables are output to `bin/`, libraries to `obj/`.

## Running Applications

Each application reads parameters from a file (default: `par/<appname>.par`). Template parameter files are in `par_example/`.

```bash
./bin/emissivity                              # uses par/emissivity.par
./bin/emissivity --parfile custom.par         # use custom parameter file
```

Parameter file format: `name = value` (one per line).

## Running Tests

```bash
cmake --build . --target gsl_test
./bin/gsl_test
```

No comprehensive test framework is integrated; `gsl_test` validates the GSL dependency.

## Dependencies

- **cfitsio**: FITS file I/O (required)
- **GSL** (GNU Scientific Library): required for most tools
- **HDF5**: required for `pointsource_mapper`

System-specific library paths go in `cmake.local` (not tracked by git). Example:
```cmake
set(CFITSIO_PATH /opt/software/cfitsio/4.6.3)
```

## Architecture

### Three-layer design

**Layer 1 — Core ray tracer** ([src/raytracer/](src/raytracer/))
- `raytracer.h`: Template class integrating null geodesics through Kerr spacetime. Handles photon propagation, redshift calculations, and configurable integration precision.
- Ray source classes provide different photon injection geometries:
  - `PointSource` — point-like lamppost source
  - `ImagePlane` — 2D observing plane
  - `HEALPixPointSource` — spherical HEALPix-pixelated extended sources
  - Velocity variants (`*_vel`) for Doppler-shifted sources

**Layer 2 — Scientific utilities** ([src/include/](src/include/))
- `kerr.h`: Kerr spacetime math — event horizon, ISCO, disc velocities, metric tensor, coordinate transforms
- `par_file.h` / `par_args.h`: Parameter parsing from files and command-line overrides
- `fits_output.h`: FITS format I/O for astronomical data
- `healpix.h`: HEALPix spherical pixelization

**Layer 3 — Scientific applications** (various `src/` subdirectories)
Each application picks a ray source, runs ray traces, and computes specific observables:

| Subdirectory | Executables | Purpose |
|---|---|---|
| `emissivity/` | `emissivity` | Disc emissivity profiles |
| `imageplane/` | `imageplane_disc_image` | Ray-traced disc images (FITS) |
| `ray_paths/` | `trace_rays*` | Raw ray path output |
| `outflow/` | `outflow`, `outflow_spectrum`, `pcyg*` | Wind/jet outflow and P-Cygni profiles |
| `mapper/` | (mapper libs) | Geometric mapping, redshift/beaming |
| `source_tracer/` | (source tracer) | X-ray reverberation |
| `return_radiation/` | `disc_source_photonfrac*` | Radiation return fractions |
| `lamppost/` | `pointsource_to_disc` | Lamppost geometry |
| `healpix/` | `healpix_to_disc`, `healpix_disc_source_photonfrac` | Extended source simulations |

### Adding a new application

1. Create a subdirectory under `src/`
2. Write a `CMakeLists.txt` defining `add_executable()` and `target_link_libraries()` (link `raytracer` and any needed utility libs)
3. Add `add_subdirectory(<name>)` to `src/CMakeLists.txt`

### Parameter access pattern in application code

```cpp
par_file.get_parameter<double>("spin")       // required; throws if missing
par_args.get_parameter<double>("--spin")     // optional command-line override
```
