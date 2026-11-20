# optics-forward

A small, self-contained C++17 project that implements forward ray tracing through
even-aspheric optical surfaces on the CPU.

This project keeps the forward (light-propagation) path only:

* intersection of an incident ray with a rotationally symmetric even-aspheric
  surface `z = zpos + f(r^2)`;
* a safeguarded hybrid intersection (`IntersectionSafe`) that builds a
  sign-changing bracket and alternates Newton steps with bisection, falling back
  to preconditioned Newton when a bracket cannot be established;
* refraction at the hit point (Snell's law, total internal reflection reporting);
* normals at the hit point;
* CPU tracing of one bundle or several batches of rays, optionally in parallel
  with OpenMP;
* the two forward-search strategies used before Newton iterations:
  interval search and ray marching;
* optional Q^(con) (Jacobi) polynomial basis via `CreatePolySurface(true)`
  (the caller must convert the coefficients with `ConvertToJacobi()`), polynomial
  re-basing and float-precision re-based coefficients;
* iteration statistics (`GetStats()`) for evaluating starting-point strategies;
* forward tracing regression tests against the reference data.

## Layout

```
include/
  even_asphere.h   packed surface/solver parameter definitions
  isurface.h       forward-only surface interface
  surface.h        ray/surface intersection and CPU tracing base class
  poly_surface.h   even-aspheric polynomial surface
  utils.h          lens JSON loading, ray generation, test I/O helpers
  tests.h          forward tracing test declarations and timing helper
src/
  surface.cpp      intersection, Newton/ray-march starting point, CPU kernel
  surface_host.cpp cross-section helpers and debug printing
  poly_surface.cpp even-aspheric surface function and derivatives
  utils.cpp        JSON/data/PNG/EXR helpers
  tests.cpp        forward tracing regression tests
  main.cpp         test runner entry point
ref_data/          reference hit distances and rays used by the tests
103.json           reference lens system
Double Gauss 28 degree field.json
ext/               LiteMath, nlohmann/json, tinyexr, stb
```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

OpenMP is optional. If it is found, ray loops are parallelized; otherwise the
same code runs serially (the `#pragma omp` directives are ignored).

## Run

Run from the project root, because the tests load `103.json`,
`Double Gauss 28 degree field.json` and `ref_data/*` using relative paths:

```bash
./build/optics-forward --output ./output
```

Options:

```
--output <dir>          output directory (default: ./output)
--raymarch              use ray marching to locate Newton's start point
--max-iterations <n>    max Newton iterations (default: 10)
--precision <eps>       intersection tolerance (default: 1e-7)
--help
```

The process returns `0` only if every test passes.

## Tests

The binary runs the forward tracing regression tests:

| Test | Checks |
| --- | --- |
| `test_default_surface` | 512x512 ray bundle against `ref_data/test_surface_t_hit.txt` |
| `test_debug_surface` | 100x100 custom rays and a reference validity mask |
| `test_double_gauss_singleray` | single ray through a multi-surface Double Gauss system |
| `test_nonsquare_ray_batch` | 128x96 non-square bundle, 7 surfaces of `103.json` |
| `test_two_batches` | two 64x64 batches at different field angles |
| `test_batches_sequential` | all batches at once vs. batch-by-batch tracing |
| `test_intersect_behind` | a hit behind the ray origin is reported correctly |
| `raymarch_matches_interval` | ray-marching and interval-search starting points converge to the same hits |

Run them through CTest as well:

```bash
cd build && ctest --output-on-failure
```

CTest sets the working directory to the project root automatically.

Expected summary:

```
All tests PASSED!
```

## API sketch

```cpp
#include "poly_surface.h"

auto surface = CreatePolySurface(); // for Jacobi: CreatePolySurface(true)
                                    // and ConvertToJacobi(surface_params) first

lens_sim::IntersectionInitParameters init = lens_sim::MakeIntersectionInitParameters();
init.maxIterations = 24;
init.subdivideRaySegment = 10;               // interval search; 1000000 = ray march
surface->SetInitIntersectionParams(init);

lens_sim::IntersectionParameters trace = lens_sim::MakeIntersectionParameters();
trace.hitEps = 1e-7;
trace.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;

std::vector<double> surface_params;          // packed as lens_sim::PackedParameters
double zpos = 0.0, eta_inv = 1.0;
surface->UpdateSurface(surface_params.data(), POLY_RPARAMS_MAX_PARAMS);

surface->TraceRays(trace, &zpos, surface_params.data(), POLY_RPARAMS_MAX_PARAMS, &eta_inv,
                   ray_pos, ray_dir, ray_mask,
                   hit_pos, hit_dir, hit_mask, hit_t, hit_normal, hpd,
                   n_rays_x, n_rays_y, n_batches, 1);
```

After tracing, `surface->GetStats()` returns the accumulated iteration counts
for that call.

### Surface parameter layout

`PackedParameters` contains, in order:

```
aperture, curvature, conic, rnorm, rparams[14]
```

The corresponding `surface_params` indices are `POLY_APERTURE_ID`,
`POLY_CURVATURE_ID`, `POLY_CONIC_ID`, `POLY_RNORM_ID` and
`POLY_RPARAMS_FIRST_ID`. Surfaces are stored and traced directly in physical
coordinates.

The output flags follow the convention:

* ray positions are always written;
* `OUTFLAG_STORE_RAYS` controls `hit_dir`;
* `OUTFLAG_STORE_NORMAL` controls `hit_normal`;
* `OUTFLAG_STORE_T` controls `hit_t` and `hit_mask`.
