#pragma once

#include <cstdint>

namespace lens_sim
{
  /// Maximum number of polynomial ("r-parameter") coefficients of an
  /// even aspheric surface.
  static constexpr const uint32_t kMaxRParamsNumber = 14;

  /// Parameters of one even aspheric surface.
  ///
  /// The surface is rotationally symmetric around the optical axis Oz and is
  /// described by
  ///
  ///     z = zpos + f(r^2),
  ///
  /// where f() is the conic term plus an even polynomial in the normalized
  /// radius rho = r / rnorm.
  struct PackedParameters
  {
    double aperture;  ///< aperture radius
    double curvature; ///< curvature (inverse radius)
    double conic;
    double rnorm;     ///< normalization radius
    double rparams[kMaxRParamsNumber];
  };

  struct IntersectionInitParameters
  {
    uint32_t maxIterations;
    /// If this value is >1, then the ray segment between the bounding planes
    /// of the lens is subdivided into @c subdivideRaySegment intervals; then
    /// find the FIRST where the ray goes past the surface and use its middle
    /// as the initial point for Newton iterations.
    uint32_t subdivideRaySegment;
    bool useQtype; ///< use Q^(con) (Jacobi) polynomial representation
  };

  struct IntersectionParameters
  {
    double   hitEps;     ///< tolerance of the Newton iterations
    double   minT;       ///< consider hit with t < minT invalid
    int32_t  rebasePoly; ///< re-base polynomials around the nearest center
    int32_t  useFloat;   ///< float precision for re-based polynomials
    uint32_t flags;      ///< combination of OUTFLAG_* values
  };

  static inline IntersectionInitParameters MakeIntersectionInitParameters()
  {
    IntersectionInitParameters res;
    res.maxIterations = 24;
    res.subdivideRaySegment = 0;
    res.useQtype = false;
    return res;
  }

  static inline IntersectionParameters MakeIntersectionParameters()
  {
    IntersectionParameters res;
    res.hitEps = 1e-9;
    res.minT = -1.0;
    res.rebasePoly = 0;
    res.useFloat = 0;
    res.flags = 0;
    return res;
  }

} // namespace lens_sim
