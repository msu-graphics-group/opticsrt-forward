#pragma once

#include "even_asphere.h"

#include <cstdint>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

static constexpr const uint32_t POLY_APERTURE_ID  = 0;
static constexpr const uint32_t POLY_CURVATURE_ID = 1;
static constexpr const uint32_t POLY_CONIC_ID     = 2;
static constexpr const uint32_t POLY_RNORM_ID     = 3;
static constexpr const uint32_t POLY_RPARAMS_FIRST_ID   = 4;
static constexpr const uint32_t POLY_RPARAMS_MAX_PARAMS = lens_sim::kMaxRParamsNumber;

static constexpr const uint32_t OUTFLAG_STORE_RAYS   = 1;
static constexpr const uint32_t OUTFLAG_STORE_NORMAL = 2;
static constexpr const uint32_t OUTFLAG_STORE_T      = 4;

using lens_sim::IntersectionParameters;
using lens_sim::IntersectionInitParameters;
using lens_sim::PackedParameters;

/// Data about a hit point, enough to compute the normal and the refracted ray.
struct HitPointData
{
  double dir[3]; ///< ray direction (unit vector) BEFORE the hit
  double t;      ///< ray path length from its origin to this hit point
  double f;      ///< f = f(r^2), i.e. xyz.z WITHOUT zpos
  double df;     ///< f' = dz/d(r^2)
  double d2f;    ///< f'' = d{dz/d(r^2)}/d(r^2)
};

static inline HitPointData make_hpd()
{
  HitPointData res;
  res.dir[0] = res.dir[1] = res.dir[2] = 0;
  res.t   = 0.0;
  res.f   = 0.0;
  res.df  = 0.0;
  res.d2f = 0.0;
  return res;
}

/// Cumulative work counters for the last TraceRays() call.
/// They show how much seeding (ray marching) and Newton refinement the
/// intersection algorithm needed for the traced bundle.
struct Stats
{
  uint32_t newton_iters   = 0;
  uint32_t raymarch_iters = 0;
};

/// Rotationally symmetric lens surface z = zpos + f(r^2).
/// The axis of rotation is Oz and it is directed from the first to the last
/// surface, i.e. from the "object" to the "image".
///
/// This interface is forward-tracing only: it evaluates intersections with an
/// incident ray and the refraction of the ray at the hit point.
class ISurface
{
  public:
    virtual ~ISurface() {};

    /// The number of actually used parameters in 'surface_params'.
    virtual int GetNumParam() const = 0;

    /// The value of 'used_rparams' passed to UpdateSurface().
    virtual int GetUsedRparams() const = 0;

    /// Left bounding plane, WITHOUT 'zpos'.
    virtual double GetZmin() const = 0;
    /// Right bounding plane, WITHOUT 'zpos'.
    virtual double GetZmax() const = 0;
    /// Radius where the surface ends.
    virtual double GetMaxRadius() const = 0;

    /// Iteration counters collected during the last TraceRays() call.
    virtual Stats GetStats() const { return Stats{}; }

    /// Update the lens design parameters.
    /// @param[in] surface_params - packed surface parameters, see PackedParameters
    /// @param[in] used_rparams   - number of polynomial coefficients used
    virtual void UpdateSurface(const double* surface_params, const uint32_t used_rparams) = 0;

    virtual void UpdateIntersectionParamsCPU(const IntersectionParameters& params) = 0;

    /// Save a cross-section table, i.e. z vs. r (not vs. r^2).
    virtual void SaveSagCrossSection(const std::string& path, size_t n_points, const double* surface_params) = 0;
    /// Save the slope cross-section table, i.e. dz/dr vs. r (not vs. r^2).
    virtual void SaveSagSlopeCrossSection(const std::string& path, size_t n_points, const double* surface_params) = 0;

    /// Calculate the lens section z=f(r) on a uniform mesh in 'r'.
    virtual void CalcCrossSection(size_t n_points, std::vector<double>& r_arr, std::vector<double>& z_arr) const = 0;

    virtual void SetIntersectionParams(const IntersectionParameters& params) {}
    virtual void SetInitIntersectionParams(const IntersectionInitParameters& params) {}

    virtual IntersectionParameters GetIntersectionParams() { return lens_sim::MakeIntersectionParameters(); }
    virtual IntersectionInitParameters GetInitIntersectionParams() { return lens_sim::MakeIntersectionInitParameters(); }

    /**
     * Trace several bundles, each made of many rays, through this surface.
     *
     * @param[in] params         - intersection algorithm parameters
     * @param[in] zpos           - surface position, size [1] or [n_batches]
     * @param[in] surface_params - packed surface parameters
     * @param[in] used_rparams   - number of polynomial coefficients used
     * @param[in] eta_inv        - RELATIVE refraction index eta(out)/eta(in),
     *                             size [n_batches]
     * @param[in] ray_pos        - incident ray origins, size [n_rays_x*n_rays_y*n_batches*3]
     * @param[in] ray_dir        - incident ray directions, size [n_rays_x*n_rays_y*n_batches*3]
     * @param[in] ray_mask       - incident ray mask, size [n_rays_x*n_rays_y*n_batches]
     * @param[out] hit_pos       - hit positions, size [n_rays_x*n_rays_y*n_batches*3]
     * @param[out] hit_dir       - refracted ray directions, size [n_rays_x*n_rays_y*n_batches*3]
     * @param[out] hit_mask      - valid hit mask, size [n_rays_x*n_rays_y*n_batches]
     * @param[out] hit_t         - hit distances along rays, size [n_rays_x*n_rays_y*n_batches]
     * @param[out] hit_normal    - normals at hit positions, size [n_rays_x*n_rays_y*n_batches*3]
     * @param[out] hpd           - per-hit data (may be nullptr if not required)
     */
    virtual void TraceRays(const IntersectionParameters& params,
                           const double* zpos, const double* surface_params,
                           const uint32_t used_rparams,
                           const double* eta_inv,
                           const double* ray_pos, const double* ray_dir, const uint8_t* ray_mask,
                           double* hit_pos, double* hit_dir, uint8_t* hit_mask, double* hit_t, double* hit_normal,
                           HitPointData* hpd,
                           uint32_t n_rays_x, uint32_t n_rays_y, uint32_t n_batches, uint32_t z_pos_size) = 0;

    /// Report execution time in milliseconds.
    virtual void GetExecutionTime(const char* a_funcName, float a_out[4]) {}

    virtual void CommitDeviceData() {}

    /// Debug method.
    virtual void DebugPrint(const std::string& a_what, std::ostream& a_out) {}
};

/// Factory for the even-aspheric polynomial surface implementation.
std::shared_ptr<ISurface> CreatePolySurface(bool use_jacobi = false);
