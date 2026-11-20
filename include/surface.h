#pragma once

// C headers
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cassert>

// C++ headers
#include <algorithm>
#include <array>
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "isurface.h"
#include "LiteMath.h"

typedef LiteMath::double3 Vect3d;
typedef LiteMath::double3 Point3d;

/// Enables the preliminary ray-marching search for the Newton starting point.
/// Set by default; the "subdivideRaySegment == 1000000" mode uses it.
#define NEW_START_POINT4NEWTON

/// Surface base class: forward tracing only.
///
/// A surface is z = zpos + f(r^2) where r^2 = x^2 + y^2. Derived classes
/// define f(), f' and f''. This class implements the actual ray/surface
/// intersection, refraction and the bulk CPU ray-tracing kernel.
class Surface : public ISurface
{
  public:
    Surface();

    static constexpr double MAX_DISTANCE      = 1.e5;
    static constexpr double NEWTON_STEP_BOUND = 5;
    static constexpr double EPSILON           = 1e-6;

    /// The number of actually used parameters in 'surface_params'.
    virtual int GetNumParam() const = 0;

    /// Set intersection computation parameters.
    void SetIntersectionParams(const IntersectionParameters& params) override
    {
      m_hitEps     = params.hitEps;
      m_minT       = params.minT;
      m_flags      = params.flags;
      m_useRebasePoly = params.rebasePoly;
      m_useFloat      = params.useFloat;
    }

    void SetInitIntersectionParams(const IntersectionInitParameters& params) override
    {
      m_maxIterations       = params.maxIterations;
      m_subdivideRaySegment = params.subdivideRaySegment;
      m_jacobi              = params.useQtype;
    }

    IntersectionParameters GetIntersectionParams() override
    {
      IntersectionParameters p;
      p.hitEps = m_hitEps;
      p.minT = m_minT;
      p.flags = m_flags;
      p.rebasePoly = m_useRebasePoly;
      p.useFloat = m_useFloat;
      return p;
    }

    IntersectionInitParameters GetInitIntersectionParams() override
    {
      IntersectionInitParameters p;
      p.maxIterations = m_maxIterations;
      p.subdivideRaySegment = m_subdivideRaySegment;
      p.useQtype = m_jacobi;
      return p;
    }

    void UpdateIntersectionParamsCPU(const IntersectionParameters& params) override;

    /// Save cross-section table, i.e. z vs. r (not vs. r^2!).
    virtual void SaveSagCrossSection(const std::string& path, size_t n_points, const double* surface_params) override;
    /// Save the slope cross-section table, i.e. dz/dr vs. r (not vs. r^2!).
    virtual void SaveSagSlopeCrossSection(const std::string& path, size_t n_points, const double* surface_params) override;

    /// Calculate the lens section z=f(r) on a uniform mesh in 'r'.
    virtual void CalcCrossSection(size_t n_points, std::vector<double>& r_arr, std::vector<double>& z_arr) const override;

    inline double GetAperture() const { return aperture; };

    /// The value of z (without adding zpos!) for r=0.
    inline double GetZCenter() const { return center_z; };

    virtual bool IsFlat() const { return false; }

    /// Left bounding plane, WITHOUT 'zpos'.
    virtual double GetZmin() const override { return zmin; };
    /// Right bounding plane, WITHOUT 'zpos'.
    virtual double GetZmax() const override { return zmax; };
    /// Radius where the surface ends.
    virtual double GetMaxRadius() const override { return std::min(aperture, max_radius); };

    /// Iteration counters collected during the last TraceRays() call.
    Stats GetStats() const override { return {m_newton_iterations, m_raymarch_iterations}; }

    virtual void GetParam(std::vector<double>& param) const = 0;

    /// Surface-specific ID of the given parameter.
    virtual int GetParamID(int i) const = 0;

    /// Trace rays through this surface (see ISurface::TraceRays for details).
    virtual void TraceRays(const IntersectionParameters& params,
                           const double* zpos, const double* surface_params,
                           const uint32_t used_rparams,
                           const double* eta_inv,
                           const double* ray_pos, const double* ray_dir, const uint8_t* ray_mask,
                           double* hit_pos, double* hit_dir, uint8_t* hit_mask, double* hit_t, double* hit_normal,
                           HitPointData* hpd,
                           uint32_t n_rays_x, uint32_t n_rays_y, uint32_t n_batches, uint32_t z_pos_size) override;

    void GetExecutionTime(const char* a_funcName, float a_out[4]) override
    {
      if (std::string(a_funcName) == "UpdateSurface")
        a_out[0] = update_time;
      else
        a_out[0] = trace_time;
      a_out[1] = update_time;
    };

    virtual std::string GetResourcesRootDir() { return m_resourcesDir; }

    /// @name Operations with HitPointData
    //@{
    /// Sets all data to 0.
    static void Reset(HitPointData& data);

    /// Calculates the direction of the surface normal.
    static Vect3d CalcNormal(const double r2, const double df, const Point3d& hitPos, bool unit = true);

    /// Calculates the direction of refraction. Returns 'false' for T.I.R.
    static bool CalcRefrDir(const double r2, const double df, const Vect3d& dir, const Point3d& hitPos, double eta, Vect3d& dir_out, Vect3d& n);
    //@}

  protected:
    /// Calculates and remembers the bounds for 'z' WITHOUT zpos.
    /// Usually called only within @c CalcZminZmax()!
    virtual void kernel1D_CalcZminZmax(int npt, double dr_fraction);

#ifdef NEW_START_POINT4NEWTON
    /// Calculates the bounds of f, f' and f'' for all levels of subdivision.
    void CalcBoundsOfFuncAndDeriv(int nlevels = 4, int nfirst = 4, int reduction = 6);
#endif

  public:
    /// @name Surface function; defined in derived classes
    //@{
    /// Calculates the surface function and its derivative, i.e. f(r^2) and df/d(r^2).
    /// @return true if the surface is defined in the point, false otherwise. In the latter case
    ///              the output values are calculated for r2=(max_radius)^2.
    virtual bool CalcSurfFuncAndDeriv(double r2, double& f, double& df) const = 0;

    /// Calculates f(r^2), f'(r^2) and f''(r^2).
    /// @return true if the surface is defined in the point, false otherwise. In the latter case
    ///              the output values are calculated for r2=(max_radius)^2.
    virtual bool CalcSurfFuncAndDerivs(double r2, double& f, double& df, double& d2f) const = 0;
    //@}

  protected:
    /// @name Ray tracing
    //@{
    virtual void kernel1D_TraceRays(const double* zpos, const double* eta_inv,
                                    const double* ray_pos, const double* ray_dir, const uint8_t* ray_mask,
                                    double* hit_pos, double* hit_dir, uint8_t* hit_mask, double* hit_t, double* hit_normal,
                                    HitPointData* hpd,
                                    uint32_t n_iters, uint32_t n_rays_square, uint32_t z_pos_size);

  public: // To let users call Intersection() directly
    /// Finds the intersection of the ray with the surface. Returns hit data, or
    /// a zero-filled structure with ray_valid=false on failure.
    HitPointData Intersection(double zpos, Vect3d& origin, const Vect3d& dir, bool& ray_valid) const;

    /// Same as above, additionally reports the number of Newton and ray-marching
    /// iterations spent on this ray.
    HitPointData Intersection(double zpos, Vect3d& origin, const Vect3d& dir, bool& ray_valid,
                              uint32_t& newton_iters, uint32_t& raymarch_iters) const;

    /// Safeguarded hybrid intersection: builds a sign-changing bracket from the
    /// surface z-slab or the interval-search seed and then alternates Newton
    /// steps with bisection. Falls back to Intersection() if no bracket can be
    /// established or if the refinement does not reach a tight residual.
    HitPointData IntersectionSafe(double zpos, Vect3d& origin, const Vect3d& dir, bool& ray_valid,
                                  uint32_t& newton_iters, uint32_t& raymarch_iters) const;
  protected:

    /// Move the ray close to the surface where its origin obeys r <= aperture.
    /// @return natural parameter 't' (the distance from the ray origin to its new point)
    /// or -1 if failed (i.e. all the available variants are outside of max_radius)
    bool MoveRayCloseToLens(double zpos, const Point3d& origin, const Vect3d& direction, Point3d& origin_new, double& t0,
                            Point3d& nearest, uint32_t& raymarch_iters) const;

#ifdef NEW_START_POINT4NEWTON
    /// When the ray enters the lens domain cylinder of max_radius between the bounding planes.
    double RayEntersLensDomain(double zpos, const Point3d& origin, const Vect3d& direction) const;

    /// Finds the low bound of 't' for the intersection interval.
    int ApproxIntersection(const Vect3d& origin,
                           const Vect3d& direction,
                           const double zpos,
                           int& iter,
                           const double threshold1,
                           const double threshold2,
                           double& t) const;
#endif
    //@}

    void PrintError(const std::string& str);

  protected:
    /// Whether the coefficients relate to Q^{con}_n(rho^2), otherwise to powers of rho^2.
    bool m_jacobi = false;

    /// Define how the intersection is calculated, with which tolerance etc.
    double   m_hitEps;
    double   m_minT;
    uint32_t m_flags;
    int32_t  m_useRebasePoly = 0; ///< re-base polynomials (PolySurface only)
    int32_t  m_useFloat      = 0; ///< float precision for re-based polynomials

    uint32_t m_maxIterations;
    uint32_t m_subdivideRaySegment;

    /// Aperture, i.e. the radius past which the surface is cut.
    double aperture;

    /// @name Slaves calculated from surface parameters
    //@{
    /// Max allowed radius for the surface.
    double max_radius;

    /// z-coordinate of the center (z for x=y=0), WITHOUT zpos.
    double center_z;

    /// Bounds of z, WITHOUT zpos. Safe margins are applied.
    float zmin, zmax;
    //@}

#ifdef NEW_START_POINT4NEWTON
    /// min/max of f'' over cells subdividing [0, r^2_max], per subdivision level.
    std::vector<std::vector<double> > d2f_max_arr_arr, d2f_min_arr_arr;

    /// min/max of f' over cells subdividing [0, r^2_max], per subdivision level.
    std::vector<std::vector<double> > df_max_arr_arr, df_min_arr_arr;

    /// min/max of f over cells subdividing [0, r^2_max], per subdivision level.
    std::vector<std::vector<double> > f_max_arr_arr, f_min_arr_arr;
#endif

    std::string m_resourcesDir = "data";
    float update_time = 0.0f;
    float trace_time  = 0.0f;

    /// Accumulated during TraceRays(); see GetStats().
    uint32_t m_newton_iterations = 0;
    uint32_t m_raymarch_iterations = 0;
};
