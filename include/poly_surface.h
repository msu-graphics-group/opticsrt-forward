#pragma once
#include "surface.h"

/// Make it possible to re-base polynomials; the actual use is controlled by m_useRebasePoly.
#define USE_REBASED_POLYNOMIAL

/// Converts polynomial coefficients from the power-series form to the
/// Q^(con) (Jacobi) form, in place.
bool ConvertToJacobi(std::vector<double>& surface_params);

/// Even aspheric surface:
///
///     z = zpos + conic_term + rho^4 * (p[0] + p[1]*rho^2 + ...),
///
/// where rho = r / rnorm. Up to kMaxRParamsNumber polynomial coefficients are
/// supported. This is the forward-tracing-only CPU implementation.
class PolySurface : public Surface
{
  public:
    PolySurface() : Surface()
    {
      // f(r=0) is 0 for this formula.
      center_z = 0;
      m_jacobi = false;
    };

    PolySurface(bool use_jacobi) : Surface()
    {
      center_z = 0;
      m_jacobi = use_jacobi;
    };

    virtual ~PolySurface() = default;

    bool UsesJacobi() const { return m_jacobi; }

    void CommitDeviceData() override {};

    virtual void GetParam(std::vector<double>& param) const override;

    /// Surface-specific ID of the given parameter.
    virtual int GetParamID(int i) const override;

    /// The number of actually used parameters in 'surface_params'.
    virtual int GetNumParam() const override;

    /// The value of 'used_rparams' passed to UpdateSurface().
    virtual int GetUsedRparams() const override;

    void UpdateSurface(const double* surface_params, const uint32_t used_rparams) override;

    void TraceRays(const IntersectionParameters& inters_params,
                   const double* zpos,
                   const double* surface_params,
                   const uint32_t used_rparams,
                   const double* eta_inv,
                   const double* ray_pos,
                   const double* ray_dir,
                   const uint8_t* ray_mask,
                   double* hit_pos,
                   double* hit_dir,
                   uint8_t* hit_mask,
                   double* hit_t,
                   double* hit_normal,
                   HitPointData* hpd,
                   uint32_t n_rays_x, uint32_t n_rays_y, uint32_t n_batches, uint32_t z_pos_size) override;

    bool IsFlat() const override;

  protected:
    virtual void kernel1D_UpdateInterParams(const IntersectionParameters& inters_params, int a_dummyIters);
    virtual void kernel1D_UpdateSurface(const double* surface_params, uint32_t used_rparams);

#ifdef USE_REBASED_POLYNOMIAL
    /// Prepare coefficients of the "re-based" polynomial.
    virtual void kernel1D_CalcRebased(double rho_max = 1.1, int n_centres = 11);
#endif

    /// @name Surface function
    //@{
    /// Calculates f(r^2) and df/d(r^2).
    /// @return true if the surface is defined in the point, false otherwise.
    virtual bool CalcSurfFuncAndDeriv(double r2, double& f, double& df) const override;

    /// Calculates f(r^2), f'(r^2) and f''(r^2).
    /// @return true if the surface is defined in the point, false otherwise.
    virtual bool CalcSurfFuncAndDerivs(double r2, double& f, double& df, double& d2f) const override;

    /// Calculates max(radius^2) for which the surface function is defined.
    double CalcMaxR2() const;
    //@}

    void DebugPrint(const std::string& a_what, std::ostream& a_out) override;

  protected:
    double m_rnormInv;
    bool   m_isFlat;
    PackedParameters surf_data;

    /// The actual number of polynomials used, <= kMaxRParamsNumber.
    int32_t npoly;

#ifdef USE_REBASED_POLYNOMIAL
    /// Coefficients of the polynomial expanded around rho_base (outer index).
    std::vector<std::vector<double> > coeffs_rebased;
    std::vector<std::vector<float> >  coeffs_rebased4;

    /// Step of rho_base^2 between expansion centers.
    double rho2_base_step;
    float  rho2_base_step4;
#endif
};
