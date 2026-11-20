#include "poly_surface.h"

#include <cfloat>
#include <iostream>
#include <memory>
#include <ostream>

using LiteMath::SQR;

using std::min;
using std::max;
using std::isnan;
using std::abs;

std::shared_ptr<ISurface> CreatePolySurface(bool use_jacobi) // Please do not use unique_ptr
{                                             // https://stackoverflow.com/questions/36752678/clang-returning-stdunique-ptr-with-type-conversion
  return std::make_shared<PolySurface>(use_jacobi);
}


/// Calculates coefficients for expansion of the Jacobi polynomials in power series
///              P^{0,4}_n(z)=SUM coeff[n][m]*z^m
/*static*/ bool CalcJacobiCoeffs(int nmax, std::vector< std::vector<double> > &coeff)
  {
  coeff.resize(nmax + 1);

  int n, m;

  for (n = 0; n < nmax + 1; n++)
    {
    coeff[n].resize(nmax + 1);
    std::fill(coeff[n].begin(), coeff[n].end(), 0.0);
    }

  coeff[0][0] = 1;

  coeff[1][0] = -2;
  coeff[1][1] =  3;

  coeff[2][0] = 1;
  coeff[2][1] = -7;
  coeff[2][2] = 7;

  for (n = 3; n <= nmax; n++)
    {
    const double rev_a = 1.0 / (n * (n + 4) * (n + 1));
    const double b = (n + 2) * (n + 1);
    //const double c = (2 * n + 3) * (b * xi - 4);
    const double d = (2 * n + 3) * b;
    const double e = (n - 1) * (n + 3) * (n + 2);
    const double f = 4 * (2 * n + 3);

    coeff[n][0] = rev_a * (-f * coeff[n - 1][0] - e * coeff[n - 2][0]);
    for (m = 1; m <= n - 2; m++)
      coeff[n][m] = rev_a * (d * coeff[n - 1][m - 1] - f * coeff[n - 1][m] - e * coeff[n - 2][m]);
    coeff[n][n - 1] = rev_a * (d * coeff[n - 1][n - 2] - f * coeff[n - 1][n - 1]);
    coeff[n][n]     = rev_a * (d * coeff[n - 1][n - 1]);
    }

  return true;
  }


/// Calculates coefficients for expansion of Q^{con}_n polynomials in power series
///              Q^{con}_n(z)=SUM coeff[n][m]*z^m
/*static*/ bool CalcQconCoeffs(int nmax, std::vector< std::vector<double> >& coeff)
  {
  std::vector< std::vector<double> > jacobi_coeff;
  if (!CalcJacobiCoeffs(nmax, jacobi_coeff))
    {
    assert(false);
    return false;
    }

  coeff.resize(nmax + 1);
  for (int n = 0; n < nmax + 1; n++)
    {
    coeff[n].resize(nmax + 1);
    std::fill(coeff[n].begin(), coeff[n].end(), 0.0);
    }

  for (int n = 0; n <= nmax; n++)
    {
    double power = 1; // = 2^k
    for (int k = 0; k <= n; k++)
      {
      double sum = jacobi_coeff[n][k], prod = 1;
      for (int l = 1; l <= n - k; l++)
        {
        prod *= (-1.0) * double(l + k) / double(l);
        sum += jacobi_coeff[n][l + k] * prod;
        }
      coeff[n][k] = power * sum;
      power *= 2;
      }
    }

  return true;
  }

/// Converts representation from simple power series to Q^{con}_n polynomials
/*static*/ bool ConvPowerSeries2Qcon(const std::vector<double> &A, std::vector<double>& C_con)
  {
  const int nmax = (int)A.size() - 1;

  std::vector< std::vector<double> > T;
  if (!CalcQconCoeffs(nmax, T))
    {
    assert(false);
    return false;
    }

  C_con.resize(nmax + 1);

  for (int n = nmax; n >= 0; n--)
    {
    double sum = 0;
    for (int k = n + 1; k <= nmax; k++)
      sum += T[k][n] * C_con[k];

    C_con[n] = (A[n] - sum) / T[n][n];
    }

  return true;
  }

/// Converts representation from Q^{con}_n polynomials to simple power series
/*static*/ bool ConvQcon2PowerSeries(const std::vector<double>& C_con, std::vector<double>& A)
  {
  const int nmax = (int)C_con.size() - 1;

  std::vector< std::vector<double> > T;
  if (!CalcQconCoeffs(nmax, T))
    {
    assert(false);
    return false;
    }

  A.resize(nmax + 1);

  for (int m = 0; m <= nmax; m++)
    {
    double sum = 0;
    for (int n = m; n <= nmax; n++)
      sum += T[n][m] * C_con[n];

    A[m] = sum;
    }

  return true;
  }

/*static*/ bool ConvertToJacobi(std::vector<double>&surface_params)
  {
  std::vector<double> coeffs, jacobi_coeffs;

  coeffs.resize(POLY_RPARAMS_MAX_PARAMS);
  for (int j = 0; j < POLY_RPARAMS_MAX_PARAMS; j++)
    coeffs[j] = surface_params[POLY_RPARAMS_FIRST_ID + j];

  ConvPowerSeries2Qcon(coeffs, jacobi_coeffs);

  for (int j = 0; j < POLY_RPARAMS_MAX_PARAMS; j++)
    surface_params[POLY_RPARAMS_FIRST_ID + j] = jacobi_coeffs[j];

  return true;
  }


bool PolySurface::IsFlat() const {return m_isFlat; };

/// Calculates max(radius^2) for which the surface function is defined
double PolySurface::CalcMaxR2() const
{
  const double conic_factor = (1 + surf_data.conic) * SQR(surf_data.curvature);
  const double max_r2 = surf_data.conic > -1.0 ? 1.0 / conic_factor : double(FLT_MAX);
  return max_r2;
}


void PolySurface::kernel1D_UpdateInterParams(const IntersectionParameters& inters_params, int a_dummyIters)
{
  for (int i = 0; i < a_dummyIters; i++)
  {
    m_hitEps     = inters_params.hitEps;
    m_minT       = inters_params.minT;
    m_flags      = inters_params.flags;

    m_useRebasePoly = inters_params.rebasePoly;
    m_useFloat      = inters_params.useFloat;
  }
}

void PolySurface::kernel1D_UpdateSurface(const double* surface_params, uint32_t used_rparams)
{
  surf_data.conic     = surface_params[POLY_CONIC_ID];
  surf_data.curvature = surface_params[POLY_CURVATURE_ID];
  surf_data.rnorm     = surface_params[POLY_RNORM_ID];
  surf_data.aperture  = surface_params[POLY_APERTURE_ID];

  aperture    = surface_params[POLY_APERTURE_ID];
  m_rnormInv  = 1.0 / surface_params[POLY_RNORM_ID];
  npoly       = int(used_rparams);

  max_radius = sqrt(CalcMaxR2());
  center_z   = 0;

  bool flat = abs(surf_data.curvature) < 1e-10; // slicer eats expression otherwise
  m_isFlat = flat;

  for (uint32_t i = 0; i < used_rparams; i++)
  {
    surf_data.rparams[i] = surface_params[POLY_RPARAMS_FIRST_ID + i];
    m_isFlat = m_isFlat && (abs(surf_data.rparams[i]) < 1e-10);
  }
}

#ifdef USE_REBASED_POLYNOMIAL
/// Prepare coefficients of the "re-based" polynomial
void PolySurface::kernel1D_CalcRebased(double rho_max, int n_centres)
  {
  coeffs_rebased.resize(n_centres);
  coeffs_rebased4.resize(n_centres);
  rho2_base_step  = rho_max / (n_centres - 1);
  rho2_base_step4 = (float)rho2_base_step;

  for (int i = 0; i < n_centres; i++)
    {
    coeffs_rebased[i].resize(npoly);
    coeffs_rebased4[i].resize(npoly);

    // Central point of the (Taylor) expansion
    const double rho2_base = i * rho2_base_step;

    for (int k = 0; k < npoly; k++)
      {
      double sum = surf_data.rparams[k], prod = 1;
      for (int l = 1; l < npoly - k; l++)
        {
        prod *= rho2_base * double(l + k) / double(l);
        sum += surf_data.rparams[l + k] * prod;
        }
      coeffs_rebased[i][k] = sum;
      coeffs_rebased4[i][k] = (float)sum;
      }
    }
  }
#endif

/// Calculates the surface function and its derivative, i.e. f(r^2) and df/d(r^2)
/// @return true if the surface is defined in the point, false otherwise. In the latter case
///              the output values are calculated for r2=(max_radius)^2.
/// @attention The surface section is z=f(x^2+y^2)+zpos, 'zpos' is not used in this method!
bool PolySurface::CalcSurfFuncAndDeriv(double r2, double& f, double& df) const
{
#if 1
  const double r2_max = SQR(min(max_radius, aperture));
#else
  const double r2_max = SQR(max_radius);
#endif
  double r2_excess = 0;
  bool clipped = false;
  if (r2 > r2_max)
  {
    r2_excess = r2 - r2_max;
    r2 = r2_max - double(EPSILON);
    clipped = true;
  }

  const double conic_factor = (1 + surf_data.conic) * SQR(surf_data.curvature);
  const double root = sqrt(1 - conic_factor * r2);
  const double conic = surf_data.curvature / (1 + root);

  f  = r2 * conic;
  df = 0.5 * surf_data.curvature / root;

  const double rho2_factor = SQR(m_rnormInv);
  const double rho2 = r2 * rho2_factor;

  // The polynomial term in f() is rho^2 * (c[0]*rho^2 + c[1]*rho^4 + c[2]*rho^6 + ...)
#ifdef USE_REBASED_POLYNOMIAL
  if (m_useRebasePoly)
    {
    // use re-based polynomials (powers not of rho^2 but of rho^2-rho0^2 where rho0 is the center=base point)
    const int nbases = (int)coeffs_rebased.size();
    if (m_useFloat)
      {
      // float precision for rebased polynomials

      // The nearest rho2_base
      float cr = float(rho2) / rho2_base_step4;
      int ir = int(cr);
      if (abs(cr - (ir + 1)) < abs(cr - ir))
        ir++;

      if (ir < 0)
        ir = 0;
      else if (ir > nbases - 1)
        ir = nbases - 1;

      const float rho2_base = ir * rho2_base_step4;
      const float delta = float(rho2) - rho2_base;

      const std::vector<float>& coeffs = coeffs_rebased4[ir];

      // Calculate P=c[0] + c[1]*rho^2 + c[2]*rho^4 + ... using re-based polynomial:
      // P = b[0]+b[1]*delta+... where b[] are the "rebased coefficients" and delta=rho^2-rho_base^2
      float pp = 1, poly = 0, dpoly = 0;
      for (uint32_t j = 2; j < npoly; j++)
        {
        const float tmp = coeffs[j] * pp;
        poly += tmp;
        dpoly += j * tmp;
        pp *= delta;
        }
      poly = coeffs[0] + delta * (coeffs[1] + poly * delta);
      dpoly = coeffs[1] + dpoly * delta;

      // Now transform to poly*rho^4 and its derivative (w.r.t. rho^2)
      dpoly = float(rho2) * (2 * poly + float(rho2) * dpoly);
      poly *= float(rho2 * rho2);

      f += poly;
      df += dpoly * float(rho2_factor); // d()/d(rho^2) --> d()/d(r^2)
      }
    else
      {
      // 'double' precision for rebased for polynomials

      // The nearest rho2_base
      double cr = rho2 / rho2_base_step;
      int ir = int(cr);
      if (abs(cr - (ir + 1)) < abs(cr - ir))
        ir++;

      if (ir < 0)
        ir = 0;
      else if (ir > nbases - 1)
        ir = nbases - 1;

      const double rho2_base = ir * rho2_base_step;
      const double delta = rho2 - rho2_base;

      const std::vector<double>& coeffs = coeffs_rebased[ir];

      // Calculate P=c[0] + c[1]*rho^2 + c[2]*rho^4 + ... using re-based polynomial:
      // P = b[0]+b[1]*delta+... where b[] are the "rebased coefficients" and delta=rho^2-rho_base^2

      double pp = 1, poly = 0, dpoly = 0;
      for (uint32_t j = 2; j < npoly; j++)
        {
        const double tmp = coeffs[j] * pp;
        poly += tmp;
        dpoly += j * tmp;
        pp *= delta;
        }
      poly = coeffs[0] + delta * (coeffs[1] + poly * delta);
      dpoly = coeffs[1] + dpoly * delta;

      // Now transform to poly*rho^4 and its derivative (w.r.t. rho^2)
      dpoly = rho2 * (2 * poly + rho2 * dpoly);
      poly *= rho2 * rho2;

      f += poly;
      df += dpoly * rho2_factor; // d()/d(rho^2) --> d()/d(r^2)
      }
    }
  else
#endif
    {
    // standard polynomials
    if (m_jacobi)
      {
      // Use rho^4*Q^{con}_n(rho^2)

      const double xi = 2 * rho2 - 1;

      // Starting point for the recursion of the Jacobi polynomials P^{0,4}_n(xi) where xi=2*rho^2-1

      double P_prevprev = 1;           // P^{0,4}_{n-2}(xi)
      double P_prev     = -2 + 3 * xi; // P^{0,4}_{n-1}(xi)
      double P;                        // P^{0,4}_{n}(xi)

      double dP_prevprev = 0;           // P^{0,4}_{n-2}'(xi)
      double dP_prev     = 3;           // P^{0,4}_{n-1}'(xi)
      double dP;                        // P^{0,4}_{n}'(xi)

      double sum = 0, dsum = 0;
      sum  = surf_data.rparams[0] * P_prevprev  + surf_data.rparams[1] * P_prev;
      dsum = surf_data.rparams[0] * dP_prevprev + surf_data.rparams[1] * dP_prev;
      for (uint32_t n = 2; n < npoly; n++)
        {
        const double rev_a = 1.0 / (n * (n + 4) * (n + 1));
        const double b = (n + 2) * (n + 1);
        const double c = (2 * n + 3) * (b * xi - 4);
        const double d = (2 * n + 3) * b;
        const double e = (n - 1) * (n + 3) * (n + 2);

        P  = rev_a * (c * P_prev - e * P_prevprev);
        dP = rev_a * (d * P_prev + c * dP_prev - e * dP_prevprev);

        const double &coeff = surf_data.rparams[n];

        sum  += coeff * P;
        dsum += coeff * dP;

        // Promote to n+1
        P_prevprev  = P_prev;
        P_prev      = P;
        dP_prevprev = dP_prev;
        dP_prev     = dP;
        }
      f  += rho2 * rho2 * sum;
      df += 2 * rho2 * (sum + rho2 * dsum) * rho2_factor;
      }
    else
      {
      // simple power series
      double pp = rho2_factor * rho2, sum = 0;
      for (uint32_t j = 0; j < npoly; j++)
        {
        const double tmp = surf_data.rparams[j] * pp;
        sum += tmp;
        df += (2 + j) * tmp;
        pp *= rho2;
        }
      f += sum * r2;
      }
    }

  if (clipped)
  {
    // continue f(r^2) by f(r2_max) + f'(r2_max) * (r^2-r2_max)
    f += df * r2_excess;
  }

  return !clipped;
}

/// Calculates the surface function and its derivativeS, i.e. f(r^2), f'(r^2) and f''(r^2)
/// @return true if the surface is defined in the point, false otherwise. In the latter case
///              the output values are calculated for r2=(max_radius)^2.
/// @attention The surface section is z=f(x^2+y^2)+zpos, 'zpos' is not used in this method!
bool PolySurface::CalcSurfFuncAndDerivs(double r2, double& f, double& df, double& d2f) const
  {
#if 1
  const double r2_max = SQR(min(max_radius, aperture));
#else
  const double r2_max = SQR(max_radius);
#endif
  double r2_excess = 0;
  bool clipped = false;
  if (r2 > r2_max)
    {
    r2_excess = r2 - r2_max;
    r2 = r2_max - double(EPSILON);
    clipped = true;
    }

  const double conic_factor = (1 + surf_data.conic) * SQR(surf_data.curvature);
  const double s = 1 - conic_factor * r2;
  const double root = sqrt(s);
  const double conic = surf_data.curvature / (1 + root);

  f = r2 * conic;
  df = 0.5 * surf_data.curvature / root;
  d2f = 0.25 * (1 + surf_data.conic) * surf_data.curvature * surf_data.curvature * surf_data.curvature
    / (s * root);

  const double rho2_factor = SQR(m_rnormInv);
  const double rho2 = r2 * rho2_factor;

#ifdef USE_REBASED_POLYNOMIAL
  if (m_useRebasePoly)
    {
    // use re-based polynomials (powers not of rho^2 but of rho^2-rho0^2 where rho0 is the center=base point)
    const int nbases = (int)coeffs_rebased.size();

    // The nearest rho2_base
    double cr = rho2 / rho2_base_step;
    int ir = int(cr);
    if (abs(cr - (ir + 1)) < abs(cr - ir))
      ir++;

    if (ir < 0)
      ir = 0;
    else if (ir > nbases - 1)
      ir = nbases - 1;

    const double rho2_base = ir * rho2_base_step;
    const double delta = rho2 - rho2_base;

    const std::vector<double>& coeffs = coeffs_rebased[ir];

    // Calculate P=c[0] + c[1]*rho^2 + c[2]*rho^4 + ... using re-based polynomial:
    // P = b[0]+b[1]*delta+... where b[] are the "rebased coefficients" and delta=rho^2-rho_base^2
    double pp = 1, poly = 0, dpoly = 0, d2poly = 0;
    for (uint32_t j = 2; j < npoly; j++)
      {
      const double tmp = coeffs[j] * pp;
      poly += tmp;
      dpoly += j * tmp;
      d2poly += (j - 1) * j * tmp;
      pp *= delta;
      }
    poly = coeffs[0] + delta * (coeffs[1] + poly * delta);
    dpoly = coeffs[1] + dpoly * delta;

    // Now transform to poly*rho^4 and its derivative (w.r.t. rho^2)
    d2poly = 2 * poly + rho2 * (4 * dpoly + rho2 * d2poly);
    dpoly = rho2 * (2 * poly + rho2 * dpoly);
    poly *= rho2 * rho2;

    f += poly;
    df += dpoly * rho2_factor;         // d()/d(rho^2) --> d()/d(r^2)
    d2f += d2poly * rho2_factor * rho2_factor; // d()/d(rho^2) --> d()/d(r^2)
    }
  else
#endif
    {
    // standard polynomials
    if (m_jacobi)
      {
      // Use rho^4*Q^{con}_n(rho^2)

      const double xi = 2 * rho2 - 1;

      // Starting point for the recursion of the Jacobi polynomials P^{0,4}_n(xi) where xi=2*rho^2-1

      double P_prevprev = 1;       // P^{0,4}_{n-2}(xi)
      double P_prev = -2 + 3 * xi; // P^{0,4}_{n-1}(xi)
      double P;                    // P^{0,4}_{n}(xi)

      double dP_prevprev = 0;      // P^{0,4}_{n-2}'(xi)
      double dP_prev = 3;          // P^{0,4}_{n-1}'(xi)
      double dP;                   // P^{0,4}_{n}'(xi)

      double d2P_prevprev = 0;     // P^{0,4}_{n-2}"(xi)
      double d2P_prev = 0;         // P^{0,4}_{n-1}"(xi)
      double d2P;                  // P^{0,4}_{n}"(xi)

      double sum = 0, dsum = 0, d2sum = 0;

      sum   = surf_data.rparams[0] * P_prevprev   + surf_data.rparams[1] * P_prev;
      dsum  = surf_data.rparams[0] * dP_prevprev  + surf_data.rparams[1] * dP_prev;
      d2sum = surf_data.rparams[0] * d2P_prevprev + surf_data.rparams[1] * d2P_prev;
      for (uint32_t n = 2; n < npoly; n++)
        {
        const double rev_a = 1.0 / (n * (n + 4) * (n + 1));
        const double b = (n + 2) * (n + 1);
        const double c = (2 * n + 3) * (b * xi - 4);
        const double d = (2 * n + 3) * b;
        const double e = (n - 1) * (n + 3) * (n + 2);

        P = rev_a * (c * P_prev - e * P_prevprev);
        dP = rev_a * (d * P_prev + c * dP_prev - e * dP_prevprev);
        d2P = rev_a * (2 * d * dP_prev + c * d2P_prev - e * d2P_prevprev);

        const double& coeff = surf_data.rparams[n];

        sum   += coeff * P;
        dsum  += coeff * dP;
        d2sum += coeff * d2P;

        // Promote to n+1
        P_prevprev = P_prev;
        P_prev = P;
        dP_prevprev = dP_prev;
        dP_prev = dP;
        d2P_prevprev = d2P_prev;
        d2P_prev = d2P;
        }
      f   += rho2 * rho2 * sum;
      df  += 2 * rho2 * (sum + rho2 * dsum) * rho2_factor;
      d2f += 2 * (sum + 2 * rho2 * (2 * dsum + rho2 * d2sum)) * rho2_factor * rho2_factor;
      }
    else
      {
      // simple power series
      double pp = rho2_factor * rho2, sum_f = 0, sum_d2f = 0;
      // The polynomial term in f() is rho^2 * (c[0]*rho^2 + c[1]*rho^4 + c[2]*rho^6 + ...)
      for (uint32_t j = 0; j < npoly; j++)
        {
        const double tmp = surf_data.rparams[j] * pp;
        sum_f += tmp;
        df += (2 + j) * tmp;
        sum_d2f += (1 + j) * (2 + j) * tmp;
        pp *= rho2;
        }
      f += sum_f * r2;
      if (r2 == 0)
        d2f += 2.0 * surf_data.rparams[0] * rho2_factor * rho2_factor;
      else
        d2f += sum_d2f / r2;
      }
    }

  if (clipped)
  {
    // continue f(r^2) by f(r2_max) + f'(r2_max) * (r^2-r2_max)
    f += df * r2_excess;
  }

  return !clipped;
}

void PolySurface::GetParam(std::vector<double>& param) const
  {
  param.resize(GetNumParam(), 0.0);
  std::fill(param.begin(), param.end(), 0.0);

  param[POLY_APERTURE_ID]  = aperture;
  param[POLY_CONIC_ID]     = surf_data.conic;
  param[POLY_CURVATURE_ID] = surf_data.curvature;
  param[POLY_RNORM_ID]     = surf_data.rnorm;

  for (uint32_t i = 0; i < npoly; i++)
    param[POLY_RPARAMS_FIRST_ID + i] = surf_data.rparams[i];
  }

/// Surface-specific ID of the given parameter
int PolySurface::GetParamID(int i) const
  {
  return i;
  }


/// The number of actually used parameters in 'surface_params'.
int PolySurface::GetNumParam() const
  {
  return POLY_RPARAMS_FIRST_ID + npoly;
  }

/// The value of 'used_rparams' passed to UpdateSurface()
int PolySurface::GetUsedRparams() const
  {
  return npoly;
  }


void PolySurface::UpdateSurface(const double* surface_params, const uint32_t used_rparams)
{
  m_newton_iterations = 0;
  m_raymarch_iterations = 0;
  auto start = std::chrono::high_resolution_clock::now();

  const uint32_t npt    = 101; // resolution to calculate zmin, zmax
  const uint32_t iters1 = min(used_rparams, POLY_RPARAMS_MAX_PARAMS);

  kernel1D_UpdateSurface(surface_params, iters1);

#ifdef USE_REBASED_POLYNOMIAL
  if (m_useRebasePoly)
    kernel1D_CalcRebased(1.1, 11);
#endif

  kernel1D_CalcZminZmax(npt, 1.0 / (npt - 1));

#ifdef NEW_START_POINT4NEWTON
  // Calculate bounds of f'() over the system of intervals. At each level, [0,rmax]
  // is subdivided into n intervals of equal length and min(f'), max(f') over each
  // is calculated. On the next level, the number of intervals increases 'reduction'
  // times (and length decreases).
  const int nlevels = 1, reduction = 4, nfirst = 4;

  auto start1 = std::chrono::steady_clock::now();
  CalcBoundsOfFuncAndDeriv(nlevels, nfirst, reduction);
  auto current1 = std::chrono::steady_clock::now();

  std::chrono::duration<double> elapsed1 = current1 - start1;
  ///  printf("CalcBoundsOfDeriv() took %10.3f sec\n", elapsed1.count());
#endif

  update_time = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()) * 0.001f;
}

void PolySurface::TraceRays(const IntersectionParameters& inters_params,
                            const double* zpos, const double* surface_params,
                            const uint32_t used_rparams,
                            const double* eta_inv,
                            const double* ray_pos, const double* ray_dir, const uint8_t* ray_mask,
                            double* hit_pos, double* hit_dir, uint8_t* hit_mask, double* hit_t, double* hit_normal,
                            HitPointData* hpd,
                            uint32_t n_rays_x, uint32_t n_rays_y, uint32_t n_batches, uint32_t z_pos_size)
{
  auto start = std::chrono::high_resolution_clock::now();

  m_newton_iterations = 0;
  m_raymarch_iterations = 0;
  kernel1D_UpdateInterParams(inters_params, 1);
  kernel1D_TraceRays(zpos, eta_inv, ray_pos, ray_dir, ray_mask, hit_pos, hit_dir, hit_mask, hit_t, hit_normal, hpd, n_rays_x*n_rays_y*n_batches, n_rays_x*n_rays_y, z_pos_size);

  trace_time = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()) * 0.001f;
}
