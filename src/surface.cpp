#include "surface.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using LiteMath::double2;
using LiteMath::make_double3;
using LiteMath::SQR;
using LiteMath::dot;
using LiteMath::length;

using std::min;
using std::max;
using std::isnan;
using std::abs;

/// Aperture is treated as a hard cut; the surface function is continued
/// linearly past it (see CalcSurfFuncAndDeriv in poly_surface.cpp).
const double excess_of_aperture = 1.00;

/// New raymarching which uses quadratic approx for f(t) while linear for r^2(t)
#define NEWRAYMARCH

///#define PRINT_RAYMARCH_NITER


#ifdef NEW_START_POINT4NEWTON
/// Calculates boundaries of derivative f'()
/// @param[in] df_min_arr - array of min(f') over uniform cells of size dr2
/// @param[in] df_max_arr - array of max(f') over uniform cells of size dr2
/// @param[in] dr2        - cell size (delta of r^2)
/// @param[in] b          - origin.x * v.x + origin.y * v.y, i.e. ray's d(r^2)/dt at t=0
/// @param[in] r2         - current value of r^2
/// @param[out] df_min    - low  bound of f'() in [r2_left, r2_right]
/// @param[out] df_max    - high bound of f'() in [r2_left, r2_right]
/// @param[out] r2_left   - low  bound of r^2 for which the estimate of f'() holds
/// @param[out] r2_right  - high bound of r^2 for which the estimate of f'() holds
int CalcDerivBounds(const std::vector<double>& df_min_arr,
                    const std::vector<double>& df_max_arr,
                    const double dr2, const double b,
                    double r2,
                    double& df_min, double& df_max,
                    double& r2_left, double& r2_right)
  {
  const int ni = (int)df_min_arr.size();

  int i = min(int(r2 / dr2), ni - 1);
  r2_left = dr2 * i;
  r2_right = dr2 * (i + 1);

  // Bounds of f'() in this interval of r^2
  df_min = df_min_arr[i];
  df_max = df_max_arr[i];

  if (r2 > r2_right - dr2 * 0.001 && b > Surface::EPSILON) // d(r^2)/dt at t0 is =b
    {
    i = min(i + 1, ni - 1);
    // Bounds of f'() in the TWO intervals of r^2
    df_min = min(df_min, df_min_arr[i]);
    df_max = max(df_max, df_max_arr[i]);

    r2_right = dr2 * (i + 1);
    }
  if (r2 < r2_left + dr2 * 0.001 && b < -Surface::EPSILON) // d(r^2)/dt at t0 is =b
    {
    i = max(i - 1, 0);
    // Bounds of f'() in the TWO intervals of r^2
    df_min = min(df_min, df_min_arr[i]);
    df_max = max(df_max, df_max_arr[i]);

    r2_left = dr2 * i;
    }
  return i;
  }

/// Calculates boundaries of f()
/// @param[in] f_min_arr  - array of min(f) over uniform cells of size dr2
/// @param[in] f_max_arr  - array of max(f) over uniform cells of size dr2
/// @param[in] dr2        - cell size (delta of r^2)
/// @param[in] b          - origin.x * v.x + origin.y * v.y, i.e. ray's d(r^2)/dt at t=0
/// @param[in] r2         - current value of r^2
/// @param[out] f_min     - low  bound of f() in [r2_left, r2_right]
/// @param[out] f_max     - high bound of f() in [r2_left, r2_right]
/// @param[out] r2_left   - low  bound of r^2 for which the estimate of f'() holds
/// @param[out] r2_right  - high bound of r^2 for which the estimate of f'() holds
double2 CalcFuncBounds(const std::vector<double>& f_min_arr,
                       const std::vector<double>& f_max_arr,
                       const double dr2, const double b,
                       const double2& r2_range,
                       double& r2_left, double& r2_right)
  {
  const int ni = (int)f_min_arr.size();

  double c_first = r2_range[0] / dr2;
  double c_second = r2_range[1] / dr2;
  int i_first = min(int(c_first), ni - 1);
  int i_second = min(int(c_second), ni - 1);

  if (c_second > 0.999 && b > Surface::EPSILON && i_second < ni - 1) // d(r^2)/dt at t0 is =b
    i_second++;
  if (c_first < 0.001 && b < -Surface::EPSILON && i_first > 0) // d(r^2)/dt at t0 is =b
    i_first--;


  r2_left = dr2 * i_first;
  r2_right = dr2 * (i_second + 1);

  double2 f_range(+1e10, -1e10);
  for (int i = i_first; i <= i_second; i++)
    {
    f_range[0] = min(f_range[0], f_min_arr[i]);
    f_range[1] = max(f_range[1], f_max_arr[i]);
    }

  return f_range;
  }

double Table1(double b_minus_s, double st_aster, double a)
  {
  if (st_aster <= 0)
    {
    double discr = b_minus_s * b_minus_s - 4 * a * st_aster;
    if (discr < 0)
      return -1;
    else
      {
      double d = sqrt(discr);
      return (-b_minus_s + d/*sqrt(discr)*/) / (2 * a);
      }
    }
  else
    {
    if (b_minus_s > 0)
      return -1;

    double discr = b_minus_s * b_minus_s - 4 * a * st_aster;
    if (discr < 0)
      return -1;
    else
      return (-b_minus_s - sqrt(discr)) / (2 * a);
    }
  }

/// Constraint = ray path length until r^2(t) leaves [r2_left,r2_right]
/// @param[in] a    = (1-v.z^2) where v is the direction of the ray
/// @param[in] b    = origin.x * v.x + origin.y * v.y
/// @param[in] low  = low   bound of r^2-r0^2 = r2_left  - r2, <= 0
/// @param[in] high = upper bound of r^2-r0^2 = r2_right - r2, >= 0
double Constraint(const double a,
                  const double b,
                  const double low,
                  const double high)
  {
  if (fabs(a) < 1e-12)
    {
    assert(fabs(b) < 1e-10);
    return 1e10;
    }

  if (b <= 0)
    {
    // The root of a*t^2+b*t = low such that t<=t_extr
    double discr = b * b + 4 * a * low;
    if (discr >= 0) // if misses the low bound, constraint is by the high bound!
      {
      double t = (-b - sqrt(discr)) / (2 * a);
      assert(t >= 0 && (a * t + b) * t <= high);
      return t;
      }
    }

  // The root of a*t^2+b*t = high such that t>=0
  double discr = b * b + 4 * a * high;
  if (discr >= 0)
    {
    double t = max((-b + sqrt(discr)) / (2 * a), 0.0);
    assert(t >= 0 && (a * t + b) * t >= low);
    return t;
    }
  else if (discr > -1e-12)
    {
    double t = max(-b / (2 * a), 0.0);
    assert(t >= 0 && (a * t + b) * t >= low);
    return t;
    }
  else
    {
    assert(false);
    return -1;
    }
  }

/// Whether a*t^2+b*t is inside [low,high] if for t=0 it is inside
/// @param[in] a    = (1-v.z^2) where v is the direction of the ray
/// @param[in] b    = origin.x * v.x + origin.y * v.y
/// @param[in] low  = low   bound of r^2-r0^2 = r2_left  - r2, <= 0
/// @param[in] high = upper bound of r^2-r0^2 = r2_right - r2, >= 0
int InsideBounds(const double t,
                 const double a,
                 const double b,
                 const double low,
                 const double high)
  {
  const double parabola = (a * t + b) * t;
  if (b <= 0)
    {
    // t=0 is to the left of the parabola's minimum. Check the lower bound
    if (t * (2 * a) <= -b) // t <= t_extremum
      {
      if (parabola <= low - 1e-10)
        return -1;
      else
        return 0;
      }
    else
      {
      if (-b * b < (4 * a) * (low - 1e-10)) // minimum is below the low bound?
        return -1;
      }
    }

  // Parabola is ascending, so it is enough to check the upper bound
  if (parabola > high + 1e-10)
    return +1;

  return 0;
  }

/// The range of the ray's r^2(t)=a*t^2+b*t+c when 't' changes in the given range
double2 CalcR2Range(const double a, const double b, const double c, const double2& t_range)
  {
  const double t_extr = -b / (2 * a);

  double2 r2_range;
  if (t_range[1] <= t_extr || t_range[0] >= t_extr)
    {
    // monotone
    r2_range[0] = c + (a * t_range[0] + b) * t_range[0];
    r2_range[1] = c + (a * t_range[1] + b) * t_range[1];

    if (r2_range[1] < r2_range[0])
      std::swap(r2_range[0], r2_range[1]);
    }
  else
    {
    // extremum (always minimum!) is inside
    r2_range[0] = c + 0.5 * b * t_extr;
    if (t_extr - t_range[0] >= t_range[1] - t_extr)
      r2_range[1] = c + (a * t_range[0] + b) * t_range[0];
    else
      r2_range[1] = c + (a * t_range[1] + b) * t_range[1];
    }
  return r2_range;
  }

/// @param[in] slope   = slope of the red line (=v.z/f')
/// @param[in] t_aster = the value of 't' where the red line is 0 (=F0/f').
/// @param[in] low     = low   bound of r^2-r0^2 = r2_left  - r2, <= 0
/// @param[in] high    = upper bound of r^2-r0^2 = r2_right - r2, >= 0
/// @param[in] a    = (1-v.z^2) where v is the direction of the ray
/// @param[in] b    = origin.x * v.x + origin.y * v.y
double MinFreePathLength(const double2& slopes,
                         const double t_aster,
                         const double low,
                         const double high,
                         const double a,
                         const double b)
  {
  // Calculate intersection of parabola and red line

  double t;
  double t1 = Table1(b - slopes[0], slopes[0] * t_aster, a);
  double t2 = Table1(b - slopes[1], slopes[1] * t_aster, a);

  t = t1;
  if (t1 < 0 || (t2 < t1 && t2 >= 0 && t2 < 1e5))
    t = t2;

  if (t >= 0)
    {
    const int res = InsideBounds(t, a, b, low, high);
    if (res == 0)
      return t; // ray intersection is BEFORE r^2(t) leaves allowed bounds
    }

  // The ray leaves the allowed bounds first
  t = Constraint(a, b, low, high);
  const int rc = InsideBounds(t, a, b, low, high);
  if (rc != 0)
    assert(false);

  return t;
  }


/// Calculates free ray path from current ray position (at t=0).
/// @param[in] df_min_arr - array of min(f') over uniform cells of size dr2
/// @param[in] df_max_arr - array of max(f') over uniform cells of size dr2
/// @param[in] r2max      - upper bound of the interval of r^2 which @c df_min_arr
///                         had been calculated
/// @param[in] direction  - ray direction
/// @param[in] r2         - current value of r^2
/// @param[in] a          - (1-direction.z^2)
/// @param[in] b          - origin.x * v.x + origin.y * v.y, i.e. ray's d(r^2)/dt at t=0
/// @param[in] F_t0       - the value of the function F(t=0). The function
///                         F(t)=z_pos+f(r^2(t)) - z(t) where r^2(t)=x^2(t)+y^2(t)
///                         and (x(t),y(t),z(t)) is the ray position at time t.
///                         F=0 when the ray intersects the surface.
double RayStep(const std::vector<double>& df_min_arr,
               const std::vector<double>& df_max_arr,
               const double r2max,
               const Vect3d& direction,
               const double r2,
               const double a,
               const double b,
               const double F_t0)
  {
  const int ni = (int)df_max_arr.size();
  const double dr2 = r2max / ni;

  double df_min, df_max, r2_left, r2_right;


  // In the interval r2_left<=r^2<=r2_right, f' is between df_min and df_max
  // If r2 is close to the cell boundary, the interval r2_left<=r^2<=r2_right
  // can include TWO cells
  /*const int i = */CalcDerivBounds(df_min_arr, df_max_arr, dr2, b, r2, df_min, df_max, r2_left, r2_right);

  const double t_aster = F_t0 / direction.z;
  const double2 slopes(direction.z / df_min, direction.z / df_max);

  const double t0 = MinFreePathLength(slopes, t_aster, r2_left - r2, r2_right - r2, a, b);
  return t0;
  }

/// Calculates free ray path from current ray position (at t=0).
/// @param[in] df_min_arr_arr - arrays of min(f') over uniform cells of size dr2,
///                             for ALL LEVELS of sudvivision (level = outer array index)
/// @param[in] df_max_arr_arr - arrays of max(f') over uniform cells of size dr2
///                             for ALL LEVELS of sudvivision (level = outer array index)
/// @param[in] r2max      - upper bound of the interval of r^2 which @c df_min_arr
///                         had been calculated
/// @param[in] direction  - ray direction
/// @param[in] r2         - current value of r^2
/// @param[in] a          - (1-direction.z^2)
/// @param[in] b          - origin.x * v.x + origin.y * v.y, i.e. ray's d(r^2)/dt at t=0
/// @param[in] F_t0       - the value of the function F(t=0). The function
///                         F(t)=z_pos+f(r^2(t)) - z(t) where r^2(t)=x^2(t)+y^2(t)
///                         and (x(t),y(t),z(t)) is the ray position at time t.
///                         F=0 when the ray intersects the surface.
double RayStep(const std::vector<std::vector<double> >& df_min_arr_arr,
               const std::vector<std::vector<double> >& df_max_arr_arr,
               const double r2max,
               const Vect3d& direction,
               const double r2,
               const double a,
               const double b,
               const double F_t0)
  {
  double t0 = 0, t0_;

  const int nlevels = (int)df_min_arr_arr.size();

  for (int level = 0; level < nlevels; level++)
    {
    t0_ = RayStep(df_min_arr_arr[level], df_max_arr_arr[level], r2max,
      direction, r2, a, b, F_t0);

    if (t0_ > t0)
      {
      t0 = t0_;
#if 1
      if (t0 > 0.02)///5e-2)
        break;
#endif
      }
    }
  return t0;
  }
#endif

#ifdef NEWRAYMARCH
/// Finds the smallest positive root of the quadratic eq. ax^2+bx+c=0
/// @return The smallest positive root or infinity if such a root is absent
double Quadratic(double a, double b, double c)
  {
  if (a == 0)
    return -c / b;

  double d = b * b - 4 * a * c;
  if (d < -1e-10)
    return 1e10;
  else if (d <= 0)
    return -b / (2 * a);

  d = sqrt(d);

  const double s = 1 / (2 * a);
  double t1 = s * (-b - d), t2 = s * (-b + d);

  if (t2 < t1)
    std::swap(t1, t2);

  if (t1 >= -1e-10)
    return std::max(0.0, t1);
  else if (t2 >= -1e-10)
    return std::max(0.0, t2);
  else
    return 1e10;
  }


/// Calculates free ray path from current ray position (at t=0).
/// @param[in] r2         - current value of r^2
/// @param[in] a          - (1-direction.z^2)
/// @param[in] b          - origin.x * v.x + origin.y * v.y, i.e. ray's d(r^2)/dt at t=0
/// @param[in] F0         - the value of the function F(t=0). The function
///                         F(t)=z_pos+f(r^2(t)) - z(t) where r^2(t)=x^2(t)+y^2(t)
///                         and (x(t),y(t),z(t)) is the ray position at time t.
///                         F=0 when the ray intersects the surface.
/// @param[in] df         - the value of f'(r^2); derivative is in r^2 as usual
/// @param[in] d2f_min_arr - array of min(f") over uniform cells of size dr2
/// @param[in] d2f_max_arr - array of max(f") over uniform cells of size dr2
/// @param[in] r2max      - upper bound of the interval of r^2 which @c df_min_arr
///                         had been calculated
double RayStep2(const double r2,
                const double vz,
                const double a,
                const double b,
                const double F0,
                const double df,
                const std::vector<double>& d2f_min_arr,
                const std::vector<double>& d2f_max_arr,
                const double r2max)
  {
  const int ni = (int)d2f_max_arr.size();
  const double dr2 = r2max / ni;

  double d2f_min, d2f_max, r2_left, r2_right;

  // In the interval r2_left<=r^2<=r2_right, f" is between d2f_min and d2f_max
  // If r2 is close to the cell boundary, the interval r2_left<=r^2<=r2_right
  // can include TWO cells
  /*const int i = */CalcDerivBounds(d2f_min_arr, d2f_max_arr, dr2, b, r2, d2f_min, d2f_max, r2_left, r2_right);

  // Constraint tc: when the ray first leaves the range r2_left<=r^2<=r2_right
  const double tc = Constraint(a, b, r2_left - r2, r2_right - r2);

  // If the signs of the slope of the two lines bounding our parabola are DIFFERENT,
  // then the squared parabola can be 0, i.e. the min slope for it is 0!
  const double beta_min = b * (b + a * tc) <= 0 ? 0 : std::min(SQR(b), SQR(b + a * tc));
  const double beta_max = std::max(SQR(b), SQR(b + a * tc));

  const double gamma_min = 0.5 * d2f_min * (d2f_min <= 0 ? beta_max : beta_min);
  const double gamma_max = 0.5 * d2f_max * (d2f_max <= 0 ? beta_min : beta_max);

  double t1 = Quadratic(a * df + gamma_min, b * df - vz, F0);
  double t2 = Quadratic(a * df + gamma_max, b * df - vz, F0);

  return std::min(tc, std::min(t1, t2));
  }

/// Calculates free ray path from current ray position (at t=0).
/// @param[in] r2         - current value of r^2
/// @param[in] a          - (1-direction.z^2)
/// @param[in] b          - origin.x * v.x + origin.y * v.y, i.e. ray's d(r^2)/dt at t=0
/// @param[in] F0         - the value of the function F(t=0). The function
///                         F(t)=z_pos+f(r^2(t)) - z(t) where r^2(t)=x^2(t)+y^2(t)
///                         and (x(t),y(t),z(t)) is the ray position at time t.
///                         F=0 when the ray intersects the surface.
/// @param[in] df         - the value of f'(r^2); derivative is in r^2 as usual
/// @param[in] d2f_min_arr_arr - arrays of min(f") over uniform cells of size dr2,
///                              for ALL LEVELS of sudvivision (level = outer array index)
/// @param[in] d2f_max_arr_arr - arrays of max(f") over uniform cells of size dr2
///                              for ALL LEVELS of sudvivision (level = outer array index)
/// @param[in] r2max      - upper bound of the interval of r^2 which @c df_min_arr
///                         had been calculated
double RayStep2(const double r2,
                const double vz,
                const double a,
                const double b,
                const double F0,
                const double df,
                const std::vector<std::vector<double> >& d2f_min_arr_arr,
                const std::vector<std::vector<double> >& d2f_max_arr_arr,
                const double r2max)
  {
  double t0 = 0, t0_;

  const int nlevels = (int)d2f_min_arr_arr.size();

  for (int level = 0; level < nlevels; level++)
    {
    t0_ = RayStep2(r2, vz, a, b, F0, df, d2f_min_arr_arr[level], d2f_max_arr_arr[level], r2max);

    if (t0_ > t0)
      {
      t0 = t0_;
#if 1
      if (t0 > 0.02)///5e-2)
        break;
#endif
      }
    }
  return t0;
  }
#endif

Surface::Surface()
{
  aperture = max_radius = zmin = zmax = 0;
  center_z = 0;
  }

HitPointData Surface::Intersection(double zpos, Vect3d& origin, const Vect3d& dir, bool& ray_valid) const
{
  uint32_t newton_iters = 0;
  uint32_t raymarch_iters = 0;
  return Intersection(zpos, origin, dir, ray_valid, newton_iters, raymarch_iters);
}

HitPointData Surface::Intersection(double zpos, Vect3d& origin, const Vect3d& dir, bool& ray_valid,
                                   uint32_t& newton_iters, uint32_t& raymarch_iters) const
{
  HitPointData hpd;
  Reset(hpd); //  For a premature return
  newton_iters = 0;
  raymarch_iters = 0;

  // Move the ray close to the lens
  double t0;
  Point3d origin_new;

  Point3d nearest;
  bool success = MoveRayCloseToLens(zpos, origin, dir, origin_new, t0, nearest, raymarch_iters);
  Point3d guess = origin_new;
  Point3d origin_initial = origin;
  if (!success)
  {
    // Either of the available variants failed (usually the point is always outside of max_radius)
    ray_valid = false;
    origin = make_double3(10, 10, 0);
    hpd.t = 0.0;
    return hpd;
  }

  uint32_t it = 0;
  double t = 0;
  double df = 0, ft = MAX_DISTANCE;
  double delta = 1e10;

  const double a = 2 * (origin_new.x * dir.x + origin_new.y * dir.y);
  const double b = 2 * (dir.x * dir.x + dir.y * dir.y);
  while (abs(delta) > 0.5 * m_hitEps && it < m_maxIterations)
  {
    const Point3d o_new = origin_new + dir * t;
    const double r2 = o_new.x * o_new.x + o_new.y * o_new.y;

    bool good = CalcSurfFuncAndDeriv(r2, ft, df); // if calculation had been done for the CLIPPED r2 to make the formula defined
    // Deviation from the surface F(t):
    ft = (ft + zpos) - o_new.z;
    // F'(t):
    const double dt_ft = df * (a + b * t) - dir.z;

#if 0
    if (ft < -NEWTON_STEP_BOUND * dt_ft) // F/F' < -NEWTON_STEP_BOUND
      delta = -NEWTON_STEP_BOUND;
    else if (ft > NEWTON_STEP_BOUND * dt_ft) // F/F' > NEWTON_STEP_BOUND
      delta = NEWTON_STEP_BOUND;
    else
#endif
      delta = ft / dt_ft;

    t -= delta;
    it++;
  }

  Point3d o_new = origin_new + dir * t;
  double r2 = o_new.x * o_new.x + o_new.y * o_new.y; // =r^2

#if 1 // Ose experiment
  bool valid = r2 <= SQR(excess_of_aperture * aperture) && r2 <= SQR(excess_of_aperture * max_radius);
  if (r2 > SQR(excess_of_aperture * max_radius))
#else
  bool valid = r2 <= SQR(aperture) && r2 <= SQR(max_radius);
  if (r2 > SQR(max_radius))
#endif
  {
    ray_valid = false;
    // For compatibility with tests, hit point must be 0 for invalid rays
    origin.x = origin.y = origin.z = 0;
    newton_iters = 0;
    raymarch_iters = 0;
  }
  else
  {
    CalcSurfFuncAndDerivs(r2, hpd.f, hpd.df, hpd.d2f);
    ft = (hpd.f + zpos) - o_new.z;
    ray_valid = valid && (abs(ft) < m_hitEps) && (t + t0 > m_minT);
    newton_iters = it;
    origin = o_new;

    hpd.dir[0] = dir.x;
    hpd.dir[1] = dir.y;
    hpd.dir[2] = dir.z;

    // hpd.xyz[0] = origin.x;
    // hpd.xyz[1] = origin.y;
    // hpd.xyz[2] = origin.z;

    hpd.t   = t + t0;
    // hpd.r2  = r2;
///    hpd.d2f = Calc2ndDeriv(r2);
  }

  hpd.t = t + t0;

  // debug/visualization output of initial guess vs final intersection
  #if 0
  {
    std::ofstream f("found_points_vs_initial_guess.txt", std::ios_base::app);

    if (valid && abs(origin.y) < aperture )
      f << guess.x << " " << guess.y << " " << guess.z << " - " << origin.x << " " << origin.y << " " << origin.z << " " << std::endl;
    else
      f << guess.x << " " << guess.y << " " << guess.z << " - " << origin.x << " " << origin.y << " " << 1000 << " " << std::endl;

    f.close();
  }

  {
    std::ofstream f("ray_origin_ray_dir.txt", std::ios_base::app);

    if (valid)
      f << origin_initial.x << " " << origin_initial.y << " " << origin_initial.z << " - " << dir.x << " " << dir.y << " " << dir.z << " " << std::endl;

    f.close();
  }
  #endif

  return hpd;
}

HitPointData Surface::IntersectionSafe(double zpos, Vect3d& origin, const Vect3d& dir, bool& ray_valid,
                                       uint32_t& newton_iters, uint32_t& raymarch_iters) const
{
  HitPointData hpd;
  Reset(hpd);
  ray_valid = false;
  newton_iters = 0;
  raymarch_iters = 0;

  const double ox = origin.x, oy = origin.y, oz = origin.z;
  const double dx = dir.x,    dy = dir.y,    dz = dir.z;

  // r^2(t) = A*t^2 + B*t + C
  const double A = dx*dx + dy*dy;
  const double B = 2.0 * (ox*dx + oy*dy);
  const double C = ox*ox + oy*oy;

  const double r2_eval_max = SQR(excess_of_aperture * std::min(aperture, max_radius));
  const double r2_ap_max   = SQR(excess_of_aperture * aperture);
  const double r2_max_max  = SQR(excess_of_aperture * max_radius);

  // F(t) = z_surface(r^2(t)) + zpos - z_ray(t), and its derivative.
  auto evalF = [&](double t, double& ft, double& dt_ft) -> bool
  {
    const double r2 = A*t*t + B*t + C;
    if (!(r2 >= 0.0) || r2 > r2_eval_max)
      return false;

    double sag = 0.0, dsag = 0.0;
    if (!CalcSurfFuncAndDeriv(r2, sag, dsag))
      return false;

    ft    = (sag + zpos) - (oz + dz*t);
    dt_ft = dsag * (2.0*A*t + B) - dz;
    return std::isfinite(ft) && std::isfinite(dt_ft);
  };

  // If the safeguarded bracket cannot be established, or the final residual
  // is not tight enough, use the proven preconditioned Newton routine.
  auto fallback = [&]() -> HitPointData
  {
    origin.x = ox;
    origin.y = oy;
    origin.z = oz;
    const uint32_t safe_newton = newton_iters;
    const uint32_t safe_raymarch = raymarch_iters;
    uint32_t fallback_newton = 0;
    uint32_t fallback_raymarch = 0;
    HitPointData res = Intersection(zpos, origin, dir, ray_valid, fallback_newton, fallback_raymarch);
    newton_iters = safe_newton + fallback_newton;
    raymarch_iters = safe_raymarch + fallback_raymarch;
    return res;
  };

  auto finalize = [&](double t) -> HitPointData
  {
    const double r2 = A*t*t + B*t + C;
    if (!(r2 >= 0.0) || r2 > r2_max_max)
    {
      ray_valid = false;
      origin.x = origin.y = origin.z = 0.0;
      newton_iters = 0;
      raymarch_iters = 0;
      return hpd;
    }

    CalcSurfFuncAndDerivs(r2, hpd.f, hpd.df, hpd.d2f);
    const double ft_check = (hpd.f + zpos) - (oz + dz*t);

    ray_valid = (r2 <= r2_ap_max) &&
                (r2 <= r2_max_max) &&
                (std::abs(ft_check) < m_hitEps) &&
                (t > m_minT);

    origin.x = ox + dx*t;
    origin.y = oy + dy*t;
    origin.z = oz + dz*t;

    hpd.dir[0] = dx;
    hpd.dir[1] = dy;
    hpd.dir[2] = dz;
    hpd.t = t;
    return hpd;
  };

  // Reuse the usual preconditioning to obtain a good initial guess t_seed.
  double t_seed = m_minT;
  const bool interval_search = (m_subdivideRaySegment < 1000000) && (m_subdivideRaySegment > 1);
  Point3d origin_new;
  Point3d nearest;
  {
    double t0 = 0.0;
    if (MoveRayCloseToLens(zpos, origin, dir, origin_new, t0, nearest, raymarch_iters))
      t_seed = t0;
  }

  // For interval search the seed is already the midpoint of the bracket where
  // F changes sign; the residual gives its magnitude. Otherwise use the z-slab.
  double delta_z = 0.0;
  if (interval_search)
  {
    double ft = 0.0, df = 0.0;
    if (evalF(t_seed, ft, df))
      delta_z = std::fabs(ft);
  }

  double ta = m_minT;
  double tb = std::numeric_limits<double>::infinity();

  if (interval_search && t_seed > 0)
  {
    ta = t_seed;
    tb = t_seed + delta_z;
  }
  else
  {
    const double zlo = zpos + zmin;
    const double zhi = zpos + zmax;

    if (std::abs(dz) > 1e-12)
    {
      const double t1 = (zlo - oz) / dz;
      const double t2 = (zhi - oz) / dz;
      ta = std::max(ta, std::min(t1, t2));
      tb = std::min(tb, std::max(t1, t2));
    }
    else
    {
      if (oz < zlo || oz > zhi)
        return fallback();
    }
  }

  // Intersect the candidate interval with the radial domain where the surface
  // function can be evaluated.
  if (A > 1e-24)
  {
    const double disc = B*B - 4.0*A*(C - r2_eval_max);
    if (disc < 0.0)
      return fallback();

    const double sq = std::sqrt(disc);
    const double tr1 = (-B - sq) / (2.0*A);
    const double tr2 = (-B + sq) / (2.0*A);

    ta = std::max(ta, std::min(tr1, tr2));
    tb = std::min(tb, std::max(tr1, tr2));
  }
  else
  {
    if (C > r2_eval_max)
      return fallback();
  }

  if (!(tb > ta))
    return fallback();

  // Verify that the bracket actually changes the sign of F.
  double fa = 0.0, dfa = 0.0, fb = 0.0, dfb = 0.0;
  const bool ok_a = evalF(ta, fa, dfa);
  const bool ok_b = evalF(tb, fb, dfb);

  if (!ok_a || !ok_b || (fa * fb > 0.0))
    return fallback();

  if (fa == 0.0)
  {
    HitPointData res = finalize(ta);
    return ray_valid ? res : fallback();
  }
  if (fb == 0.0)
  {
    HitPointData res = finalize(tb);
    return ray_valid ? res : fallback();
  }

  double t = std::clamp(t_seed, ta, tb);
  double f = 0.0, df = 0.0;

  if (!evalF(t, f, df))
  {
    t = 0.5 * (ta + tb);
    if (!evalF(t, f, df))
      return fallback();
  }

  // Keep the current iterate inside the sign-changing interval.
  if (fa * f <= 0.0)
  {
    tb = t;
    fb = f;
  }
  else
  {
    ta = t;
    fa = f;
  }

  for (uint32_t it = 0; it < m_maxIterations; ++it)
  {
    newton_iters = it + 1;

    if (std::abs(f) <= m_hitEps)
      break;

    // Prefer Newton whenever it stays inside the sign-changing bracket;
    // otherwise fall back to bisection.
    double t_next = 0.5 * (ta + tb);
    if (std::abs(df) > 1e-30)
    {
      const double t_newton = t - f / df;
      if (t_newton > ta && t_newton < tb)
        t_next = t_newton;
    }

    if (std::abs(t_next - t) <= 0.5 * m_hitEps)
    {
      t = t_next;
      evalF(t, f, df);
      break;
    }

    t = t_next;

    if (!evalF(t, f, df))
    {
      t = 0.5 * (ta + tb);
      if (!evalF(t, f, df))
        return fallback();
    }

    if (fa * f <= 0.0)
    {
      tb = t;
      fb = f;
    }
    else
    {
      ta = t;
      fa = f;
    }

    if ((tb - ta) <= 0.5 * m_hitEps && std::abs(f) <= m_hitEps)
      break;
  }

  HitPointData result = finalize(t);
  return ray_valid ? result : fallback();
}

void Surface::kernel1D_TraceRays(const double* in_zpos, const double* in_eta_inv,
                                 const double* in_ray_pos, const double* in_ray_dir, const uint8_t* in_ray_mask,
                                 double* out_hit_pos, double* out_hit_dir, uint8_t* out_hit_mask, double* out_hit_t, double* out_hit_normal,
                                 HitPointData* in_hpd,
                                 uint32_t n_iters, uint32_t n_rays_square, uint32_t z_pos_size)
{
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (uint32_t tid = 0; tid < n_iters; ++tid)
  {
    if (in_ray_mask[tid] > 0)
    {
      const uint32_t batchId = tid / n_rays_square;
      const double   zpos    = z_pos_size == 1 ? in_zpos[0] : in_zpos[batchId];

      const Point3d input_origin = make_double3(in_ray_pos[tid*3 + 0], in_ray_pos[tid*3 + 1], in_ray_pos[tid*3 + 2]);
      const Vect3d  ray_direction = make_double3(in_ray_dir[tid*3 + 0], in_ray_dir[tid*3 + 1], in_ray_dir[tid*3 + 2]);
      Point3d ray_origin = input_origin;

      bool ray_valid = true;
      uint32_t newton_iters = 0;
      uint32_t raymarch_iters = 0;
      HitPointData hpd = IntersectionSafe(zpos, ray_origin, ray_direction, ray_valid,
                                          newton_iters, raymarch_iters); // write new ray pos to 'ray_origin'

      // The safeguarded method may fail to build a bracket, or may stop at a
      // residual that is not tight enough for the caller. In that case fall back
      // to the proven preconditioned Newton path.
      const double residual = ray_valid ? std::abs(hpd.f + zpos - ray_origin.z) : 0.0;
      if (!ray_valid || residual > 0.5 * m_hitEps)
      {
        ray_origin = input_origin;
        uint32_t ni2 = 0, ri2 = 0;
        hpd = Intersection(zpos, ray_origin, ray_direction, ray_valid, ni2, ri2);
        newton_iters += ni2;
        raymarch_iters += ri2;
      }

      #pragma omp atomic update
      m_newton_iterations += newton_iters;

      #pragma omp atomic update
      m_raymarch_iterations += raymarch_iters;

      Vect3d out_normal        = make_double3(0,0,0);
      Vect3d out_ray_direction = make_double3(0,0,0);
      if (ray_valid)
      {
        const double eta = 1.0 / in_eta_inv[batchId];
        const double r2 = ray_origin.x * ray_origin.x + ray_origin.y * ray_origin.y;
        const Vect3d hpd_dir = make_double3(hpd.dir[0], hpd.dir[1], hpd.dir[2]);
        ray_valid = CalcRefrDir(r2, hpd.df, hpd_dir, ray_origin, eta, out_ray_direction, out_normal);
      }

      if (in_hpd != nullptr)
        in_hpd[tid] = hpd;

      if((m_flags & OUTFLAG_STORE_T) != 0)
      {
        out_hit_t   [tid] = ray_valid ? hpd.t : -1;
        out_hit_mask[tid] = uint8_t(ray_valid);
      }

      out_hit_pos   [tid*3 + 0] = ray_origin.x;
      out_hit_pos   [tid*3 + 1] = ray_origin.y;
      out_hit_pos   [tid*3 + 2] = ray_origin.z;

      if((m_flags & OUTFLAG_STORE_RAYS) != 0)
      {
        out_hit_dir   [tid*3 + 0] = out_ray_direction.x;
        out_hit_dir   [tid*3 + 1] = out_ray_direction.y;
        out_hit_dir   [tid*3 + 2] = out_ray_direction.z;
      }

      if((m_flags & OUTFLAG_STORE_NORMAL) != 0)
      {
        out_hit_normal[tid*3 + 0] = out_normal.x;
        out_hit_normal[tid*3 + 1] = out_normal.y;
        out_hit_normal[tid*3 + 2] = out_normal.z;
      }
    }
    else
    {
      if((m_flags & OUTFLAG_STORE_T) != 0)
      {
        out_hit_t   [tid] = 1.0;
        out_hit_mask[tid] = 0;
      }

      out_hit_pos   [tid*3 + 0] = 1.0;
      out_hit_pos   [tid*3 + 1] = 1.0;
      out_hit_pos   [tid*3 + 2] = 1.0;
      if((m_flags & OUTFLAG_STORE_RAYS) != 0)
      {
        out_hit_dir   [tid*3 + 0] = 1.0;
        out_hit_dir   [tid*3 + 1] = 1.0;
        out_hit_dir   [tid*3 + 2] = 1.0;
      }

      if((m_flags & OUTFLAG_STORE_NORMAL) != 0)
      {
        out_hit_normal[tid*3 + 0] = 1.0;
        out_hit_normal[tid*3 + 1] = 1.0;
        out_hit_normal[tid*3 + 2] = 1.0;
      }
    }
  }
}

void Surface::TraceRays(const IntersectionParameters& inters_params,
                        const double* zpos, const double* surface_params,
                        const uint32_t used_rparams,
                        const double* eta_inv,
                        const double* ray_pos, const double* ray_dir, const uint8_t* ray_mask,
                        double* hit_pos, double* hit_dir, uint8_t* hit_mask, double* hit_t, double* hit_normal,
                        HitPointData* hpd,
                        uint32_t n_rays_x, uint32_t n_rays_y, uint32_t n_batches, uint32_t z_pos_size)
{
  ///  SetIntersectionParams(inters_params);
  // if(n_rays_x*n_rays_y*n_batches > m_hpd.capacity())
  // {
  //   PrintError("Not enough memory for hit point data");
  //   return;
  // }
  auto start = std::chrono::high_resolution_clock::now();
  SetIntersectionParams(inters_params); // TODO: move this inside kernel
  UpdateSurface(surface_params, used_rparams);
  m_newton_iterations = 0;
  m_raymarch_iterations = 0;
  kernel1D_TraceRays(zpos, eta_inv, ray_pos, ray_dir, ray_mask, hit_pos, hit_dir, hit_mask, hit_t, hit_normal, hpd, n_rays_x*n_rays_y*n_batches, n_rays_x*n_rays_y, z_pos_size);
  trace_time = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()) * 0.001f;
}

#ifdef NEW_START_POINT4NEWTON
/// Finds the low bound of 't' intersection the interval
/// @param[in] origin     - ray position for t=0
/// @param[in] direction  - ray direction
/// @param[in] zpos       - surface shift as usual
/// @param[in,out] iter   - On INPUT, allowed (max) number the number until complete, usually 100
///                         On EXIT, the last iterations actually done
/// @param[in] threshold1 - stop when the "free path" from the CURRENT ray point is <= threshold1
///                         usually = 0.001
/// @param[in] threshold2 - stop when F for the CURRENT ray point is <= threshold2
///                         usually = 0.001
/// @param[out] t0        - low bound of free path.
/// @return  0 succeded,
///         -1 failure in calculation,
///         -2 if intersection is (seems?) impossible,
///         -3 stop criteria had not been achieved diring allowed 'niter',
///            the 't0' <b>make sence and can be used</b>
///         -4 algorithm faced a controversy
int Surface::ApproxIntersection(const Vect3d& origin,
                                const Vect3d& direction,
                                const double zpos,
                                int& iter,
                                const double threshold1,
                                const double threshold2,
                                double& t0) const
  {
  const double r2max = SQR(min(aperture, max_radius));
  const double a = SQR(direction.x) + SQR(direction.y);
  const double niter = iter;

  if (fabs(a) < 1e-12)
    {
    // Ray parallel to the optical axis Oz
    double r2, f, df;

    // r^2(t)
    r2 = origin.x * origin.x + origin.y * origin.y;
    CalcSurfFuncAndDeriv(r2, f, df);

    t0 = ((zpos + f) - origin.z) / direction.z;
    return 0;
    }

  double r2, f, df, dt, F_t0 = 1e99;
  Vect3d origin_new = origin;

  t0 = 0;
  for (iter = 0; iter < niter; iter++)
    {
    const double b = 2 * (origin_new.x * direction.x + origin_new.y * direction.y);

    // r^2(t0)
    r2 = origin_new.x * origin_new.x + origin_new.y * origin_new.y;

    if (std::isnan(r2))
      {
      assert(false);
      return -1;
      }

    if (r2 > r2max + 1e-10 && b > 1e-10)  // reached the boundary of allowed r => no intersection inside domain
      return -2;

    CalcSurfFuncAndDeriv(r2, f, df);

    // The deviation at t0:
    const double F_t0_ = (f + zpos) - origin_new.z;

    if (fabs(F_t0_) < 1e-14/*1e-11*/) // aimed at flat surfaces
      return EXIT_SUCCESS;

    if (fabs(F_t0_) < threshold2/* && fabs(F_t0_) < fabs(F_t0)*/)
      {
      // We are near the intersection. Use Newton:
      const double dt_F = b * df - direction.z; // = F'(0):
      const double delta = -F_t0_ / dt_F;
      if (fabs(delta) < threshold1)
        {
        if (t0 + delta > 1e-5 || (delta < 1e-10 && t0 + delta > 1e-12))
          {
          t0 += delta;
          //        origin_new = origin + t0 * direction;
          return EXIT_SUCCESS;
          }
        }
      }

    if (iter > 0 && (F_t0_ * F_t0 < 0/* || fabs(F_t0_) < EPSILON*/)) // changed sign => intersection!
      {
      // A sign change is a valid bracket even if the residual is still
      // larger than EPSILON: Newton from this bracket will converge.
      return EXIT_SUCCESS;
      }
    F_t0 = F_t0_;

#ifdef NEWRAYMARCH
    dt = RayStep2(r2, direction.z, a, b, F_t0, df, d2f_min_arr_arr, d2f_max_arr_arr, r2max);
#else
    dt = RayStep(df_min_arr_arr, df_max_arr_arr, r2max, direction, r2, a, b, F_t0);
#endif

    if (dt == 1e10 || dt < 0)
      return -2; // never intersects

    t0 += dt;
    origin_new = origin + t0 * direction;

#if 0
    if (t0 < 1.0)
      {
      const double r2_first = SQR(origin_new.x + (t0 - dt) * direction.x)
        + SQR(origin_new.y + (t0 - dt) * direction.y);
      const double r2_second = SQR(origin_new.x + t0 * direction.x)
        + SQR(origin_new.y + t0 * direction.y);

      if (r2_first != r2_second && fabs(r2_first - r2_second) < threshold * (r2_right - r2_left))
        return EXIT_SUCCESS;
      if (iter > (1 * niter / 2) && fabs(r2 - r2_first) < threshold * (r2_right - r2_left))
        return EXIT_SUCCESS;
      }
#endif
    }
  return -3;
  }


/// When the ray enters the lens domain cylinder of max_radius between the bounding planes
double Surface::RayEntersLensDomain(double zpos, const Point3d& origin, const Vect3d& direction) const
  {
  const double localAperture = aperture;   // VisitMemberExpr bug work around
  const double localmax_radius = max_radius; // VisitMemberExpr bug work around
  const double r2max = SQR(min(localAperture, localmax_radius));
  double t0 = -1, r2;

  if (direction.z >= 0)
    t0 = (zmin + zpos - origin.z) / direction.z; // The FIRST entry
  else
    t0 = (zmax + zpos - origin.z) / direction.z; // The FIRST entry

  t0 = max(0.0, t0);

  const Point3d origin_new = origin + t0 * direction;
  // r^2(t0)
  r2 = origin_new.x * origin_new.x + origin_new.y * origin_new.y;

  if (r2 > r2max || t0 < 0)
    {
    // Instead, find when r(t) = max allowed radius
    const double a = SQR(direction.x) + SQR(direction.y);
    const double b = 2 * (origin.x * direction.x + origin.y * direction.y);
    const double c = SQR(origin.x) + SQR(origin.y) - r2max;

    double dt1, discr;

    // a*t0^2 + b*t0 + c = 0
    double D = b * b - 4 * a * c;
    if (fabs(D) < 1e-12)
      {
      if (a == 0 && b == 0)
        return -1;
      dt1 = -b / (2 * a);
      }
    else
      {
      dt1 = 1e10;
      if (D >= 0)
        {
        double dt_a, dt_b;
        discr = sqrt(D);
        dt_a = (-b - discr) / (2 * a);
        dt_b = (-b + discr) / (2 * a);

        // smallest POSITIVE root
        if (dt_a < 0)
          dt1 = dt_b;
        else
          dt1 = dt_a;
        }
      }
    if (dt1 < t0 || dt1 == 1e10)
      {
      t0 = -1;
      return -1;
      }
    t0 = dt1;
    }
  return t0;
  }
#endif

/// Move the ray close to the lens surface where its origin obeys r<=aperture & r<= max_radius
/// @return natural parameter 't' (the distance from the ray origin to its new point)
/// or -1 if failed (i.e. all the available variants are outside of the max_radius)
/// @param zpos --- Shift of the lens surface along Oz (the value to be added to f(r^2))
bool Surface::MoveRayCloseToLens(double zpos, const Point3d& origin, const Vect3d& direction, Point3d& origin_new, double &t0,
                                 Point3d& nearest, uint32_t& raymarch_iters) const
{
  const double localAperture   = aperture;   // VisitMemberExpr bug work around
  const double localmax_radius = max_radius; // VisitMemberExpr bug work around
#if 1 // Ose experiment
  const double r2max = SQR(excess_of_aperture * min(localAperture, localmax_radius));
#else
  const double r2max = SQR(min(localAperture, localmax_radius));
#endif
  double r2;

  // The method to be used:
  const bool simple = (m_subdivideRaySegment <= 1) || IsFlat();
  const bool ray_marching = m_subdivideRaySegment == 1000000;
  const bool interval_search = !simple && !ray_marching;
  (void)nearest; // reserved for alternative starting-point strategies

  const double a = SQR(direction.x) + SQR(direction.y);

  // When the ray FIRST enters the domain where CalcFuncAndDeriv() works
  if (fabs(a) < 1e-12)
    {
    // Ray parallel to the optical axis Oz
    double r2, f, df;

    // r^2(t)
    r2 = origin.x * origin.x + origin.y * origin.y;
    CalcSurfFuncAndDeriv(r2, f, df);

    t0 = ((zpos + f) - origin.z) / direction.z;
    origin_new = origin + t0 * direction;
    return t0;
    }

#ifdef NEW_START_POINT4NEWTON
  if (ray_marching)
    {
    // When the ray FIRST enters the domain where CalcFuncAndDeriv() works
    t0 = RayEntersLensDomain(zpos, origin, direction);
    if (t0 < 0)
      return false;

    origin_new = origin + t0 * direction;

#if 0 // preliminary bracketing
    const double t_min = t0;
    const double b = 2 * (origin.x * direction.x + origin.y * direction.y);

    const std::vector<double>& f_min_arr = f_min_arr_arr[f_min_arr_arr.size() - 2];
    const std::vector<double>& f_max_arr = f_max_arr_arr[f_min_arr_arr.size() - 2];
    const int ni = (int)f_max_arr.size();
    const double dr2 = r2max / ni;
    const double r2_0 = origin.x * origin.x + origin.y * origin.y;

    double2 t_range, r2_range, f_range(zmin, zmax);
    double r2_left, r2_right, z_left = zmin, z_right = zmax, delta, delta_prev;
    int iter;

    if (fabs(a) > 1e-12)
      {
      for (iter = 0; iter < 20; iter++)
        {
        if (direction.z >= 0)
          {
          t_range[0] = (f_range[0] + zpos - origin.z) / direction.z; // The FIRST entry
          t_range[1] = (f_range[1] + zpos - origin.z) / direction.z;
          }
        else
          {
          t_range[0] = (f_range[1] + zpos - origin.z) / direction.z; // The FIRST entry
          t_range[1] = (f_range[0] + zpos - origin.z) / direction.z;
          }

        // Do not allow t which makes the point OUTSIDE the domain
        t_range[0] = max(t_range[0], t_min);
        t_range[1] = max(t_range[1], t_min);

        if (t_range[0] < 0)
          break;///t_range[0] = 0;

        delta = t_range[1] - t_range[0];
        if (iter > 0 && delta >= 0.95 * delta_prev)
          break;
        delta_prev = delta;

        r2_range = CalcR2Range(a, b, r2_0, t_range);
        if (r2_range[1] = r2max)
          r2_range[1] = r2max;

        f_range = CalcFuncBounds(f_min_arr, f_max_arr, dr2, b, r2_range, r2_left, r2_right);
        }

      if (t_range[0] > 0) // If ray origin is INSIDE DOMAIN, this time can be <0!
        t0 = t_range[0];

      origin_new = origin + t0 * direction;
      }
#endif


    // Finds the lower bound of 't' for intersection, counted from t0. It is 'dt'
    double dt;
    int niter = 100;
    int rc = ApproxIntersection(origin_new, direction, zpos, niter, 1e-3, 1e-3, dt);
    raymarch_iters += uint32_t(niter);

#ifdef PRINT_RAYMARCH_NITER
    printf("rc=%3i, %10i iterations\n", rc, niter);
#endif

#if 0
    if (rc != -3 && niter > 50)
      {
      origin_new = origin + (t0 + dt) * direction;
      r2 = origin_new.x * origin_new.x + origin_new.y * origin_new.y;
      printf("niter=%10i, rc=%3i, (r/rmax)^2=%10.5f\n", niter, rc, r2 / r2max);
      }
#endif

    if (rc < 0)
      {
      if (rc != -3) // -3 means that ApproxIntersection() did not reach acuracy. But it it did give smth...
        return false;
      }

    if (dt > 0 && dt < 1e10)
      {
      t0 += dt;
      origin_new = origin + t0 * direction;
      }

    return true;
    }
#endif

  if (simple)
    {
#if 0
    // When the ray FIRST enters the domain where CalcFuncAndDeriv() works
    t0 = RayEntersLensDomain(zpos, origin, direction);
#else
    // Move the ray to the middle plane of the lens
    t0 = ((center_z + zpos) - origin.z) / direction.z;
    origin_new = origin + t0 * direction;

    r2 = origin_new.x * origin_new.x + origin_new.y * origin_new.y;
    if (r2 > r2max)
    {
      // Most likely the hit point will be outside the domain where sag() is defined...

      // Use another plane for the preliminary intersection
      t0 = ((zmin + zpos) - origin.z) / direction.z;
      origin_new = origin + t0 * direction;

      r2 = origin_new.x * origin_new.x + origin_new.y * origin_new.y;
      if (r2 > r2max)
      {
        t0 = ((zmax + zpos) - origin.z) / direction.z;
        origin_new = origin + t0 * direction;
        r2 = origin_new.x * origin_new.x + origin_new.y * origin_new.y;

        if (r2 > r2max)
          return false; // Neither of the two bounding planes helped. Failure...
      }
    }
#endif
  }
  else if (interval_search)
  {
    // Split the ray segment between the bounding plane into nj segments,
    // find the one where the ray goes past the surface and use its middle
    // as the starting point for Newton
    const int nj = m_subdivideRaySegment;

    const double t_min = max(0.0, min(((zmin + zpos) - origin.z) / direction.z,
                                      ((zmax + zpos) - origin.z) / direction.z));
    const double t_max = max(((zmin + zpos) - origin.z) / direction.z,
                             ((zmax + zpos) - origin.z) / direction.z);
    const double dt = (t_max - t_min) / (nj - 1);

    double delta_z, delta_z_prev = 0; //  = z_surface - z_ray
    int j0 = -1000000;
    t0 = -1;

    for (int j = 0; j < nj; j++)
    {
      const double tt = t_min + j * dt;
      const Point3d o_new = origin + direction * tt;
      const double r2 = o_new.x * o_new.x + o_new.y * o_new.y;

      double ft, df;
      // When r is outside the aperture the function is evaluated for the
      // clipped r^2. Its value is still usable for locating the sign change:
      // discarding such samples makes grazing rays miss their intersection.
      CalcSurfFuncAndDeriv(r2, ft, df);
      delta_z = (ft + zpos) - o_new.z; // = z_surface - z_ray

      if (j > -nj && delta_z * delta_z_prev < 0)
        {
        // the sign changes between j-1 and j
        if (j0 == -1000000) // the FIRST intersection!!!
          {
          j0 = j;
          break;
          }
        }
      delta_z_prev = delta_z;
    }
    if (j0 != -1000000)
      {
      t0 = t_min + (j0 - 0.5) * dt;
      origin_new = origin + direction * t0;
      }
    else
      return false;
  }

  return true;
}

/// Calculates and remembers the bounds for 'z' WITHOUT zpos
/// Usually called only within @c CalcZminZmax()!
///
/// Before the call, @c aperture and @c max_radius must be initialized
///
/// @param[in] ntab - the number of values of 'r' between 0 and max(max_radius,aperture)
///                   to calculate the bounding @c zmin, @c zmax from.
void Surface::kernel1D_CalcZminZmax(int npt, double dr_fraction)
{
  zmin = double(FLT_MAX);
  zmax = double(-FLT_MAX);

  double f_prev, df_prev;
  bool good = CalcSurfFuncAndDeriv(0.0, f_prev, df_prev);
  if (good && !isnan(f_prev))
    {
    zmin = min(zmin, float(f_prev));
    zmax = max(zmax, float(f_prev));
    }

  for (int j = 1; j < npt; j++)
  {
    const double localAperture   = aperture;                        // VisitMemberExpr bug work around
    const double localmax_radius = max_radius;                      // VisitMemberExpr bug work around
    const double dr = min(localAperture, localmax_radius) * dr_fraction;  // VisitMemberExpr bug work around
    const double r = (j == npt - 1 ? min(localAperture, localmax_radius) : j * dr);
    double f, df;
    bool good = CalcSurfFuncAndDeriv(r * r, f, df);
    if (good && !isnan(f))
    {
      zmin = min(zmin, float(f));
      zmax = max(zmax, float(f));

      if (df * df_prev < 0)
        {
///        printf("j=%10i %10.5f %10.5f    %10.5f %10.5f\n", j, df, df_prev, f, f_prev);

        // In this interval min(f) and/or max(f) can be INSIDE and differ from the boundary values
        // Approximate f(r^2) with a cubic function of r^2-r_prev^2, using two boundary values of f and f'
        // and calculate its extrema. Parabola is a3*x^2+a2*x^2+a1*x+a0, a0=f_prev
        double r2_prev = SQR((j - 1) * dr);
        double x_max = SQR(r) - r2_prev;
        const double c1 = (df - df_prev) / x_max;
        const double c2 = (f - f_prev - df_prev * x_max) / (x_max * x_max);

        const double a0 = f_prev;
        const double a1 = df_prev;
        const double a2 = 3 * c2 - c1;
        const double a3 = (c1 - 2 * c2) / x_max;

        double d = a2 * a2 - 3 * a1 * a3;
        if (d >= 0)
          {
          d = sqrt(d);
          // points of extrema:
          const double x1 = (-a2 - d) / (3 * a3), x2 = (-a2 + d) / (3 * a3);

          if (x1 >= 0 && x1 <= x_max)
            {
            double f1, df1;
            bool good = CalcSurfFuncAndDeriv(x1 + r2_prev, f1, df1);
            if (good && !isnan(f1))
              {
///              printf("x1=%10.5f, f1=%10.5f\n", x1, f1);
              zmin = min(zmin, float(f1));
              zmax = max(zmax, float(f1));
              }
            }

          if (x2 >= 0 && x2 <= x_max)
            {
            double f2, df2;
            bool good = CalcSurfFuncAndDeriv(x2 + r2_prev, f2, df2);
            if (good && !isnan(f2))
              {
///              printf("x2=%10.5f, f2=%10.5f\n", x2, f2);
              zmin = min(zmin, float(f2));
              zmax = max(zmax, float(f2));
              }
            }
          }
        }
    }
    f_prev  = f;
    df_prev = df;
  }
///  printf("zmin=%10.5f zmax=%10.5f\n", zmin, zmax);


  // These bounds were calculated approximately because of finite number of tested points
  // Apply some safe margins to cure:
  zmin = zmin - 0.001f * (zmax - zmin);
  zmax = zmax + 0.001f * (zmax - zmin);
}

void Surface::UpdateIntersectionParamsCPU(const IntersectionParameters& params)
{
  m_hitEps = params.hitEps;
  m_minT = params.minT;
  m_flags = params.flags;

  m_useRebasePoly = params.rebasePoly;
  m_useFloat = params.useFloat;
}



void Surface::Reset(HitPointData &data)
{
  data.dir[0] = data.dir[1] = data.dir[2] = 0;
  // data.xyz[0] = data.xyz[1] = data.xyz[2] = 0;
  data.t = data.df = data.d2f = 0;
}

/// Calculates direction of the surface normal

Vect3d Surface::CalcNormal(const double r2, const double df, const Point3d& hitPos, bool unit)
{
  Vect3d n = make_double3(-2 * hitPos.x * df, -2 * hitPos.y * df, 1);
  if (unit)
    n *= 1 / sqrt(1 + r2 * SQR(2 * df));
  return n;
}

bool Surface::CalcRefrDir(const double r2, const double df, const Vect3d& dir, const Point3d& hitPos,  double eta, Vect3d& dir_out, Vect3d& n)
{
//n = CalcNormal(data, true); // missed "CalcNormal" bug work around
{
  n = make_double3(-2 * hitPos.x * df, -2 * hitPos.y * df, 1);
  n *= 1 / sqrt(1 + r2 * SQR(2 * df));
}

double cos_in = dot(n, dir);

const double tmp = cos_in * cos_in + eta * eta - 1;
if (tmp < 0)
  return false;

// Expression below assumes orientation of the normal such that (n,v0) >= 0, ensure this:
if (cos_in < 0)
  {
  n = -n;
  cos_in = -cos_in;
  }
dir_out = (sqrt(tmp) - cos_in) * n + dir;
dir_out *= 1 / eta;

#if 0 // TEST
assert(std::abs(length(dir_out) - 1) < EPSILON);
double sin_in = sqrt(1 - cos_in * cos_in);
double sin_out = sqrt(1 - SQR(dot(n, dir_out)));
assert(std::abs(sin_out - sin_in / eta) < EPSILON);
///  printf("%f %f %f %f\n", eta, sin_in, sin_out, sin_in / eta);
#endif
return true;
}



void Surface::CalcBoundsOfFuncAndDeriv(int nlevels, int nfirst, int reduction)
  {
  const double r2_max = SQR(min(aperture, max_radius));

#ifdef NEWRAYMARCH
  d2f_min_arr_arr.resize(nlevels);
  d2f_max_arr_arr.resize(nlevels);
#endif
  df_min_arr_arr.resize(nlevels);
  df_max_arr_arr.resize(nlevels);
  f_min_arr_arr.resize(nlevels);
  f_max_arr_arr.resize(nlevels);

  double f, df, d2f;

  int n = nfirst;
  for (int level = 0; level < nlevels; level++)
    {
#ifdef NEWRAYMARCH
    d2f_min_arr_arr[level].resize(n);
    d2f_max_arr_arr[level].resize(n);
#endif
    df_min_arr_arr[level].resize(n);
    df_max_arr_arr[level].resize(n);
    f_min_arr_arr[level].resize(n);
    f_max_arr_arr[level].resize(n);
    n *= reduction;
    }

  n = (int)df_min_arr_arr[nlevels - 1].size();
  const double dr2 = r2_max / n;
  const int npoints = 20;

  for (int i = 0; i < n; i++)
    {
    const double r2_left = dr2 * i;
    const double r2_right = dr2 * (i + 1);

    const double dr2_ = (r2_right - r2_left) / (npoints - 1);

#ifdef NEWRAYMARCH
    double d2f_min = +1e10, d2f_max = -1e10;
#endif
    double df_min = +1e10, df_max = -1e10, f_min = +1e10, f_max = -1e10;
    for (int j = 0; j < npoints; j++)
      {
      const double r2 = r2_left + dr2_ * j;
#ifdef NEWRAYMARCH
      const bool good = CalcSurfFuncAndDerivs(r2, f, df, d2f);
#else
      const bool good = CalcSurfFuncAndDeriv(r2, f, df);
#endif
      if (good)
        {
#ifdef NEWRAYMARCH
        d2f_min = min(d2f, d2f_min);
        d2f_max = max(d2f, d2f_max);
#endif
        df_min = min(df, df_min);
        df_max = max(df, df_max);

        f_min = min(f, f_min);
        f_max = max(f, f_max);
        }
      }
#ifdef NEWRAYMARCH
    d2f_min_arr_arr[nlevels - 1][i] = d2f_min;
    d2f_max_arr_arr[nlevels - 1][i] = d2f_max;
#endif
    df_min_arr_arr[nlevels - 1][i] = df_min;
    df_max_arr_arr[nlevels - 1][i] = df_max;
    f_min_arr_arr[nlevels - 1][i] = f_min;
    f_max_arr_arr[nlevels - 1][i] = f_max;
    }

  for (int level = nlevels - 1; level >= 1; level--)
    {
    n = (int)df_min_arr_arr[level - 1].size();
    for (int i = 0; i < n; i++)
      {
#ifdef NEWRAYMARCH
      double d2f_min = +1e10, d2f_max = -1e10;
#endif
      double df_min = +1e10, df_max = -1e10, f_min = +1e10, f_max = -1e10;
      for (int j = 0; j < reduction; j++)
        {
#ifdef NEWRAYMARCH
        d2f_min = min(d2f_min, d2f_min_arr_arr[level][i * reduction + j]);
        d2f_max = max(d2f_max, d2f_max_arr_arr[level][i * reduction + j]);
#endif
        df_min = min(df_min, df_min_arr_arr[level][i * reduction + j]);
        df_max = max(df_max, df_max_arr_arr[level][i * reduction + j]);
        f_min = min(f_min, f_min_arr_arr[level][i * reduction + j]);
        f_max = max(f_max, f_max_arr_arr[level][i * reduction + j]);
        }
#ifdef NEWRAYMARCH
      d2f_min_arr_arr[level - 1][i] = d2f_min;
      d2f_max_arr_arr[level - 1][i] = d2f_max;
#endif
      df_min_arr_arr[level - 1][i] = df_min;
      df_max_arr_arr[level - 1][i] = df_max;
      f_min_arr_arr[level - 1][i] = f_min;
      f_max_arr_arr[level - 1][i] = f_max;
      }
    }
  }
