#include "surface.h"
#include <iostream>

/// Save cross-section table i.e. z vs. r (not vs r^2!)
/// @attention w/o 'zpos'!
void Surface::SaveSagCrossSection(const std::string& path, size_t n_points, const double* surface_params)
  {
  FILE* fd;
  fd = fopen(path.data(), "w");

  const double rmax = std::min(aperture, max_radius);
  const double dr = 2 * rmax / n_points;
  double r, f, df;
  for (r = -rmax; r <= rmax + EPSILON; r += dr)
    {
    CalcSurfFuncAndDeriv(r * r, f, df);
    fprintf(fd, "r = %10.5f, z = %10.5f\n", r, f);
    }
  fclose(fd);
  }

/// Save the slope cross-section table i.e. dz/dr z vs. r (not vs r^2!)
void Surface::SaveSagSlopeCrossSection(const std::string& path, size_t n_points, const double* surface_params)
  {
  FILE* fd;
  fd = fopen(path.data(), "w");

  const double rmax = std::min(aperture, max_radius);
  const double dr = 2 * rmax / n_points;
  double r, f, df;
  for (r = -rmax; r <= rmax + EPSILON; r += dr)
    {
    CalcSurfFuncAndDeriv(r * r, f, df);
    fprintf(fd, "r = %10.5f, dz/dr = %10.5f\n", r, 2 * r * df);
    }
  fclose(fd);
  }

/// Calculate lens section z=f(r) for uniform mesh in 'r'
void Surface::CalcCrossSection(size_t n_points, std::vector<double>& r_arr, std::vector<double>& z_arr) const
  {
  r_arr.resize(n_points);
  z_arr.resize(n_points);

  const double rmax = std::min(max_radius, aperture);
  const double dr = 2 * rmax / (n_points - 1);
  double r, f, df;
  for (size_t i = 0; i < n_points; i++)
    {
    r = -rmax + dr * i;
    CalcSurfFuncAndDeriv(r * r, f, df);
    r_arr[i] = r;
    z_arr[i] = f;
    }
  }

void Surface::PrintError(const std::string& str)
{
  std::cout << str << std::endl;
}

#include "poly_surface.h"

void PolySurface::DebugPrint(const std::string& a_what, std::ostream& a_out)
{
  if(a_what == "surf_data")
  {
    a_out << "surf_data: {" << std::endl; 
    a_out << "surf_data.conic    : " << surf_data.conic     << std::endl;
    a_out << "surf_data.curvature: " << surf_data.curvature << std::endl;
    a_out << "surf_data.rnorm    : " << surf_data.rnorm     << std::endl;
    a_out << "surf_data.aperture : " << surf_data.aperture  << std::endl;
    a_out << "aperture           : " << aperture            << std::endl;
    a_out << "npoly              : " << npoly               << std::endl;
    
    for (uint32_t j = 0; j < npoly; j++) 
    {
      a_out << "surf_data.rparams[" << j << "]: " << surf_data.rparams[j] << std::endl;
    }
    
    // CalcSlaveValues() : 
    
    // max(r) for which the surface expression is well-defined
    a_out << "max_radius: " << max_radius << std::endl;

    // This is f(r=0), which for the formula used is 0
    a_out << "center_z: " << center_z << std::endl;

    a_out << "zmin: " << zmin << std::endl;
    a_out << "zmax: " << zmax << std::endl;

    a_out << "}" << std::endl;
  }

  if(a_what == "m_params")
  {
    a_out << "m_params: {" << std::endl; 
    a_out << "m_params.EPSILON           = " << Surface::EPSILON << std::endl;
    a_out << "m_params.hitEps            = " << m_hitEps << std::endl;
    a_out << "m_params.MAX_DISTANCE      = " << Surface::MAX_DISTANCE << std::endl;
    a_out << "m_params.NEWTON_STEP_BOUND = " << Surface::NEWTON_STEP_BOUND << std::endl;
    a_out << "} " << std::endl; 
  }

}
