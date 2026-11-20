#pragma once
#include "LiteMath.h"

using Float  = double;
using Float3 = LiteMath::double3;

#include <cassert>
#include <string>
#include <cstdint>
#include <vector>
#include <array>

#include "isurface.h"

/// Description of one surface read from a lens JSON file.
struct SurfaceParams
{
  uint32_t id = 0;
  Float aperture_radius = 2.457944025403378;
  Float zdist = 0.0;
  Float curvature = 1.0 / 3.58880536446178;
  Float conic = 0.37196;
  Float rnorm = 2.457943930331;
  Float eta = 1.54;
  Float abbe_num = 55.6;

  // By default use fixed eta and abbe_num; if use_sellmeier == true, compute
  // eta and abbe_num from sellmeier_coeffs for the requested wavelength.
  bool use_sellmeier = false;

  std::string material_name = "N-SK2";

  // Sellmeier coefficients (B) and resonance wavelengths squared (C), interleaved.
  std::vector<Float> sellmeier_coeffs = {
    1.28189012, 0.0072719164, 0.257738258,
    0.0242823527, 0.96818604, 110.377773,
    0.0, 0.0
  };

  // The 2nd degree coefficient is assumed to always be 0.
  std::vector<Float> poly_coeffs = {
    0.0,
    -0.04607460173359,
    -0.01600407708108,
    -0.04692170620822,
    -0.0638401104121,
    0.1282613807002,
    -0.1323516385156,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };
};

struct Statistics {
    float mean;
    float stddev;
    float variance;

    Statistics() : mean(0.0f), stddev(0.0f), variance(0.0f) {}
};

Statistics computeStatistics(const std::vector<float>& values);

// F (486.1 nm) - Blue hydrogen line (Fraunhofer F-line)
// d (587.6 nm) - Yellow helium line (Fraunhofer d-line)
// C (656.3 nm) - Red hydrogen line (Fraunhofer C-line)
static const std::array<Float, 3> wave_FdC = {486.1, 587.6, 656.3};

static inline std::vector<Float> linspace(Float start, Float stop, size_t num)
{
  assert(num > 0);

  std::vector<Float> res(num, 0.0);
  if(num == 1)
  {
    res[0] = (stop + start) * 0.5;
  }

  Float step = (stop - start) / Float(num - 1);
  for(size_t idx = 0; idx < num; ++idx)
  {
    res[idx] = start + step * Float(idx);
  }

  return res;
}

bool load_surfaces_from_json(const std::string& path, std::vector<SurfaceParams>& surfaces, bool skip_plane_surfaces = true,
                             bool add_image_plane = false);

void convert_surface(const SurfaceParams& params, std::vector<double>& surface, double& zpos, double& eta,
                     const Float wavelength = 555.0);

bool save_data_to_exr(const std::vector<double> &data, int32_t w, int32_t h, const std::string &path);

void save_t_hit_txt(const std::vector<Float>& t_hit, int32_t w, int32_t h, const std::string &path);
void save_mask_to_png(const std::string& path, int width, int height, int channels, const std::vector<uint8_t>& data);

void save_points_txt(const std::vector<Float3> &pts, const std::vector<int32_t> &valid, int32_t w, int32_t h, const std::string &path, bool delimit_rows = false);

std::vector<uint8_t> ray_mask_from_png(const std::string& path);

std::vector<Float3> load_points_txt(int32_t w, int32_t h, const std::string &path);
std::vector<Float> load_t_hit_txt(int32_t w, int32_t h, const std::string &path);

bool load_rays(const std::string &path, std::vector<Float3> &out_ray_pos, std::vector<Float3> &out_ray_dir);

bool save_t_hit_image(const std::string& experiment_name, std::vector<Float> t_hit, const std::vector<uint8_t> &t_valid,
                      Float obj_dist, int32_t w, int32_t h, bool propagate_ray = false, Float propagate_z = 0.0);

bool compare_hit_points(const std::vector<Float3>& t_hit_ref, const std::vector<Float3>& t_hit_tst,
                        const std::string& test_name, Float threshold, bool skip_zeros, bool verbose);

bool compare_t_hit_values(const std::vector<Float>& t_hit_ref, const std::vector<Float>& t_hit_tst,
                          const std::string& test_name, Float threshold, bool verbose);

void get_eta_and_abbe_number(const SurfaceParams& surf, const std::vector<Float>& wavelengths,
                             std::vector<Float>& eta, Float& abbe_num);

void generate_rays(uint32_t num_rays_x, uint32_t num_rays_y, Float obj_dist, Float field_angle, Float entrance_pupil_radius,
                   std::vector<Float3> &ray_origin, std::vector<Float3> &ray_dir);

void generate_rays_batches(uint32_t num_rays_x, uint32_t num_rays_y, std::vector<Float> field_angles,
                           Float obj_dist, Float entrance_pupil_radius,
                           std::vector<Float3> &ray_origin, std::vector<Float3> &ray_dir);
