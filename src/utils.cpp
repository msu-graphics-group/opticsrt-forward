#include "tinyexr.h"

#include <sstream>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "LiteMath.h"

#include "utils.h"

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "stb_image.h"
#include "stb_image_write.h"

using LiteMath::normalize;
using LiteMath::length;
using LiteMath::dot;

bool surface_from_json(const json& j, SurfaceParams& surf)
{
  auto found = j.find("id");
  if(found != j.end())
  {
    found->get_to(surf.id);
  }

  found = j.find("curv_rad");
  if(found != j.end())
  {
    Float r = 1.0;
    if(found->is_string())
    {
      r = std::numeric_limits<Float>::infinity();
    }
    else 
    {
      found->get_to(r);
    }
    
    surf.curvature = 1.0 / r;
  }
  else 
  {
    return false;
  }

  found = j.find("k");
  if(found != j.end())
  {
    found->get_to(surf.conic);
  }
  else 
  {
    return false;
  }

  found = j.find("half_of_aperture");
  if(found != j.end())
  {
    found->get_to(surf.aperture_radius);
  }
  else 
  {
    return false;
  }

  found = j.find("r_norm");
  if(found != j.end())
  {
    found->get_to(surf.rnorm);
  }
  else 
  {
    return false;
  }

  found = j.find("aspheric_coeffs");
  if(found != j.end())
  {
    auto sz = found->size();
    surf.poly_coeffs.resize(sz);
    for(int i = 0; i < sz; ++i)
      surf.poly_coeffs[i] = j["aspheric_coeffs"][i].get<Float>();
  }
  else 
  {
    return false;
  }

  found = j.find("n_and_abbe_num");
  if(found != j.end())
  {
    auto name_node = found->find("name");
    if(name_node != found->end())
    {
      auto name_str = name_node->get<std::string>();
      if(name_str == "air")
      {
        surf.eta = 1.0;
        surf.abbe_num = 1.0;
      }
      else
      {
        auto coeff_node = found->find("coeff");
        if(coeff_node != found->end())
        {
          surf.use_sellmeier = true;
          auto sz = coeff_node->size();
          surf.sellmeier_coeffs.resize(sz);
          for(int i = 0; i < sz; ++i)
            surf.sellmeier_coeffs[i] = coeff_node->at(i).get<Float>();
        }
      }      
    }
    else 
    {
      auto str = found->get<std::string>();
      auto pos = str.find_first_of(',');
      surf.eta = std::stod(str.substr(0, pos));
      surf.abbe_num = std::stod(str.substr(pos + 1));
    }
    
  }
  else 
  {
    return false;
  }

  return true;
}

bool load_surfaces_from_json(const std::string& path, std::vector<SurfaceParams>& surfaces, bool skip_plane_surfaces, bool add_image_plane) 
{
  std::ifstream in(path);
  if (in.bad())
  {
    std::cout << "Error opening lens json file: " << path << std::endl;
    return false;
  }

  try
  {
    json j = json::parse(in);
    // size_t total_lenses = j.size() / 2;
    Float z_pos = 0.0;
    for(int i = 0; i < j.size(); ++i) 
    {
      if(j[i]["curv_rad"].is_string())
      {
        auto str = j[i]["curv_rad"].get<std::string>();
        if(str != "inf" || skip_plane_surfaces)
          continue;
      }
        
      SurfaceParams surf;
      surf.zdist = z_pos;
      surface_from_json(j[i], surf);
      z_pos += j[i]["thickness"].get<Float>();
      surfaces.push_back(std::move(surf));
    }

    if(add_image_plane)
    {
      SurfaceParams image;
      image.zdist = z_pos;
      image.conic = 0.0;
      image.aperture_radius = 20.0;
      image.rnorm = 20.0;
      image.curvature = 0.0;
      image.eta = 1.0;
      image.id = surfaces.size();
      image.material_name = "sensor";
      image.abbe_num = 0;
      
      std::fill(image.poly_coeffs.begin(), image.poly_coeffs.end(), 0.0);
      std::fill(image.sellmeier_coeffs.begin(), image.sellmeier_coeffs.end(), 0.0);

      surfaces.push_back(std::move(image));
    }
  }
  catch (json::exception& ex)
  {
    std::cerr << "lens .json (" << path << ") parse error : " << ex.what() << '\n';
  }
  catch (std::exception& ex)
  {
    std::cerr << "lens .json (" << path << ") exception : " << ex.what() << '\n';
  }

  return true;
}

bool save_data_to_exr(const std::vector<float> &data, int32_t w, int32_t h, const std::string &path)
{
  EXRHeader header;
  InitEXRHeader(&header);

  EXRImage image;
  InitEXRImage(&image);

  float* image_ptr[1];
  image_ptr[0] = const_cast<float*>(data.data()); 

  image.num_channels = 1;
  image.images = (unsigned char**)image_ptr;
  image.width = w;
  image.height = h;

  header.num_channels = 1;
  header.channels = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);

  strncpy(header.channels[0].name, "Y", 255); 
  header.channels[0].name[strlen("Y")] = '\0';

  header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
  header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
  for (int i = 0; i < header.num_channels; i++)
  {
    header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
    header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored in .EXR
  }

  const char* err;
  int ret = SaveEXRImageToFile(&image, &header, path.c_str(), &err);
  if (ret != TINYEXR_SUCCESS)
  {
    std::cerr << "Save EXR error: " << err << std::endl;
    return false;
  }
  std::cout << "Saved exr file: " << path.c_str() << std::endl;

  
  free(header.requested_pixel_types);
  free(header.pixel_types);
  free(header.channels);

  return true;
}

bool save_data_to_exr(const std::vector<double> &data, int32_t w, int32_t h, const std::string &path)
{
  std::vector<float> conv(data.begin(), data.end());
  return save_data_to_exr(conv, w, h, path);
}


std::vector<uint8_t> ray_mask_from_png(const std::string& path)
{
  int w, h, c;
  stbi_info(path.c_str(), &w, &h, &c);
  
  if(w <= 0 || h <= 0 || c <= 0)
  {
    std::cerr << "ray_mask_from_png: bad image " << path << std::endl;
    return {};
  }

  unsigned char *pixels = stbi_load(path.c_str(), &w, &h, &c, 1);

  if(c != 1)
  {
    std::cerr << "ray_mask_from_png: image has more than 1 channel " << path << std::endl;
    return {};
  }

  std::vector<uint8_t> result(w * h);
  memcpy(result.data(), pixels, result.size());

  stbi_image_free(pixels);

  return result;
}


void save_mask_to_png(const std::string& path, int width, int height, int channels, const std::vector<uint8_t>& data)
{
  std::string outfile = path + ".png";
  std::vector<uint8_t> copy(data);
  for(auto &a: copy)
    a *= 255;
  stbi_write_png(outfile.c_str(), width, height, channels, copy.data(), width * channels);
}

// save hit points in row major order
void save_points_txt(const std::vector<Float3> &pts, const std::vector<int32_t> &valid, int32_t w, int32_t h, const std::string &path, bool delimit_rows)
{
  std::ofstream f(path);

  for(int32_t y = 0; y < h; ++y)
  {
    if(delimit_rows)
      f << " --- " << y << " --- " << std::endl;
    for(int32_t x = 0; x < w; ++x)
    {
      f << std::setprecision(16) 
        << pts[x + y * w].x * valid[x + y * w] << " "
        << pts[x + y * w].y * valid[x + y * w] << " "
        << pts[x + y * w].z * valid[x + y * w] << std::endl;
    }
  }
  f.close();
}

std::vector<Float3> load_points_txt(int32_t w, int32_t h, const std::string &path)
{
  std::vector<Float3> pts;
  pts.reserve(size_t(w * h));
  std::ifstream f(path);
  if(!f.is_open())
  {
    std::cerr << "Failed to open intersection points file: "<< path << std::endl;
    return pts;
  }
  {
    for(int32_t y = 0; y < h; ++y)
    {
      for(int32_t x = 0; x < w; ++x)
      {
        Float3 pt;
        f >> pt.x >> pt.y >> pt.z;
        pts.push_back(pt);
      }
    }
    f.close();
  }

  return pts;
}

bool load_rays(const std::string &path, std::vector<Float3> &out_ray_pos, std::vector<Float3> &out_ray_dir)
{
  std::ifstream f(path);
  if(!f.is_open())
  {
    std::cerr << "Failed to open rays file: "<< path << std::endl;
    return false;
  }
  
  std::string line;
  while (std::getline(f, line))
  {
    std::istringstream iss(line);

    Float3 ray_pos;
    iss >> ray_pos.x >> ray_pos.y >> ray_pos.z;

    Float3 ray_dir;
    iss >> ray_dir.x >> ray_dir.y >> ray_dir.z;
    
    out_ray_pos.push_back(ray_pos);
    out_ray_dir.push_back(ray_dir);
  }
    
  f.close();
  
  return true;
}

std::vector<Float> load_t_hit_txt(int32_t w, int32_t h, const std::string &path)
{
  std::vector<Float> t_hit;
  t_hit.reserve(size_t(w * h));
  std::ifstream f(path);
  if(!f.is_open())
  {
    std::cerr << "Failed to open t_hit file: "<< path << std::endl;
    return t_hit;
  }
  {
    for(int32_t y = 0; y < h; ++y)
    {
      for(int32_t x = 0; x < w; ++x)
      {
        Float t;
        f >> t;
        t_hit.push_back(t);
      }
    }
    f.close();
  }

  return t_hit;
}


void save_t_hit_txt(const std::vector<Float>& t_hit, int32_t w, int32_t h, const std::string &experiment_name)
{
  std::string filename = experiment_name + "_t_hit.txt";
  std::ofstream f(filename);

  for(int32_t y = 0; y < h; ++y)
  {
    for(int32_t x = 0; x < w; ++x)
    {
      f << std::setprecision(16) << t_hit[x + y * w] << std::endl;
    }
  }
  f.close();
}


bool save_t_hit_image(const std::string& experiment_name, std::vector<Float> t_hit, const std::vector<uint8_t> &t_valid, 
                      Float obj_dist, int32_t w, int32_t h, bool propagate_ray, Float propagate_z)
{
  #pragma omp parallel for
  for(uint32_t i = 0; i < t_hit.size(); ++i)
  {
    if(propagate_ray)
    {
      t_hit[i] = t_hit[i] * (obj_dist / propagate_z) - obj_dist;
    }
    else
    {
      t_hit[i] -= obj_dist;
    }
  
    t_hit[i] *= t_valid[i];
  }

  std::string filename = experiment_name + "_t_hit.exr";

  return save_data_to_exr(t_hit, w, h, filename);
}

bool compare_hit_points(const std::vector<Float3>& t_hit_ref, const std::vector<Float3>& t_hit_tst,
                        const std::string& test_name, Float threshold, bool skip_zeros, bool verbose)
{
  if(t_hit_ref.size() != t_hit_tst.size())
  {
    std::cout << "[compare_hit_points : " << test_name << "] FAILED. Reference and test points vectors have different sizes!" << std::endl;
    return false;
  }

  bool no_errors = true;
  size_t zeros = 0;
  for(size_t i = 0; i < t_hit_ref.size(); ++i)
  {
    bool point_is_not_valid = (std::abs(length(t_hit_ref[i])) < 1e-9) || (std::abs(length(t_hit_tst[i])) < 1e-9);
    bool point_is_not_valid_both = (std::abs(length(t_hit_ref[i])) < 1e-9) && (std::abs(length(t_hit_tst[i])) < 1e-9);
    
    if (point_is_not_valid && !point_is_not_valid_both)
      zeros++;

    if(skip_zeros && point_is_not_valid)
      continue;

    auto absolute_err = std::abs(t_hit_ref[i].x - t_hit_tst[i].x) + 
                        std::abs(t_hit_ref[i].y - t_hit_tst[i].y) + 
                        std::abs(t_hit_ref[i].z - t_hit_tst[i].z);

    if (absolute_err > threshold || verbose)
    {
      std::cout << "[compare_hit_points : " << test_name << "] abs. error = " << absolute_err << std::endl;
      std::cout << std::setprecision(16) << "\tREF: " << t_hit_ref[i].x << " " << t_hit_ref[i].y << " " << t_hit_ref[i].z << std::endl;
      std::cout << std::setprecision(16) << "\tTST: " << t_hit_tst[i].x << " " << t_hit_tst[i].y << " " << t_hit_tst[i].z << std::endl;
      no_errors = false;
    }
  }

  if(no_errors)
  {
    std::cout << "[compare_hit_points : " << test_name << "] PASSED"  << std::endl;
  }

  if(zeros > 0)
  {
    std::cout << "[compare_hit_points : " << test_name << "] point is not valid in ref or test: " << zeros << std::endl;
  }

  return no_errors;
}


bool compare_t_hit_values(const std::vector<Float>& t_hit_ref, const std::vector<Float>& t_hit_tst,
                          const std::string& test_name, Float threshold, bool verbose)
{
  if(t_hit_ref.size() != t_hit_tst.size())
  {
    std::cout << "[compare_t_hit_values : " << test_name << "] FAILED. Reference and test vectors have different sizes!" << std::endl;
    return false;
  }

  bool no_errors = true;
  double avg_error = 0.0;
  size_t error_val_number = 0;
  for(size_t i = 0; i < t_hit_ref.size(); ++i)
  {
    if(std::abs(t_hit_ref[i]) < threshold)
      continue;
    
    auto absolute_err = std::abs(t_hit_ref[i] - t_hit_tst[i]);

    if (absolute_err > threshold)
    {
      if(verbose)
      {
        std::cout << "[compare_t_hit_values : " << test_name << "] abs. error = " << absolute_err << std::endl;
        std::cout << std::setprecision(16) << i << "\tREF: " << t_hit_ref[i] << std::endl;
        std::cout << std::setprecision(16) << i << "\tTST: " << t_hit_tst[i] << std::endl;
      }
      no_errors = false;
      error_val_number++;
      avg_error += absolute_err;
    }
  }

  if(no_errors)
  {
    std::cout << "[compare_t_hit_values : " << test_name << "] PASSED"  << std::endl;
  }

  if(error_val_number > 0)
  {
    std::cout << "[compare_t_hit_values : " << test_name << "] total values with error > " << threshold << " :" << error_val_number << std::endl;
    std::cout << "[compare_t_hit_values : " << test_name << "] average error: " << avg_error / double(error_val_number) << std::endl;
  }

  return no_errors;
}

void generate_rays(uint32_t num_rays_x, uint32_t num_rays_y, Float obj_dist, Float field_angle, Float entrance_pupil_radius,
                   std::vector<Float3> &ray_origin, std::vector<Float3> &ray_dir)
{
  static constexpr Float TO_RADIANS = 3.14159265358979323846 / 180.0;
  ray_origin.resize(num_rays_x * num_rays_y);
  ray_dir.resize(num_rays_x * num_rays_y);

  Float3 o {0.0, obj_dist * std::tan(field_angle * TO_RADIANS), -obj_dist};

  const double dx = num_rays_x > 1 ? 2 * entrance_pupil_radius / (num_rays_x - 1) : 0;
  const double dy = num_rays_y > 1 ? 2 * entrance_pupil_radius / (num_rays_y - 1) : 0;

  #ifndef _DEBUG
  #pragma omp parallel for
  #endif
  for(uint32_t y = 0; y < num_rays_y; ++y)
  {
    for(uint32_t x = 0; x < num_rays_x; ++x)
    {
      double entrance_x = num_rays_x > 1 ? -entrance_pupil_radius + x * dx : 0;
      double entrance_y = -entrance_pupil_radius + y * dy;
      Float3 entrance_pos = {entrance_x, entrance_y, 0 };

      Float3 d = entrance_pos - o;
      d = normalize(d);

      ray_dir[x + y * num_rays_x] = d;
      ray_origin[x + y * num_rays_x] = o;
    }
  }
}


void generate_rays_batches(uint32_t num_rays_x, uint32_t num_rays_y, std::vector<Float> field_angles, 
                           Float obj_dist, Float entrance_pupil_radius,
                           std::vector<Float3> &ray_origin, std::vector<Float3> &ray_dir)
{
  static constexpr Float TO_RADIANS = 3.14159265358979323846 / 180.0;
  const size_t num_batches = field_angles.size();
  const size_t batch_size = num_rays_x * num_rays_y;

  ray_origin.resize(batch_size * num_batches);
  ray_dir.resize(batch_size * num_batches);

  auto entrance_x = linspace(-entrance_pupil_radius, entrance_pupil_radius, num_rays_x);
  auto entrance_y = linspace(-entrance_pupil_radius, entrance_pupil_radius, num_rays_y);

  for(uint32_t i = 0; i < num_batches; ++i)
  {
    const auto angle = field_angles[i];
    Float3 o {0.0, obj_dist * std::tan(angle * TO_RADIANS), -obj_dist};
  
    #ifndef _DEBUG
    #pragma omp parallel for
    #endif
    for(uint32_t y = 0; y < num_rays_y; ++y)
    {
      for(uint32_t x = 0; x < num_rays_x; ++x)
      {
        Float3 entrance_pos = {entrance_x[x], entrance_y[y], 0};
        
        Float3 d = entrance_pos - o;
        d = normalize(d);

        ray_dir[x + y * num_rays_x + i * batch_size] = d;
        ray_origin[x + y * num_rays_x + i * batch_size] = o;
      }
    }
  }
}

void get_eta_and_abbe_number(const SurfaceParams& surf, const std::vector<Float>& wavelengths, 
                             std::vector<Float>& eta, Float& abbe_num)
{ 
  if(eta.size() != wavelengths.size())
  {
    eta.resize(wavelengths.size());
  }
  std::fill(eta.begin(), eta.end(), 0.0);
    
  // Process each input wavelength
  for (size_t wl_idx = 0; wl_idx < wavelengths.size(); ++wl_idx) 
  {
      double rn_squared_minus_1 = 0.0;
      double wavelength_sq = wavelengths[wl_idx] * wavelengths[wl_idx] * 1.0e-6; // to micrometers
      
      for (size_t i = 0; i < surf.sellmeier_coeffs.size(); i += 2) 
      {
          rn_squared_minus_1 += surf.sellmeier_coeffs[i] * wavelength_sq / (wavelength_sq - surf.sellmeier_coeffs[i + 1]);
      }
      
      eta[wl_idx] = std::sqrt(rn_squared_minus_1 + 1.0);
  }
  
  // Calculate refractive indices at standard F, d, C wavelengths for Abbe number
  std::vector<double> rn_FdC(3, 0.0);
  for (size_t band = 0; band < 3; ++band) 
  {
    double rn_squared_minus_1 = 0.0;
    double wavelength_sq = wave_FdC[band] * wave_FdC[band] * 1.0e-6;
    
    for (size_t i = 0; i < surf.sellmeier_coeffs.size(); i += 2)  
    {
      rn_squared_minus_1 += surf.sellmeier_coeffs[i] * wavelength_sq / (wavelength_sq - surf.sellmeier_coeffs[i + 1]);
    }
    
    rn_FdC[band] = std::sqrt(rn_squared_minus_1 + 1.0);
  }
  
  // Abbe number = (n_d - 1)/(n_F - n_C)
  abbe_num = (rn_FdC[1] - 1.0) / (rn_FdC[0] - rn_FdC[2]);
}


void convert_surface(const SurfaceParams& params, std::vector<double>& surface, double& zpos, double& eta,
                     const Float wavelength)
{
  if(params.use_sellmeier)
  {
    std::vector<Float> eta_vec {0.0};
    Float abbe_num = 0.0;
    get_eta_and_abbe_number(params, { wavelength }, eta_vec, abbe_num);
    eta = eta_vec[0];
  }
  else 
  {
    eta  = params.eta;
  }

  zpos = params.zdist;

  surface.resize(sizeof(PackedParameters)/sizeof(double));
  surface[POLY_APERTURE_ID] = params.aperture_radius;
  surface[POLY_CURVATURE_ID] = params.curvature;
  surface[POLY_CONIC_ID] = params.conic;
  surface[POLY_RNORM_ID] = params.rnorm;
  for(uint32_t i = 0; i < POLY_RPARAMS_MAX_PARAMS; ++i)
    surface[POLY_RPARAMS_FIRST_ID + i] = params.poly_coeffs[i + 1];
}

Statistics computeStatistics(const std::vector<float>& values) 
{
  Statistics stats;
  
  if (values.empty()) {
      return stats;
  }
  
  // Calculate mean
  float sum = std::accumulate(values.begin(), values.end(), 0.0f);
  stats.mean = sum / values.size();
  
  // Calculate variance and standard deviation
  float squared_diff_sum = 0.0f;
  for (float val : values) {
      squared_diff_sum += (val - stats.mean) * (val - stats.mean);
  }
  
  stats.variance = squared_diff_sum / values.size();
  stats.stddev = std::sqrt(stats.variance);
  
  return stats;
}
