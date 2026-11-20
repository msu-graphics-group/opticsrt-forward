#include "even_asphere.h"
#include "utils.h"

#include <vector>
#include "isurface.h"
#include "surface.h"
#include "poly_surface.h"

#include "LiteMath.h"

#include "tests.h"

#include <sstream>
#include <functional>
#include <iostream>
#include <fstream>

using LiteMath::float4;
using LiteMath::double3;
using LiteMath::length;
using LiteMath::dot;
using LiteMath::normalize;

extern bool ConvertToJacobi(std::vector<double>& surface_params);

static constexpr bool CHECK_REF = false;

template<typename T>
void CompareBuffer(const std::string& a_name, const std::vector<T>& a_data1,  const std::vector<T>& a_data2)
{
  if(a_data1.size() != a_data2.size()) {
    std::cout << "CompareBuffer("<< a_name.c_str() << "); different size: " << a_data1.size() << " != " << a_data2.size();
    return;
  }

  for(size_t i=0;i<a_data1.size();i++) {
    if(std::abs(a_data1[i] - a_data2[i]) > 1e-10) {
      std::cout << "CompareBuffer("<< a_name.c_str() << "); different value at(" << i << "): " << a_data1[i] << " != " << a_data2[i] << std::endl;
      return;
    }
  }
}

template<>
void CompareBuffer<uint8_t>(const std::string& a_name, const std::vector<uint8_t>& a_data1,  const std::vector<uint8_t>& a_data2)
{
  if(a_data1.size() != a_data2.size()) {
    std::cout << "CompareBuffer("<< a_name.c_str() << "); different size: " << a_data1.size() << " != " << a_data2.size();
    return;
  }

  for(size_t i=0;i<a_data1.size();i++) {
    if(std::abs(a_data1[i] - a_data2[i]) > 1e-10) {
      std::cout << "CompareBuffer("<< a_name.c_str() << "); different value at(" << i << "): " << int(a_data1[i]) << " != " << int(a_data2[i]) << std::endl;
      return;
    }
  }
}

template<>
void CompareBuffer<double3>(const std::string& a_name, const std::vector<double3>& a_data1,  const std::vector<double3>& a_data2)
{
  if(a_data1.size() != a_data2.size()) {
    std::cout << "CompareBuffer("<< a_name.c_str() << "); different size: " << a_data1.size() << " != " << a_data2.size();
    return;
  }

  for(size_t i=0;i<a_data1.size();i++) {
    auto diff = a_data1[i] - a_data2[i];
    if(std::abs(dot(diff,diff)) > 1e-20) {
      std::cout << "CompareBuffer("<< a_name.c_str() << "); different value at(" << i << "): (" << a_data1[i].x << "," << a_data1[i].y << "," << a_data1[i].z << ") != ";
      std::cout <<                                                                          "(" << a_data2[i].x << "," << a_data2[i].y << "," << a_data2[i].z << ")" << std::endl;
      return;
    }
  }
}

void output_performance_stats(const std::filesystem::path& out_dir, Timings& timings, uint32_t num_surfaces)
{
  float div = static_cast<float>(num_surfaces);

  std::stringstream perf_ss;
  perf_ss << "Forward tracing average execution times : " << std::endl;
  perf_ss << "\tTraceRays exec  = "        << timings.traceTimings.x / div << " ms " << std::endl;
  perf_ss << "\tData copy to device  = "   << timings.traceTimings.y / div << " ms " << std::endl;
  perf_ss << "\tData copy from device  = " << timings.traceTimings.w / div << " ms " << std::endl;

  std::cout << perf_ss.str();

  std::filesystem::create_directories(out_dir);
  auto out_stream = std::ofstream((out_dir / "avg_performance.txt").string(), std::ios_base::trunc);
  if(!out_stream.bad())
  {
    out_stream << perf_ss.str();
  }
}

struct TestFunc
{
  bool(*func)(const std::shared_ptr<ISurface>&, const std::shared_ptr<ISurface>&, const std::filesystem::path&);
  std::string name;

  TestFunc(bool(*a_func)(const std::shared_ptr<ISurface>&, const std::shared_ptr<ISurface>&, const std::filesystem::path&), std::string a_name)
  {
    func = a_func;
    name = a_name;
  }
};

bool run_all_tests(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir)
{
  bool passed = true;
  std::vector<TestFunc> tests =
  {
    { &test_default_surface,         "test_default_surface" },
    { &test_debug_surface,           "test_debug_surface" },
    { &test_double_gauss_singleray,  "test_double_gauss_singleray" },
    { &test_nonsquare_ray_batch,     "test_nonsquare_ray_batch" },
    { &test_two_batches,             "test_two_batches"},
    { &test_batches_sequential,      "test_batches_sequential"},
    { &test_intersect_behind,        "test_intersect_behind"},
    { &test_raymarch_matches_interval, "test_raymarch_matches_interval"},
  };

  std::vector<std::string> failed_tests;
  for(const auto& test: tests)
  {
    auto res = test.func(pSim, pSimRef, out_dir);
    if(!res)
      failed_tests.push_back(test.name);

    passed = passed && res;
  }

  if(passed)
  {
    std::cout << "--------------------" << std::endl << "All tests PASSED!" << std::endl;
  }
  else
  {
    std::cout << "--------------------" << std::endl << "Some tests FAILED:" << std::endl;
    for(const auto& test_name: failed_tests)
    {
      std::cout << "\t" << test_name << std::endl;
    }
  }

  return passed;
}


bool test_default_surface(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string experiment_name = "default_surface";
  auto save_path = out_dir;
  save_path.append(experiment_name);

  std::cout << "TEST: " << experiment_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.hitEps = 1e-7;
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;
  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 10;
  algorithm_params_init.subdivideRaySegment = 10;

  uint32_t num_rays = 512;
  double obj_dist = 1000.0;
  double field_angle = 5.0;

  SurfaceParams default_surface;
  // trace generated rays
  {
    std::vector<double3>  ray_origin;
    std::vector<double3>  ray_dir;

    generate_rays(num_rays, num_rays, obj_dist, field_angle,
      default_surface.aperture_radius, ray_origin, ray_dir);

    size_t total_rays = num_rays * num_rays;
    std::vector<uint8_t> ray_mask(total_rays, 1);

    std::vector<double> zpos = {0.0};
    std::vector<double> eta = {1.0};
    std::vector<double> surface_params;

    convert_surface(default_surface, surface_params, zpos[0], eta[0]);

    const PolySurface* pSimPolySurf = dynamic_cast<const PolySurface*>(pSim.get());
    if (pSimPolySurf != NULL)
      {
      if (pSimPolySurf->UsesJacobi())
        ConvertToJacobi(surface_params);
      }

    // The number of surface parameters used, for each surface
    uint32_t n_surface_params = POLY_RPARAMS_MAX_PARAMS;

    pSim->SetInitIntersectionParams(algorithm_params_init);
    pSim->UpdateSurface(surface_params.data(), n_surface_params);
    pSim->CommitDeviceData();

    std::vector<Float3> out_hit_pos(total_rays, {0., 0., 0.});
    std::vector<Float3> out_hit_dir(total_rays, {0., 0., 0.});
    std::vector<Float3> out_hit_normal(total_rays, {0., 0., 0.});
    std::vector<Float>   out_t_hit(total_rays, 0.);
    std::vector<uint8_t> out_t_valid(total_rays, 1);
    std::vector<HitPointData> hpd(total_rays, make_hpd());

    uint32_t BATCHES_NUMBER = 1;


    auto timings = measure_performance(pSim, {algorithm_params, zpos.data(), surface_params.data(),
                                              n_surface_params, eta.data(),
                                              (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                                              (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                                              hpd.data(),
                                              512u, 512u, BATCHES_NUMBER, uint32_t(zpos.size())});

    auto ref_t_hit_points = load_t_hit_txt(512, 512, "ref_data/test_surface_t_hit.txt");
    passed = compare_t_hit_values(ref_t_hit_points, out_t_hit, "default surface t_hit", algorithm_params.hitEps, true);

    save_t_hit_image(save_path.string(), out_t_hit, out_t_valid, obj_dist, num_rays, num_rays);

    output_performance_stats(out_dir, timings, 1);
  }

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}

bool test_nonsquare_ray_batch(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string test_name = "nonsquare_ray_batch_103";
  std::filesystem::path test_data = "ref_data/nonsquare_ray_batch_103";

  std::cout << "TEST: " << test_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;
  algorithm_params.useFloat = true;
  algorithm_params.rebasePoly = true;
  algorithm_params.hitEps = algorithm_params.useFloat ? 1e-6 : 2e-9;

  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 24;
  algorithm_params_init.subdivideRaySegment = 1;

  pSim->SetInitIntersectionParams(algorithm_params_init);
  pSim->CommitDeviceData();

  uint32_t num_rays_x = 128;
  uint32_t num_rays_y = 96;

  std::vector<SurfaceParams> surfaces;
  load_surfaces_from_json("103.json", surfaces, true, false);

  Float obj_dist = 1000.;
  Float field_angle = 5.0;

  std::vector<uint8_t> ray_mask(num_rays_x * num_rays_y, 1);
  std::vector<Float3>  out_hit_pos(num_rays_x * num_rays_y, {0., 0., 0.});
  std::vector<Float3>  out_hit_dir(num_rays_x * num_rays_y, {0., 0., 0.});
  std::vector<Float3>  out_hit_normal(num_rays_x * num_rays_y, {0., 0., 0.});
  std::vector<Float>   out_t_hit(num_rays_x * num_rays_y, 0.);
  std::vector<uint8_t> out_t_valid(num_rays_x * num_rays_y, 1);

  std::vector<HitPointData> hpd(num_rays_x * num_rays_y, make_hpd());

  std::vector<double> eta_prev = { 1.0 }; // BEFORE the surface
  for(uint32_t surface_idx = 0; surface_idx < surfaces.size(); ++surface_idx)
  {
    SurfaceParams surf = surfaces[surface_idx];

    std::cout << "INTERSECTING SURFACE: " << surf.id << std::endl;

    std::stringstream ss;
    ss << test_name << "_" << surf.id << "_" << num_rays_x << "X" << num_rays_y << "rays_" << field_angle << "deg";
    std::string experiment_name = ss.str();

    auto save_path = out_dir;
    save_path.append(experiment_name);

    std::vector<Float3>  ray_origin;
    std::vector<Float3>  ray_dir;

    generate_rays(num_rays_x, num_rays_y, obj_dist, field_angle,
      surf.aperture_radius, ray_origin, ray_dir);

    std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
    std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
    std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
    std::fill(out_t_hit.begin(), out_t_hit.end(), 0.0);
    std::fill(out_t_valid.begin(), out_t_valid.end(), 1);

    std::vector<double> zpos = {0.0};
    std::vector<double> eta = { 1.0 }, eta_inv = { 1.0 };
    std::vector<double> surface_params;
    convert_surface(surf, surface_params, zpos[0], eta[0]);

    eta_inv[0] = 1;

    // The number of polynomials parameters used, for each surface
    uint32_t used_raparms = POLY_RPARAMS_MAX_PARAMS;

    pSim->SetIntersectionParams(algorithm_params);
    pSim->UpdateSurface(surface_params.data(), used_raparms);

    pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta_inv.data(),
                    (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                    (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                    hpd.data(),
                    num_rays_x, num_rays_y, 1, zpos.size());

    if(std::filesystem::is_directory(test_data))
    {
      std::filesystem::path ref_path = test_data;

      ref_path.append("103_" + std::to_string(surf.id) + "_t_hit.txt");
      if(std::filesystem::is_regular_file(ref_path))
      {
        auto ref_t_hit_points = load_t_hit_txt(num_rays_x, num_rays_y, ref_path.string());
        passed = compare_t_hit_values(ref_t_hit_points, out_t_hit, ref_path.string() + "_TEST", algorithm_params.hitEps, false);
      }
      // save_t_hit_txt(out_t_hit, args.num_rays, args.num_rays, save_path.string());
    }

    {
      save_t_hit_image(save_path.string(), out_t_hit, out_t_valid, obj_dist, num_rays_x, num_rays_y);
    }

    save_mask_to_png(save_path.string(), num_rays_x, num_rays_y, 1, out_t_valid);

    eta_prev = eta;
  }

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}


bool test_two_batches(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string test_name = "two_batches_103";

  std::filesystem::path test_data_1 = "ref_data/two_batches_103_10deg";
  std::filesystem::path test_data_2 = "ref_data/two_batches_103_15deg";
  Float field_angle_1 = 10.0;
  Float field_angle_2 = 15.0;


  std::cout << "TEST: " << test_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;
  algorithm_params.useFloat = true;
  algorithm_params.rebasePoly = true;
  algorithm_params.hitEps = algorithm_params.useFloat ? 1e-6 : 2e-9;

  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 24;
  algorithm_params_init.subdivideRaySegment = 1;

  pSim->SetInitIntersectionParams(algorithm_params_init);
  pSim->CommitDeviceData();

  uint32_t num_rays_x = 64;
  uint32_t num_rays_y = 64;
  uint32_t num_batches = 2;

  std::vector<SurfaceParams> surfaces;
  load_surfaces_from_json("103.json", surfaces, true, false);

  Float obj_dist = 1000.;

  std::vector<uint8_t> ray_mask(num_rays_x * num_rays_y * num_batches, 1);
  std::vector<Float3>  out_hit_pos(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
  std::vector<Float3>  out_hit_dir(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
  std::vector<Float3>  out_hit_normal(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
  std::vector<Float>   out_t_hit(num_rays_x * num_rays_y * num_batches, 0.);
  std::vector<uint8_t> out_t_valid(num_rays_x * num_rays_y * num_batches, 1);
  std::vector<HitPointData> hpd(num_rays_x * num_rays_y * num_batches, make_hpd());

  std::vector<double> eta_prev = { 1.0 }; // BEFORE the surface
  for(uint32_t surface_idx = 0; surface_idx < surfaces.size(); ++surface_idx)
  {
    SurfaceParams surf = surfaces[surface_idx];

    std::cout << "INTERSECTING SURFACE: " << surf.id << std::endl;

    std::stringstream ss;
    ss << test_name << "_" << surf.id << "_" << num_rays_x << "X" << num_rays_y << "rays_"
       << field_angle_1 << "deg" << "_" << field_angle_2 << "deg";
    std::string experiment_name = ss.str();

    auto save_path = out_dir;
    save_path.append(experiment_name);

    std::vector<Float3>  ray_origin_1;
    std::vector<Float3>  ray_dir_1;
    std::vector<Float3>  ray_origin_2;
    std::vector<Float3>  ray_dir_2;

    generate_rays(num_rays_x, num_rays_y, obj_dist, field_angle_1,
      surf.aperture_radius, ray_origin_1, ray_dir_1);

    generate_rays(num_rays_x, num_rays_y, obj_dist, field_angle_2,
      surf.aperture_radius, ray_origin_2, ray_dir_2);

    ray_origin_1.reserve(ray_origin_1.size() + ray_origin_2.size());
    ray_origin_1.insert(ray_origin_1.end(), ray_origin_2.begin(), ray_origin_2.end());

    ray_dir_1.reserve(ray_dir_1.size() + ray_dir_2.size());
    ray_dir_1.insert(ray_dir_1.end(), ray_dir_2.begin(), ray_dir_2.end());

    std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
    std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
    std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
    std::fill(out_t_hit.begin(), out_t_hit.end(), 0.0);
    std::fill(out_t_valid.begin(), out_t_valid.end(), 1);

    std::vector<double> zpos = {0.0 , 0.0 };
    std::vector<double> eta = { 1.0, 1.0 }, eta_inv = { 1.0, 1.0 };
    std::vector<double> surface_params;
    convert_surface(surf, surface_params, zpos[0], eta[0]);
    zpos[1] = zpos[0];


    // The number of polynomials parameters used, for each surface
    uint32_t used_raparms = POLY_RPARAMS_MAX_PARAMS;

    pSim->UpdateSurface(surface_params.data(), used_raparms);

    pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta_inv.data(),
                    (const double*)ray_origin_1.data(), (const double*)ray_dir_1.data(), ray_mask.data(),
                    (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                    hpd.data(),
                    num_rays_x, num_rays_y, num_batches, zpos.size());

    if(std::filesystem::is_directory(test_data_1))
    {
      std::filesystem::path ref_path = test_data_1;

      ref_path.append("103_" + std::to_string(surf.id) + "_t_hit.txt");
      if(std::filesystem::is_regular_file(ref_path))
      {
        std::vector<Float> out_t_hit_1(out_t_hit.begin(), out_t_hit.begin() + num_rays_x * num_rays_y);
        auto ref_t_hit_points = load_t_hit_txt(num_rays_x, num_rays_y, ref_path.string());
        passed = compare_t_hit_values(ref_t_hit_points, out_t_hit_1, ref_path.string(), algorithm_params.hitEps, false);
      }
      // save_t_hit_txt(out_t_hit, args.num_rays, args.num_rays, save_path.string());
    }

    if(std::filesystem::is_directory(test_data_2))
    {
      std::filesystem::path ref_path = test_data_2;

      ref_path.append("103_" + std::to_string(surf.id) + "_t_hit.txt");
      if(std::filesystem::is_regular_file(ref_path))
      {
        std::vector<Float> out_t_hit_2(out_t_hit.begin() + num_rays_x * num_rays_y, out_t_hit.end());
        auto ref_t_hit_points = load_t_hit_txt(num_rays_x, num_rays_y, ref_path.string());
        passed = passed && compare_t_hit_values(ref_t_hit_points, out_t_hit_2, ref_path.string(), algorithm_params.hitEps, false);
      }
      // save_t_hit_txt(out_t_hit, args.num_rays, args.num_rays, save_path.string());
    }

    // {
    //   save_t_hit_image(save_path.string(), out_t_hit, out_t_valid, obj_dist, num_rays_x, num_rays_y);
    // }

    // save_mask_to_png(save_path.string(), num_rays_x, num_rays_y, 1, out_t_valid);

  }

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}

bool test_batches_sequential(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string test_name = "batches_sequential";

  std::vector<Float> field_angles = {0.0, 5.0, 10.0, 20.0};

  std::cout << "TEST: " << test_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.hitEps = 2e-9;
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;

  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 24;
  algorithm_params_init.subdivideRaySegment = 1;

  pSim->SetInitIntersectionParams(algorithm_params_init);
  pSim->CommitDeviceData();

  const uint32_t num_rays_x  = 64;
  const uint32_t num_rays_y  = 32;
  const size_t   num_batches = field_angles.size();
  const Float    obj_dist    = 1000.;
  const uint32_t used_raparms = POLY_RPARAMS_MAX_PARAMS;

  std::vector<SurfaceParams> surfaces;
  load_surfaces_from_json("103.json", surfaces, false, false);

  // {
  //   std::vector<Float3>  ray_origin;
  //   std::vector<Float3>  ray_dir;
  //   generate_rays_batches(num_rays_x, num_rays_y, field_angles, obj_dist, surfaces[0].aperture_radius, ray_origin, ray_dir);

  //   std::vector<Float3>  ray_origin_seq;
  //   std::vector<Float3>  ray_dir_seq;
  //   for(int i = 0; i < field_angles.size(); ++i)
  //   {
  //     const Float angle = field_angles[i];
  //     std::vector<Float3>  ray_origin_tmp;
  //     std::vector<Float3>  ray_dir_tmp;
  //     generate_rays(num_rays_x, num_rays_y, obj_dist, angle, surfaces[0].aperture_radius,
  //                   ray_origin_tmp, ray_dir_tmp);
  //     // std::copy(ray_origin_tmp.begin(), ray_origin_tmp.end(), ray_origin_seq.begin() + i * num_rays_x * num_rays_y);
  //     // std::copy(ray_dir_tmp.begin(), ray_dir_tmp.end(), ray_dir_seq.begin() + i * num_rays_x * num_rays_y);
  //     std::copy(ray_dir_tmp.begin(), ray_dir_tmp.end(), std::back_inserter(ray_dir_seq));
  //     std::copy(ray_origin_tmp.begin(), ray_origin_tmp.end(), std::back_inserter(ray_origin_seq));
  //   }
  //   bool a = areVectorsAlmostEqual(ray_dir, ray_dir_seq, 1e-8);
  //   bool b = areVectorsAlmostEqual(ray_origin, ray_origin_seq, 1e-8);

  //   std::cout << "Generated rays are " << (a == b) << std::endl;
  // }


  std::vector<Float>  out_t_hit(num_rays_x * num_rays_y * num_batches, 0.);

  std::cout << "Intersecting surfaces by all batches at once..." << std::endl;
  // trace all batches
  {
    std::vector<uint8_t> ray_mask(num_rays_x * num_rays_y * num_batches, 1);
    std::vector<Float3>  out_hit_dir(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
    std::vector<Float3>  out_hit_normal(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
    std::vector<uint8_t> out_t_valid(num_rays_x * num_rays_y * num_batches, 1);
    std::vector<Float3>  out_hit_pos(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
    std::vector<HitPointData> hpd(num_rays_x * num_rays_y * num_batches, make_hpd());
    std::vector<Float3>  ray_origin;
    std::vector<Float3>  ray_dir;

    std::vector<double> eta_prev(num_batches, 1.0); // BEFORE the surface
    for(uint32_t surface_idx = 0; surface_idx < surfaces.size(); ++surface_idx)
    {
        SurfaceParams surf = surfaces[surface_idx];

      std::stringstream ss;
      ss << test_name << "_" << surf.id << "_" << num_rays_x << "X" << num_rays_y << "rays_"
         << field_angles.size() << "batches";
      std::string experiment_name = ss.str();

      auto save_path = out_dir;
      save_path.append(experiment_name);

      if(surface_idx == 0)
      {
        generate_rays_batches(num_rays_x, num_rays_y, field_angles, obj_dist, surf.aperture_radius, ray_origin, ray_dir);
      }
      else
      {
        std::copy(out_hit_pos.begin(), out_hit_pos.end(), ray_origin.begin());
        std::copy(out_hit_dir.begin(), out_hit_dir.end(), ray_dir.begin());
        std::copy(out_t_valid.begin(), out_t_valid.end(), ray_mask.begin());
      }

      std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
      std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
      std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
      std::fill(out_t_hit.begin(), out_t_hit.end(), 0.0);
      std::fill(out_t_valid.begin(), out_t_valid.end(), 1);


      std::vector<double> zpos = {0.0 /*, 0.0 */};
      std::vector<double> eta(num_batches, 1.0);
      std::vector<double> eta_inv(num_batches, 1.0);
      std::vector<double> surface_params;

      double eta_surf = 0.0;
      convert_surface(surf, surface_params, zpos[0], eta_surf);

      std::fill(eta.begin(), eta.end(), eta_surf);

      // zpos[1] = zpos[0];

      for(int e = 0; e < eta.size(); ++e)
      {
        eta_inv[e] = eta_prev[e] / eta[e];
      }

      pSim->UpdateSurface(surface_params.data(), used_raparms);

      pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta_inv.data(),
                      (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                      (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                      hpd.data(), num_rays_x, num_rays_y, num_batches, zpos.size());

      eta_prev = eta;
    }
  }

  std::vector<Float>  out_t_hit_seq(num_rays_x * num_rays_y * num_batches, 0.);
  // std::vector<Float3> out_hit_pos_seq(num_rays_x * num_rays_y * num_batches, {0., 0., 0.});
  std::vector<Float3> out_path_batch_1_ray_21_seq(surfaces.size(), {0., 0., 0.});
  // trace batch by batch
  for(int i = 0; i < num_batches; ++i)
  {
    Float field_angle = field_angles[i];

    std::vector<uint8_t> ray_mask(num_rays_x * num_rays_y * 1, 1);
    std::vector<Float3>  out_hit_dir(num_rays_x * num_rays_y * 1, {0., 0., 0.});
    std::vector<Float3>  out_hit_normal(num_rays_x * num_rays_y * 1, {0., 0., 0.});
    std::vector<uint8_t> out_t_valid(num_rays_x * num_rays_y * 1, 1);

    std::vector<Float3>  out_hit_pos(num_rays_x * num_rays_y * 1, {0., 0., 0.});
    std::vector<Float>   out_t_hit_tmp(num_rays_x * num_rays_y, 0.);
    std::vector<HitPointData> hpd(num_rays_x * num_rays_y, make_hpd());
    std::vector<Float3>  ray_origin;
    std::vector<Float3>  ray_dir;

    std::cout << "Intersecting surfaces batch by batch: " << i << ", angle = " << field_angle << std::endl;
    std::vector<double> eta_prev = { 1.0 }; // BEFORE the surface
    for(uint32_t surface_idx = 0; surface_idx < surfaces.size(); ++surface_idx)
    {
        SurfaceParams surf = surfaces[surface_idx];

      std::stringstream ss;
      ss << test_name << "_" << surf.id << "_" << num_rays_x << "X" << num_rays_y << "rays_"
         << field_angle << "deg";
      std::string experiment_name = ss.str();

      auto save_path = out_dir;
      save_path.append(experiment_name);

      if(surface_idx == 0)
      {
        generate_rays(num_rays_x, num_rays_y, obj_dist, field_angle, surf.aperture_radius,
                   ray_origin, ray_dir);
      }
      else
      {
        std::copy(out_hit_pos.begin(), out_hit_pos.end(), ray_origin.begin());
        std::copy(out_hit_dir.begin(), out_hit_dir.end(), ray_dir.begin());
        std::copy(out_t_valid.begin(), out_t_valid.end(), ray_mask.begin());
      }

      std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
      std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
      std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
      std::fill(out_t_hit_tmp.begin(), out_t_hit_tmp.end(), 0.0);
      std::fill(out_t_valid.begin(), out_t_valid.end(), 1);


      std::vector<double> zpos = {0.0 /*, 0.0 */};
      std::vector<double> eta = { 1.0 }, eta_inv = { 1.0 };
      std::vector<double> surface_params;
      convert_surface(surf, surface_params, zpos[0], eta[0]);
      // zpos[1] = zpos[0];

      eta_inv[0] = eta_prev[0] / eta[0];

      pSim->UpdateSurface(surface_params.data(), used_raparms);

      pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta_inv.data(),
                      (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                      (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit_tmp.data(),
                      (double*)out_hit_normal.data(), hpd.data(), num_rays_x, num_rays_y, 1, zpos.size());
      eta_prev = eta;
    }

    std::copy(out_t_hit_tmp.begin(), out_t_hit_tmp.end(), out_t_hit_seq.begin() + i * num_rays_x * num_rays_y);
  }

  passed = compare_t_hit_values(out_t_hit_seq, out_t_hit, test_name, algorithm_params.hitEps, true);

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}


bool test_debug_surface(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string experiment_name = "debug_surface";
  auto save_path = out_dir;
  save_path.append(experiment_name);

  std::cout << "TEST: " << experiment_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;
  algorithm_params.useFloat   = true;
  algorithm_params.rebasePoly = !true;
  algorithm_params.hitEps = algorithm_params.useFloat ? 1.1e-6 : 3e-9;


  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 24;
  algorithm_params_init.subdivideRaySegment = 10;

  pSim->SetInitIntersectionParams(algorithm_params_init);
  pSim->CommitDeviceData();

  uint32_t num_rays = 100;

  // trace generated rays
  {
    std::vector<double3>  ray_origin = load_points_txt(num_rays, num_rays, "ref_data/debug_ray_pos.txt");
    std::vector<double3>  ray_dir = load_points_txt(num_rays, num_rays, "ref_data/debug_ray_dir.txt");

    size_t total_rays = num_rays * num_rays;
    std::vector<uint8_t> ray_mask(total_rays, 1);

    std::vector<double> zpos = {0.0};
    std::vector<double> eta = {0.5988};
    std::vector<double> surface_params = {6.3800e+00, 0.17967849548140979, -1.2116704598129759,
                                          6.552198202737241, -74.855935601961348, 863.47157089691007, -7975.9241319766843,
                                          55309.020087541918, -284929.32741243445, 1085886.4212789999, -3062225.1028069998,
                                          6386243.4709989997, -9791775.9918930009, 10874632.94327, -8496288.3587129992,
                                          4422713.6429850003,-1375665.938911, 193280.54386899999};

    const PolySurface* pSimPolySurf = dynamic_cast<const PolySurface*>(pSim.get());
    if (pSimPolySurf != NULL)
      {
      if (pSimPolySurf->UsesJacobi())
        ConvertToJacobi(surface_params);
      }

    // convert_surface(surf, surface_params, zpos[0], eta[0]);

    uint32_t n_surface_params = POLY_RPARAMS_MAX_PARAMS;

    std::vector<Float3> out_hit_pos(total_rays, {0., 0., 0.});
    std::vector<Float3> out_hit_dir(total_rays, {0., 0., 0.});
    std::vector<Float3> out_hit_normal(total_rays, {0., 0., 0.});
    std::vector<Float>   out_t_hit(total_rays, 0.);
    std::vector<uint8_t> out_t_valid(total_rays, 1);
    std::vector<HitPointData> hpd(total_rays, make_hpd());

    pSim->SetIntersectionParams(algorithm_params);
    pSim->UpdateSurface(surface_params.data(), n_surface_params);
    pSim->CommitDeviceData();

    //pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(),
    //                n_surface_params,  eta.data(),
    //                (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
    //                (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
    //                num_rays, 1u);

    auto timings = measure_performance(pSim, {algorithm_params, zpos.data(), surface_params.data(),
                                             n_surface_params, eta.data(),
                                             (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                                             (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                                             hpd.data(),
                                             num_rays, num_rays, 1u, uint32_t(zpos.size())});

    if(CHECK_REF)
    {
      auto out_hit_pos2 = out_hit_pos;
      auto out_hit_dir2 = out_hit_dir;
      auto out_t_valid2 = out_t_valid;
      auto out_hit_normal2 = out_hit_normal;

      std::fill(out_hit_pos.begin(), out_hit_pos.end(),       double3(0.0,0.0,0.0));
      std::fill(out_hit_dir.begin(), out_hit_dir.end(),       double3(0.0,0.0,0.0));
      std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3(0.0,0.0,0.0));
      std::fill(out_t_valid.begin(), out_t_valid.end(), 0);

      pSimRef->UpdateSurface(surface_params.data(), n_surface_params);
      pSimRef->CommitDeviceData();

      pSimRef->TraceRays(algorithm_params, zpos.data(), surface_params.data(),
                         n_surface_params,  eta.data(),
                         (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                         (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                         hpd.data(),
                         num_rays, num_rays, 1u, zpos.size());

      CompareBuffer("out_hit_pos", out_hit_pos, out_hit_pos2);
      CompareBuffer("out_hit_dir", out_hit_dir, out_hit_dir2);
      CompareBuffer("out_t_valid", out_t_valid, out_t_valid2);
      CompareBuffer("out_hit_normal", out_hit_normal, out_hit_normal2);
    }


    auto ref_mask = ray_mask_from_png("ref_data/debug_ref.png");

    uint32_t total_nonequal = 0u;
    for(uint32_t i = 0; i < num_rays * num_rays; ++i)
    {
      if(ref_mask[i] != 0 && out_t_valid[i] == 0)
      {
        std::cout << "[test_debug_surface]: point valid in reference but not valid in test, point index: " << i << std::endl;
        total_nonequal += 1;
      }
    }

    auto ref_t_hit = load_t_hit_txt(num_rays, num_rays, "ref_data/debug_surface_t_hit.txt");
#if 1 // ose
    for (uint32_t i = 0; i < num_rays * num_rays; ++i)
      {
      if (!ref_mask[i])
        ref_t_hit[i] = -1;
      }
#endif
    passed = compare_t_hit_values(ref_t_hit, out_t_hit, "debug surface t_hit", algorithm_params.hitEps, false);

    if(total_nonequal > 0)
    {
      std::cout << "[test_debug_surface]: FAILED!" << std::endl;
      std::cout << "\tTotal points valid in reference and not valid in test: " << total_nonequal << std::endl;
    }
    else
    {
      std::cout << "[test_debug_surface]: PASSED!" << std::endl;
    }


    save_t_hit_image(save_path.string(), out_t_hit, out_t_valid, -ray_origin[0].z, num_rays, num_rays);

    save_mask_to_png(save_path.string(), num_rays, num_rays, 1, out_t_valid);
  }

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}


bool test_intersect_behind(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string experiment_name = "intersect_behind";
  auto save_path = out_dir;
  save_path.append(experiment_name);

  std::cout << "TEST: " << experiment_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.hitEps = 3e-9;
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;
  algorithm_params.minT   = -0.5001;

  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 24;
  algorithm_params_init.subdivideRaySegment = 1;

  pSim->SetInitIntersectionParams(algorithm_params_init);
  pSim->CommitDeviceData();

  uint32_t num_rays = 1;

  std::vector<Float3>  ray_origin = {Float3(0.0, 0.0, 0.5)};
  std::vector<Float3>  ray_dir = {Float3(0.0, 0.0, 1.0)};

  std::vector<uint8_t> ray_mask(num_rays * num_rays, 1);
  std::vector<Float3>  out_hit_pos(num_rays * num_rays, {0., 0., 0.});
  std::vector<Float3>  out_hit_dir(num_rays * num_rays, {0., 0., 0.});
  std::vector<Float3>  out_hit_normal(num_rays * num_rays, {0., 0., 0.});
  std::vector<Float>   out_t_hit(num_rays * num_rays, 0.);
  std::vector<uint8_t> out_t_valid(num_rays * num_rays, 1);
  std::vector<HitPointData> hpd(num_rays * num_rays, make_hpd());

  SurfaceParams surf;
  surf.aperture_radius = 10;
  surf.conic = 0.0;
  surf.curvature = 0.0;
  surf.rnorm = 10.0;
  std::fill(surf.poly_coeffs.begin(), surf.poly_coeffs.end(), 0.0);

  {
    std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
    std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
    std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
    std::fill(out_t_hit.begin(), out_t_hit.end(), 0.0);
    std::fill(out_t_valid.begin(), out_t_valid.end(), 1);

    std::vector<double> zpos = {0.0};
    std::vector<double> eta = { 1.0 };
    std::vector<double> surface_params;

    convert_surface(surf, surface_params, zpos[0], eta[0]);
    // std::cout << "eta = " << eta[0] << std::endl;

    uint32_t used_raparms = POLY_RPARAMS_MAX_PARAMS;

    pSim->UpdateSurface(surface_params.data(), used_raparms);

    pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta.data(),
                    (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                    (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                    hpd.data(), 1, 1, 1, zpos.size());
  }

  passed = out_t_valid[0] && out_t_hit[0] < 0.0;

  // ray_origin[0].y = 1.0;
  // ray_origin[0].z = 0.6;
  ray_dir[0] = normalize(Float3(0.0, 1.0, 1.0));
  {
    std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
    std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
    std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
    std::fill(out_t_hit.begin(), out_t_hit.end(), 0.0);
    std::fill(out_t_valid.begin(), out_t_valid.end(), 1);

    std::vector<double> zpos = {0.0};
    std::vector<double> eta = { 1.0 };
    std::vector<double> surface_params;

    convert_surface(surf, surface_params, zpos[0], eta[0]);
    // std::cout << "eta = " << eta[0] << std::endl;

    uint32_t used_raparms = POLY_RPARAMS_MAX_PARAMS;

    pSim->UpdateSurface(surface_params.data(), used_raparms);

    pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta.data(),
                    (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                    (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                    hpd.data(), 1, 1, 1, zpos.size());
  }

  passed = passed && (!out_t_valid[0]) && (out_t_hit[0] < 0.0);

  // std::vector<Float3> ref_t_hit_point = {Float3(0.000000000000000E+00, 9.930611122611568E-01, 1.329884193450000E+01)};
  // passed = compare_hit_points(ref_t_hit_point, out_hit_pos, experiment_name + "_TEST", algorithm_params.hitEps, false, false);

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}


bool test_double_gauss_singleray(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path& out_dir)
{
  bool passed = false;
  std::string experiment_name = "double_gauss_singleray";

  std::cout << "TEST: " << experiment_name << std::endl;

  IntersectionParameters     algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.hitEps = 2e-9;
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;


  IntersectionInitParameters algorithm_params_init = lens_sim::MakeIntersectionInitParameters();
  algorithm_params_init.maxIterations = 24;
  algorithm_params_init.subdivideRaySegment = 1;

  pSim->SetInitIntersectionParams(algorithm_params_init);
  pSim->CommitDeviceData();

  uint32_t num_rays = 1;

  std::vector<SurfaceParams> surfaces;
  load_surfaces_from_json("Double Gauss 28 degree field.json", surfaces, false, true);

  std::vector<Float3>  ray_origin = {Float3(0.000000000000000E+00, 0.000000000000000E+00, -1.000000000000000E+00)};
  std::vector<Float3>  ray_dir = {Float3(0.000000000000000E+00, 9.950371902099893E-02, 9.950371902099893E-01)};

  std::vector<uint8_t> ray_mask(num_rays * num_rays, 1);
  std::vector<Float3>  out_hit_pos(num_rays * num_rays, {0., 0., 0.});
  std::vector<Float3>  out_hit_dir(num_rays * num_rays, {0., 0., 0.});
  std::vector<Float3>  out_hit_normal(num_rays * num_rays, {0., 0., 0.});
  std::vector<Float>   out_t_hit(num_rays * num_rays, 0.);
  std::vector<uint8_t> out_t_valid(num_rays * num_rays, 1);
  std::vector<HitPointData> hpd(num_rays * num_rays, make_hpd());

  std::vector<double> eta_prev = { 1.0 }; // BEFORE the surface
  for(uint32_t surface_idx = 0; surface_idx < surfaces.size(); ++surface_idx)
  {
    const auto &surf = surfaces[surface_idx];
    // std::cout << "INTERSECTING SURFACE: " << surf.id << std::endl;

    if(surface_idx != 0)
    {
      std::copy(out_hit_pos.begin(), out_hit_pos.end(), ray_origin.begin());
      std::copy(out_hit_dir.begin(), out_hit_dir.end(), ray_dir.begin());
      std::copy(out_t_valid.begin(), out_t_valid.end(), ray_mask.begin());
    }

    std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
    std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
    std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
    std::fill(out_t_hit.begin(), out_t_hit.end(), 0.0);
    std::fill(out_t_valid.begin(), out_t_valid.end(), 1);

    std::vector<double> zpos = {0.0};
    std::vector<double> eta = { 1.0 }, eta_inv = { 1.0 };
    std::vector<double> surface_params;

    convert_surface(surf, surface_params, zpos[0], eta[0]);
    // std::cout << "eta = " << eta[0] << std::endl;

    eta_inv[0] = eta_prev[0] / eta[0];

    uint32_t used_raparms = POLY_RPARAMS_MAX_PARAMS;

    pSim->UpdateSurface(surface_params.data(), used_raparms);

    pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), used_raparms, eta_inv.data(),
                    (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                    (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), out_t_valid.data(), out_t_hit.data(), (double*)out_hit_normal.data(),
                    hpd.data(),
                    1, 1, 1, zpos.size());

    eta_prev = eta;
  }

  std::vector<Float3> ref_t_hit_point = {Float3(0.000000000000000E+00, 9.930611122611568E-01, 1.329884193450000E+01)};
  passed = compare_hit_points(ref_t_hit_point, out_hit_pos, experiment_name + "_TEST", algorithm_params.hitEps, false, false);

  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}

bool test_raymarch_matches_interval(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path&)
{
  const std::string test_name = "raymarch_matches_interval";
  std::cout << "TEST: " << test_name << std::endl;

  const uint32_t num_rays_x = 64;
  const uint32_t num_rays_y = 64;
  const uint32_t total_rays = num_rays_x * num_rays_y;
  const float    obj_dist   = 1000.0f;
  const float    field_angle = 5.0f;

  SurfaceParams surf;
  std::vector<Float3> ray_origin, ray_dir;
  generate_rays(num_rays_x, num_rays_y, obj_dist, field_angle, surf.aperture_radius,
                ray_origin, ray_dir);

  std::vector<uint8_t> ray_mask(total_rays, 1);
  std::vector<Float3>  out_hit_pos(total_rays, {0., 0., 0.});
  std::vector<Float3>  out_hit_dir(total_rays, {0., 0., 0.});
  std::vector<Float3>  out_hit_normal(total_rays, {0., 0., 0.});
  std::vector<HitPointData> hpd(total_rays, make_hpd());

  std::vector<double> zpos = {0.0};
  std::vector<double> eta = {1.0};
  std::vector<double> surface_params;
  convert_surface(surf, surface_params, zpos[0], eta[0]);

  IntersectionParameters algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.hitEps = 1e-7;
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;

  IntersectionInitParameters interval_init = lens_sim::MakeIntersectionInitParameters();
  interval_init.maxIterations = 24;
  interval_init.subdivideRaySegment = 10;

  IntersectionInitParameters raymarch_init = interval_init;
  raymarch_init.subdivideRaySegment = 1000000;

  std::vector<Float>   t_interval(total_rays, 0.0f);
  std::vector<uint8_t> valid_interval(total_rays, 0);
  std::vector<Float>   t_raymarch(total_rays, 0.0f);
  std::vector<uint8_t> valid_raymarch(total_rays, 0);

  auto trace = [&](const IntersectionInitParameters& init,
                   std::vector<Float>& t_hit, std::vector<uint8_t>& valid)
  {
    pSim->SetInitIntersectionParams(init);
    pSim->UpdateSurface(surface_params.data(), POLY_RPARAMS_MAX_PARAMS);
    std::fill(out_hit_pos.begin(), out_hit_pos.end(), double3{0., 0., 0.});
    std::fill(out_hit_dir.begin(), out_hit_dir.end(), double3{0., 0., 0.});
    std::fill(out_hit_normal.begin(), out_hit_normal.end(), double3{0., 0., 0.});
    std::fill(t_hit.begin(), t_hit.end(), 0.0f);
    std::fill(valid.begin(), valid.end(), 0);
    pSim->TraceRays(algorithm_params, zpos.data(), surface_params.data(), POLY_RPARAMS_MAX_PARAMS, eta.data(),
                    (const double*)ray_origin.data(), (const double*)ray_dir.data(), ray_mask.data(),
                    (double*)out_hit_pos.data(), (double*)out_hit_dir.data(), valid.data(), t_hit.data(),
                    (double*)out_hit_normal.data(), hpd.data(),
                    num_rays_x, num_rays_y, 1, uint32_t(zpos.size()));
  };

  trace(interval_init, t_interval, valid_interval);
  trace(raymarch_init,  t_raymarch, valid_raymarch);
  const Stats raymarch_stats = pSim->GetStats();

  size_t valid_count = 0;
  size_t mismatched_validity = 0;
  size_t mismatched_t = 0;
  Float max_t_diff = 0.0f;

  for(uint32_t i = 0; i < total_rays; ++i)
  {
    if(valid_interval[i] != valid_raymarch[i])
    {
      mismatched_validity++;
      continue;
    }
    if(!valid_interval[i])
      continue;

    valid_count++;
    const Float diff = std::abs(t_interval[i] - t_raymarch[i]);
    max_t_diff = std::max(max_t_diff, diff);
    if(diff > 1e-6f)
      mismatched_t++;
  }

  const bool passed = (valid_count > 0) &&
                      (mismatched_validity == 0) &&
                      (mismatched_t == 0) &&
                      (raymarch_stats.raymarch_iters > 0);

  std::cout << "  valid rays: " << valid_count << " / " << total_rays << std::endl;
  std::cout << "  validity mismatches: " << mismatched_validity << std::endl;
  std::cout << "  t_hit mismatches > 1e-6: " << mismatched_t
            << " (max diff " << max_t_diff << ")" << std::endl;
  std::cout << "  raymarch iterations: " << raymarch_stats.raymarch_iters << std::endl;
  std::cout << (passed ? "[raymarch_matches_interval]: PASSED!" : "[raymarch_matches_interval]: FAILED!") << std::endl;
  std::cout << "********************************************" << std::endl << std::endl;

  return passed;
}
