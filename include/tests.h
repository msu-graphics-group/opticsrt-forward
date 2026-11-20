#pragma once
#include "LiteMath.h"

#include <vector>
#include <filesystem>
#include <memory>
#include <iostream>

#include "isurface.h"
#include "utils.h"

struct Timings
{
  LiteMath::float4 traceTimings;
};

struct RunArguments
{
  lens_sim::IntersectionParameters params;
  const double* zpos;
  const double* surface_params;
  const uint32_t used_rparams;
  const double* eta_inv;
  const double* ray_pos;
  const double* ray_dir;
  const uint8_t* ray_mask;
  double* hit_pos;
  double* hit_dir;
  uint8_t* hit_mask;
  double* hit_t;
  double* hit_normal;
  HitPointData* hpd;
  uint32_t n_rays_x;
  uint32_t n_rays_y;
  uint32_t n_batches;
  uint32_t z_pos_size;
};

static Timings measure_performance(const std::shared_ptr<ISurface>& pSim, RunArguments args, bool print = false)
{
  constexpr int WARMUP_RUNS = 3;
  constexpr int RUNS_NUM    = 3;

  pSim->SetIntersectionParams(args.params);

  for(int i = 0; i < WARMUP_RUNS; ++i)
  {
    pSim->TraceRays(args.params,
                    args.zpos, args.surface_params,
                    args.used_rparams,
                    args.eta_inv,
                    args.ray_pos,  args.ray_dir,  args.ray_mask,
                    args.hit_pos,  args.hit_dir,  args.hit_mask, args.hit_t, args.hit_normal, args.hpd,
                    args.n_rays_x, args.n_rays_y, args.n_batches, args.z_pos_size);
  }

  Timings res {{0.0f, 0.0f, 0.0f, 0.0f}};

  for(int i = 0; i < RUNS_NUM; ++i)
  {
    pSim->TraceRays(args.params,
                    args.zpos, args.surface_params,
                    args.used_rparams,
                    args.eta_inv,
                    args.ray_pos,  args.ray_dir,  args.ray_mask,
                    args.hit_pos,  args.hit_dir,  args.hit_mask, args.hit_t, args.hit_normal, args.hpd,
                    args.n_rays_x, args.n_rays_y, args.n_batches, args.z_pos_size);

    float timings[4] = {0,0,0,0};
    pSim->GetExecutionTime("TraceRays", timings);
    res.traceTimings.x += timings[0];
    res.traceTimings.y += timings[1];
    res.traceTimings.z += timings[2];
    res.traceTimings.w += timings[3];
  }

  res.traceTimings /= float(RUNS_NUM);

  if(print)
  {
    std::cout << "per surface average time (whole cmd buffer) : " << std::endl;
    std::cout << "TraceRays(exec)  = " << res.traceTimings.x << " ms " << std::endl;
    std::cout << "TraceRays(copy)  = " << res.traceTimings.y + res.traceTimings.z << " ms " << std::endl;
    std::cout << "TraceRays(ovrh)  = " << res.traceTimings.w << " ms " << std::endl;
  }

  return res;
}

void output_performance_stats(const std::filesystem::path& out_dir, Timings& timings, uint32_t num_surfaces);

bool test_default_surface(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);
bool test_debug_surface(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);
bool test_double_gauss_singleray(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);
bool test_nonsquare_ray_batch(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);
bool test_two_batches(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path& out_dir);
bool test_batches_sequential(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>&, const std::filesystem::path& out_dir);
bool test_intersect_behind(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);
bool test_raymarch_matches_interval(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);

bool run_all_tests(const std::shared_ptr<ISurface>& pSim, const std::shared_ptr<ISurface>& pSimRef, const std::filesystem::path& out_dir);
