#include "even_asphere.h"
#include "isurface.h"
#include "poly_surface.h"
#include "tests.h"
#include "utils.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

static void PrintUsage(const char* exe)
{
  std::cout
    << "Usage: " << exe << " [options]\n"
    << "  --output <dir>          Directory for output files (default: ./output)\n"
    << "  --raymarch              Use ray marching to locate Newton's start point\n"
    << "  --max-iterations <n>    Max Newton iterations (default: 10)\n"
    << "  --precision <eps>       Intersection tolerance (default: 1e-7)\n"
    << "  --help                  Show this message\n";
}

int main(int argc, const char** argv)
{
  std::filesystem::path out_dir = "./output";
  bool use_ray_marching = false;
  uint32_t max_iterations = 10;
  double precision = 1e-7;

  for(int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if(arg == "--help" || arg == "-h")
    {
      PrintUsage(argv[0]);
      return 0;
    }
    else if(arg == "--output" && i + 1 < argc)
    {
      out_dir = argv[++i];
    }
    else if(arg == "--raymarch")
    {
      use_ray_marching = true;
    }
    else if(arg == "--max-iterations" && i + 1 < argc)
    {
      max_iterations = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
    }
    else if(arg == "--precision" && i + 1 < argc)
    {
      precision = std::strtod(argv[++i], nullptr);
    }
    else
    {
      std::cerr << "Unknown or incomplete option: " << arg << std::endl;
      PrintUsage(argv[0]);
      return 2;
    }
  }

  std::filesystem::create_directories(out_dir);

  std::cout << "USING CPU (forward tracing only)" << std::endl;

  // Two independent instances are deliberately created, so that a test can
  // compare the implementation against itself without shared state.
  std::shared_ptr<ISurface> pSim    = CreatePolySurface();
  std::shared_ptr<ISurface> pSimRef = CreatePolySurface();

  IntersectionParameters algorithm_params = lens_sim::MakeIntersectionParameters();
  algorithm_params.hitEps = precision;
  algorithm_params.flags = OUTFLAG_STORE_RAYS | OUTFLAG_STORE_NORMAL | OUTFLAG_STORE_T;

  IntersectionInitParameters init_params = lens_sim::MakeIntersectionInitParameters();
  init_params.maxIterations = max_iterations;
  // 0 invokes the standard method; 1000000 invokes ray marching.
  init_params.subdivideRaySegment = use_ray_marching ? 1000000 : 0;
  init_params.useQtype = false;

  pSim->SetInitIntersectionParams(init_params);
  pSimRef->SetInitIntersectionParams(init_params);
  pSim->SetIntersectionParams(algorithm_params);
  pSimRef->SetIntersectionParams(algorithm_params);
  pSim->CommitDeviceData();
  pSimRef->CommitDeviceData();

  const bool passed = run_all_tests(pSim, pSimRef, out_dir);
  return passed ? 0 : 1;
}
