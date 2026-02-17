#include <iostream>
#include <stdexcept>

#include "debug.h"
#include "ecg.h"

int main (int argc, char* argv[])
{
  debug::set_verbose_mode (argc, argv);

  if (argc < 3) {
    std::cerr
      << "ERROR: expected at least 2 arguments: input file and output file\n\n"
      << "  usage: ecg_hrv infile outfile [ V_threshold [ refractory_time ] ]\n\n";
    return 1;
  }

  double V_threshold = 0.5;
  double t_refract = 0.3;

  // use optional parameters if provided on command-line:
  if (argc > 3) {
    V_threshold = std::stod (argv[3]);
    if (argc > 4)
      t_refract = std::stod (argv[4]);
  }

  try {
    auto ecg = load_ecg (argv[1]);

    const auto peaks = detect_r_peaks (ecg, V_threshold, t_refract);
    const auto RR_intervals = compute_RR_intervals (peaks);
    const auto rmssd = compute_RMSSD (RR_intervals);
    const auto min_mean_max = compute_min_mean_max (RR_intervals);

    write_report (argv[2], argv[1], min_mean_max, rmssd);

  }
  catch (std::exception& err) {
    std::cerr << "ERROR: " << err.what() << " - aborting\n";
    return 1;
  }
  catch (...) {
    std::cerr << "ERROR: exception of unknown type was thrown - aborting\n";
    return 1;

  }

  return 0;
}
