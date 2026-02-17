#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cmath>

#include <termviz.h>

#include "debug.h"
#include "ecg.h"


ECG load_ecg (const std::string& filename)
{
  std::ifstream infile (filename);
  if (!infile)
    throw std::runtime_error ("failed to open input file \"" + filename + "\"");

  // First entry should be the marker "ECG" as verify this is an ECG file:
  std::string check;
  infile >> check;
  if (check != "ECG")
    throw std::runtime_error ("file \"" + filename + "\" does not contain ECG data in expected format");

  // Next entry should be the sampling interval in milliseconds:
  double sampling_interval;
  infile >> sampling_interval;

  // subsequent entries correspond to the recorded voltage in millivolts:
  std::vector<double> voltage;
  double val;
  while (infile >> val)
    voltage.push_back (val);

  if (voltage.empty())
    throw std::runtime_error ("no usable data found in input file \"" + filename + "\"");


  // If running in debug/verbose mode, log relevant information, and plot the ECG trace for inspection:
  if (debug::verbose_mode) {
    debug::log (std::format ("ECG traced loaded, sampling rate {} ms, {} samples, total duration {} seconds\n",
          sampling_interval, voltage.size(), voltage.size() * sampling_interval / 1000.0));
    std::vector<double> t (voltage.size());
    for (int n = 0; n < t.size(); n++)
      t[n] = n*sampling_interval / 1000.0;
    termviz::figure(1024, 256).plot (t, voltage);
  }

  return { sampling_interval, voltage };
}



std::vector<double> detect_r_peaks (const ECG& ecg, const double V_threshold, const double t_refract_seconds)
{
  std::vector<double> peaks;
  const double t_refract = t_refract_seconds * 1000.0;

  for (int n = 1; n < ecg.V.size()-1; n++) {
    if (ecg.V[n] > V_threshold && ecg.V[n] > ecg.V[n-1] && ecg.V[n] > ecg.V[n+1]) {
      double t_peak = n * ecg.sampling_interval;
      if (peaks.empty() || t_peak > peaks.back()+t_refract)
        peaks.push_back (t_peak);
    }
  }

  // If running in debug/verbose mode, log relevant information for inspection:
  if (debug::verbose_mode) {
    std::string times;
    for (const auto x : peaks)
      times += std::to_string (x) + " ";
    debug::log (std::format ("detected {} R peaks at times: {}", peaks.size(), times));
  }

  return peaks;
}




std::vector<double> compute_RR_intervals (const std::vector<double> r_peaks)
{
  std::vector<double> intervals;
  for (int n = 1; n < r_peaks.size(); n++)
    intervals.push_back (r_peaks[n] - r_peaks[n-1]);

  // If running in debug/verbose mode, log relevant information for inspection:
  if (debug::verbose_mode) {
    std::string mesg = "RR intervals: ";
    for (const auto x : intervals)
      mesg += std::to_string(x) + " ";
    debug::log (mesg);
  }

  return intervals;
}



double compute_RMSSD (const std::vector<double>& RR_intervals)
{
  double sum = 0.0;
  for (int n = 1; n < RR_intervals.size(); n++) {
    double diff = RR_intervals[n] - RR_intervals[n-1];
    sum += diff*diff;
  }

  double rmssd = std::sqrt (sum / (RR_intervals.size() - 1.0));

  // If running in debug/verbose mode, log relevant information for inspection:
  debug::log ("computed RMSSD: " + std::to_string (rmssd));

  return rmssd;
}




std::array<double,3> compute_min_mean_max (const std::vector<double>& RR_intervals)
{
  double sum = 0.0;
  double min = RR_intervals[0];
  double max = RR_intervals[0];

  for (const auto t : RR_intervals) {
    min = std::min (min, t);
    max = std::max (max, t);
    sum += t;
  }
  double mean = sum / RR_intervals.size();

  // If running in debug/verbose mode, log relevant information for inspection:
  debug::log (std::format ("computed min/mean/max: {} {} {}", min, mean, max));

  return { min, mean, max };
}




std::string stress_level (double rmssd)
{
  if (rmssd > 50.0) return "low";
  if (rmssd >= 20.0) return "moderate";
  return "high";
}




void write_report (const std::string& report_filename, const std::string& ecg_filename,
    const std::array<double,3>& min_mean_max, double rmssd)
{
  std::ofstream outfile (report_filename);
  if (!outfile)
    throw std::runtime_error ("failed to open output file \"" + report_filename + "\"");

  outfile
    << "data file: \"" << ecg_filename << "\"\n"
    << "min/mean/max RR interval: " << min_mean_max[0] << " " << min_mean_max[1] << " " << min_mean_max [2] << "\n"
    << "RMSSD: " << rmssd << "\n"
    << "stress level: " << stress_level (rmssd) << "\n";
}


