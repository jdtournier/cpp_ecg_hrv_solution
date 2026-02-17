#pragma once

#include <array>
#include <string>
#include <vector>

// A simple structure to hold both the sampling interval and recorded voltages:
struct ECG {
  double sampling_interval;   // in milliseconds
  std::vector<double> V;      // in millivolts
};

// load the ECG data contained in the file specified, and return the
// information in the form of an "ECG" structure:
ECG load_ecg (const std::string& filename);

// Detect R peaks in the ECG trace using the proposed algorithm, given the
// threshold and refractory time specified. This returns a vector of times (in
// milliseconds) for each detected peak, in the order encountered in the ECG trace.
std::vector<double> detect_r_peaks (const ECG& ecg, const double V_threshold, const double t_refract_seconds);

// convert list of R peak times into list of consecutive RR intervals in
// milliseconds:
std::vector<double> compute_RR_intervals (const std::vector<double> r_peaks);

// Compute root mean sum of squared differences (RMSSD) metric:
double compute_RMSSD (const std::vector<double>& RR_intervals);

// Compute minimum, mean & maximum value of the RR interval.
// This returns an array of 3 values, corresponding to the min, mean and max
// value, in that order.
std::array<double,3> compute_min_mean_max (const std::vector<double>& RR_intervals);

// return stress level as a string based on RMSSD:
std::string stress_level (double rmssd);

// Write out a report to the file specified, based on the computed parameters:
void write_report (const std::string& report_filename, const std::string& ecg_filename,
    const std::array<double,3>& min_mean_max, double rmssd);

