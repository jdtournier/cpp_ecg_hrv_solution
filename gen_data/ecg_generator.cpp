// ============================================================
// ECG Synthetic Data Generator
// Generates realistic ECG CSV files for 3 stress levels:
//   - Low Stress    (high HRV, RMSSD > 50 ms)
//   - Moderate Stress (medium HRV, RMSSD 20-50 ms)
//   - High Stress   (low HRV, RMSSD < 20 ms)
//
// Usage:
//   g++ -std=c++17 -Wall -o ecg_gen ecg_generator.cpp -lm
//   ./ecg_gen [duration_sec] [sampling_freq_Hz]
//
// Defaults: 10 seconds at 500 Hz
// Output:   ecg_low_stress.txt
//           ecg_moderate_stress.txt
//           ecg_high_stress.txt
//
// File format:
//   Line 1:  "ECG"                  <- magic header (identifies valid ECG file)
//   Line 2:  sampling interval (ms) <- e.g. "2.0" for 500 Hz
//   Line 3+: voltage values (mV)    <- one reading per line
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iomanip>
#include <cstdlib>

// ── Structs ──────────────────────────────────────────────

// Parameters that define a single ECG waveform component
// modeled as a Gaussian bump: amplitude * exp(-(t-center)^2 / (2*width^2))
struct WaveComponent {
    double amplitude;   // Peak height (mV), negative for Q and S waves
    double center;      // Time offset from R-peak center (seconds)
    double width;       // Gaussian width (seconds), controls sharpness
};

// Full set of parameters for one stress profile
struct ECGProfile {
    std::string name;           // e.g. "low_stress"
    std::string label;          // e.g. "Low Stress"
    double base_hr_bpm;         // Mean heart rate (beats per minute)
    double rr_variability_ms;   // Std deviation of R-R intervals (ms)
                                //   High value -> high HRV -> low stress
                                //   Low value  -> low HRV  -> high stress
    double baseline_wander_amp; // Amplitude of slow baseline drift (mV)
    double noise_std;           // Gaussian noise standard deviation (mV)
    double r_peak_amplitude;    // R-peak height (mV)

    // ECG morphology components (P, Q, R, S, T waves)
    WaveComponent p_wave;
    WaveComponent q_wave;
    WaveComponent r_wave;
    WaveComponent s_wave;
    WaveComponent t_wave;
};

// ── Function Declarations ────────────────────────────────

// Creates the 3 stress-level profiles with clinically motivated parameters
std::vector<ECGProfile> createProfiles();

// Generates R-peak times with realistic variability for a given profile
std::vector<double> generateRPeakTimes(const ECGProfile& profile,
                                        double duration_sec,
                                        std::mt19937& rng);

// Computes the ECG voltage at a given time based on nearby R-peaks
double computeECGVoltage(double t,
                         const std::vector<double>& peak_times,
                         const ECGProfile& profile,
                         std::mt19937& rng);

// Generates and writes one complete ECG file
void generateECGFile(const ECGProfile& profile,
                     double duration_sec,
                     double sampling_freq,
                     unsigned int seed);

// Prints a summary of the generated R-R intervals for verification
void printProfileSummary(const std::string& label,
                         const std::vector<double>& peak_times);


// ── Main ─────────────────────────────────────────────────

int main(int argc, char* argv[]) {

    double duration_sec   = 10.0;   // Default: 10 seconds
    double sampling_freq  = 500.0;  // Default: 500 Hz

    if (argc >= 2) {
        duration_sec = std::atof(argv[1]);
        if (duration_sec < 2.0 || duration_sec > 300.0) {
            std::cerr << "Error: Duration must be between 2 and 300 seconds." << std::endl;
            return 1;
        }
    }
    if (argc >= 3) {
        sampling_freq = std::atof(argv[2]);
        if (sampling_freq < 100.0 || sampling_freq > 2000.0) {
            std::cerr << "Error: Sampling frequency must be between 100 and 2000 Hz." << std::endl;
            return 1;
        }
    }

    std::cout << "ECG Synthetic Data Generator" << std::endl;
    std::cout << "============================" << std::endl;
    std::cout << "  Duration       : " << duration_sec << " seconds" << std::endl;
    std::cout << "  Sampling freq  : " << sampling_freq << " Hz" << std::endl;
    std::cout << "  Samples/file   : " << static_cast<int>(duration_sec * sampling_freq) << std::endl;
    std::cout << std::endl;

    std::vector<ECGProfile> profiles = createProfiles();

    unsigned int base_seed = 42;  // Fixed seed for reproducibility

    for (size_t i = 0; i < profiles.size(); i++) {
        std::cout << "Generating: " << profiles[i].label << " ..." << std::endl;
        generateECGFile(profiles[i], duration_sec, sampling_freq, base_seed + i * 100);
    }

    std::cout << "\nDone! Generated 3 ECG files." << std::endl;
    return 0;
}


// ── Profile Definitions ──────────────────────────────────

std::vector<ECGProfile> createProfiles() {
    std::vector<ECGProfile> profiles;

    // ────────────────────────────────────────────────────
    // LOW STRESS (Relaxed / Athlete at rest)
    //   - Lower resting HR (~62 bpm)
    //   - HIGH R-R variability (std ~80 ms) -> RMSSD > 50
    //   - Clean signal, slight baseline wander
    // ────────────────────────────────────────────────────
    ECGProfile low;
    low.name                = "low_stress";
    low.label               = "Low Stress (High HRV)";
    low.base_hr_bpm         = 62.0;
    low.rr_variability_ms   = 80.0;     // High variability
    low.baseline_wander_amp = 0.04;
    low.noise_std           = 0.02;
    low.r_peak_amplitude    = 1.2;
    low.p_wave  = {  0.15, -0.18,  0.035 };
    low.q_wave  = { -0.12, -0.025, 0.008 };
    low.r_wave  = {  1.20,  0.0,   0.005 };  // Sharp, tall R-peak
    low.s_wave  = { -0.18,  0.030, 0.008 };
    low.t_wave  = {  0.25,  0.20,  0.045 };  // Broad, normal T-wave
    profiles.push_back(low);

    // ────────────────────────────────────────────────────
    // MODERATE STRESS (Normal daily activity / mild stress)
    //   - Normal resting HR (~75 bpm)
    //   - MODERATE R-R variability (std ~35 ms) -> RMSSD 20-50
    //   - Some noise
    // ────────────────────────────────────────────────────
    ECGProfile moderate;
    moderate.name                = "moderate_stress";
    moderate.label               = "Moderate Stress (Medium HRV)";
    moderate.base_hr_bpm         = 75.0;
    moderate.rr_variability_ms   = 25.0;   // Moderate variability
    moderate.baseline_wander_amp = 0.05;
    moderate.noise_std           = 0.03;
    moderate.r_peak_amplitude    = 1.1;
    moderate.p_wave  = {  0.13, -0.17,  0.032 };
    moderate.q_wave  = { -0.10, -0.022, 0.007 };
    moderate.r_wave  = {  1.10,  0.0,   0.005 };
    moderate.s_wave  = { -0.20,  0.028, 0.008 };
    moderate.t_wave  = {  0.22,  0.19,  0.040 };
    profiles.push_back(moderate);

    // ────────────────────────────────────────────────────
    // HIGH STRESS (Anxiety / sympathetic dominance)
    //   - Elevated HR (~95 bpm)
    //   - LOW R-R variability (std ~10 ms) -> RMSSD < 20
    //   - More noise, slight morphology changes
    // ────────────────────────────────────────────────────
    ECGProfile high;
    high.name                = "high_stress";
    high.label               = "High Stress (Low HRV)";
    high.base_hr_bpm         = 95.0;
    high.rr_variability_ms   = 10.0;     // Very low variability
    high.baseline_wander_amp = 0.07;
    high.noise_std           = 0.04;
    high.r_peak_amplitude    = 1.0;
    high.p_wave  = {  0.10, -0.15,  0.028 };  // Slightly smaller P
    high.q_wave  = { -0.08, -0.020, 0.007 };
    high.r_wave  = {  1.00,  0.0,   0.004 };  // Slightly sharper
    high.s_wave  = { -0.22,  0.025, 0.007 };  // Slightly deeper S
    high.t_wave  = {  0.18,  0.16,  0.035 };  // Flatter T-wave
    profiles.push_back(high);

    return profiles;
}


// ── R-Peak Time Generation ───────────────────────────────

std::vector<double> generateRPeakTimes(const ECGProfile& profile,
                                        double duration_sec,
                                        std::mt19937& rng) {
    std::vector<double> peak_times;

    double mean_rr_sec = 60.0 / profile.base_hr_bpm;
    double rr_std_sec  = profile.rr_variability_ms / 1000.0;

    // Normal distribution for R-R intervals
    std::normal_distribution<double> rr_dist(mean_rr_sec, rr_std_sec);

    // Start with a small offset so the first beat isn't at t=0
    double current_time = 0.2 + std::abs(rr_dist(rng)) * 0.3;
    double min_rr = 0.33;   // Minimum R-R: ~180 bpm physiological limit
    double max_rr = 2.0;    // Maximum R-R: ~30 bpm

    while (current_time < duration_sec - 0.3) {
        peak_times.push_back(current_time);

        // Generate next R-R interval with variability
        double rr = rr_dist(rng);

        // Clamp to physiological limits
        if (rr < min_rr) rr = min_rr;
        if (rr > max_rr) rr = max_rr;

        current_time += rr;
    }

    return peak_times;
}


// ── ECG Voltage Computation ──────────────────────────────

double computeECGVoltage(double t,
                         const std::vector<double>& peak_times,
                         const ECGProfile& profile,
                         std::mt19937& rng) {
    double voltage = 0.0;

    // Baseline wander (slow sinusoidal drift)
    voltage += profile.baseline_wander_amp * std::sin(2.0 * M_PI * 0.15 * t);
    voltage += profile.baseline_wander_amp * 0.5 * std::sin(2.0 * M_PI * 0.08 * t + 0.7);

    // Gaussian noise
    std::normal_distribution<double> noise(0.0, profile.noise_std);
    voltage += noise(rng);

    // Add contribution from each heartbeat
    for (const double& pt : peak_times) {
        double dt = t - pt;

        // Only compute if we're close enough to this beat (optimization)
        if (std::abs(dt) > 0.5) continue;

        // P-wave
        double p_dt = dt - profile.p_wave.center;
        voltage += profile.p_wave.amplitude *
                   std::exp(-(p_dt * p_dt) / (2.0 * profile.p_wave.width * profile.p_wave.width));

        // Q-wave
        double q_dt = dt - profile.q_wave.center;
        voltage += profile.q_wave.amplitude *
                   std::exp(-(q_dt * q_dt) / (2.0 * profile.q_wave.width * profile.q_wave.width));

        // R-wave (the main peak)
        double r_dt = dt - profile.r_wave.center;
        voltage += profile.r_wave.amplitude *
                   std::exp(-(r_dt * r_dt) / (2.0 * profile.r_wave.width * profile.r_wave.width));

        // S-wave
        double s_dt = dt - profile.s_wave.center;
        voltage += profile.s_wave.amplitude *
                   std::exp(-(s_dt * s_dt) / (2.0 * profile.s_wave.width * profile.s_wave.width));

        // T-wave
        double t_dt = dt - profile.t_wave.center;
        voltage += profile.t_wave.amplitude *
                   std::exp(-(t_dt * t_dt) / (2.0 * profile.t_wave.width * profile.t_wave.width));
    }

    return voltage;
}


// ── File Generation ──────────────────────────────────────

void generateECGFile(const ECGProfile& profile,
                     double duration_sec,
                     double sampling_freq,
                     unsigned int seed) {

    std::mt19937 rng(seed);

    // Step 1: Generate R-peak times with HRV
    std::vector<double> peak_times = generateRPeakTimes(profile, duration_sec, rng);

    // Step 2: Generate the ECG signal
    int num_samples = static_cast<int>(duration_sec * sampling_freq);
    double sampling_interval_ms = 1000.0 / sampling_freq;  // e.g. 2.0 ms at 500 Hz

    std::string filename = "ecg_" + profile.name + ".txt";

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not create file '" << filename << "'" << std::endl;
        return;
    }

    // Line 1: Magic header to identify valid ECG files
    outfile << "ECG" << std::endl;

    // Line 2: Sampling time interval in milliseconds
    outfile << std::fixed << std::setprecision(1) << sampling_interval_ms << std::endl;

    // Lines 3+: Voltage recordings in mV, one per line
    for (int i = 0; i < num_samples; i++) {
        double t = static_cast<double>(i) / sampling_freq;
        double v = computeECGVoltage(t, peak_times, profile, rng);

        outfile << std::fixed << std::setprecision(4) << v << std::endl;
    }

    outfile.close();

    // Print summary
    std::cout << "  -> " << filename << " (" << num_samples << " samples, "
              << peak_times.size() << " heartbeats, dt="
              << sampling_interval_ms << " ms)" << std::endl;

    printProfileSummary(profile.label, peak_times);
}


// ── Profile Summary ──────────────────────────────────────

void printProfileSummary(const std::string& label,
                         const std::vector<double>& peak_times) {

    if (peak_times.size() < 2) {
        std::cout << "     (Too few peaks to compute HRV)" << std::endl;
        return;
    }

    // Compute R-R intervals
    std::vector<double> rr_ms;
    for (size_t i = 1; i < peak_times.size(); i++) {
        rr_ms.push_back((peak_times[i] - peak_times[i - 1]) * 1000.0);
    }

    // Mean RR
    double sum_rr = 0.0;
    for (double rr : rr_ms) sum_rr += rr;
    double mean_rr = sum_rr / rr_ms.size();

    // SDNN
    double sum_sq = 0.0;
    for (double rr : rr_ms) {
        double d = rr - mean_rr;
        sum_sq += d * d;
    }
    double sdnn = std::sqrt(sum_sq / rr_ms.size());

    // RMSSD
    double sum_succ = 0.0;
    for (size_t i = 1; i < rr_ms.size(); i++) {
        double d = rr_ms[i] - rr_ms[i - 1];
        sum_succ += d * d;
    }
    double rmssd = std::sqrt(sum_succ / (rr_ms.size() - 1));

    // Mean HR
    double sum_bpm = 0.0;
    double min_bpm = 999.0, max_bpm = 0.0;
    for (double rr : rr_ms) {
        double bpm = 60000.0 / rr;
        sum_bpm += bpm;
        if (bpm < min_bpm) min_bpm = bpm;
        if (bpm > max_bpm) max_bpm = bpm;
    }
    double mean_bpm = sum_bpm / rr_ms.size();

    // Stress classification
    std::string stress;
    if (rmssd > 50.0)      stress = "Low Stress";
    else if (rmssd >= 20.0) stress = "Moderate Stress";
    else                     stress = "High Stress";

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "     [" << label << "]" << std::endl;
    std::cout << "     Mean RR: " << mean_rr << " ms | "
              << "SDNN: " << sdnn << " ms | "
              << "RMSSD: " << rmssd << " ms" << std::endl;
    std::cout << "     HR: " << mean_bpm << " bpm ("
              << min_bpm << " - " << max_bpm << ") | "
              << "Classification: " << stress << std::endl;
    std::cout << std::endl;
}
