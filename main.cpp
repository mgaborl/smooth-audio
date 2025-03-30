#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
#include <sndfile.hh>
#include <algorithm>

// Constants
size_t WINDOW_SIZE = 1024;  // Window size for RMS calculation
constexpr float MIN_GAIN = 0.1f;      // Minimum gain factor
constexpr float MAX_GAIN = 10.0f;     // Maximum gain factor
constexpr float HIGH_PERCENTILE = 0.90f; // 90th percentile (target loudness)ercentile of RMS values
constexpr float LOW_PERCENTILE = 0.20f; // 90th percentile (target loudness)ercentile of RMS values
float TAIL_DROP_THRESHOLD = 0.3f; // Relative RMS drop defining the start of the tail
constexpr float TAIL_RMS_DROP_FACTOR = 0.5f; // Relative RMS drop defining the start of the tail
constexpr float TAIL_REDUCTION_FACTOR = .9f; // How much to reduce gain in the tail
float ATTACK_DELAY_MS = 0.0f; // Skip the first n milliseconds


// Function to compute RMS of a block of samples
float compute_rms(const std::vector<float>& samples, size_t start, size_t window_size) {
    float sum_squares = 0.0f;
    size_t count = 0;

    for (size_t i = start; i < start + window_size && i < samples.size(); ++i) {
        sum_squares += samples[i] * samples[i];
        ++count;
    }

    return (count > 0) ? std::sqrt(sum_squares / count) : 0.0f;
}

// Function to compute RMS values across the audio
std::vector<float> compute_rms_values(const std::vector<float>& samples) {
    size_t num_windows = (samples.size() + WINDOW_SIZE - 1) / WINDOW_SIZE;
    std::vector<float> rms_values(num_windows);

    for (size_t i = 0; i < num_windows; ++i) {
        size_t start = i * WINDOW_SIZE;
        rms_values[i] = compute_rms(samples, start, WINDOW_SIZE);
    }

    return rms_values;
}
std::pair<float, float> compute_rms_percentiles(const std::vector<float>& rms_values) {
    if (rms_values.empty()) return {0.1f, 0.1f};

    std::vector<float> sorted_rms = rms_values;
    std::sort(sorted_rms.begin(), sorted_rms.end());

    size_t low_index = static_cast<size_t>(LOW_PERCENTILE * sorted_rms.size());
    size_t high_index = static_cast<size_t>(HIGH_PERCENTILE * sorted_rms.size());

    return {sorted_rms[low_index], sorted_rms[high_index]};
}
size_t detect_tail_start(const std::vector<float>& samples) {
    float total_rms = compute_rms(samples, 0, samples.size());

    for (size_t i = 0; i < samples.size(); i += WINDOW_SIZE) {
        float remaining_rms = compute_rms(samples, i, samples.size() - i);
        if (remaining_rms < total_rms * TAIL_RMS_DROP_FACTOR) {
            return i;  // Tail starts here
        }
    }
    return samples.size();
}

// Compute a smoothed gain envelope with tail-end reduction
std::vector<float> compute_smoothed_gain_envelope(
    const std::vector<float>& samples, 
    float min_rms, 
    float target_rms, 
    size_t tail_start,
    size_t sample_rate) 
{
    size_t num_windows = (samples.size() + WINDOW_SIZE - 1) / WINDOW_SIZE;
    std::cout << "Number of windows: " << num_windows << "\n";
    std::vector<float> gain_values(num_windows);
    std::vector<float> gain_envelope(samples.size(), 1.0f);

    // Compute raw gain per window
    for (size_t i = 0; i < num_windows; ++i) {
        size_t start = i * WINDOW_SIZE;
        float rms = compute_rms(samples, start, WINDOW_SIZE);

        float gain = (rms < min_rms) ? 1.0f : std::clamp(target_rms / rms, MIN_GAIN, MAX_GAIN);

        // Reduce gain progressively in the tail
        if (start >= tail_start || start < ATTACK_DELAY_MS*sample_rate/1000) {
            gain = 1; // Avoid full silence
        }

        gain_values[i] = gain;
    }

    // Apply smoothing via linear interpolation
    for (size_t i = 0; i < num_windows - 1; ++i) {
        float g1 = gain_values[i];
        float g2 = gain_values[i + 1];

        for (size_t j = 0; j < WINDOW_SIZE && (i * WINDOW_SIZE + j) < samples.size(); ++j) {
            float t = static_cast<float>(j) / WINDOW_SIZE;  // Interpolation factor
            gain_envelope[i * WINDOW_SIZE + j] = (1.0f - t) * g1 + t * g2;
        }
    }

    // Fill in remaining samples at the end
    for (size_t i = (num_windows - 1) * WINDOW_SIZE; i < samples.size(); ++i) {
        gain_envelope[i] = gain_values[num_windows - 1];
    }

    return gain_envelope;
}
// Function to apply gain envelope to audio samples
void apply_gain(std::vector<float>& samples, const std::vector<float>& gain_envelope) {
    for (size_t i = 0; i < samples.size(); ++i) {
        samples[i] *= gain_envelope[i];
    }
}

// Function to process WAV file
bool process_wav(const std::string& input_file, const std::string& output_file) {
    SF_INFO sfinfo;
    SndfileHandle f(input_file);

    std::vector<float> samples(f.frames() * f.channels());
    f.readf(samples.data(), f.frames());

    // Compute RMS values and determine adaptive target RMS
    std::vector<float> rms_values = compute_rms_values(samples);

    // Compute and apply gain envelope
    auto [low_rms, target_rms] = compute_rms_percentiles(rms_values);
    size_t tail_start = detect_tail_start(samples);
    std::cout << "Skipping boost below RMS: " << low_rms << "\n";
    std::cout << "Tail detected at sample index: " << tail_start << " out of " << samples.size() << "\n";
    std::cout << "Adaptive Target RMS: " << target_rms << "\n";
    std::vector<float> gain_envelope = compute_smoothed_gain_envelope(samples, low_rms,  target_rms, tail_start, f.samplerate());
    std::cout << "Number of points in the envelope: " << gain_envelope.size() << "\n";
    apply_gain(samples, gain_envelope);

    // Write output file
    SndfileHandle outfile(output_file, SFM_WRITE, f.format(), f.channels(), f.samplerate());

    outfile.write(samples.data(), f.frames()*f.channels());

    std::cout << "Processing complete! Output saved to: " << output_file << "\n";
    return true;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_wav> <output_wav> <window_size=1024> <tail_drop_threshold_percentage=30> <attack_delay_ms=0>\n";
        if (argc == 2 && strcmp(argv[1], "--help") * strcmp(argv[1], "-h") == 0){
            std::cout << "Tail drop threshold is defined as (remaining_section_volume)*100/(overall_volume) - the lower it is, the more conservative tail detection gets\n";
            std::cout << "Attack delay skips the first N milliseconds of the sample (applies no boost)\n";
        }
        return 1;
    }
    if (argc >=4 ) 
        WINDOW_SIZE = std::atoi(argv[3]);
    if (argc >=5 ) 
        TAIL_DROP_THRESHOLD = static_cast<float>(std::atoi(argv[4]))/100;
    if (argc >=6 ) 
        ATTACK_DELAY_MS = std::atoi(argv[5]);
    std::cout << "Window size set to: " << WINDOW_SIZE << "\n";

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    if (!process_wav(input_file, output_file)) {
        return 1;
    }


    return 0;
}
