#pragma once
#include <vector>
namespace smooth_audio{
    
    namespace detail{
        static constexpr const double SMOOTHING_FACTOR=0.01;  // Smoothing factor for gradual gain change
        static constexpr const size_t WINDOW_SIZE=1024;  // Window size for local RMS computation
        static constexpr const double NOISE_THRESHOLD_FACTOR=0.005;  // Minimum RMS required to apply gain
        static constexpr const double NOISE_THRESHOLD=0.005;
        std::vector<float> apply_k_weighting(std::vector<float>const& samples, int sample_rate);
        double compute_rms_segment(const std::vector<float> &buffer, size_t start, size_t end);
        size_t detect_fade_start(const std::vector<float> &buffer);
        void apply_adaptive_gain(std::vector<float> &buffer);
        double compute_rms(const std::vector<float> &buffer);
        void apply_final_normalization(std::vector<float> &buffer, double original_rms);
    }
    void gain_rider(std::vector<float>& audio_buffer, int num_channels, float target_peak, float gain_smoothing_factor = 0.1f);
    float calculate_lufs_with_gating(const std::vector<float>& samples, int sampleRate, size_t block_size);
    void boost_sample(std::vector<float>& buffer);
}
