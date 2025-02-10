#include <vector>
#include <cmath>
#include <sndfile.h>
#include <iostream>
#include "smooth_audio.h"

namespace smooth_audio {
    namespace detail {
        float calculate_rms(const std::vector<float>& samples) {
            float sum = 0.0;
            for (float sample : samples) {
                sum += sample * sample;
            }
            return std::sqrt(sum / samples.size());
        }
        
        std::vector<float> apply_k_weighting(std::vector<float>const& samples, int sample_rate) {
            // High-pass filter constants (from ITU-R BS.1770)
            const float alpha = exp(-2.0 * M_PI * 38.0 / sample_rate);
            constexpr float a_weight = 1.0;
            float prev_sample = 0.0;
            std::vector<float> ret;
            ret.reserve(samples.size());
        
            for (auto&& i: samples) {
                float current_sample = i;
                ret.push_back(a_weight * (current_sample - alpha * prev_sample));
                prev_sample = current_sample;
            }
            return ret;
        }
        
        
        double compute_rms_segment(const std::vector<float> &buffer, size_t start, size_t end) {
            double sum = 0.0;
            for (size_t i = start; i < end; ++i) {
                sum += buffer[i] * buffer[i];
            }
            return sqrt(sum / (end - start));
        }
        void apply_adaptive_gain(std::vector<float> &buffer) {
            double gain = 1.0;  // Persistent gain across segments
        
            for (size_t i = 0; i < buffer.size(); i += WINDOW_SIZE) {
                size_t end = std::min(i + WINDOW_SIZE, buffer.size());
                double rms = compute_rms_segment(buffer, i, end);
        
                // Skip gain adjustment for segments below noise threshold
                if (rms < NOISE_THRESHOLD) {
                    continue;
                }
        
                double target_gain = (rms > 0) ? (TARGET_RMS / rms) : 1.0;
        
                // Gradually transition gain within the segment
                for (size_t j = i; j < end; ++j) {
                    gain += (target_gain - gain) * SMOOTHING_FACTOR;
                    buffer[j] *= gain;
                }
            }
        }
        double compute_rms(const std::vector<float> &buffer) {
            double sum = 0.0;
            for (float sample : buffer) {
                sum += sample * sample;
            }
            return sqrt(sum / buffer.size());
        }
        
        // Apply final gain correction to match original loudness
        void apply_final_normalization(std::vector<float> &buffer, double original_rms) {
            double new_rms = compute_rms(buffer);
            if (new_rms > 0) {
                double final_gain = original_rms / new_rms;
                for (float &sample : buffer) {
                    sample *= final_gain;
                }
            }
        }
    }
    float calculate_lufs_with_gating(const std::vector<float>& samples, int sampleRate, size_t block_size) {
        constexpr float GATING_THRESHOLD_DB = -70.0; 
        const size_t hop_size = block_size / 4;                         // 75% overlap
    
        // Apply K-weighting
        std::vector<float> weighted_samples = samples;
        detail::apply_k_weighting(weighted_samples, sampleRate);
    
        // Calculate short-term loudness for each block
        std::vector<float> block_rms;
        for (size_t i = 0; i + block_size <= weighted_samples.size(); i += hop_size) {
            std::vector<float> block(
                weighted_samples.begin() + i, 
                weighted_samples.begin() + i + block_size
            );
            block_rms.push_back(detail::calculate_rms(block));
        }
    
        // Convert RMS to LUFS and find the loudest block
        std::vector<float> block_lufs(block_rms.size());
        float max_lufs = -INFINITY;
        for (size_t i = 0; i < block_rms.size(); ++i) {
            block_lufs[i] = -0.691 + 10.0 * log10(block_rms[i] * block_rms[i]);
            max_lufs = std::max(max_lufs, block_lufs[i]);
        }
    
        // Apply gating: Include only blocks above the threshold
        const float gating_threshold = max_lufs + GATING_THRESHOLD_DB;
        std::vector<float> gated_blocks;
        for (size_t i = 0; i < block_lufs.size(); ++i) {
            if (block_lufs[i] >= gating_threshold) {
                gated_blocks.push_back(block_rms[i] * block_rms[i]); // Use squared RMS values
            }
        }
    
        if (gated_blocks.empty()) {
            std::cerr << "Warning: All blocks were gated out!\n";
            return -INFINITY; // Return a very low value if no blocks pass gating
        }
    
        // Calculate integrated LUFS from gated blocks
        float mean_square = 0.0;
        for (float block_power : gated_blocks) {
            mean_square += block_power;
        }
        mean_square /= gated_blocks.size();
    
        return -0.691 + 10.0 * log10(mean_square);
    }
    void boost_sample(std::vector<float>& buffer){
        double original_rms = detail::compute_rms(buffer);
        detail::apply_adaptive_gain(buffer);
        detail::apply_final_normalization(buffer, original_rms);
    }
}
