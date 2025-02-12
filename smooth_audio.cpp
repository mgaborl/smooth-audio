#include <vector>
#include <cmath>
#include <sndfile.h>
#include <iostream>
#include "smooth_audio.h"

namespace smooth_audio {
    namespace detail {
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
        // Detects where the sample naturally fades out
        size_t detect_fade_start(const std::vector<float> &buffer, float fade_threshold) {
            size_t total_samples = buffer.size();
            size_t window_size = WINDOW_SIZE;
            double prev_rms = compute_rms_segment(buffer, total_samples - window_size, total_samples);
            size_t fade_start = total_samples;
            
            // Scan backward to find where RMS starts decreasing monotonically until below threshold
            for (size_t i = total_samples - window_size; i > 0; i -= window_size) {
                double rms = compute_rms_segment(buffer, i, i + window_size);
                double dynamic_noise_threshold = std::max(rms * 0.1, NOISE_THRESHOLD);
                if (rms > prev_rms) {
                    fade_start = i + window_size;
                    break;
                }
                if (rms < dynamic_noise_threshold) {
                    fade_start = i;
                    break;
                }
                prev_rms = rms;
            }
            return fade_start;
        }
        
        // Apply adaptive gain correction with noise gating and fade preservation
        void apply_adaptive_gain(std::vector<float> &buffer) {
            float target_rms = detail::compute_rms(buffer);
            double persistent_gain = 1.0;
            size_t fade_start = detect_fade_start(buffer, target_rms*0.2);
            
            for (size_t i = 0; i < buffer.size(); i += WINDOW_SIZE) {
                size_t end = std::min(i + WINDOW_SIZE, buffer.size());
                double rms = compute_rms_segment(buffer, i, end);
                
                double dynamic_noise_threshold = std::max(rms * 0.1, NOISE_THRESHOLD);
                double target_gain = (rms > dynamic_noise_threshold) ? (target_rms / rms) : 1.0;
        
                for (size_t j = i; j < end; ++j) {
                    double fade_progress = (double)(j - fade_start) / (buffer.size() - fade_start);
                    // Gradually reduce gain adjustments after fade start
                    double fade_factor = (j > fade_start) ? exp(-5.0 * fade_progress) : 1.0;
                    persistent_gain += (target_gain - persistent_gain) * SMOOTHING_FACTOR * fade_factor;
                    buffer[j] *= persistent_gain;
                }
            }
        }
        void apply_gain(std::vector<float>& buffer, float gain) {
            for (float& sample : buffer) {
                sample *= gain;
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
            double total_rms = 0.0;
            size_t valid_samples = 0;
            
            for (size_t i = 0; i < buffer.size(); i += WINDOW_SIZE) {
                size_t end = std::min(i + WINDOW_SIZE, buffer.size());
                double rms = compute_rms_segment(buffer, i, end);
                double dynamic_noise_threshold = std::max(rms * 0.1, NOISE_THRESHOLD);
                
                if (rms >= dynamic_noise_threshold) {
                    total_rms += rms;
                    valid_samples++;
                }
            }
            
            if (valid_samples > 0) {
                double avg_rms = total_rms / valid_samples;
                double normalization_gain = new_rms / avg_rms;
                
                double rms = compute_rms(buffer);
                double dynamic_noise_threshold = std::max(rms * 0.1, NOISE_THRESHOLD);
                for (size_t i = 0; i < buffer.size(); ++i) {
                    if (compute_rms_segment(buffer, i, i + 1) >= dynamic_noise_threshold)
                        buffer[i] *= normalization_gain;
                }
            }
        }
    } // namespace  detail
    void gain_rider(std::vector<float>& audio_buffer, int num_channels, float target_peak, float gain_smoothing_factor) {
        // Define a smoothing factor for gain changes
        static float smoothed_gain = 1.0f;  // Initial gain (no adjustment)
        
        // Iterate over the entire audio buffer
        for (int i = 0; i < audio_buffer.size() / num_channels; ++i) {
            // Prepare the current buffer for each frame (if stereo or multi-channel)
            std::vector<float> frame_buffer(num_channels);
            for (int ch = 0; ch < num_channels; ++ch) {
                frame_buffer[ch] = audio_buffer[i * num_channels + ch];
            }
    
            // Calculate the peak level of the current signal
            float current_rms = detail::compute_rms(frame_buffer);
    
            if (current_rms < 1e-6f) {
                continue;
            }

            // Calculate the gain required to bring the current peak to the target level
            float target_gain = target_peak / current_rms;  // Avoid division by zero
    
            // Apply smoothing to the gain adjustment
            smoothed_gain += (target_gain - smoothed_gain) * gain_smoothing_factor;
            float gain_change = smoothed_gain - target_gain;
            if (std::abs(gain_change) > 0.2)
                smoothed_gain = target_gain + (gain_change > 0 ? 0.2 : -0.2);
    
            // Apply the smoothed gain to the frame buffer
            detail::apply_gain(frame_buffer, smoothed_gain);
    
            // Update the main audio buffer with the adjusted frame
            for (int ch = 0; ch < num_channels; ++ch) {
                audio_buffer[i * num_channels + ch] = frame_buffer[ch];
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
            block_rms.push_back(detail::compute_rms(block));
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
        std::cout << original_rms << " original RMS\n";
        gain_rider(buffer, 1, original_rms, 0.01);
        detail::apply_final_normalization(buffer, original_rms);
    }
}
