#include <sndfile.hh>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

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

// Calculate LUFS
float calculate_lufs(std::vector<float> const& samples, int sample_rate) {
    // Apply K-weighting
    auto weighted_samples = apply_k_weighting(samples, sample_rate);

    // Calculate mean squared value
    float sum_of_squares = 0.0;
    for (auto sample : weighted_samples) 
        sum_of_squares += sample * sample;
    
    float mean_square = sum_of_squares / weighted_samples.size();

    // Convert to LUFS (logarithmic scale)
    return -0.691 + 10.0 * log10(mean_square);
}
// Calculate LUFS with gating
float calculate_lufs_with_gating(const std::vector<float>& samples, int sampleRate, size_t block_size) {
    //std::cout << "Block size: " << block_size << "\n";
    constexpr float GATING_THRESHOLD_DB = -70.0; 
    // Parameters for short-term blocks (400ms with 75% overlap)
    // const size_t block_size = static_cast<size_t>(0.4 * sampleRate); // 400ms
    const size_t hop_size = block_size / 4;                         // 75% overlap

    // Apply K-weighting
    std::vector<float> weighted_samples = samples;
    apply_k_weighting(weighted_samples, sampleRate);

    // Calculate short-term loudness for each block
    std::vector<float> block_rms;
    for (size_t i = 0; i + block_size <= weighted_samples.size(); i += hop_size) {
        std::vector<float> block(
            weighted_samples.begin() + i, 
            weighted_samples.begin() + i + block_size
        );
        block_rms.push_back(calculate_rms(block));
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

std::vector<float> get_segment_volumes(std::vector<float> const& samples, int sample_rate, int segment_size){
    auto begin = std::begin(samples);
    auto end = std::end(samples);
    std::vector<float> ret;
    ret.reserve(samples.size() / segment_size);
    while(end - begin > segment_size){
        std::vector<float> tmp(begin, begin + segment_size);
        ret.push_back(calculate_lufs_with_gating(tmp, sample_rate, segment_size));
        // std::cout << (calculate_lufs_with_gating(tmp, sample_rate, segment_size)) << "\n";
        // std::cout << segment_size << "\n";

        begin += segment_size;
    }
    return ret;
}
// linear for now
struct Spline{
    std::vector<float> y_values;
    std::vector<float> x_values;
    float operator()(float dt) const{
        if (dt > 1)
            throw std::runtime_error("Spline dt is too large. Correct range: [0, 1)");
        if (y_values.empty())
            return 0.0;
        if (y_values.size() == 1)
            return y_values[0];
        float x_length = x_values.back() - x_values.front();
        float x = x_values.front()+ dt * x_length;
        auto x_after = std::upper_bound(
            std::begin(x_values),
            std::end(x_values),
            x
        );
        auto x_before = x_after - 1;
        auto y_idx = x_before - x_values.begin();
        float y_before = y_values[y_idx];
        float y_after = y_values[y_idx + 1];

        auto local_delta = (x - *x_before)/x_length;
        return (1-local_delta) * y_before + local_delta * y_after;
    };
    void add_point(float x, float y){
        // std::cout << "x: " << x << " y: " << y<< "\n";
        if (!x_values.empty() && x < x_values.back())
            throw std::runtime_error("x values are out of order");
        x_values.push_back(x);
        y_values.push_back(y);
    }
    void set_point(int idx, float x, float y){
        if(idx >= x_values.size())
            throw std::runtime_error("Index out of bounds");
        x_values[idx] = x;
        y_values[idx] = y;
    }
    std::pair<float, float> get_point(int idx) const{
        if(idx >= x_values.size())
            throw std::runtime_error("Index out of bounds");
        return {x_values.at(idx), y_values.at(idx)};
    }
    size_t size() const {
        return x_values.size();
    }
    std::string to_string() const {
        std::stringstream ret;        
        for(int i = 0; i < x_values.size(); ++i)
            ret << '(' << x_values[i] << ", " << y_values[i]<< "), ";
        return ret.str();
    }
};
bool vol_in_range(float vol, float avg, float up_d_max, float low_d_max){
        auto delta = vol - avg;
        // std::cout << vol << " " << avg << " " << up_d_max << " " << low_d_max << " " << delta << "\n";
        return delta < up_d_max && delta > low_d_max;
}
Spline interp_segments_in_range(std::vector<float> const& segments, float avg_volume, float max_pos_delta, float max_neg_delta){
    Spline ret;
    for (int i = 0; i < segments.size(); ++i){
        if (vol_in_range(segments[i], avg_volume, max_pos_delta, max_neg_delta))
            ret.add_point(i, segments[i]);
        else if (i == 0) // force the start to be at the start of the recording
            ret.add_point(0, -70);
        else if (i == segments.size() - 1) // force a last point
            ret.add_point(segments.size() - 1, -70);
    }
    return ret;
};

Spline interp_peaks(std::vector<float> const& segments, float avg_volume){
    Spline ret;
    ret.add_point(0, -70);

    // We store peaks (first derivative changes its sign from + to -)
    for (int i = 1; i < segments.size()-1; ++i){
        // prev d
        auto prev = segments[i] - segments[i-1];
        auto next = segments[i+1] - segments[i];
        if(prev > 0 && next < 0 && segments[i] > avg_volume)
            ret.add_point(i, segments[i]);
    }
    // We do not have a segment of complete silence and nothing else
    if(segments.size() > 2){
        auto [first_x, first_y] = ret.get_point(1);
        auto [last_x, last_y] = ret.get_point(ret.size() - 2);
        ret.set_point(0, first_x, first_y);
        ret.add_point(segments.size() - 1, last_y);
    }
    return ret;
}

Spline abs_peak(std::vector<float> const& segments, float threshold){
    Spline ret;
    auto max_y = *std::max(std::begin(segments), std::end(segments));
    auto min_y = *std::min(std::begin(segments), std::end(segments));
    auto old_dy_range = fabs(max_y - min_y);

    std::vector<float> tmp(std::begin(segments), std::end(segments));
    std::sort(std::begin(tmp), std::end(tmp), std::greater<float>());
    float tenth_percentile = std::accumulate(
        std::begin(tmp), 
        std::begin(tmp) + tmp.size()/10, 
        0
    ) / (tmp.size() / 10.0);
    auto new_dy_range = fabs(max_y - tenth_percentile);
    auto scaling_factor = new_dy_range/old_dy_range;
    // std::cout << tenth_percentile << "\n";

    // std::cout << "Max volume: " << max_y << " threshold: " << threshold << "\n";
    for(int i = 0; i < segments.size(); ++i)
        if(segments[i] < threshold)
            ret.add_point(i, segments[i]);
        else {
            float orig_dist = fabs(segments[i] - min_y)/old_dy_range;
            float new_volume = orig_dist * scaling_factor + tenth_percentile;
            std::cout << "Old volume: " << segments[i] << " New volume: " << new_volume << "\n";
            ret.add_point(i, new_volume);
        }
    return ret;
}

// TODO: analytical solution
std::vector<float> numeric_approx(std::vector<float> sample, float target_lufs, int sample_rate, int segment_size){
    auto current_lufs = calculate_lufs_with_gating(sample, sample_rate, segment_size);
    float tmp = 1.0;
    while(fabs(target_lufs - current_lufs) > 0.01){
        float mul = target_lufs > current_lufs ? 1.0001 : 0.9999;
        auto skip = pow(mul, fabs(target_lufs - current_lufs)*(1/0.0009));
        for(auto& i : sample)
            i *= skip;
        current_lufs = calculate_lufs_with_gating(sample, sample_rate, segment_size);
        // std::cout << current_lufs<< " " << target_lufs << "\n";
        tmp *= mul;
    }
    return sample;
}

std::vector<float> boost(
    std::vector<float> const& samples, 
    std::vector<bool> const& mask, 
    Spline const& spline, 
    int sample_rate,
    int segment_size
){

    //std::cout << "Segment size (BOOST): " << segment_size << "\n";
    std::vector<float> ret(std::begin(samples), std::end(samples));
    auto segment_count = mask.size();
    for(int i = 0; i < segment_count; ++i){
        if (!mask[i])
            continue;
        // //std::cout << calculate_lufs_with_gating(samples, sample_rate) << " gated lufs\n";
        // //std::cout << calculate_lufs_with_gating(ret, sample_rate) << " ret gated lufs\n";
        size_t sample_offset = i * segment_size;
        auto sample_start = std::begin(samples) + sample_offset;
        auto sample_end = (std::end(samples) - sample_start > segment_size) ? sample_start + segment_size : std::end(samples);
        auto current_lufs = calculate_lufs_with_gating(std::vector<float>(sample_start, sample_end), sample_rate, segment_size);
        auto target_lufs = spline(static_cast<float>(i)/(segment_count - 1));
        auto last_idx = std::min(sample_offset + segment_size, samples.size()-1);
        // std::cout << "Segment: " << i << "\nTarget LUFS: " << target_lufs << "\nCurrent LUFS: " << current_lufs << "\n";
        std::vector<float> ret_segment(
                std::begin(ret) + sample_offset, 
                std::begin(ret) + last_idx
        );
        auto res = numeric_approx(ret_segment, target_lufs, sample_rate, segment_size);
        for(int j = sample_offset; j < last_idx; ++j)
            ret[j] = res[j - sample_offset]; //std::pow(10, current_lufs - target_lufs);
        /* std::cout << "New LUFS: " << calculate_lufs_with_gating(
            std::vector<float>(
                std::begin(ret) + sample_offset, 
                std::begin(ret) + last_idx
            ), 
            sample_rate, 
            segment_size
        ) << "\n";*/
    }
    return ret;
}

std::vector<float> smooth_clip(std::vector<float> const& samples, int sample_rate, int segment_size){
    //std::cout << "Segment size (SMOOTH): " << segment_size << "\n";
    constexpr float max_upper_delta = 20;
    constexpr float max_lower_delta = -2.0;

    auto y_values = get_segment_volumes(samples, sample_rate, segment_size);
    auto avg_volume = calculate_lufs_with_gating(samples, sample_rate, segment_size);
    float relative_threshold = avg_volume - 10;
    // auto spline = interp_segments_in_range(
    //     y_values, 
    //     avg_volume, 
    //     max_upper_delta,
    //     max_lower_delta
    // );
    // auto spline = interp_peaks(y_values, avg_volume);
    auto spline = abs_peak(y_values, relative_threshold);
    // std::cout << spline.to_string() << "\n";
    // for(int i = 0; i < spline.y_values.size(); ++i)
    //     std::cout << "("<<spline.x_values[i]<< ", " << spline.y_values[i]<<")\n";
    std::vector<bool> mask(y_values.size());
    std::transform(
        std::begin(y_values), 
        std::end(y_values), 
        std::begin(mask), 
        [&](auto&& i){
            return i > relative_threshold; //|| vol_in_range(i, avg_volume, max_upper_delta, max_lower_delta);
        });
    std::cout << "Fraction of segments touched: " << 1/(static_cast<float>(mask.size()) / std::count_if(std::begin(mask) , std::end(mask), [](auto&& val) {return !val;})) << "\n";
    return boost(samples, mask, spline, sample_rate, segment_size);
};

int main(int argc, char** argv) {
    std::vector<std::string> args{argv, argv+argc};
    if (args.size() < 3){
        std::cout << "Usage: " << args[0] << " [PATH_TO_AUDIO_FILE] [PATH_TO_OUTPUT]\n";
        return 0;
    }
    SndfileHandle f(args[1]);
    std::vector<float> frames(f.frames());
    f.readf(frames.data(), f.frames());
    int sample_rate = f.samplerate();
    // auto overall_loudness = calculate_lufs_with_gating(frames, sample_rate);
    // initial_tests(frames, f.samplerate());
    auto segment_size = sample_rate / 10;
    auto overall_loudness = calculate_lufs_with_gating(frames, sample_rate, segment_size);
    // We calculate loudness 
    auto result = smooth_clip(frames, sample_rate, segment_size);
    std::cout << "Writing " << args[2] << "\n";
    SndfileHandle result_file(args[2], SFM_WRITE, f.format(), f.channels(), sample_rate);
    result_file.write(result.data(), result.size());
    // initial_tests(result, sample_rate);
    auto result_loudness = calculate_lufs_with_gating(result, sample_rate, segment_size);
    std::cout << "Original loudness: " << overall_loudness << "\n";
    std::cout << "Result loudness: " << result_loudness << "\n";
    std::cout << f.channels() << "\n";
}
