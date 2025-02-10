#include <sndfile.hh>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "smooth_audio.h"

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
    auto segment_size = sample_rate / 100;
    smooth_audio::boost_sample(frames);
    std::cout << "Writing " << args[2] << "\n";
    SndfileHandle result_file(args[2], SFM_WRITE, f.format(), f.channels(), sample_rate);
    result_file.write(frames.data(), frames.size());
    auto original_loudness = smooth_audio::calculate_lufs_with_gating(frames, sample_rate, segment_size);
    auto result_loudness = smooth_audio::calculate_lufs_with_gating(frames, sample_rate, segment_size);
    std::cout << "Original loudness: " << original_loudness << "\n";
    std::cout << "Result loudness: " << result_loudness << "\n";
}
