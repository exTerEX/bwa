#pragma once

#include "bwa/fm_index.hpp"

#include <filesystem>

namespace bwa {

void save_index(const FMIndex& index, const std::filesystem::path& output_path);

FMIndex load_index(const std::filesystem::path& input_path);

}  // namespace bwa
