#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace bwa {

struct SequenceRecord {
    std::string name;
    std::string sequence;
};

std::vector<SequenceRecord> load_fasta_fastq(const std::filesystem::path& input_path);

std::string concatenate_reference(const std::vector<SequenceRecord>& records, char spacer = 'N');

}  // namespace bwa
