#include "bwa/fasta.hpp"

#include "bwa/dna.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string trim_cr(std::string line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    return line;
}

bool looks_like_fastq_record(const std::vector<std::string>& lines, std::size_t index) {
    if (index + 2U >= lines.size()) {
        return false;
    }
    return !lines[index].empty() && lines[index][0] == '@' &&
           !lines[index + 2U].empty() && lines[index + 2U][0] == '+';
}

}  // namespace

namespace bwa {

std::vector<SequenceRecord> load_fasta_fastq(const std::filesystem::path& input_path) {
    std::ifstream in(input_path);
    if (!in) {
        throw std::runtime_error("Failed to open input file: " + input_path.string());
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(trim_cr(line));
    }

    std::vector<SequenceRecord> records;
    std::size_t i = 0U;
    while (i < lines.size()) {
        if (lines[i].empty()) {
            ++i;
            continue;
        }

        if (lines[i][0] == '>') {
            const std::string name = lines[i].substr(1);
            ++i;

            std::string sequence;
            while (i < lines.size()) {
                if (lines[i].empty()) {
                    ++i;
                    continue;
                }
                if (lines[i][0] == '>' || looks_like_fastq_record(lines, i)) {
                    break;
                }
                sequence += dna::normalize_sequence(lines[i]);
                ++i;
            }

            if (sequence.empty()) {
                throw std::runtime_error("Empty FASTA sequence for record: " + name);
            }
            records.push_back(SequenceRecord{name, sequence});
            continue;
        }

        if (looks_like_fastq_record(lines, i)) {
            const std::string name = lines[i].substr(1);
            const std::string sequence = dna::normalize_sequence(lines[i + 1U]);
            i += 3U;

            std::size_t quality_chars = 0U;
            while (i < lines.size() && quality_chars < sequence.size()) {
                quality_chars += lines[i].size();
                ++i;
            }

            if (quality_chars < sequence.size()) {
                throw std::runtime_error("Malformed FASTQ quality section for record: " + name);
            }
            if (sequence.empty()) {
                throw std::runtime_error("Empty FASTQ sequence for record: " + name);
            }
            records.push_back(SequenceRecord{name, sequence});
            continue;
        }

        throw std::runtime_error("Unrecognized record at line " + std::to_string(i + 1U) + " in " + input_path.string());
    }

    if (records.empty()) {
        throw std::runtime_error("No FASTA/FASTQ records found in: " + input_path.string());
    }

    return records;
}

std::string concatenate_reference(const std::vector<SequenceRecord>& records, char spacer) {
    if (records.empty()) {
        throw std::runtime_error("No records provided for reference concatenation");
    }

    std::string reference;
    const char normalized_spacer = dna::decode_base(dna::encode_base(spacer));

    for (const SequenceRecord& record : records) {
        const std::string cleaned = dna::normalize_sequence(record.sequence);
        if (cleaned.empty()) {
            continue;
        }
        if (!reference.empty()) {
            reference.push_back(normalized_spacer);
        }
        reference += cleaned;
    }

    if (reference.empty()) {
        throw std::runtime_error("All reference records are empty after normalization");
    }

    return reference;
}

}  // namespace bwa
