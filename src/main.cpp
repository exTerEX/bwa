#include "bwa/fasta.hpp"
#include "bwa/fm_index.hpp"
#include "bwa/index_io.hpp"

#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct IndexOptions {
    std::filesystem::path reference_path;
    std::filesystem::path output_path = "reference.bwaidx";
    std::size_t sa_sample_rate = 1U;
};

struct AlignOptions {
    std::filesystem::path index_path;
    std::filesystem::path reads_path;
    std::size_t max_edits = 0U;
    std::size_t max_hits = 20U;
    bwa::InexactSearchMode search_mode = bwa::InexactSearchMode::Auto;
    std::size_t seed_candidate_limit = 4096U;
};

void print_usage() {
    std::cout << "Usage:\n"
              << "  bwa index --reference ref.fasta --output ref.bwaidx\n"
              << "  bwa align --index ref.bwaidx --reads reads.fasta [options]\n\n"
              << "Index options:\n"
              << "  --sa-sample-rate <int>  Store every k-th SA value (default: 1, full SA)\n\n"
              << "Align options:\n"
              << "  --max-edits <int>   Max allowed edit operations (default: 0, exact mode)\n"
              << "  --max-hits <int>    Maximum hits per read to report (default: 20)\n"
              << "  --search-mode <auto|backtracking|seed>   Inexact search strategy (default: auto)\n"
              << "  --seed-candidates <int>   Max seed-and-extend candidate starts (default: 4096)\n"
              << "  --help              Show this message\n";
}

std::size_t parse_size(std::string_view value, const char* flag) {
    try {
        return static_cast<std::size_t>(std::stoull(std::string(value)));
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid value for ") + flag + ": " + std::string(value));
    }
}

std::string_view require_value(int argc, char** argv, int& i, const char* flag) {
    if (i + 1 >= argc) {
        throw std::runtime_error(std::string("Missing value for ") + flag);
    }
    ++i;
    return argv[i];
}

bwa::InexactSearchMode parse_search_mode(std::string_view value) {
    if (value == "auto") {
        return bwa::InexactSearchMode::Auto;
    }
    if (value == "backtracking") {
        return bwa::InexactSearchMode::Backtracking;
    }
    if (value == "seed" || value == "seed-and-extend") {
        return bwa::InexactSearchMode::SeedAndExtend;
    }
    throw std::runtime_error("Invalid value for --search-mode: " + std::string(value));
}

std::string_view search_mode_name(bwa::InexactSearchMode mode) {
    switch (mode) {
        case bwa::InexactSearchMode::Auto:
            return "auto";
        case bwa::InexactSearchMode::Backtracking:
            return "backtracking";
        case bwa::InexactSearchMode::SeedAndExtend:
            return "seed";
    }
    return "auto";
}

IndexOptions parse_index_options(int argc, char** argv) {
    IndexOptions options;

    for (int i = 2; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--reference") {
            options.reference_path = std::string(require_value(argc, argv, i, "--reference"));
            continue;
        }
        if (arg == "--output") {
            options.output_path = std::string(require_value(argc, argv, i, "--output"));
            continue;
        }
        if (arg == "--sa-sample-rate") {
            options.sa_sample_rate = parse_size(require_value(argc, argv, i, "--sa-sample-rate"), "--sa-sample-rate");
            continue;
        }
        if (arg == "--help") {
            print_usage();
            std::exit(0);
        }

        throw std::runtime_error("Unknown index argument: " + std::string(arg));
    }

    if (options.reference_path.empty()) {
        throw std::runtime_error("Missing required --reference argument");
    }

    return options;
}

AlignOptions parse_align_options(int argc, char** argv) {
    AlignOptions options;

    for (int i = 2; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--index") {
            options.index_path = std::string(require_value(argc, argv, i, "--index"));
            continue;
        }
        if (arg == "--reads") {
            options.reads_path = std::string(require_value(argc, argv, i, "--reads"));
            continue;
        }
        if (arg == "--max-edits") {
            options.max_edits = parse_size(require_value(argc, argv, i, "--max-edits"), "--max-edits");
            continue;
        }
        if (arg == "--max-hits") {
            options.max_hits = parse_size(require_value(argc, argv, i, "--max-hits"), "--max-hits");
            continue;
        }
        if (arg == "--search-mode") {
            options.search_mode = parse_search_mode(require_value(argc, argv, i, "--search-mode"));
            continue;
        }
        if (arg == "--seed-candidates") {
            options.seed_candidate_limit =
                parse_size(require_value(argc, argv, i, "--seed-candidates"), "--seed-candidates");
            continue;
        }
        if (arg == "--help") {
            print_usage();
            std::exit(0);
        }

        throw std::runtime_error("Unknown align argument: " + std::string(arg));
    }

    if (options.index_path.empty()) {
        throw std::runtime_error("Missing required --index argument");
    }
    if (options.reads_path.empty()) {
        throw std::runtime_error("Missing required --reads argument");
    }

    return options;
}

std::string sanitize_name(std::string name) {
    for (char& ch : name) {
        if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
            ch = '_';
        }
    }
    return name;
}

std::string format_positions(const std::vector<uint64_t>& positions) {
    if (positions.empty()) {
        return "-";
    }

    std::string joined;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (i > 0U) {
            joined.push_back(',');
        }
        joined += std::to_string(positions[i]);
    }
    return joined;
}

std::string format_matches(const std::vector<bwa::Match>& matches) {
    if (matches.empty()) {
        return "-";
    }

    std::string joined;
    for (std::size_t i = 0; i < matches.size(); ++i) {
        if (i > 0U) {
            joined.push_back(',');
        }
        joined += std::to_string(matches[i].position);
        joined.push_back(':');
        joined += std::to_string(matches[i].edits);
    }
    return joined;
}

int run_index(const IndexOptions& options) {
    const auto reference_records = bwa::load_fasta_fastq(options.reference_path);
    const std::string reference = bwa::concatenate_reference(reference_records, 'N');

    bwa::FMIndex index;
    index.set_sa_sample_rate(options.sa_sample_rate);
    index.build_from_reference(reference);

    bwa::save_index(index, options.output_path);

    std::cout << "command=index"
              << " reference_records=" << reference_records.size()
              << " reference_bases=" << index.reference_length()
              << " sa_entries=" << index.suffix_array_entry_count()
              << " sa_sample_rate=" << index.sa_sample_rate()
              << " sampled_sa=" << (index.using_sampled_sa_only() ? 1 : 0)
              << " output=" << options.output_path << '\n';

    return 0;
}

int run_align(const AlignOptions& options) {
    const bwa::FMIndex index = bwa::load_index(options.index_path);
    const auto reads = bwa::load_fasta_fastq(options.reads_path);

    std::cout << "command=align"
              << " reads=" << reads.size()
              << " reference_bases=" << index.reference_length()
              << " max_edits=" << options.max_edits
              << " max_hits=" << options.max_hits
              << " search_mode=" << search_mode_name(options.search_mode)
              << " seed_candidates=" << options.seed_candidate_limit
              << '\n';

    for (const auto& read : reads) {
        const std::string name = sanitize_name(read.name);

        if (options.max_edits == 0U) {
            const bwa::SearchInterval interval = index.exact_interval(read.sequence);
            const auto positions = index.locate(interval, options.max_hits);

            std::cout << "read=" << name
                      << " len=" << read.sequence.size()
                      << " hits=" << positions.size()
                      << " positions=" << format_positions(positions)
                      << '\n';
        } else {
            const auto matches = index.inexact_search(
                read.sequence,
                options.max_edits,
                options.max_hits,
                options.search_mode,
                options.seed_candidate_limit);
            const std::size_t best_edits = matches.empty() ? (options.max_edits + 1U) : matches.front().edits;

            std::cout << "read=" << name
                      << " len=" << read.sequence.size()
                      << " hits=" << matches.size()
                      << " best_edits=" << best_edits
                      << " matches=" << format_matches(matches)
                      << '\n';
        }
    }

    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            print_usage();
            return 1;
        }

        const std::string_view command = argv[1];
        if (command == "index") {
            const IndexOptions options = parse_index_options(argc, argv);
            return run_index(options);
        }
        if (command == "align") {
            const AlignOptions options = parse_align_options(argc, argv);
            return run_align(options);
        }
        if (command == "--help" || command == "help") {
            print_usage();
            return 0;
        }

        throw std::runtime_error("Unknown command: " + std::string(command));
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
