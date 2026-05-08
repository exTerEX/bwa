#include "bwa/index_io.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace {

constexpr std::array<char, 8> kMagicV1 = {'B', 'W', 'A', 'I', 'D', 'X', '1', '\0'};
constexpr std::array<char, 8> kMagicV2 = {'B', 'W', 'A', 'I', 'D', 'X', '2', '\0'};

template <typename T>
void write_scalar(std::ofstream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), static_cast<std::streamsize>(sizeof(T)));
    if (!out) {
        throw std::runtime_error("Failed while writing index file");
    }
}

template <typename T>
T read_scalar(std::ifstream& in) {
    T value{};
    in.read(reinterpret_cast<char*>(&value), static_cast<std::streamsize>(sizeof(T)));
    if (!in) {
        throw std::runtime_error("Failed while reading index file");
    }
    return value;
}

void write_bytes(std::ofstream& out, const void* data, std::size_t bytes) {
    if (bytes == 0U) {
        return;
    }
    out.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(bytes));
    if (!out) {
        throw std::runtime_error("Failed while writing index payload");
    }
}

void read_bytes(std::ifstream& in, void* data, std::size_t bytes) {
    if (bytes == 0U) {
        return;
    }
    in.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(bytes));
    if (!in) {
        throw std::runtime_error("Failed while reading index payload");
    }
}

}  // namespace

namespace bwa {

void save_index(const FMIndex& index, const std::filesystem::path& output_path) {
    std::ofstream out(output_path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open output index file: " + output_path.string());
    }

    if (index.sa_sample_rows().size() != index.sa_samples().size()) {
        throw std::runtime_error("Invalid index state: sampled SA row/value count mismatch");
    }

    write_bytes(out, kMagicV2.data(), kMagicV2.size());

    const uint64_t text_size = static_cast<uint64_t>(index.text_codes().size());
    const uint64_t bwt_size = static_cast<uint64_t>(index.bwt_codes().size());
    const uint64_t full_sa_size = static_cast<uint64_t>(index.suffix_array().size());
    const uint64_t sa_sample_rate = static_cast<uint64_t>(index.sa_sample_rate());
    const uint64_t sample_count = static_cast<uint64_t>(index.sa_sample_rows().size());

    write_scalar<uint64_t>(out, text_size);
    write_scalar<uint64_t>(out, bwt_size);
    write_scalar<uint64_t>(out, full_sa_size);
    write_scalar<uint64_t>(out, sa_sample_rate);
    write_scalar<uint64_t>(out, sample_count);

    write_bytes(out, index.text_codes().data(), static_cast<std::size_t>(text_size));
    write_bytes(out, index.bwt_codes().data(), static_cast<std::size_t>(bwt_size));
    write_bytes(out, index.suffix_array().data(), static_cast<std::size_t>(full_sa_size) * sizeof(uint64_t));
    write_bytes(out, index.sa_sample_rows().data(), static_cast<std::size_t>(sample_count) * sizeof(uint64_t));
    write_bytes(out, index.sa_samples().data(), static_cast<std::size_t>(sample_count) * sizeof(uint64_t));
}

FMIndex load_index(const std::filesystem::path& input_path) {
    std::ifstream in(input_path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open index file: " + input_path.string());
    }

    std::array<char, kMagicV1.size()> magic{};
    read_bytes(in, magic.data(), magic.size());

    const bool is_v1 = (std::memcmp(magic.data(), kMagicV1.data(), kMagicV1.size()) == 0);
    const bool is_v2 = (std::memcmp(magic.data(), kMagicV2.data(), kMagicV2.size()) == 0);
    if (!is_v1 && !is_v2) {
        throw std::runtime_error("Invalid index file magic: " + input_path.string());
    }

    const uint64_t text_size = read_scalar<uint64_t>(in);
    const uint64_t bwt_size = read_scalar<uint64_t>(in);

    if (is_v1) {
        const uint64_t sa_size = read_scalar<uint64_t>(in);

        std::vector<uint8_t> text(static_cast<std::size_t>(text_size), 0U);
        std::vector<uint8_t> bwt(static_cast<std::size_t>(bwt_size), 0U);
        std::vector<uint64_t> sa(static_cast<std::size_t>(sa_size), 0U);

        read_bytes(in, text.data(), text.size());
        read_bytes(in, bwt.data(), bwt.size());
        read_bytes(in, sa.data(), sa.size() * sizeof(uint64_t));

        FMIndex index;
        index.load_prebuilt(std::move(text), std::move(bwt), std::move(sa));
        return index;
    }

    const uint64_t full_sa_size = read_scalar<uint64_t>(in);
    const uint64_t sa_sample_rate = read_scalar<uint64_t>(in);
    const uint64_t sample_count = read_scalar<uint64_t>(in);

    std::vector<uint8_t> text(static_cast<std::size_t>(text_size), 0U);
    std::vector<uint8_t> bwt(static_cast<std::size_t>(bwt_size), 0U);
    std::vector<uint64_t> full_sa(static_cast<std::size_t>(full_sa_size), 0U);
    std::vector<uint64_t> sample_rows(static_cast<std::size_t>(sample_count), 0U);
    std::vector<uint64_t> sample_values(static_cast<std::size_t>(sample_count), 0U);

    read_bytes(in, text.data(), text.size());
    read_bytes(in, bwt.data(), bwt.size());
    read_bytes(in, full_sa.data(), full_sa.size() * sizeof(uint64_t));
    read_bytes(in, sample_rows.data(), sample_rows.size() * sizeof(uint64_t));
    read_bytes(in, sample_values.data(), sample_values.size() * sizeof(uint64_t));

    FMIndex index;
    index.load_prebuilt(
        std::move(text),
        std::move(bwt),
        std::move(full_sa),
        static_cast<std::size_t>(sa_sample_rate),
        std::move(sample_rows),
        std::move(sample_values));
    return index;
}

}  // namespace bwa
