#include "bwa/bit_vector.hpp"
#include "bwa/fm_index.hpp"
#include "bwa/index_io.hpp"
#include "bwa/sais.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

#define CHECK_OR_FAIL(condition, message)      \
    do {                                        \
        if (!(condition)) {                     \
            std::cerr << "[FAIL] " << message \
                      << "\n";               \
            return false;                       \
        }                                       \
    } while (false)

bool test_bit_vector_rank_select() {
    bwa::BitVector bits(128U);
    bits.set(1U);
    bits.set(3U);
    bits.set(64U);
    bits.set(127U);
    bits.build_rank();

    CHECK_OR_FAIL(bits.rank1(0U) == 0U, "rank at 0 should be zero");
    CHECK_OR_FAIL(bits.rank1(2U) == 1U, "rank should count first set bit");
    CHECK_OR_FAIL(bits.rank1(4U) == 2U, "rank should count two early set bits");
    CHECK_OR_FAIL(bits.rank1(65U) == 3U, "rank should include position 64");
    CHECK_OR_FAIL(bits.rank1(128U) == 4U, "rank at end should count all ones");

    const auto s0 = bits.select1(0U);
    const auto s1 = bits.select1(1U);
    const auto s2 = bits.select1(2U);
    const auto s3 = bits.select1(3U);

    CHECK_OR_FAIL(s0.has_value() && *s0 == 1U, "select(0) should return first set bit");
    CHECK_OR_FAIL(s1.has_value() && *s1 == 3U, "select(1) should return second set bit");
    CHECK_OR_FAIL(s2.has_value() && *s2 == 64U, "select(2) should return third set bit");
    CHECK_OR_FAIL(s3.has_value() && *s3 == 127U, "select(3) should return fourth set bit");

    return true;
}

bool test_sais_suffix_array() {
    const std::vector<uint8_t> text = {'b', 'a', 'n', 'a', 'n', 'a', 0U};
    const std::vector<int32_t> sa = bwa::sais::suffix_array_from_bytes(text);
    const std::vector<int32_t> expected = {6, 5, 3, 1, 0, 4, 2};
    CHECK_OR_FAIL(sa == expected, "SA-IS output mismatch on banana$");
    return true;
}

bool test_fm_index_exact_search() {
    bwa::FMIndex index;
    index.build_from_reference("GATTACA");

    const bwa::SearchInterval interval = index.exact_interval("TACA");
    const auto positions = index.locate(interval, 10U);

    CHECK_OR_FAIL(positions.size() == 1U, "exact interval should locate one position");
    CHECK_OR_FAIL(positions[0] == 3U, "TACA should start at offset 3");

    const bwa::SearchInterval absent = index.exact_interval("CCCC");
    CHECK_OR_FAIL(absent.empty(), "absent k-mer should return empty interval");

    return true;
}

bool has_match(const std::vector<bwa::Match>& matches, uint64_t position, std::size_t max_edits) {
    for (const auto& match : matches) {
        if (match.position == position && match.edits <= max_edits) {
            return true;
        }
    }
    return false;
}

bool test_fm_index_inexact_search() {
    bwa::FMIndex index;
    index.build_from_reference("GATTACA");

    const auto mismatch = index.inexact_search("GACTACA", 1U, 10U);
    CHECK_OR_FAIL(has_match(mismatch, 0U, 1U), "one-mismatch read should map with <=1 edit");

    const auto deletion = index.inexact_search("GATACA", 1U, 10U);
    CHECK_OR_FAIL(has_match(deletion, 0U, 1U), "one-deletion read should map with <=1 edit");

    const auto insertion = index.inexact_search("GATTTACA", 1U, 10U);
    CHECK_OR_FAIL(has_match(insertion, 0U, 1U), "one-insertion read should map with <=1 edit");

    return true;
}

std::size_t best_edit_for_position(const std::vector<bwa::Match>& matches, uint64_t position, std::size_t fallback) {
    std::size_t best = fallback;
    for (const auto& match : matches) {
        if (match.position == position && match.edits < best) {
            best = match.edits;
        }
    }
    return best;
}

bool test_seed_and_extend_mode() {
    const std::string reference = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
    bwa::FMIndex index;
    index.build_from_reference(reference);

    const std::string query = "ACGTTGCATGTAGCATGATGCATGAGAGCT";

    const auto backtracking = index.inexact_search(
        query,
        1U,
        20U,
        bwa::InexactSearchMode::Backtracking,
        4096U);

    const auto seeded = index.inexact_search(
        query,
        1U,
        20U,
        bwa::InexactSearchMode::SeedAndExtend,
        4096U);

    CHECK_OR_FAIL(has_match(backtracking, 0U, 1U), "backtracking mode should map query at position 0");
    CHECK_OR_FAIL(has_match(seeded, 0U, 1U), "seed-and-extend mode should map query at position 0");

    const std::size_t backtracking_best = best_edit_for_position(backtracking, 0U, 99U);
    const std::size_t seeded_best = best_edit_for_position(seeded, 0U, 99U);
    CHECK_OR_FAIL(backtracking_best == seeded_best, "seed mode should match backtracking best edit at true locus");

    return true;
}

bool test_sampled_sa_locate_and_roundtrip() {
    const std::string reference = "TTTACGGTACCGTTAACCGGTTACGGTACCGTTAACCGGAAA";

    bwa::FMIndex full_index;
    full_index.build_from_reference(reference);

    bwa::FMIndex sampled_index;
    sampled_index.set_sa_sample_rate(4U);
    sampled_index.build_from_reference(reference);

    CHECK_OR_FAIL(sampled_index.using_sampled_sa_only(), "sampled index should use sampled SA-only mode");
    CHECK_OR_FAIL(sampled_index.suffix_array().empty(), "sampled index should not keep full SA in memory");
    CHECK_OR_FAIL(!sampled_index.sa_sample_rows().empty(), "sampled index should retain sampled rows");

    const bwa::SearchInterval full_interval = full_index.exact_interval("ACCGTT");
    const bwa::SearchInterval sampled_interval = sampled_index.exact_interval("ACCGTT");

    const auto full_positions = full_index.locate(full_interval, 50U);
    const auto sampled_positions = sampled_index.locate(sampled_interval, 50U);
    CHECK_OR_FAIL(full_positions == sampled_positions, "sampled SA locate should match full SA locate results");

    const std::filesystem::path tmp = std::filesystem::temp_directory_path() / "bwa_sampled_roundtrip.bwaidx";
    bwa::save_index(sampled_index, tmp);
    bwa::FMIndex loaded = bwa::load_index(tmp);
    std::filesystem::remove(tmp);

    CHECK_OR_FAIL(loaded.using_sampled_sa_only(), "loaded sampled index should remain in sampled SA-only mode");

    const auto loaded_positions = loaded.locate(loaded.exact_interval("ACCGTT"), 50U);
    CHECK_OR_FAIL(loaded_positions == full_positions, "round-tripped sampled index should preserve locate results");

    const auto seeded_matches = loaded.inexact_search(
        "ACCGTTA",
        1U,
        20U,
        bwa::InexactSearchMode::SeedAndExtend,
        4096U);
    CHECK_OR_FAIL(!seeded_matches.empty(), "round-tripped sampled index should support seed-and-extend matching");

    return true;
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, bool (*)()>> tests = {
        {"bit-vector rank/select", test_bit_vector_rank_select},
        {"SA-IS suffix array", test_sais_suffix_array},
        {"FM-index exact search", test_fm_index_exact_search},
        {"FM-index inexact search", test_fm_index_inexact_search},
        {"seed-and-extend mode", test_seed_and_extend_mode},
        {"sampled SA locate and roundtrip", test_sampled_sa_locate_and_roundtrip},
    };

    bool all_ok = true;
    for (const auto& test : tests) {
        const bool ok = test.second();
        if (ok) {
            std::cout << "[PASS] " << test.first << "\n";
        } else {
            all_ok = false;
        }
    }

    return all_ok ? 0 : 1;
}
