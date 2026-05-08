#pragma once

#include "bwa/bit_vector.hpp"
#include "bwa/dna.hpp"
#include "bwa/sais.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace bwa {

enum class InexactSearchMode {
    Auto,
    Backtracking,
    SeedAndExtend,
};

struct SearchInterval {
    std::size_t begin = 0;
    std::size_t end = 0;

    [[nodiscard]] bool empty() const {
        return begin >= end;
    }

    [[nodiscard]] std::size_t size() const {
        return empty() ? 0U : (end - begin);
    }
};

struct Match {
    uint64_t position = 0;
    std::size_t edits = 0;
};

class FMIndex {
  public:
    void set_sa_sample_rate(std::size_t sample_rate) {
        sa_sample_rate_ = std::max<std::size_t>(1U, sample_rate);
    }

    [[nodiscard]] std::size_t sa_sample_rate() const {
        return sa_sample_rate_;
    }

    [[nodiscard]] bool using_sampled_sa_only() const {
        return sampled_sa_only_;
    }

    void build_from_reference(std::string_view reference) {
        const std::string normalized = dna::normalize_sequence(reference);
        text_codes_ = dna::encode_sequence(normalized);
        text_codes_.push_back(dna::kSentinel);

        const std::vector<int32_t> sa_i32 = sais::suffix_array_from_bytes(text_codes_);
        sa_.assign(sa_i32.begin(), sa_i32.end());

        if (!is_valid_sa_range(sa_, text_codes_.size())) {
            sa_ = build_suffix_array_fallback(text_codes_);
        }

        bwt_codes_.assign(sa_.size(), dna::kSentinel);
        for (std::size_t i = 0; i < sa_.size(); ++i) {
            const std::size_t suffix_pos = static_cast<std::size_t>(sa_[i]);
            const std::size_t prev = (suffix_pos == 0U) ? (text_codes_.size() - 1U) : (suffix_pos - 1U);
            bwt_codes_[i] = text_codes_[prev];
        }

        if (!is_lf_consistent(sa_, bwt_codes_)) {
            sa_ = build_suffix_array_fallback(text_codes_);
            bwt_codes_.assign(sa_.size(), dna::kSentinel);
            for (std::size_t i = 0; i < sa_.size(); ++i) {
                const std::size_t suffix_pos = static_cast<std::size_t>(sa_[i]);
                const std::size_t prev = (suffix_pos == 0U) ? (text_codes_.size() - 1U) : (suffix_pos - 1U);
                bwt_codes_[i] = text_codes_[prev];
            }
        }

        if (sa_sample_rate_ > 1U) {
            build_sa_samples_from_full();
        } else {
            sa_sample_rows_.clear();
            sa_samples_.clear();
            sa_sample_marks_.reset(bwt_codes_.size());
            sampled_sa_only_ = false;
        }

        rebuild_aux();
    }

    void load_prebuilt(
        std::vector<uint8_t> text_codes,
        std::vector<uint8_t> bwt_codes,
        std::vector<uint64_t> suffix_array,
        std::size_t sa_sample_rate = 1U,
        std::vector<uint64_t> sa_sample_rows = {},
        std::vector<uint64_t> sa_samples = {}) {
        text_codes_ = std::move(text_codes);
        bwt_codes_ = std::move(bwt_codes);
        sa_ = std::move(suffix_array);
        sa_sample_rate_ = std::max<std::size_t>(1U, sa_sample_rate);
        sa_sample_rows_ = std::move(sa_sample_rows);
        sa_samples_ = std::move(sa_samples);

        if (text_codes_.empty() || text_codes_.back() != dna::kSentinel) {
            throw std::runtime_error("Invalid index data: missing sentinel in text");
        }

        if (text_codes_.size() != bwt_codes_.size()) {
            throw std::runtime_error("Invalid index data: inconsistent text/BWT sizes");
        }

        if (!sa_.empty() && text_codes_.size() != sa_.size()) {
            throw std::runtime_error("Invalid index data: inconsistent full SA size");
        }

        if (!sa_sample_rows_.empty()) {
            build_sa_sample_marks_from_rows();
            sampled_sa_only_ = sa_.empty();
        } else if (sa_sample_rate_ > 1U && !sa_.empty()) {
            build_sa_samples_from_full();
        } else {
            sa_sample_marks_.reset(bwt_codes_.size());
            sa_samples_.clear();
            sampled_sa_only_ = false;
        }

        if (sa_.empty() && sa_sample_rows_.empty()) {
            throw std::runtime_error("Invalid index data: missing both full SA and sampled SA");
        }

        rebuild_aux();
    }

    [[nodiscard]] SearchInterval exact_interval(std::string_view pattern) const {
        if (bwt_codes_.empty()) {
            return {0U, 0U};
        }

        const std::vector<uint8_t> encoded = dna::encode_sequence(pattern);
        return exact_interval_encoded(encoded.data(), encoded.data() + encoded.size());
    }

    [[nodiscard]] std::vector<uint64_t> locate(
        SearchInterval interval,
        std::size_t max_hits = std::numeric_limits<std::size_t>::max()) const {
        std::vector<uint64_t> positions;
        if (interval.empty() || max_hits == 0U) {
            return positions;
        }

        const std::size_t max_pos = reference_length();
        positions.reserve(std::min(interval.size(), max_hits));

        for (std::size_t i = interval.begin; i < interval.end && positions.size() < max_hits; ++i) {
            const uint64_t pos = resolve_suffix_row(i);
            if (pos < max_pos) {
                positions.push_back(pos);
            }
        }

        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
        return positions;
    }

    [[nodiscard]] std::vector<Match> inexact_search(
        std::string_view pattern,
        std::size_t max_edits,
        std::size_t max_hits = 50U,
        InexactSearchMode mode = InexactSearchMode::Auto,
        std::size_t seed_candidate_limit = 4096U) const {
        std::vector<Match> matches;
        if ((sa_.empty() && sa_sample_rows_.empty()) || max_hits == 0U) {
            return matches;
        }

        const std::vector<uint8_t> query = dna::encode_sequence(pattern);
        if (query.empty()) {
            return matches;
        }

        if (max_edits == 0U) {
            const SearchInterval interval = exact_interval_encoded(query.data(), query.data() + query.size());
            const auto positions = locate(interval, max_hits);
            matches.reserve(positions.size());
            for (const uint64_t pos : positions) {
                matches.push_back(Match{pos, 0U});
            }
            return matches;
        }

        const bool should_try_seed =
            (mode == InexactSearchMode::SeedAndExtend) ||
            (mode == InexactSearchMode::Auto &&
             query.size() >= std::max<std::size_t>(24U, (max_edits + 1U) * 6U));

        if (should_try_seed) {
            matches = inexact_search_seed_extend(query, max_edits, max_hits, seed_candidate_limit);
            if (!matches.empty() || mode == InexactSearchMode::SeedAndExtend) {
                return matches;
            }
        }

        return inexact_search_backtracking(query, max_edits, max_hits);
    }

    [[nodiscard]] std::size_t reference_length() const {
        return text_codes_.empty() ? 0U : (text_codes_.size() - 1U);
    }

    [[nodiscard]] std::size_t suffix_array_entry_count() const {
        return bwt_codes_.size();
    }

    [[nodiscard]] std::string reference_string() const {
        if (text_codes_.empty()) {
            return {};
        }
        std::vector<uint8_t> without_sentinel(text_codes_.begin(), text_codes_.end() - 1);
        return dna::decode_sequence(without_sentinel, false);
    }

    [[nodiscard]] const std::vector<uint8_t>& text_codes() const {
        return text_codes_;
    }

    [[nodiscard]] const std::vector<uint8_t>& bwt_codes() const {
        return bwt_codes_;
    }

    [[nodiscard]] const std::vector<uint64_t>& suffix_array() const {
        return sa_;
    }

    [[nodiscard]] const std::vector<uint64_t>& sa_sample_rows() const {
        return sa_sample_rows_;
    }

    [[nodiscard]] const std::vector<uint64_t>& sa_samples() const {
        return sa_samples_;
    }

  private:
    [[nodiscard]] SearchInterval exact_interval_encoded(const uint8_t* begin, const uint8_t* end) const {
        if (bwt_codes_.empty()) {
            return {0U, 0U};
        }
        if ((begin == nullptr) || (end == nullptr) || (begin > end)) {
            return {0U, 0U};
        }

        SearchInterval interval{0U, bwt_codes_.size()};
        for (const uint8_t* it = end; it != begin;) {
            --it;
            interval = backward_extend(*it, interval);
            if (interval.empty()) {
                break;
            }
        }

        return interval;
    }

    [[nodiscard]] static bool is_valid_sa_range(const std::vector<uint64_t>& sa, std::size_t text_size) {
        if (sa.size() != text_size) {
            return false;
        }

        std::vector<uint8_t> seen(text_size, 0U);
        for (const uint64_t value : sa) {
            if (value >= static_cast<uint64_t>(text_size)) {
                return false;
            }
            const std::size_t idx = static_cast<std::size_t>(value);
            if (seen[idx] != 0U) {
                return false;
            }
            seen[idx] = 1U;
        }

        return true;
    }

    [[nodiscard]] static bool is_lf_consistent(const std::vector<uint64_t>& sa, const std::vector<uint8_t>& bwt) {
        if (sa.size() != bwt.size() || sa.empty()) {
            return false;
        }

        std::array<uint64_t, 256> counts{};
        counts.fill(0U);
        for (const uint8_t symbol : bwt) {
            ++counts[symbol];
        }

        std::array<uint64_t, 256> c_table{};
        c_table.fill(0U);

        uint64_t running = 0U;
        for (std::size_t i = 0; i < counts.size(); ++i) {
            c_table[i] = running;
            running += counts[i];
        }

        std::array<uint64_t, 256> seen{};
        seen.fill(0U);

        for (std::size_t row = 0; row < sa.size(); ++row) {
            const uint8_t symbol = bwt[row];
            const std::size_t lf_row = static_cast<std::size_t>(c_table[symbol] + seen[symbol]);
            const uint64_t expected = (sa[row] == 0U) ? static_cast<uint64_t>(sa.size() - 1U) : (sa[row] - 1U);

            if (lf_row >= sa.size() || sa[lf_row] != expected) {
                return false;
            }

            ++seen[symbol];
        }

        return true;
    }

    [[nodiscard]] static std::vector<uint64_t> build_suffix_array_fallback(const std::vector<uint8_t>& text) {
        std::vector<uint64_t> sa(text.size(), 0U);
        std::iota(sa.begin(), sa.end(), uint64_t{0});

        std::sort(sa.begin(), sa.end(), [&](uint64_t lhs, uint64_t rhs) {
            if (lhs == rhs) {
                return false;
            }

            const std::size_t n = text.size();
            std::size_t i = static_cast<std::size_t>(lhs);
            std::size_t j = static_cast<std::size_t>(rhs);

            while (i < n && j < n) {
                if (text[i] != text[j]) {
                    return text[i] < text[j];
                }
                ++i;
                ++j;
            }

            return i == n;
        });

        return sa;
    }

    [[nodiscard]] std::vector<Match> inexact_search_backtracking(
        const std::vector<uint8_t>& query,
        std::size_t max_edits,
        std::size_t max_hits) const {
        std::vector<Match> matches;
        if (search_symbols_.empty()) {
            return matches;
        }

        struct StateKey {
            int32_t pattern_index = -1;
            uint64_t left = 0;
            uint64_t right = 0;

            bool operator==(const StateKey& other) const {
                return pattern_index == other.pattern_index && left == other.left && right == other.right;
            }
        };

        struct StateHasher {
            std::size_t operator()(const StateKey& key) const {
                std::size_t h = static_cast<std::size_t>(key.left * 11400714819323198485ull);
                h ^= static_cast<std::size_t>(key.right + 0x9e3779b97f4a7c15ull + (h << 6U) + (h >> 2U));
                h ^= static_cast<std::size_t>(static_cast<uint64_t>(key.pattern_index) +
                                              0x517cc1b727220a95ull + (h << 6U) + (h >> 2U));
                return h;
            }
        };

        std::unordered_map<StateKey, std::size_t, StateHasher> best_state_edits;
        std::unordered_map<uint64_t, std::size_t> best_hit_edits;

        const SearchInterval full{0U, bwt_codes_.size()};

        const std::size_t locate_limit =
            (max_hits > std::numeric_limits<std::size_t>::max() / 4U) ? max_hits : (max_hits * 4U);

        std::function<void(int32_t, SearchInterval, std::size_t)> dfs;
        dfs = [&](int32_t pattern_index, SearchInterval interval, std::size_t edits) {
            if (edits > max_edits || interval.empty()) {
                return;
            }

            const StateKey key{pattern_index, static_cast<uint64_t>(interval.begin), static_cast<uint64_t>(interval.end)};
            const auto state_it = best_state_edits.find(key);
            if (state_it != best_state_edits.end() && state_it->second <= edits) {
                return;
            }
            best_state_edits[key] = edits;

            if (pattern_index < 0) {
                const auto located = locate(interval, locate_limit);
                for (const uint64_t pos : located) {
                    const auto hit_it = best_hit_edits.find(pos);
                    if (hit_it == best_hit_edits.end() || edits < hit_it->second) {
                        best_hit_edits[pos] = edits;
                    }
                }
                return;
            }

            const uint8_t expected = query[static_cast<std::size_t>(pattern_index)];

            for (const uint8_t symbol : search_symbols_) {
                const SearchInterval next = backward_extend(symbol, interval);
                if (next.empty()) {
                    continue;
                }
                const std::size_t next_edits = edits + ((symbol == expected) ? 0U : 1U);
                if (next_edits <= max_edits) {
                    dfs(pattern_index - 1, next, next_edits);
                }
            }

            if (edits < max_edits) {
                dfs(pattern_index - 1, interval, edits + 1U);

                for (const uint8_t symbol : search_symbols_) {
                    const SearchInterval next = backward_extend(symbol, interval);
                    if (!next.empty()) {
                        dfs(pattern_index, next, edits + 1U);
                    }
                }
            }
        };

        dfs(static_cast<int32_t>(query.size()) - 1, full, 0U);

        matches.reserve(best_hit_edits.size());
        for (const auto& [position, edits] : best_hit_edits) {
            matches.push_back(Match{position, edits});
        }

        std::sort(matches.begin(), matches.end(), [](const Match& lhs, const Match& rhs) {
            if (lhs.edits != rhs.edits) {
                return lhs.edits < rhs.edits;
            }
            return lhs.position < rhs.position;
        });

        if (matches.size() > max_hits) {
            matches.resize(max_hits);
        }

        return matches;
    }

    [[nodiscard]] std::vector<Match> inexact_search_seed_extend(
        const std::vector<uint8_t>& query,
        std::size_t max_edits,
        std::size_t max_hits,
        std::size_t seed_candidate_limit) const {
        std::vector<Match> matches;
        if (search_symbols_.empty() || query.empty() || seed_candidate_limit == 0U) {
            return matches;
        }

        const std::size_t part_count = max_edits + 1U;
        if (part_count == 0U || query.size() < part_count) {
            return matches;
        }

        struct Seed {
            std::size_t offset = 0U;
            std::size_t length = 0U;
        };

        std::vector<Seed> seeds;
        seeds.reserve(part_count);

        std::size_t cursor = 0U;
        const std::size_t base_len = query.size() / part_count;
        const std::size_t rem = query.size() % part_count;
        for (std::size_t i = 0; i < part_count; ++i) {
            const std::size_t len = base_len + ((i < rem) ? 1U : 0U);
            if (len == 0U) {
                continue;
            }
            seeds.push_back(Seed{cursor, len});
            cursor += len;
        }

        if (seeds.empty()) {
            return matches;
        }

        const std::size_t per_seed_locate_limit = std::min<std::size_t>(
            seed_candidate_limit,
            std::max<std::size_t>(max_hits * 8U, 128U));

        std::unordered_set<uint64_t> candidate_starts;
        candidate_starts.reserve(seed_candidate_limit);

        bool saturated = false;
        for (const Seed& seed : seeds) {
            const SearchInterval interval =
                exact_interval_encoded(query.data() + seed.offset, query.data() + seed.offset + seed.length);
            if (interval.empty()) {
                continue;
            }

            const auto seed_positions = locate(interval, per_seed_locate_limit);
            for (const uint64_t seed_pos : seed_positions) {
                const int64_t base_start = static_cast<int64_t>(seed_pos) - static_cast<int64_t>(seed.offset);

                for (int64_t shift = -static_cast<int64_t>(max_edits);
                     shift <= static_cast<int64_t>(max_edits);
                     ++shift) {
                    const int64_t candidate = base_start + shift;
                    if (candidate < 0) {
                        continue;
                    }

                    const std::size_t candidate_start = static_cast<std::size_t>(candidate);
                    if (candidate_start >= reference_length()) {
                        continue;
                    }

                    candidate_starts.insert(static_cast<uint64_t>(candidate_start));
                    if (candidate_starts.size() >= seed_candidate_limit) {
                        saturated = true;
                        break;
                    }
                }

                if (saturated) {
                    break;
                }
            }

            if (saturated) {
                break;
            }
        }

        std::unordered_map<uint64_t, std::size_t> best_hit_edits;
        best_hit_edits.reserve(candidate_starts.size());

        for (const uint64_t start : candidate_starts) {
            const std::size_t edits =
                bounded_edit_distance_at_start(query, static_cast<std::size_t>(start), max_edits);
            if (edits > max_edits) {
                continue;
            }

            const auto it = best_hit_edits.find(start);
            if (it == best_hit_edits.end() || edits < it->second) {
                best_hit_edits[start] = edits;
            }
        }

        matches.reserve(best_hit_edits.size());
        for (const auto& [position, edits] : best_hit_edits) {
            matches.push_back(Match{position, edits});
        }

        std::sort(matches.begin(), matches.end(), [](const Match& lhs, const Match& rhs) {
            if (lhs.edits != rhs.edits) {
                return lhs.edits < rhs.edits;
            }
            return lhs.position < rhs.position;
        });

        if (matches.size() > max_hits) {
            matches.resize(max_hits);
        }

        return matches;
    }

    [[nodiscard]] std::size_t bounded_edit_distance_at_start(
        const std::vector<uint8_t>& query,
        std::size_t reference_start,
        std::size_t max_edits) const {
        const std::size_t ref_len = reference_length();
        if (reference_start >= ref_len) {
            return max_edits + 1U;
        }

        const std::size_t query_len = query.size();
        const std::size_t max_ref_span = std::min<std::size_t>(ref_len - reference_start, query_len + max_edits);
        const std::size_t min_ref_span = (query_len > max_edits) ? (query_len - max_edits) : 0U;
        if (max_ref_span < min_ref_span) {
            return max_edits + 1U;
        }

        const std::size_t inf = max_edits + 1U;
        std::vector<std::size_t> prev(max_ref_span + 1U, inf);
        std::vector<std::size_t> curr(max_ref_span + 1U, inf);

        for (std::size_t j = 0U; j <= max_ref_span && j <= max_edits; ++j) {
            prev[j] = j;
        }

        for (std::size_t i = 1U; i <= query_len; ++i) {
            std::fill(curr.begin(), curr.end(), inf);

            const std::size_t j_start = (i > max_edits) ? (i - max_edits) : 0U;
            const std::size_t j_end = std::min<std::size_t>(max_ref_span, i + max_edits);

            if (j_start == 0U) {
                curr[0] = i;
            }

            std::size_t row_min = inf;
            for (std::size_t j = j_start; j <= j_end; ++j) {
                std::size_t best = curr[j];

                if (prev[j] < inf) {
                    best = std::min(best, prev[j] + 1U);
                }
                if (j > 0U && curr[j - 1U] < inf) {
                    best = std::min(best, curr[j - 1U] + 1U);
                }
                if (j > 0U && prev[j - 1U] < inf) {
                    const uint8_t ref_base = text_codes_[reference_start + (j - 1U)];
                    const std::size_t sub_cost = (query[i - 1U] == ref_base) ? 0U : 1U;
                    best = std::min(best, prev[j - 1U] + sub_cost);
                }

                curr[j] = best;
                row_min = std::min(row_min, best);
            }

            if (row_min > max_edits) {
                return max_edits + 1U;
            }

            prev.swap(curr);
        }

        const std::size_t lower = (query_len > max_edits) ? (query_len - max_edits) : 0U;
        const std::size_t upper = std::min<std::size_t>(max_ref_span, query_len + max_edits);

        std::size_t best = inf;
        for (std::size_t j = lower; j <= upper; ++j) {
            best = std::min(best, prev[j]);
        }

        return (best <= max_edits) ? best : (max_edits + 1U);
    }

    [[nodiscard]] uint64_t resolve_suffix_row(std::size_t row) const {
        if (!sa_.empty()) {
            return sa_[row];
        }

        if (sa_sample_rows_.empty() || sa_samples_.empty()) {
            throw std::runtime_error("No suffix array data available for locate");
        }

        const std::size_t n = bwt_codes_.size();
        std::size_t current = row;
        std::size_t steps = 0U;

        while (steps <= n) {
            if (sa_sample_marks_.test(current)) {
                const std::size_t rank = sa_sample_marks_.rank1(current + 1U);
                const uint64_t sampled_value = sa_samples_[rank - 1U];
                return (sampled_value + static_cast<uint64_t>(steps)) % static_cast<uint64_t>(n);
            }

            const uint8_t symbol = bwt_codes_[current];
            current = static_cast<std::size_t>(c_table_[symbol]) +
                      static_cast<std::size_t>(rank_symbol(symbol, current));
            ++steps;
        }

        throw std::runtime_error("LF walk exceeded text length while resolving sampled SA");
    }

    void build_sa_samples_from_full() {
        if (sa_.empty()) {
            throw std::runtime_error("Cannot build SA samples without a full suffix array");
        }

        sa_sample_rows_.clear();
        sa_samples_.clear();

        const std::size_t n = sa_.size();
        sa_sample_marks_.reset(n);

        const std::size_t reserve_hint = (n / sa_sample_rate_) + 2U;
        sa_sample_rows_.reserve(reserve_hint);
        sa_samples_.reserve(reserve_hint);

        for (std::size_t row = 0U; row < n; ++row) {
            if (sa_[row] % static_cast<uint64_t>(sa_sample_rate_) == 0U) {
                sa_sample_rows_.push_back(static_cast<uint64_t>(row));
                sa_samples_.push_back(sa_[row]);
                sa_sample_marks_.set(row);
            }
        }

        if (sa_sample_rows_.empty()) {
            throw std::runtime_error("SA sampling produced no samples");
        }

        sa_sample_marks_.build_rank();
        sampled_sa_only_ = true;
        sa_.clear();
    }

    void build_sa_sample_marks_from_rows() {
        if (sa_sample_rows_.size() != sa_samples_.size()) {
            throw std::runtime_error("Invalid sampled SA data: row/value count mismatch");
        }

        const std::size_t n = bwt_codes_.size();
        sa_sample_marks_.reset(n);

        uint64_t prev_row = 0U;
        bool first = true;
        for (const uint64_t row : sa_sample_rows_) {
            if (row >= static_cast<uint64_t>(n)) {
                throw std::runtime_error("Invalid sampled SA row outside index bounds");
            }
            if (!first && row <= prev_row) {
                throw std::runtime_error("Invalid sampled SA rows: expected strict ascending order");
            }

            sa_sample_marks_.set(static_cast<std::size_t>(row));
            prev_row = row;
            first = false;
        }

        sa_sample_marks_.build_rank();
    }

    [[nodiscard]] uint64_t rank_symbol(uint8_t symbol, std::size_t pos) const {
        if (pos > bwt_codes_.size()) {
            pos = bwt_codes_.size();
        }

        const int16_t occ_index = symbol_to_occ_[symbol];
        if (occ_index < 0) {
            return 0U;
        }

        return static_cast<uint64_t>(occ_[static_cast<std::size_t>(occ_index)].rank1(pos));
    }

    [[nodiscard]] SearchInterval backward_extend(uint8_t symbol, SearchInterval interval) const {
        if (interval.empty()) {
            return {0U, 0U};
        }

        const std::size_t left = static_cast<std::size_t>(c_table_[symbol]) +
                                 static_cast<std::size_t>(rank_symbol(symbol, interval.begin));
        const std::size_t right = static_cast<std::size_t>(c_table_[symbol]) +
                                  static_cast<std::size_t>(rank_symbol(symbol, interval.end));

        if (left >= right) {
            return {0U, 0U};
        }

        return {left, right};
    }

    void rebuild_aux() {
        if (text_codes_.size() != bwt_codes_.size()) {
            throw std::runtime_error("Index components have mismatched text/BWT lengths");
        }
        if (!sa_.empty() && text_codes_.size() != sa_.size()) {
            throw std::runtime_error("Index components have mismatched full SA length");
        }
        if (sa_.empty() && sa_sample_rows_.empty()) {
            throw std::runtime_error("Index has no suffix array representation");
        }
        if (!sa_sample_rows_.empty() && sa_sample_marks_.size() != bwt_codes_.size()) {
            throw std::runtime_error("Sampled SA bitvector has invalid length");
        }

        std::array<uint64_t, 256> counts{};
        counts.fill(0U);
        for (const uint8_t symbol : bwt_codes_) {
            ++counts[symbol];
        }

        uint64_t running = 0U;
        for (std::size_t i = 0; i < counts.size(); ++i) {
            c_table_[i] = running;
            running += counts[i];
        }
        c_table_[256] = running;

        symbol_to_occ_.fill(-1);
        symbols_.clear();
        search_symbols_.clear();
        occ_.clear();

        symbols_.reserve(8);
        search_symbols_.reserve(8);

        for (std::size_t symbol = 0; symbol < counts.size(); ++symbol) {
            if (counts[symbol] == 0U) {
                continue;
            }

            symbol_to_occ_[symbol] = static_cast<int16_t>(symbols_.size());
            symbols_.push_back(static_cast<uint8_t>(symbol));
            occ_.emplace_back(bwt_codes_.size());

            if (symbol != dna::kSentinel) {
                search_symbols_.push_back(static_cast<uint8_t>(symbol));
            }
        }

        for (std::size_t pos = 0; pos < bwt_codes_.size(); ++pos) {
            const uint8_t symbol = bwt_codes_[pos];
            const int16_t occ_index = symbol_to_occ_[symbol];
            if (occ_index >= 0) {
                occ_[static_cast<std::size_t>(occ_index)].set(pos);
            }
        }

        for (BitVector& bit_vector : occ_) {
            bit_vector.build_rank();
        }
    }

    std::vector<uint8_t> text_codes_;
    std::vector<uint8_t> bwt_codes_;
    std::vector<uint64_t> sa_;
    std::size_t sa_sample_rate_ = 1U;
    bool sampled_sa_only_ = false;
    std::vector<uint64_t> sa_sample_rows_;
    std::vector<uint64_t> sa_samples_;
    BitVector sa_sample_marks_;

    std::array<uint64_t, 257> c_table_{};
    std::array<int16_t, 256> symbol_to_occ_{};

    std::vector<uint8_t> symbols_;
    std::vector<uint8_t> search_symbols_;
    std::vector<BitVector> occ_;
};

}  // namespace bwa
