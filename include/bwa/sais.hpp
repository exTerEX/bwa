#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace bwa::sais {

namespace detail {

inline bool is_lms(const std::vector<uint8_t>& is_s_type, int32_t i) {
    return i > 0 && is_s_type[static_cast<std::size_t>(i)] != 0 &&
           is_s_type[static_cast<std::size_t>(i - 1)] == 0;
}

inline std::vector<int32_t> bucket_sizes(const std::vector<int32_t>& s, int32_t alphabet_size) {
    std::vector<int32_t> sizes(static_cast<std::size_t>(alphabet_size), 0);
    for (const int32_t symbol : s) {
        ++sizes[static_cast<std::size_t>(symbol)];
    }
    return sizes;
}

inline std::vector<int32_t> bucket_heads(const std::vector<int32_t>& sizes) {
    std::vector<int32_t> heads(sizes.size(), 0);
    int32_t sum = 0;
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        heads[i] = sum;
        sum += sizes[i];
    }
    return heads;
}

inline std::vector<int32_t> bucket_tails(const std::vector<int32_t>& sizes) {
    std::vector<int32_t> tails(sizes.size(), 0);
    int32_t sum = 0;
    for (std::size_t i = 0; i < sizes.size(); ++i) {
        sum += sizes[i];
        tails[i] = sum;
    }
    return tails;
}

inline bool lms_substring_equal(
    const std::vector<int32_t>& s,
    const std::vector<uint8_t>& is_s_type,
    int32_t lhs,
    int32_t rhs) {
    if (lhs == rhs) {
        return true;
    }

    const int32_t n = static_cast<int32_t>(s.size());
    int32_t offset = 0;
    while (true) {
        const int32_t li = lhs + offset;
        const int32_t ri = rhs + offset;
        if (li >= n || ri >= n) {
            return li == ri;
        }

        if (s[static_cast<std::size_t>(li)] != s[static_cast<std::size_t>(ri)] ||
            is_s_type[static_cast<std::size_t>(li)] != is_s_type[static_cast<std::size_t>(ri)]) {
            return false;
        }

        const bool lhs_is_lms = is_lms(is_s_type, li);
        const bool rhs_is_lms = is_lms(is_s_type, ri);
        if (lhs_is_lms && rhs_is_lms) {
            return true;
        }
        if (lhs_is_lms != rhs_is_lms) {
            return false;
        }

        ++offset;
    }
}

inline std::vector<int32_t> sais_impl(const std::vector<int32_t>& s, int32_t alphabet_size) {
    const int32_t n = static_cast<int32_t>(s.size());
    std::vector<int32_t> sa(static_cast<std::size_t>(n), -1);

    if (n == 0) {
        return sa;
    }
    if (n == 1) {
        sa[0] = 0;
        return sa;
    }

    std::vector<uint8_t> is_s_type(static_cast<std::size_t>(n), 0);
    is_s_type[static_cast<std::size_t>(n - 1)] = 1;
    for (int32_t i = n - 2; i >= 0; --i) {
        const int32_t curr = s[static_cast<std::size_t>(i)];
        const int32_t next = s[static_cast<std::size_t>(i + 1)];
        if (curr < next) {
            is_s_type[static_cast<std::size_t>(i)] = 1;
        } else if (curr == next) {
            is_s_type[static_cast<std::size_t>(i)] = is_s_type[static_cast<std::size_t>(i + 1)];
        }
    }

    const std::vector<int32_t> sizes = bucket_sizes(s, alphabet_size);

    auto induce_sort = [&](const std::vector<int32_t>& ordered_lms) {
        std::fill(sa.begin(), sa.end(), -1);

        auto tails = bucket_tails(sizes);
        for (int32_t i = static_cast<int32_t>(ordered_lms.size()) - 1; i >= 0; --i) {
            const int32_t p = ordered_lms[static_cast<std::size_t>(i)];
            const int32_t c = s[static_cast<std::size_t>(p)];
            sa[static_cast<std::size_t>(--tails[static_cast<std::size_t>(c)])] = p;
        }

        auto heads = bucket_heads(sizes);
        for (int32_t i = 0; i < n; ++i) {
            const int32_t p = sa[static_cast<std::size_t>(i)];
            if (p <= 0) {
                continue;
            }
            const int32_t prev = p - 1;
            if (is_s_type[static_cast<std::size_t>(prev)] == 0) {
                const int32_t c = s[static_cast<std::size_t>(prev)];
                sa[static_cast<std::size_t>(heads[static_cast<std::size_t>(c)]++)] = prev;
            }
        }

        tails = bucket_tails(sizes);
        for (int32_t i = n - 1; i >= 0; --i) {
            const int32_t p = sa[static_cast<std::size_t>(i)];
            if (p <= 0) {
                continue;
            }
            const int32_t prev = p - 1;
            if (is_s_type[static_cast<std::size_t>(prev)] != 0) {
                const int32_t c = s[static_cast<std::size_t>(prev)];
                sa[static_cast<std::size_t>(--tails[static_cast<std::size_t>(c)])] = prev;
            }
        }
    };

    std::vector<int32_t> lms_positions;
    lms_positions.reserve(static_cast<std::size_t>(n));
    for (int32_t i = 1; i < n; ++i) {
        if (is_lms(is_s_type, i)) {
            lms_positions.push_back(i);
        }
    }

    induce_sort(lms_positions);

    std::vector<int32_t> sorted_lms;
    sorted_lms.reserve(lms_positions.size());
    for (const int32_t p : sa) {
        if (p >= 0 && is_lms(is_s_type, p)) {
            sorted_lms.push_back(p);
        }
    }

    std::vector<int32_t> lms_name(static_cast<std::size_t>(n), -1);
    int32_t current_name = 0;
    int32_t previous_lms = -1;
    for (const int32_t p : sorted_lms) {
        if (previous_lms >= 0 && !lms_substring_equal(s, is_s_type, previous_lms, p)) {
            ++current_name;
        }
        lms_name[static_cast<std::size_t>(p)] = current_name;
        previous_lms = p;
    }

    const int32_t name_count = current_name + 1;

    std::vector<int32_t> reduced;
    reduced.reserve(lms_positions.size());
    for (const int32_t p : lms_positions) {
        reduced.push_back(lms_name[static_cast<std::size_t>(p)]);
    }

    std::vector<int32_t> reduced_sa;
    if (name_count == static_cast<int32_t>(reduced.size())) {
        reduced_sa.assign(reduced.size(), -1);
        for (int32_t i = 0; i < static_cast<int32_t>(reduced.size()); ++i) {
            reduced_sa[static_cast<std::size_t>(reduced[static_cast<std::size_t>(i)])] = i;
        }
    } else {
        reduced_sa = sais_impl(reduced, name_count);
    }

    std::vector<int32_t> ordered_lms;
    ordered_lms.reserve(lms_positions.size());
    for (const int32_t idx : reduced_sa) {
        ordered_lms.push_back(lms_positions[static_cast<std::size_t>(idx)]);
    }

    induce_sort(ordered_lms);
    return sa;
}

}  // namespace detail

inline std::vector<int32_t> build_suffix_array(const std::vector<int32_t>& symbols, int32_t alphabet_size) {
    if (symbols.empty()) {
        return {};
    }
    if (alphabet_size <= 0) {
        throw std::runtime_error("alphabet_size must be positive");
    }
    if (symbols.back() != 0) {
        throw std::runtime_error("input symbols must be sentinel-terminated with 0");
    }
    for (const int32_t value : symbols) {
        if (value < 0 || value >= alphabet_size) {
            throw std::runtime_error("input symbol outside alphabet range");
        }
    }

    return detail::sais_impl(symbols, alphabet_size);
}

inline std::vector<int32_t> suffix_array_from_bytes(const std::vector<uint8_t>& text) {
    if (text.empty()) {
        return {};
    }

    int32_t max_symbol = 0;
    std::vector<int32_t> symbols;
    symbols.reserve(text.size());
    for (const uint8_t value : text) {
        symbols.push_back(static_cast<int32_t>(value));
        if (static_cast<int32_t>(value) > max_symbol) {
            max_symbol = static_cast<int32_t>(value);
        }
    }

    return build_suffix_array(symbols, max_symbol + 1);
}

}  // namespace bwa::sais
