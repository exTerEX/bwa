#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace bwa {

class BitVector {
  public:
    BitVector() = default;

    explicit BitVector(std::size_t bit_count) {
        reset(bit_count);
    }

    void reset(std::size_t bit_count) {
        size_bits_ = bit_count;
        words_.assign((bit_count + 63U) / 64U, 0U);
        superblock_rank_.clear();
        block_rank_.clear();
        ones_ = 0;
    }

    std::size_t size() const {
        return size_bits_;
    }

    void set(std::size_t pos) {
        assert(pos < size_bits_);
        words_[pos / 64U] |= (uint64_t{1} << (pos % 64U));
    }

    bool test(std::size_t pos) const {
        assert(pos < size_bits_);
        return (words_[pos / 64U] >> (pos % 64U)) & 1U;
    }

    void build_rank() {
        const std::size_t word_count = words_.size();
        const std::size_t superblock_count = (word_count + kWordsPerSuperblock - 1U) / kWordsPerSuperblock;

        superblock_rank_.assign(superblock_count + 1U, 0U);
        block_rank_.assign(word_count, 0U);

        std::size_t total = 0;
        for (std::size_t i = 0; i < word_count; ++i) {
            const std::size_t superblock = i / kWordsPerSuperblock;
            if (i % kWordsPerSuperblock == 0U) {
                superblock_rank_[superblock] = static_cast<uint64_t>(total);
            }
            const std::size_t rank_before_superblock = static_cast<std::size_t>(superblock_rank_[superblock]);
            block_rank_[i] = static_cast<uint16_t>(total - rank_before_superblock);
            total += static_cast<std::size_t>(std::popcount(words_[i]));
        }

        if (!superblock_rank_.empty()) {
            superblock_rank_.back() = static_cast<uint64_t>(total);
        }
        ones_ = total;
    }

    std::size_t rank1(std::size_t pos) const {
        if (pos > size_bits_) {
            pos = size_bits_;
        }
        if (pos == 0U || words_.empty()) {
            return 0U;
        }

        const std::size_t word = pos / 64U;
        const std::size_t offset = pos % 64U;

        if (word >= words_.size()) {
            return ones_;
        }

        const std::size_t superblock = word / kWordsPerSuperblock;
        std::size_t rank = static_cast<std::size_t>(superblock_rank_[superblock]);
        rank += static_cast<std::size_t>(block_rank_[word]);

        if (offset > 0U) {
            const uint64_t mask = (uint64_t{1} << offset) - 1U;
            rank += static_cast<std::size_t>(std::popcount(words_[word] & mask));
        }

        return rank;
    }

    std::optional<std::size_t> select1(std::size_t kth_zero_indexed) const {
        if (kth_zero_indexed >= ones_) {
            return std::nullopt;
        }

        const std::size_t target = kth_zero_indexed + 1U;

        std::size_t lo = 0U;
        std::size_t hi = superblock_rank_.size() - 1U;
        while (lo + 1U < hi) {
            const std::size_t mid = lo + (hi - lo) / 2U;
            if (static_cast<std::size_t>(superblock_rank_[mid]) < target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        const std::size_t start_word = lo * kWordsPerSuperblock;
        const std::size_t end_word = std::min(start_word + kWordsPerSuperblock, words_.size());

        for (std::size_t word = start_word; word < end_word; ++word) {
            const std::size_t rank_before_word = static_cast<std::size_t>(superblock_rank_[lo]) +
                                                 static_cast<std::size_t>(block_rank_[word]);
            const uint64_t bits = words_[word];
            const std::size_t pop = static_cast<std::size_t>(std::popcount(bits));
            if (rank_before_word + pop < target) {
                continue;
            }

            const std::size_t need = target - rank_before_word;
            const std::size_t bit_index = select_in_word(bits, need);
            const std::size_t absolute = word * 64U + bit_index;
            if (absolute < size_bits_) {
                return absolute;
            }
            return std::nullopt;
        }

        return std::nullopt;
    }

    std::size_t ones() const {
        return ones_;
    }

  private:
    static constexpr std::size_t kWordsPerSuperblock = 8U;

    static std::size_t select_in_word(uint64_t bits, std::size_t kth_one_1_based) {
        std::size_t remaining = kth_one_1_based;
        uint64_t value = bits;
        while (remaining > 1U) {
            value &= (value - 1U);
            --remaining;
        }
        return static_cast<std::size_t>(std::countr_zero(value));
    }

    std::size_t size_bits_ = 0U;
    std::vector<uint64_t> words_;
    std::vector<uint64_t> superblock_rank_;
    std::vector<uint16_t> block_rank_;
    std::size_t ones_ = 0U;
};

}  // namespace bwa
