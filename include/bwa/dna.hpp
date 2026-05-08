#pragma once

#include <cctype>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace bwa::dna {

constexpr uint8_t kSentinel = 0;
constexpr uint8_t kA = 1;
constexpr uint8_t kC = 2;
constexpr uint8_t kG = 3;
constexpr uint8_t kT = 4;
constexpr uint8_t kN = 5;

inline uint8_t encode_base(char base) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A':
            return kA;
        case 'C':
            return kC;
        case 'G':
            return kG;
        case 'T':
            return kT;
        default:
            return kN;
    }
}

inline char decode_base(uint8_t code) {
    switch (code) {
        case kSentinel:
            return '$';
        case kA:
            return 'A';
        case kC:
            return 'C';
        case kG:
            return 'G';
        case kT:
            return 'T';
        default:
            return 'N';
    }
}

inline std::vector<uint8_t> encode_sequence(std::string_view sequence) {
    std::vector<uint8_t> encoded;
    encoded.reserve(sequence.size());
    for (const char ch : sequence) {
        if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
            continue;
        }
        encoded.push_back(encode_base(ch));
    }
    return encoded;
}

inline std::string decode_sequence(const std::vector<uint8_t>& encoded, bool stop_at_sentinel = false) {
    std::string sequence;
    sequence.reserve(encoded.size());
    for (const uint8_t code : encoded) {
        if (stop_at_sentinel && code == kSentinel) {
            break;
        }
        sequence.push_back(decode_base(code));
    }
    return sequence;
}

inline std::string normalize_sequence(std::string_view sequence) {
    std::string normalized;
    normalized.reserve(sequence.size());
    for (const char ch : sequence) {
        if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
            continue;
        }
        normalized.push_back(decode_base(encode_base(ch)));
    }
    return normalized;
}

}  // namespace bwa::dna
