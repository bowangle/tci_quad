#pragma once

#include <iostream>
#include <string>
#include <algorithm>
#include <limits>

namespace util {

using i128 = __int128;

// ===============================
// to_string
// ===============================
inline std::string to_string(i128 x) {
    if (x == 0) return "0";

    bool neg = x < 0;
    if (neg) x = -x;

    std::string s;

    while (x > 0) {
        s.push_back('0' + (x % 10));
        x /= 10;
    }

    if (neg) s.push_back('-');

    std::reverse(s.begin(), s.end());
    return s;
}


// ===============================
// string -> i128
// ===============================
inline i128 parse_i128(const std::string& s) {
    i128 x = 0;
    bool neg = false;

    size_t i = 0;
    if (s[0] == '-') {
        neg = true;
        i = 1;
    }

    for (; i < s.size(); i++) {
        x = x * 10 + (s[i] - '0');
    }

    return neg ? -x : x;
}

inline std::istream& operator>>(std::istream& is, i128& x) {
    std::string s;
    is >> s;
    x = parse_i128(s);
    return is;
}

// ===============================
// helpers
// ===============================
inline i128 abs(i128 x) {
    return x < 0 ? -x : x;
}

inline i128 pow10(int exp) {
    i128 r = 1;
    for (int i = 0; i < exp; i++) r *= 10;
    return r;
}

} // namespace util

// ===============================
// ostream support
// ===============================
inline std::ostream& operator<<(std::ostream& os, util::i128 x) {
    return os << util::to_string(x);
}
