// SPDX-License-Identifier: MIT
// Copyright (c) 2025 matthieu.jeannin


#ifndef MINDEX_H
#define MINDEX_H

#include <vector>
#include <ostream>
#include <initializer_list>
#include <utility>

struct MultiIndex {
    std::vector<int> data;

    MultiIndex() = default;
    explicit MultiIndex(size_t nBits)
        : data(nBits, 0) {}

    explicit MultiIndex(std::vector<int> v) : data(std::move(v)) {}
    MultiIndex(std::initializer_list<int> v) : data(v) {}

    MultiIndex operator+(const MultiIndex& other) const {
        MultiIndex out;
        out.data.reserve(data.size() + other.data.size());
        out.data.insert(out.data.end(), data.begin(), data.end());
        out.data.insert(out.data.end(), other.data.begin(), other.data.end());
        return out;
    }

    bool operator==(const MultiIndex& other) const noexcept {
        return data == other.data;
    }
};

// -------- stream output --------

inline std::ostream& operator<<(std::ostream& os, const MultiIndex& mi) {
    os << '[';
    for (std::size_t i = 0; i < mi.data.size(); ++i) {
        os << mi.data[i];
        if (i + 1 < mi.data.size())
            os << ',';
    }
    os << ']';
    return os;
}

#endif // MINDEX_H
