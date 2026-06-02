// SPDX-License-Identifier: MIT
// Copyright (c) 2025 matthieu.jeannin

#pragma once

#include <Eigen/Dense>
#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include <iostream>

template <typename Scalar>
class LazyMat {
public:
    using Index = int;
    using Func  = std::function<Scalar(Index, Index)>;

    LazyMat(Index rows, Index cols, Func f)
        : n_(rows), m_(cols), f_(std::move(f)) {}

    Index rows() const { return n_; }
    Index cols() const { return m_; }

    /* -------------------------
       Scalar access
       ------------------------- */
    Scalar operator()(Index i, Index j) const {
        check_bounds(i, j);
        return get_entry(i, j);
    }

    /* -------------------------
       Full row / column
       ------------------------- */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
    row(Index i) const {
        return row(i, 0, m_);
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
    col(Index j) const {
        return col(j, 0, n_);
    }

    /* -------------------------
       Partial row / column
       ------------------------- */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
    row(Index i, Index j0, Index j1) const {
        check_row(i);
        check_col_range(j0, j1);

        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(j1 - j0);
        for (Index j = j0; j < j1; ++j)
            v(j - j0) = get_entry(i, j);
        return v;
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
    col(Index j, Index i0, Index i1) const {
        check_col(j);
        check_row_range(i0, i1);

        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(i1 - i0);
        for (Index i = i0; i < i1; ++i)
            v(i - i0) = get_entry(i, j);
        return v;
    }

    /* -------------------------
       2D block
       ------------------------- */
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    block(Index i0, Index i1, Index j0, Index j1) const {
        check_row_range(i0, i1);
        check_col_range(j0, j1);

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
            A(i1 - i0, j1 - j0);

        for (Index i = i0; i < i1; ++i)
            for (Index j = j0; j < j1; ++j)
                A(i - i0, j - j0) = get_entry(i, j);

        return A;
    }

    /* -------------------------
       Full materialization
       ------------------------- */
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    eval() const {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
            A(n_, m_);

        for (Index i = 0; i < n_; ++i)
            for (Index j = 0; j < m_; ++j)
                A(i, j) = get_entry(i, j);

        return A;
    }

    std::size_t cached_entries() const {
        return cache_.size();
    }

private:
    Index n_, m_;
    Func f_;

    mutable std::unordered_map<std::uint64_t, Scalar> cache_;

    static std::uint64_t key(Index i, Index j) {
        return (static_cast<std::uint64_t>(i) << 32)
             | static_cast<std::uint32_t>(j);
    }

    Scalar get_entry(Index i, Index j) const {
        auto k = key(i, j);
        auto it = cache_.find(k);
        if (it != cache_.end())
            return it->second;

        Scalar v = f_(i, j);
        cache_.emplace(k, v);
        return v;
    }

    void check_bounds(Index i, Index j) const {
        check_row(i);
        check_col(j);
    }

    void check_row(Index i) const {
        if (i < 0 || i >= n_)
            throw std::out_of_range("Row index out of range");
    }

    void check_col(Index j) const {
        if (j < 0 || j >= m_)
            throw std::out_of_range("Column index out of range");
    }

    void check_row_range(Index i0, Index i1) const {
        if (i0 < 0 || i1 > n_ || i0 > i1)
            throw std::out_of_range("Invalid row range");
    }

    void check_col_range(Index j0, Index j1) const {
        if (j0 < 0 || j1 > m_ || j0 > j1)
            throw std::out_of_range("Invalid column range");
    }
};
