// SPDX-License-Identifier: MIT
// Copyright (c) 2025 matthieu.jeannin


#ifndef PRRRLU_H
#define PRRRLU_H

#include <eigen3/Eigen/Dense>
#include <vector>
#include <utility>
#include <optional>
#include <algorithm>
#include <stdexcept>
#include "lazymat.h"

template <typename Scalar>
class PRRLU {
public:
    using Index      = int;
    using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

    using Vec   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Mat   = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Pivot = std::pair<Index, Index>;

    // M is a LazyMat, ini_piv is (i, j)
    PRRLU(LazyMat<Scalar> M, Pivot ini_piv)
        : M_(std::move(M)),
          n_(M_.rows()),
          m_(M_.cols())
    {
        const Index i = ini_piv.first;
        const Index j = ini_piv.second;

        pivot_.push_back({i, j});
        pivot_rows_.push_back(i);
        pivot_cols_.push_back(j);

        const Scalar d0 = M_(i, j);
        D_.push_back(d0);

        // For complex (and even real), use magnitude to decide "zero pivot"
        if (Eigen::numext::abs(d0) == RealScalar(0)) {
            throw std::runtime_error("Initial pivot is zero; cannot initialize PRRLU.");
        }

        // L0 = M[:, j] / d0
        // U0 = M[i, :] / d0
        const Vec L0 = M_.col(j) / d0;
        const Vec U0 = M_.row(i) / d0;

        // Store as rank-rows: L, U have shape (rank x n) and (rank x m)
        L_.resize(1, n_);
        U_.resize(1, m_);
        L_.row(0) = L0.transpose();
        U_.row(0) = U0.transpose();

        pivotError_.push_back(Eigen::numext::abs(d0)); // RealScalar
    }

    Index rows() const { return n_; }
    Index cols() const { return m_; }
    Index rank() const { return static_cast<Index>(D_.size()); }

    // reconstruct() = sum_k D[k] * outer(L[k], U[k])
    Mat reconstruct() const {
        Mat R = Mat::Zero(n_, m_);
        const Index r = rank();
        for (Index k = 0; k < r; ++k) {
            const Vec Lk = L_.row(k).transpose(); // (n)
            const Vec Uk = U_.row(k).transpose(); // (m)
            R.noalias() += D_[k] * (Lk * Uk.transpose());
        }
        return R;
    }

    // reconstruct_ij(i,j) = sum_k D[k] * L(k,i) * U(k,j)
    Scalar reconstruct_ij(Index i, Index j) const {
        const Index r = rank();
        Scalar s = Scalar(0);
        for (Index k = 0; k < r; ++k) {
            s += D_[k] * L_(k, i) * U_(k, j);
        }
        return s;
    }

    Scalar residual_ij(Index i, Index j) const {
        if (contains(pivot_rows_, i)) return Scalar(0);
        if (contains(pivot_cols_, j)) return Scalar(0);
        return M_(i, j) - reconstruct_ij(i, j);
    }

    // residual_full(): returns dense residual matrix with pivot rows/cols zeroed.
    Mat residual_full() const {
        Mat res = M_.eval() - reconstruct();

        // enforce structural zeros
        for (Index i : pivot_rows_) res.row(i).setZero();
        for (Index j : pivot_cols_) res.col(j).setZero();

        return res;
    }

    LazyMat<Scalar> residual_lazy() const {
        // WARNING: captures `this`. Returned LazyMat must not outlive PRRLU.
        return LazyMat<Scalar>(n_, m_, [this](Index i, Index j) -> Scalar {
            return this->residual_ij(i, j);
        });
    }

    // addPivot(ip, jp, row, col) where row is residual row at ip, col residual col at jp
    // addPivot(ip, jp, row, col) where row is residual row at ip, col residual col at jp
    Scalar addPivot(Index ip, Index jp, const Vec& row, const Vec& col) {
        const Scalar d_new = row(jp);
        if (Eigen::numext::abs(d_new) == RealScalar(0)) {
            // this pivot is useless
            return Scalar(0);
        }

        pivot_.push_back({ip, jp});
        pivot_rows_.push_back(ip);
        pivot_cols_.push_back(jp);

        const Vec u_new = row / d_new; // length m
        const Vec l_new = col / d_new; // length n

        D_.push_back(d_new);

        // Stack onto L_, U_ (new row appended)
        const Index k = rank() - 1; // new k after push
        L_.conservativeResize(k + 1, n_);
        U_.conservativeResize(k + 1, m_);
        L_.row(k) = l_new.transpose();
        U_.row(k) = u_new.transpose();

        pivotError_.push_back(Eigen::numext::abs(d_new)); // RealScalar
        return d_new;
    }

    bool iterate(bool full = true, Index nbrook = 5) {
        std::optional<SearchResult> res;

        if (full) res = full_search();
        else      res = rook_search(nbrook);

        if (!res.has_value()) return false;

        const auto& out = *res;
        return Eigen::numext::abs(addPivot(out.ip, out.jp, out.row, out.col)) > RealScalar(0);
    }

    // full_search(): compute dense residual, find argmax |residual|, return pivot info.
    struct SearchResult {
        Index ip;
        Index jp;
        Vec row; // length m
        Vec col; // length n
    };

    std::optional<SearchResult> full_search() const {
        Mat R = residual_full();

        // IMPORTANT: for complex Scalars, cwiseAbs() is RealScalar-valued
        const auto A = R.cwiseAbs(); // Matrix<RealScalar, ...>

        Index ip = 0, jp = 0;
        RealScalar best = RealScalar(0);

        for (Index i = 0; i < n_; ++i) {
            for (Index j = 0; j < m_; ++j) {
                const RealScalar v = A(i, j);
                if (v > best) {
                    best = v;
                    ip = i; jp = j;
                }
            }
        }

        if (!(best > RealScalar(0))) {
            return std::nullopt;
        }

        // row = res[ip, :], col = res[:, jp]
        Vec row(m_), col(n_);
        row = R.row(ip).transpose();
        col = R.col(jp);

        return SearchResult{ip, jp, row, col};
    }



    // rook_search(max_iter): pivoting using lazy residual
    std::optional<SearchResult> rook_search(Index max_iter = 5) const {
        auto R = residual_lazy();

        Index ip = first_non_pivot_row();
        Index jp = first_non_pivot_col();
        if (ip < 0 || jp < 0) return std::nullopt;

        for (Index it = 0; it < max_iter; ++it) {
            // Maximize over current row ip
            Index j_new = jp;
            RealScalar best_row = RealScalar(0);

            for (Index j = 0; j < m_; ++j) {
                const RealScalar v = Eigen::numext::abs(R(ip, j));
                if (v > best_row) {
                    best_row = v;
                    j_new = j;
                }
            }

            if (!(best_row > RealScalar(0))) return std::nullopt;

            // Maximize over selected column j_new
            Index i_new = ip;
            RealScalar best_col = RealScalar(0);

            for (Index i = 0; i < n_; ++i) {
                const RealScalar v = Eigen::numext::abs(R(i, j_new));
                if (v > best_col) {
                    best_col = v;
                    i_new = i;
                }
            }

            if (!(best_col > RealScalar(0))) return std::nullopt;

            const Index ip_old = ip;
            const Index jp_old = jp;

            ip = i_new;
            jp = j_new;

            // Rook condition: current row/column maximizers agree
            if (ip == ip_old && jp == jp_old) {
                break;
            }
        }

        // Return the current candidate, even if max_iter was reached.
        Vec row(m_), col(n_);
        for (Index j = 0; j < m_; ++j) row(j) = R(ip, j);
        for (Index i = 0; i < n_; ++i) col(i) = R(i, jp);

        if (!(Eigen::numext::abs(row(jp)) > RealScalar(0))) {
            return std::nullopt;
        }

        return SearchResult{ip, jp, row, col};
    }

    // enlarge: assuming M_big's top-left is old M.
    void enlarge(const LazyMat<Scalar>& M_big) {
        const Index r_old = n_;
        const Index c_old = m_;
        const Index r_new = M_big.rows();
        const Index c_new = M_big.cols();

        const Index added_rows = r_new - r_old;
        const Index added_cols = c_new - c_old;

        const Index r = rank();

        // STEP 1: extend L and U with zeros
        if (added_rows > 0) {
            Mat L_ext(r, r_new);
            L_ext.leftCols(r_old) = L_;
            L_ext.rightCols(added_rows).setZero();
            L_ = std::move(L_ext);
        }

        if (added_cols > 0) {
            Mat U_ext(r, c_new);
            U_ext.leftCols(c_old) = U_;
            U_ext.rightCols(added_cols).setZero();
            U_ = std::move(U_ext);
        }

        // store enlarged M and update dims
        M_ = M_big;
        n_ = r_new;
        m_ = c_new;

        // STEP 3: update extended parts of L[k] and U[k]
        for (Index k = 0; k < r; ++k) {
            const Index ip = pivot_[k].first;
            const Index jp = pivot_[k].second;
            const Scalar d_k = D_[k];

            if (Eigen::numext::abs(d_k) == RealScalar(0)) continue;

            if (added_rows > 0) {
                for (Index rr = r_old; rr < r_new; ++rr) {
                    Scalar res = M_(rr, jp);
                    for (Index i = 0; i < k; ++i) {
                        res -= D_[i] * L_(i, rr) * U_(i, jp);
                    }
                    L_(k, rr) = res / d_k;
                }
            }

            if (added_cols > 0) {
                for (Index cc = c_old; cc < c_new; ++cc) {
                    Scalar res = M_(ip, cc);
                    for (Index i = 0; i < k; ++i) {
                        res -= D_[i] * L_(i, ip) * U_(i, cc);
                    }
                    U_(k, cc) = res / d_k;
                }
            }
        }
    }

    // Pivot errors are REAL (magnitudes), even when Scalar is complex.
    std::vector<RealScalar> pivot_errors() const {
        return pivotError_;
    }

    RealScalar last_pivot_error() const {
        if (pivotError_.empty())
            throw std::runtime_error("No pivots have been added yet");
        return pivotError_.back();
    }

    const Mat& L() const { return L_; }
    const Mat& U() const { return U_; }
    const std::vector<Index>& pivot_rows() const { return pivot_rows_; }
    const std::vector<Index>& pivot_cols() const { return pivot_cols_; }

private:
    LazyMat<Scalar> M_;
    Index n_{0}, m_{0};

    std::vector<Pivot> pivot_;
    std::vector<Index> pivot_rows_;
    std::vector<Index> pivot_cols_;

    std::vector<Scalar> D_;
    Mat L_; // (rank x n)
    Mat U_; // (rank x m)

    std::vector<RealScalar> pivotError_; // magnitudes

    static bool contains(const std::vector<Index>& v, Index x) {
        return std::find(v.begin(), v.end(), x) != v.end();
    }

    Index first_non_pivot_row() const {
        for (Index i = 0; i < n_; ++i) {
            if (!contains(pivot_rows_, i)) return i;
        }
        return -1;
    }

    Index first_non_pivot_col() const {
        for (Index j = 0; j < m_; ++j) {
            if (!contains(pivot_cols_, j)) return j;
        }
        return -1;
    }
};

#endif // PRRRLU_H
