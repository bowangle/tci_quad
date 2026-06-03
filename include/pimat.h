// SPDX-License-Identifier: MIT
// Copyright (c) 2025 matthieu.jeannin
#pragma once

#ifndef PIMAT_H
#define PIMAT_H

#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

#include "lazymat.h"
#include "prrrlu.h"
#include "mindex.h"

template <typename Scalar>
using TTCore = std::vector<
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
>;

// ---------------- PiMat ----------------

template <typename Scalar>
class PiMat {
public:
    using Index      = int;
    using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

    using Mat      = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vec      = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using BaseFunc = std::function<Scalar(const MultiIndex&)>;
    using Core     = TTCore<Scalar>;

    PiMat(BaseFunc f,
          std::vector<MultiIndex> Iset,
          std::vector<MultiIndex> Jset,
          std::vector<MultiIndex> d_loc1,
          std::vector<MultiIndex> d_loc2,
          Index i_ini = 0,
          Index j_ini = 0)
        : f_base_(std::move(f)),
          Iset_(std::move(Iset)),
          Jset_(std::move(Jset)),
          d_loc1_(std::move(d_loc1)),
          d_loc2_(std::move(d_loc2)),
          nI_(Iset_.size()),
          nJ_(Jset_.size()),
          n1_(d_loc1_.size()),
          n2_(d_loc2_.size()),
          f_local_([this](Index i_loc, Index j_loc) -> Scalar {
              Index i_p = i_loc / static_cast<Index>(n1_);
              Index p1  = i_loc % static_cast<Index>(n1_);
              Index p2  = j_loc % static_cast<Index>(n2_);
              Index j_p = j_loc / static_cast<Index>(n2_);

              MultiIndex idx =
                  Iset_[static_cast<std::size_t>(i_p)] +
                  d_loc1_[static_cast<std::size_t>(p1)] +
                  d_loc2_[static_cast<std::size_t>(p2)] +
                  Jset_[static_cast<std::size_t>(j_p)];

              return f_base_(idx);
          }),
          M_(static_cast<Index>(nI_ * n1_),
             static_cast<Index>(nJ_ * n2_),
             f_local_)
    {
        prrlu_ = std::make_unique<PRRLU<Scalar>>(
            M_, std::make_pair(i_ini, j_ini)
        );

        auto [ip, jp] = local_to_multi(i_ini, j_ini);
        Ipiv_.push_back(ip);
        Jpiv_.push_back(jp);

        // IMPORTANT: pivot errors are RealScalar-valued (magnitudes)
        pivot_error_ = prrlu_->pivot_errors();
    }

    // ---------------- local → multi ----------------

    std::pair<MultiIndex, MultiIndex>
    local_to_multi(Index i_loc, Index j_loc) const {
        Index i_p = i_loc / static_cast<Index>(n1_);
        Index p1  = i_loc % static_cast<Index>(n1_);
        Index p2  = j_loc % static_cast<Index>(n2_);
        Index j_p = j_loc / static_cast<Index>(n2_);

        return {
            Iset_[static_cast<std::size_t>(i_p)] + d_loc1_[static_cast<std::size_t>(p1)],
            d_loc2_[static_cast<std::size_t>(p2)] + Jset_[static_cast<std::size_t>(j_p)]
        };
    }

    // ---------------- iteration ----------------

    bool iterate(bool full = false, Index nbrook = 5) {
        bool ok = prrlu_->iterate(full, nbrook);

        if (ok) {
            Index i = prrlu_->pivot_rows().back();
            Index j = prrlu_->pivot_cols().back();

            auto [ip, jp] = local_to_multi(i, j);
            Ipiv_.push_back(ip);
            Jpiv_.push_back(jp);

            pivot_error_ = prrlu_->pivot_errors();
        } else {
            pivot_error_.push_back(RealScalar(0));
        }
        return ok;
    }

    // ---------------- enlarge ----------------

    void add_I(const MultiIndex& I_new) {
        Iset_.push_back(I_new);
        ++nI_;
        rebuild_and_enlarge_();
    }

    void add_J(const MultiIndex& J_new) {
        Jset_.push_back(J_new);
        ++nJ_;
        rebuild_and_enlarge_();
    }

    // ---------------- TT cores ----------------

    Core get_TPm() const {
        const Mat& Lr = prrlu_->L();
        Index r = static_cast<Index>(Lr.rows());
        Index n = static_cast<Index>(Lr.cols());
        if (r == 0) return {};

        Mat L = Lr.transpose();

        Mat Mpiv(r, r);
        const auto& piv = prrlu_->pivot_rows();
        for (Index k = 0; k < r; ++k) {
            Mpiv.row(k) = L.row(piv[static_cast<std::size_t>(k)]);
        }

        Mat Minv = Mpiv.fullPivLu().solve(Mat::Identity(r, r));
        Mat T = L * Minv;

        Core core(static_cast<std::size_t>(n1_));
        for (Index p1 = 0; p1 < static_cast<Index>(n1_); ++p1) {
            Mat slice(static_cast<Index>(nI_), r);
            for (Index i = 0; i < static_cast<Index>(nI_); ++i) {
                slice.row(i) = T.row(i * static_cast<Index>(n1_) + p1);
            }
            core[static_cast<std::size_t>(p1)] = std::move(slice);
        }
        return core;
    }

    Core get_T_left() const {
        Index r = static_cast<Index>(Jpiv_.size());
        Core core(static_cast<std::size_t>(n1_));
        for (Index p1 = 0; p1 < static_cast<Index>(n1_); ++p1) {
            Mat slice(static_cast<Index>(nI_), r);
            for (Index i = 0; i < static_cast<Index>(nI_); ++i) {
                for (Index k = 0; k < r; ++k) {
                    slice(i, k) =
                        f_base_(Iset_[static_cast<std::size_t>(i)] +
                                d_loc1_[static_cast<std::size_t>(p1)] +
                                Jpiv_[static_cast<std::size_t>(k)]);
                }
            }
            core[static_cast<std::size_t>(p1)] = std::move(slice);
        }
        return core;
    }

    Core get_T_right() const {
        Index r = static_cast<Index>(Ipiv_.size());
        Core core(static_cast<std::size_t>(n2_));
        for (Index p2 = 0; p2 < static_cast<Index>(n2_); ++p2) {
            Mat slice(r, static_cast<Index>(nJ_));
            for (Index k = 0; k < r; ++k) {
                for (Index j = 0; j < static_cast<Index>(nJ_); ++j) {
                    slice(k, j) =
                        f_base_(Ipiv_[static_cast<std::size_t>(k)] +
                                d_loc2_[static_cast<std::size_t>(p2)] +
                                Jset_[static_cast<std::size_t>(j)]);
                }
            }
            core[static_cast<std::size_t>(p2)] = std::move(slice);
        }
        return core;
    }

    // Pivot errors are REAL magnitudes, even when Scalar is complex.
    const std::vector<RealScalar>& pivot_errors() const {
        return pivot_error_;
    }

    RealScalar last_pivot_error() const {
        if (pivot_error_.empty())
            throw std::runtime_error("PiMat: no pivot error available");
        return pivot_error_.back();
    }

    const std::vector<MultiIndex>& Ipiv() const { return Ipiv_; }
    const std::vector<MultiIndex>& Jpiv() const { return Jpiv_; }

private:
    void rebuild_and_enlarge_() {
        LazyMat<Scalar> M_big(
            static_cast<Index>(nI_ * n1_),
            static_cast<Index>(nJ_ * n2_),
            f_local_
        );
        prrlu_->enlarge(M_big);
        M_ = std::move(M_big);
    }

private:
    BaseFunc f_base_;
    std::vector<MultiIndex> Iset_, Jset_;
    std::vector<MultiIndex> d_loc1_, d_loc2_;
    std::vector<MultiIndex> Ipiv_, Jpiv_;

    std::size_t nI_, nJ_, n1_, n2_;

    std::function<Scalar(Index, Index)> f_local_;
    LazyMat<Scalar> M_;
    std::unique_ptr<PRRLU<Scalar>> prrlu_;

    // IMPORTANT: RealScalar (magnitudes), not Scalar
    std::vector<RealScalar> pivot_error_;
};

#endif // PIMAT_H
