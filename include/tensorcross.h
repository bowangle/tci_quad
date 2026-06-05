// SPDX-License-Identifier: MIT
// Copyright (c) 2025 matthieu.jeannin
#pragma once

#ifndef TCI_H
#define TCI_H

#include <vector>
#include <functional>
#include <algorithm>
#include <stdexcept>

#include "pimat.h"
#include "tensortrain.h"
#include "mindex.h"

template <typename Scalar>
class TCI {
public:
    using Index      = int;
    using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

    using BaseFunc = std::function<Scalar(const MultiIndex&)>;
    using Core     = TTCore<Scalar>;

    // --------------------------------------------------
    // constructor
    // --------------------------------------------------
    TCI(BaseFunc f,
        const std::vector<Index>& l_d,
        const MultiIndex& pivot0)
        : f_(std::move(f)),
          l_d_(l_d),
          N_(static_cast<Index>(l_d.size())),
          forward_(true)
    {
        if (pivot0.data.size() != static_cast<std::size_t>(N_))
            throw std::invalid_argument("TCI: pivot0 has wrong dimension");

        l_pi_.reserve(static_cast<std::size_t>(N_ > 0 ? (N_ - 1) : 0));

        // build PiMat chain
        for (Index site = 0; site < N_ - 1; ++site) {

            // I = pivot0[:site]
            std::vector<MultiIndex> Iset = {
                MultiIndex(
                    std::vector<int>(
                        pivot0.data.begin(),
                        pivot0.data.begin() + site
                    )
                )
            };

            // J = pivot0[site+2:]
            std::vector<MultiIndex> Jset = {
                MultiIndex(
                    std::vector<int>(
                        pivot0.data.begin() + site + 2,
                        pivot0.data.end()
                    )
                )
            };

            // local physical indices
            std::vector<MultiIndex> d_loc1, d_loc2;
            for (Index d = 0; d < l_d_[static_cast<std::size_t>(site)]; ++d)
                d_loc1.emplace_back(std::initializer_list<int>{d});
            for (Index d = 0; d < l_d_[static_cast<std::size_t>(site + 1)]; ++d)
                d_loc2.emplace_back(std::initializer_list<int>{d});

            l_pi_.emplace_back(
                f_,
                std::move(Iset),
                std::move(Jset),
                std::move(d_loc1),
                std::move(d_loc2),
                pivot0.data[static_cast<std::size_t>(site)],
                pivot0.data[static_cast<std::size_t>(site + 1)]
            );
        }

        // IMPORTANT: store a REAL-valued error measure (magnitude)
        pivotError_.push_back(Eigen::numext::abs(f_(pivot0)));
    }

    // --------------------------------------------------
    // single-site iteration
    // --------------------------------------------------
    void iterate_at(Index site) {
        bool ok = l_pi_[static_cast<std::size_t>(site)].iterate();

        if (!ok) return;

        const MultiIndex& ipiv = l_pi_[static_cast<std::size_t>(site)].Ipiv().back();
        const MultiIndex& jpiv = l_pi_[static_cast<std::size_t>(site)].Jpiv().back();

        if (site < N_ - 2)
            l_pi_[static_cast<std::size_t>(site + 1)].add_I(ipiv);

        if (site > 0)
            l_pi_[static_cast<std::size_t>(site - 1)].add_J(jpiv);
    }


    // --------------------------------------------------
    // full sweep
    // --------------------------------------------------
    void iterate() {
        std::vector<RealScalar> l_err;

        if (N_ < 2) {
            // degenerate case: no PiMat blocks, keep last error
            pivotError_.push_back(pivotError_.empty() ? RealScalar(0) : pivotError_.back());
            forward_ = !forward_;
            return;
        }

        l_err.reserve(static_cast<std::size_t>(N_ - 1));

        if (forward_) {
            for (Index site = 0; site < N_ - 1; ++site) {
                iterate_at(site);
                l_err.push_back(
                    l_pi_[static_cast<std::size_t>(site)].pivot_errors().back()
                );
            }
        } else {
            for (Index site = N_ - 2; site >= 0; --site) {
                iterate_at(site);
                l_err.push_back(
                    l_pi_[static_cast<std::size_t>(site)].pivot_errors().back()
                );
                if (site == 0) break; // avoid signed underflow
            }
        }

        forward_ = !forward_;

        const RealScalar sweep_err =
            l_err.empty()
                ? RealScalar(0)
                : *std::max_element(l_err.begin(), l_err.end());

        pivotError_.push_back(sweep_err);
    }

    // --------------------------------------------------
    // core extraction
    // --------------------------------------------------
    Core get_core_at(Index site) const {
        if (site < N_ - 1)
            return l_pi_[static_cast<std::size_t>(site)].get_TPm();
        else
            return l_pi_[static_cast<std::size_t>(site - 1)].get_T_right();
    }

    // --------------------------------------------------
    // build TensorTrain
    // --------------------------------------------------
    TensorTrain<Scalar> get_TensorTrain() const {
        std::vector<Core> cores;
        cores.reserve(static_cast<std::size_t>(N_));

        for (Index site = 0; site < N_; ++site)
            cores.push_back(get_core_at(site));

        return TensorTrain<Scalar>(std::move(cores));
    }

    // --------------------------------------------------
    // accessors
    // --------------------------------------------------
    const std::vector<RealScalar>& pivotErrors() const {
        return pivotError_;
    }

    RealScalar last_pivot_error() const {
        if (pivotError_.empty())
            throw std::runtime_error("TCI: no pivot error available");
        return pivotError_.back();
    }

private:
    BaseFunc f_;
    std::vector<Index> l_d_;
    Index N_;

    bool forward_;
    std::vector<PiMat<Scalar>> l_pi_;

    // IMPORTANT: RealScalar (magnitudes), not Scalar
    std::vector<RealScalar> pivotError_;
};

#endif // TCI_H
