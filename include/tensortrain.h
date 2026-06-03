// SPDX-License-Identifier: MIT
// Copyright (c) 2025 matthieu.jeannin
#pragma once

#ifndef TENSORTRAIN_H
#define TENSORTRAIN_H

#include <vector>
#include <stdexcept>
#include <Eigen/Dense>

// Reuse your TTCore alias
template <typename Scalar>
using TTCore = std::vector<
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
>;

template <typename Scalar>
class TensorTrain {
public:
    using Index = int;
    using Mat   = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    explicit TensorTrain(std::vector<TTCore<Scalar>> cores)
        : N_(cores.size()), cores_(std::move(cores))
    {
        if (N_ == 0)
            throw std::invalid_argument("TensorTrain: empty core list");

        // Basic consistency checks
        for (Index k = 0; k < N_; ++k) {
            if (cores_[k].empty())
                throw std::invalid_argument("TensorTrain: empty TT core");

            if (k > 0) {
                if (cores_[k-1][0].cols() != cores_[k][0].rows())
                    throw std::invalid_argument(
                        "TensorTrain: rank mismatch between cores"
                    );
            }
        }

        // boundary ranks must be 1
        if (cores_.front()[0].rows() != 1 ||
            cores_.back()[0].cols() != 1)
            throw std::invalid_argument(
                "TensorTrain: boundary TT ranks must be 1"
            );
    }

    Index order() const { return N_; }


    // ---------------- sum of all entries ----------------
    Scalar sum1() const
    {
        // Start with S^{(0)} = sum_p G^{(0)}[:,p,:]
        Mat acc = Mat::Zero(
            cores_[0][0].rows(),
            cores_[0][0].cols()
        );

        for (const Mat& Gp : cores_[0]) {
            acc += Gp;
        }

        // Successive contractions
        for (Index k = 1; k < N_; ++k) {
            Mat Sk = Mat::Zero(
                cores_[k][0].rows(),
                cores_[k][0].cols()
            );

            for (const Mat& Gp : cores_[k]) {
                Sk += Gp;
            }

            acc = acc * Sk;
        }

        // Result must be 1x1
        return acc(0, 0);
    }

    int max_bond_dimension() const
    {
        int max_r = 0;

        for (int k = 0; k < N_; ++k) {
            const auto& G0 = cores_[k][0];  // representative slice

            max_r = std::max(max_r, (int)G0.rows());
            max_r = std::max(max_r, (int)G0.cols());
        }

        return max_r;
    }

    // ---------------- evaluation ----------------
    // inds.size() == N
    Scalar eval(const std::vector<Index>& inds) const {
        if (inds.size() != static_cast<std::size_t>(N_))
            throw std::invalid_argument("TensorTrain::eval: wrong index size");

        // core = G1[:,i1,:]
        Mat core = slice_(0, inds[0]);

        // successive contractions
        for (Index k = 1; k < N_; ++k) {
            core = core * slice_(k, inds[k]);
        }

        // scalar result (1×1)
        return core(0, 0);
    }

private:
    // return G^{(k)}[:,p,:]
    const Mat& slice_(Index k, Index p) const {
        if (p < 0 || p >= static_cast<Index>(cores_[k].size()))
            throw std::out_of_range("TensorTrain: physical index out of range");
        return cores_[k][p];
    }

private:
    Index N_;
    std::vector<TTCore<Scalar>> cores_;
};

#endif // TENSORTRAIN_H


