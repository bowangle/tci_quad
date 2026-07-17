#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <complex>

#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/complex128.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "tensorcross.h"
#include "tensortrain.h"
#include "mindex.h"

// ============================================================
//  The function to learn:
//
//    f(x_0, ..., x_{N-1})  =  exp(-0.1 * sum_i x_i)
//                            * cos(0.3 * sum_i (i+1)*x_i)
//                            / (1 + 0.05 * sum_i x_i^2)
//
//  Indices x_i in {0, ..., d-1}.  d=6, N=8 -> 6^8 ~ 1.7M entries.
//  The function is smooth but nontrivial in every variable.
// ============================================================

// ----- double version -----
double f_double(const MultiIndex& mi)
{
    const int N = static_cast<int>(mi.data.size());
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for (int k = 0; k < N; ++k) {
        double xk = static_cast<double>(mi.data[k]);
        s1 += xk;
        s2 += static_cast<double>(k + 1) * xk;
        s3 += xk * xk;
    }
    return std::exp(-0.1 * s1) * std::cos(0.3 * s2) / (1.0 + 0.05 * s3);
}

// ----- quad version -----
using float128 = boost::multiprecision::float128;

float128 f_quad(const MultiIndex& mi)
{
    const int N = static_cast<int>(mi.data.size());
    float128 s1 = 0, s2 = 0, s3 = 0;
    for (int k = 0; k < N; ++k) {
        float128 xk(mi.data[k]);
        s1 += xk;
        s2 += float128(k + 1) * xk;
        s3 += xk * xk;
    }
    using boost::multiprecision::exp;
    using boost::multiprecision::cos;
    return exp(float128("-0.1") * s1)
         * cos(float128("0.3") * s2)
         / (float128(1) + float128("0.05") * s3);
}

using cquad = std::complex<float128>;
using complex128 = boost::multiprecision::complex128;
// ----- quad complex version -----
cquad f_complex128(const MultiIndex& mi)
{
    const int N = static_cast<int>(mi.data.size());

    cquad s1 = cquad(0.0), s2 = cquad(0.0), s3 = cquad(0.0);

    for (int k = 0; k < N; ++k) {
        cquad xk(mi.data[k]);

        s1 += xk;
        s2 += cquad(k + 1) * xk;
        s3 += xk * xk;
    }

    using std::exp;
    using std::cos;

    return exp(cquad(-0.1) * s1)
         * cos(cquad(0.3) * s2)
         / (cquad(1) + cquad(0.05) * s3);
}

// ============================================================
//  Helper: run TCI and report errors
// ============================================================
template <typename Scalar, typename Func>
void run_demo(const std::string& label,
              Func               f_eval,
              int                N,      // number of sites
              int                d,      // physical dimension per site
              int                n_iter, // TCI sweeps
              const std::vector<std::vector<int>>& test_points)
{
    std::cout << "\n========================================\n";
    std::cout << "  " << label << "\n";
    std::cout << "  N=" << N << "  d=" << d
              << "  sweeps=" << n_iter << "\n";
    std::cout << "========================================\n";

    std::vector<int>  l_d(N, d);
    MultiIndex pivot0(std::vector<int>(N, 0));
    TCI<Scalar> tci(f_eval, l_d, pivot0);

    for (int it = 0; it < n_iter; ++it) {
        tci.iterate();
        std::cout << "  sweep " << std::setw(3) << (it + 1)
                  << "  |pivot error| = "
                  << std::setprecision(6) << std::scientific
                  << static_cast<double>(tci.last_pivot_error())
                  << "\n";
    }

    auto tt = tci.get_TensorTrain();

    std::cout << "\n  Point-wise check (TT vs exact):\n";
    std::cout << std::setprecision(20) << std::fixed;

    for (const auto& idx : test_points) {
        MultiIndex mi(idx);
        Scalar tt_val  = tt.eval(idx);
        Scalar exact   = f_eval(mi);
        // For complex Scalar use abs; for real scalars std::abs works fine.
        double err = static_cast<double>(
            Eigen::numext::abs(tt_val - exact));

        std::cout << "  idx=";
        for (int v : idx) std::cout << v << " ";
        std::cout << "\n";
        std::cout << "    chi = " << tt.max_bond_dimension()   << "\n";
        std::cout << "    exact = " << exact   << "\n";
        std::cout << "    TT    = " << tt_val  << "\n";
        std::cout << "    |err| = " << std::setprecision(6)
                  << std::scientific << err << "\n";
    }

    std::cout << "\n  sum over all entries (TT): "
              << std::setprecision(20) << std::fixed
              << tt.sum1() << "\n";
}

int main()
{
    // ----------------------------------------------------------
    // Problem parameters
    // ----------------------------------------------------------
    const int N      = 8;   // number of tensor modes
    const int d      = 6;   // physical dimension (indices 0..5)
    const int sweeps = 60;  // TCI sweeps

    // A few test points to check pointwise accuracy
    std::vector<std::vector<int>> test_pts = {
        {0, 1, 2, 3, 4, 5, 0, 1},
        {3, 3, 3, 3, 3, 3, 3, 3},
        {5, 4, 3, 2, 1, 0, 5, 4},
        {1, 0, 4, 2, 5, 3, 0, 2},
    };

    // ----------------------------------------------------------
    // Double precision run
    // ----------------------------------------------------------
    run_demo<double>(
        "DOUBLE PRECISION",
        f_double,
        N, d, sweeps,
        test_pts
    );

    // ----------------------------------------------------------
    // Quad precision run
    // ----------------------------------------------------------
    run_demo<float128>(
        "QUAD PRECISION (float128)",
        f_quad,
        N, d, sweeps,
        test_pts
    );

    // ----------------------------------------------------------
    // cQuad precision run
    // ----------------------------------------------------------
    run_demo<cquad>(
        "Quad complex PRECISION (complex128)",
        f_complex128,
        N, d, sweeps,
        test_pts
    );

    return 0;
}
