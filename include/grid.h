#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
#include "mindex.h"

#include "int128.h"
#include <boost/multiprecision/float128.hpp>

using float128 = boost::multiprecision::float128;

// QTGrid support up to nBit=63. It might fail for 
template <typename Scalar, typename Sint>
class QTGrid {
public:
    using scalar_type = Scalar;
    using index_type = Sint;

    Scalar a;
    Scalar b;

    int nBits;
    Sint N;

    Scalar dx;
    Sint k_offset;
    Scalar x_ref;

    QTGrid(Scalar a_, Scalar b_, int nBits_,
           Sint k_offset_ = 0)
        : a(a_), b(b_), nBits(nBits_), k_offset(k_offset_), x_ref(a_)
    {
        if (nBits < 0 || nBits > 63) {
            throw std::invalid_argument("nBits must be in [0, 63]");
        }

        N = Sint(1) << nBits;
        dx = (b_ - a_) / Scalar(N);
    }

    MultiIndex coord_to_id(Scalar x) const {
        using boost::multiprecision::llround;
        using std::llround;

        Sint k = static_cast<Sint>(llround((x - a) * N/ (b-a)));

        if (k == N) k = N - 1;

        if (!(k >= 0 && k <= N - 1))
        {
            std::cout << k << "\n" << N << "\n";
            throw std::out_of_range("k is outside [0, N-1]");
        }
            
        return _bits(k);
    }

    Scalar id_to_coord(const MultiIndex& bits) const {
        Sint k = 0;
        for (int i = 0; i < nBits; ++i)
            k |= (Sint(bits.data[i]) << i);

        return _coords(k);
    }

private:
    MultiIndex _bits(Sint k) const {
        MultiIndex bits(nBits);

        for (int i = 0; i < nBits; ++i)
            bits.data[i] = (k >> i) & 1;

        return bits;
    }

    Scalar _coords(Sint k) const {
        return x_ref + Scalar(k - k_offset) * dx;
    }
};

using GridDouble = QTGrid<double, long long>;
using GridQuadFast = QTGrid<float128, long long>;
using GridQuad   = QTGrid<float128, util::i128>;