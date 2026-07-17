#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <type_traits>

#include <nlohmann/json.hpp>

// for llround in coord_to_id
#include <boost/multiprecision/float128.hpp>

#include "mindex.h"

struct dd_128;  // fwd decl for lossy JSON double conversion in save_json

using json = nlohmann::json;

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

    QTGrid(const std::string& filename)
    {
        auto [a_, b_, nBits_] = parse_json(filename);
        *this = QTGrid(a_, b_, nBits_);
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

    void save_json(const std::string& base) const
    {
        json j;
        json j_2;
        std::string suffix;
        std::string suffix_2;

        auto to_double_lossy = [](const Scalar& v) -> double {
            if constexpr (std::is_same_v<Scalar, dd_128>) {
                return v._hi();
            } else {
                return static_cast<double>(v);
            }
        };
        j["a"] = to_double_lossy(a);
        j["b"] = to_double_lossy(b);
        suffix = "_grid.json";
        j["nBits"] = nBits;
        
        std::ostringstream oss_a;
        oss_a << std::setprecision(std::numeric_limits<Scalar>::digits10 + 5) << a;
        j_2["a"] = oss_a.str();

        std::ostringstream oss_b;
        oss_b << std::setprecision(std::numeric_limits<Scalar>::digits10 + 5) << b;
        j_2["b"] = oss_b.str();
        suffix_2 = "_grid_E.json";
        j_2["nBits"] = nBits;

        std::ofstream file(base + suffix);
        if (!file.is_open())
            throw std::runtime_error("Cannot open file");
        std::ofstream file_2(base + suffix_2);
        if (!file_2.is_open())
            throw std::runtime_error("Cannot open file_2");

        file << j.dump(4);
        file_2 << j_2.dump(4);
    }

private:

    static std::tuple<Scalar, Scalar, int> parse_json(const std::string& filename)
    {
        json j;

        std::ifstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Cannot open file");

        file >> j;

        int nBits = j.at("nBits").get<int>();

        Scalar a, b;

        if (j["a"].is_string())
        {
            std::istringstream(j.at("a").get<std::string>()) >> a;
            std::istringstream(j.at("b").get<std::string>()) >> b;
        }
        else
        {
            a = Scalar(j.at("a").get<double>());
            b = Scalar(j.at("b").get<double>());
        }

        return {a, b, nBits};
    }

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