#pragma once

#include <complex>
#include <fstream>
#include <Eigen/Dense>

template<class T>
struct Cube {
    size_t n_rows, n_cols, n_slices;
    std::vector<T> data; // in this format, we want the ordering [left, physical, right]

    T& operator()(size_t i, size_t j, size_t k) {
        return data[k*n_rows*n_cols + i*n_cols + j];
    }

    const T& operator()(size_t i, size_t j, size_t k) const {
        return data[k*n_rows*n_cols + i*n_cols + j];
    }
};

template<typename cT>
void save_arma_cube_vector(std::ostream& out,
                            const std::vector<Cube<cT>>& Xs)
{
    out << Xs.size() << "\n";

    for (size_t t = 0; t < Xs.size(); ++t)
    {
        out << "ARMA_CUB_TXT_FC016\n";

        const auto& X = Xs[t];

        out << X.n_rows << " "
            << X.n_cols << " "
            << X.n_slices << "\n";

        out << std::scientific;
        out << std::setprecision(std::numeric_limits<typename Eigen::NumTraits<cT>::Real>::digits10 + 5);

        for (size_t k = 0; k < X.n_slices; ++k)
        {
            for (size_t i = 0; i < X.n_rows; ++i)
            {
                for (size_t j = 0; j < X.n_cols; ++j)
                {
                    const auto& z = X(i, j, k);

                    out << "("
                        << z.real() << ","
                        << z.imag()
                        << ") ";
                }
            }
            out << "\n";
        }
        out << "\n";
    }
}

template<typename cT>
void save_vec_cube(std::vector<Cube<cT>> cores, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file)
        throw std::runtime_error("Cannot open file");
    save_arma_cube_vector<cT>(file, cores);
}