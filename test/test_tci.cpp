#include <functional>
#include <complex>

#include "runner.h"
#include "grid.h"

#include <boost/multiprecision/float128.hpp>
using float128 = boost::multiprecision::float128;


template <typename Scalar>
std::function<std::complex<Scalar>(Scalar)> make_function_sin()
{
    using boost::multiprecision::sin;

    return [](Scalar E) -> std::complex<Scalar> {
        return std::complex<Scalar>(sin(E)+sin(Scalar("0.5")*E), sin(E)*sin(Scalar("0.01")*E));
    };
}


template <typename Scalar>
void TCI_sin(int nBit, int n_iter, bool do_save=false, const std::string& filename="", int nb_point_out=1000){
    using Complex = std::complex<Scalar>;
    bool is_lesser = true;
    int id_lead_1 = 0;
    int id_lead_2 = 1;

    Scalar a_1 = Scalar("-5");
    Scalar b_1 = Scalar("10");
    QTGrid<Scalar, long long> grid(a_1, b_1, nBit);


    TCI_param tci_param = TCI_param(grid.nBits, n_iter);

    std::function<std::complex<Scalar>(Scalar)> f_sin = make_function_sin<Scalar>();
    TCI_Runner<Scalar, long long> tci_runner(grid, tci_param, f_sin);

    std::vector<Scalar> v = {Scalar("2.05"), Scalar("2.3")};

    tci_runner.fit(Scalar("2.05"), v, true, do_save, filename, nb_point_out);
}

int main() {
    int nBit = 15;
    int n_iter = 30;
    bool do_save = true;
    std::string file_prefix = "test/test_out";
    int nb_point_out = 1000;

    TCI_sin<float128>(nBit, n_iter, do_save, file_prefix, nb_point_out);
    
    
}