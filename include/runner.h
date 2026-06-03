#pragma

#include "grid.h"
#include "tensortrain.h"
#include "tensorcross.h"
#include "mindex.h"
#include "error_calc.h"


class TCI_param{
    public:

    int nBit;
    int nb_iter;
};

// f_type can be a scalar type or a complex type:
template <typename Scalar, typename Sint>
class TCI_Runner{
    public:
    
    using Complex = std::complex<Scalar>;

    const QTGrid<Scalar, Sint> grid;
    const TCI_param tci_param;
    std::function<Complex(Scalar)> function_x;  // original function
    std::function<Complex(MultiIndex)> function_id; // function on grid
    const std::vector<int> l_d;

    
    TCI_Runner(
        QTGrid<Scalar, Sint> grid_,
        TCI_param tci_param_,
        std::function<Complex(Scalar)> function_x_
    )
    : 
        grid(grid_),
        tci_param(tci_param_),
        function_x(function_x_),
        function_id(func_to_grid(function_x_, grid_)),
        l_d(std::vector<int>(grid_.nBits, 2))
    {}

    std::function<Complex(MultiIndex)>
    func_to_grid(
        std::function<Complex(Scalar)> f,
        const QTGrid<Scalar, Sint>& qt_grid)
    {
        return [&qt_grid, f](const MultiIndex& id) -> Complex {
            return f(qt_grid.id_to_coord(id));
        };
    }

    void fit(bool verbose = true){
        if (verbose){
            std::cout << "\n========================================\n";
            std::cout << "  N=" << grid.nBits << "  d=" << 2
                    << "  sweeps=" << tci_param.nb_iter << "\n";
            std::cout << "========================================\n";
        }
        MultiIndex pivot0(std::vector<int>(grid.nBits, 0));
        TCI<Complex> tci(function_id, l_d, pivot0);
        
        for (int it = 0; it < tci_param.nb_iter; ++it) {
            tci.iterate();
            std::cout << "  sweep " << std::setw(3) << (it + 1)
                  << "  |pivot error| = "
                  << std::setprecision(6) << std::scientific
                  << static_cast<double>(tci.last_pivot_error())
                  << "\n";
            if (it%10==0){
                auto tt_temp = tci.get_TensorTrain();
                std::cout << "Max bond dim:" <<tt_temp.max_bond_dimension() << "\n";
                TTErrorOnGrid<Scalar> error = 
                    error_TT_on_grid_point(tt_temp, function_id, grid, 1000);
                std::cout << error << "\n";
            }
        }
        auto tt_temp = tci.get_TensorTrain();
        std::cout << "Max bond dim:" <<tt_temp.max_bond_dimension() << "\n";
        TTErrorOnGrid<Scalar> error = 
            error_TT_on_grid_point(tt_temp, function_id, grid, 1000);
        std::cout << error << "\n";

    }


};
