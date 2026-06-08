#pragma once

#include "grid.h"
#include "tensortrain.h"
#include "tensorcross.h"
#include "mindex.h"
#include "output_tci.h"

#include "tt_sav_load.h"


class TCI_param{
    public:

    int nBit;
    int nb_iter;

    TCI_param(int nBit_, int nb_iter_)
    :   nBit(nBit_),
        nb_iter(nb_iter_)
    {}
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

    void fit(
        const Scalar E_init,
        const std::vector<Scalar> other_E,
        bool verbose = true,
        bool do_save=false,
        const std::string& file_prefix="",
        int nb_point_res=1000
    ){
        if (verbose){
            std::cout << "\n========================================\n";
            std::cout << "  N=" << grid.nBits << "  d=" << 2
                    << "  sweeps=" << tci_param.nb_iter << "\n";
            std::cout << "========================================\n";
        }
        MultiIndex pivot0 = grid.coord_to_id(E_init);

        TCI<Complex> tci(function_id, l_d, pivot0);
        
        int nb_additional_pivot = other_E.size();

        for (int i=0; i<nb_additional_pivot; i++){
            MultiIndex pivot_temp = grid.coord_to_id(E_init);
            //tci.addpivot_all_bound(pivot_temp);
        }

        for (int it = 0; it < tci_param.nb_iter; ++it) {
            tci.iterate();
            if (verbose){
                std::cout << "  sweep " << std::setw(3) << (it + 1)
                    << "  |pivot error| = "
                    << std::setprecision(6) << std::scientific
                    << static_cast<double>(tci.last_pivot_error())
                    << "\n";
            }
        }
        TensorTrain<Complex> tt_temp = tci.get_TensorTrain();
        std::vector<Cube<Complex>> cube_cores = tt_temp.convert_to_cubes();

        save_vec_cube(cube_cores, file_prefix + ".tt");

        if (verbose){
            std::cout << "Max bond dim:" <<tt_temp.max_bond_dimension() << "\n";
        }
        
        TTErrorOnGrid<Scalar> error = 
            error_TT_on_grid_point(tt_temp, function_id, grid, nb_point_res);
        
        if (verbose){
            std::cout << error << "\n";
        }
        std::vector<Scalar> pivot_error = tci.pivotErrors();

        // save the tt and the grid
        if (do_save){
            save_TTErrorOnGrid(error, file_prefix + "_error.dat");
            grid.save_json(file_prefix);
        }
    }


};
