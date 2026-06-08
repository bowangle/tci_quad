#pragma once

#include "grid.h"
#include "tensortrain.h"
#include "tensorcross.h"
#include "mindex.h"
#include "output_tci.h"

#include "tt_sav_load.h"

#include <spdlog/spdlog.h>


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
    int counter = 0 ;
    std::shared_ptr<spdlog::logger> const logger;

    
    TCI_Runner(
        QTGrid<Scalar, Sint> grid_,
        TCI_param tci_param_,
        std::function<Complex(Scalar)> function_x_,
        std::shared_ptr<spdlog::logger> logger_ = nullptr
    )
    : 
        grid(grid_),
        tci_param(tci_param_),
        function_x(function_x_),
        function_id(func_to_grid(function_x_, grid_, logger_)),
        l_d(std::vector<int>(grid_.nBits, 2)),
        logger(logger_),
        counter(0)
    {}

    std::function<Complex(MultiIndex)>
    func_to_grid(
        std::function<Complex(Scalar)> f,
        const QTGrid<Scalar, Sint>& qt_grid,
        std::shared_ptr<spdlog::logger> logger_ = nullptr
    )
    {
        if (logger_){
                return [this, &qt_grid, f](const MultiIndex& id) -> Complex {
                    counter++;
                    return f(qt_grid.id_to_coord(id));
                };
            }
        else{
            return [&qt_grid, f](const MultiIndex& id) -> Complex {
                    return f(qt_grid.id_to_coord(id));
                };
        }
    }

    template <typename T>
    std::string to_string_stream(const T& v) {
        std::ostringstream oss;
        oss << v;
        return oss.str();
    }

    void fit(
        const Scalar E_init,
        const std::vector<Scalar> other_E,
        bool verbose = true,
        bool do_save=false,
        const std::string& file_prefix="",
        int nb_point_res=1000
    ){
        auto log = [&](const std::string& msg) {
            if (logger) {
                logger->info(msg);
            } 
            if (verbose) {
                std::cout << msg << "\n";
            }
        };

        log("========== FIT PARAMETERS ==========");

        log("E_init = " + to_string_stream(E_init));

        std::string other_E_str = "[";
        for (size_t i = 0; i < other_E.size(); ++i) {
            other_E_str += to_string_stream(other_E[i]);
            if (i + 1 < other_E.size()) other_E_str += ", ";
        }
        other_E_str += "]";

        log("other_E = " + other_E_str);

        log("verbose = " + to_string_stream(verbose ? "true" : "false"));
        log("do_save = " + to_string_stream(do_save ? "true" : "false"));

        log("file_prefix = " + file_prefix);

        log("nb_point_res = " + to_string_stream(nb_point_res));

        log("============== TCI params ==============");
        log("  N = " + std::to_string(grid.nBits));
        log("  d = 2");
        log("  sweeps = " + std::to_string(tci_param.nb_iter));
        log("========================================");

        MultiIndex pivot0 = grid.coord_to_id(E_init);

        TCI<Complex> tci(function_id, l_d, pivot0);
        
        int nb_additional_pivot = other_E.size();

        for (int i=0; i<nb_additional_pivot; i++){
            MultiIndex pivot_temp = grid.coord_to_id(E_init);
            //tci.addpivot_all_bound(pivot_temp);
        }

        for (int it = 0; it < tci_param.nb_iter; ++it) {
            tci.iterate();


            std::ostringstream oss_pivot;

            oss_pivot << std::scientific
                    << std::setprecision(std::numeric_limits<Scalar>::max_digits10);

            oss_pivot << "sweep " << (it + 1)
                    << " |pivot error| = "
                    << tci.last_pivot_error();

            log(oss_pivot.str());
        }
        TensorTrain<Complex> tt_temp = tci.get_TensorTrain();
        std::vector<Cube<Complex>> cube_cores = tt_temp.convert_to_cubes();

        save_vec_cube(cube_cores, file_prefix + ".tt");

        log("Max bond dim: " + std::to_string(tt_temp.max_bond_dimension()));
        
        log("Nb function call:" + std::to_string(counter));

        TTErrorOnGrid<Scalar> error = 
            error_TT_on_grid_point(tt_temp, function_id, grid, nb_point_res);
        
        std::ostringstream oss;
        oss << error;
        log(oss.str());

        std::vector<Scalar> pivot_error = tci.pivotErrors();

        // save the tt and the grid
        if (do_save){
            save_TTErrorOnGrid(error, file_prefix + "_error.dat");
            grid.save_json(file_prefix);
        }
    }


};
