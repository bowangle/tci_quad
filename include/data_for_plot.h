
#include "error_calc.h"


// class to save data for the plotter
template <typename Scalar>
class DatPlotter{
    TTErrorOnGrid<Scalar> data_error;
    std::vector<Scalar> pivot_error;
};