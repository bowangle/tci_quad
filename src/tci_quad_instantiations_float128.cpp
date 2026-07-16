#include "tensorcross.h"

#include <boost/multiprecision/float128.hpp>
#include "grid.h"
#include "type_int128.h"
#include "runner.h"

using float128 = boost::multiprecision::float128;
using Complex128 = std::complex<float128>;
using GridDouble = QTGrid<double, long long>;
using GridQuadFast = QTGrid<float128, long long>;
using GridQuad   = QTGrid<float128, util::i128>;

template class TCI<float128>;
template class TCI<Complex128>;
template class TCI<double>;

template class TCI_Runner<float128, long long>;

template struct TTErrorOnGrid<float128>;
