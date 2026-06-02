#include "tensorcross.h"

#include <boost/multiprecision/float128.hpp>

using float128 = boost::multiprecision::float128;
using Complex128 = std::complex<float128>;

template class TCI<float128>;
template class TCI<Complex128>;
template class TCI<double>;