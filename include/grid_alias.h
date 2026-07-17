#include <type_double_double.h>
#include <type_float128_boost.h>
#include <type_int128.h>
#include "grid.h"

using GridDouble = QTGrid<double, long long>;
using GridQuadFast = QTGrid<float128, long long>;
using GridQuad   = QTGrid<float128, util::i128>;
