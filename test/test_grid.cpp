#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include "grid.h"
#include "int128.h"

#include <boost/multiprecision/float128.hpp>
using float128 = boost::multiprecision::float128;

template <typename Grid>
void test_one_grid_roudtrip(const char* name, const int nBit) {
    std::cout << "[TEST] " << name << "nBit" << nBit << "\n";

    using Scalar = typename Grid::scalar_type;

    Grid g(-5, 10, nBit);

    Scalar x = Scalar(1.25);

    std::cout << "a="<< g.a <<"\n";
    std::cout << "b="<< g.b <<"\n";
    std::cout << "N="<< g.N <<"\n";
    std::cout << "nBits="<< g.nBits <<"\n";

    auto bits = g.coord_to_id(x);
    Scalar x2 = g.id_to_coord(bits);
    // x2 is x mapped on grid

    auto bits_2 = g.coord_to_id(x2);
    Scalar x3 = g.id_to_coord(bits_2);

    std::cout << "x = " << x <<"; x2 = " << x2 << " -> x3 = " << x3 << "\n";

    if constexpr (std::is_same_v<Scalar, double>) {
        assert(std::abs(double(x2 - x3)) < 1e-6);
    }
}

void test_grid_roundtrip(){
    test_one_grid_roudtrip<GridDouble>("double grid", 8);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 8);

    test_one_grid_roudtrip<GridDouble>("double grid", 30);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 30);

    test_one_grid_roudtrip<GridDouble>("double grid", 48);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 48);

    test_one_grid_roudtrip<GridDouble>("double grid", 53);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 53);

    test_one_grid_roudtrip<GridQuad>("float128 grid", 110);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 113);
}

template <typename Scalar, typename Sint>
void testgrid(){
    Scalar a_1 = Scalar(-5.);
    Scalar b_1 = Scalar(10);
    Sint nBit_1 = Sint(15);

    QTGrid<Scalar, Sint> grid (a_1, b_1, nBit_1);
}


int main() {
    
    
    testgrid<float128, util::i128>();
    testgrid<double, long long>();

    test_grid_roundtrip();
}

