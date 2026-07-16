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
    test_one_grid_roudtrip<GridQuadFast>("float128 fast grid", 8);

    test_one_grid_roudtrip<GridDouble>("double grid", 30);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 30);
    test_one_grid_roudtrip<GridQuadFast>("float128 fast grid", 30);

    test_one_grid_roudtrip<GridDouble>("double grid", 40);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 40);
    test_one_grid_roudtrip<GridQuadFast>("float128 fast grid", 40);

    test_one_grid_roudtrip<GridDouble>("double grid", 48);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 48);
    test_one_grid_roudtrip<GridQuadFast>("float128 fast grid", 48);

    test_one_grid_roudtrip<GridDouble>("double grid", 53);
    test_one_grid_roudtrip<GridQuad>("float128 grid", 53);
    test_one_grid_roudtrip<GridQuadFast>("float128 fast grid", 53);

    test_one_grid_roudtrip<GridQuad>("float128 grid", 63);
}

template <typename Scalar, typename Sint>
void testgrid(){
    Scalar a_1 = Scalar(-5.);
    Scalar b_1 = Scalar(10);
    Sint nBit_1 = Sint(15);

    QTGrid<Scalar, Sint> grid (a_1, b_1, nBit_1);
}

template <typename Scalar, typename Sint>
void test_save_load_roundtrip()
{
    Scalar a_1 = Scalar(-5.0);
    Scalar b_1 = Scalar(10.0);
    Sint nBit_1 = Sint(15);

    QTGrid<Scalar, Sint> grid(a_1, b_1, nBit_1);

    // ---------------- SAVE ----------------
    std::string filename = "grid_test";
    grid.save_json(filename);

    // ---------------- LOAD ----------------
    QTGrid<Scalar, Sint> grid2(filename + "_grid_E.json");

    // ---------------- CHECKS ----------------

    auto almost_equal = [](Scalar x, Scalar y) {
        using std::abs;
        using boost::multiprecision::abs;
        return abs(x - y) < Scalar(1e-12);
    };

    std::cout << "a: " << grid.a << " vs " << grid2.a << "\n";
    std::cout << "b: " << grid.b << " vs " << grid2.b << "\n";
    std::cout << "nBits: " << grid.nBits << " vs " << grid2.nBits << "\n";
    std::cout << "k_offset: " << grid.k_offset << " vs " << grid2.k_offset << "\n";
    std::cout << "x_ref: " << grid.x_ref << " vs " << grid2.x_ref << "\n";

    assert(almost_equal(grid.a, grid2.a));
    assert(almost_equal(grid.b, grid2.b));
    assert(grid.nBits == grid2.nBits);
    assert(grid.k_offset == grid2.k_offset);
    assert(almost_equal(grid.x_ref, grid2.x_ref));

    // derived consistency
    assert(grid2.N == (Sint(1) << grid2.nBits));

    Scalar dx1 = (grid.b - grid.a) / Scalar(grid.N);
    Scalar dx2 = (grid2.b - grid2.a) / Scalar(grid2.N);

    assert(almost_equal(dx1, dx2));

    std::cout << "QTGrid roundtrip test PASSED\n";
}

int main() {
    
    
    testgrid<float128, util::i128>();
    testgrid<double, long long>();

    test_save_load_roundtrip<double, long long>();
    test_save_load_roundtrip<float128, util::i128>();

    test_grid_roundtrip();
}

