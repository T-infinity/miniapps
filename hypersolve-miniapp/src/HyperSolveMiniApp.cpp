#include "Element.h"
#include "HyperSolveMiniApp.h"
#include "Ddata.h"

template<typename T, size_t NumEqns, size_t NCorners>
using CornerSolution = std::array<std::array<T, NumEqns>, NCorners>;

template<typename T>
std::array<std::array<T, 5>, 4> LoopCellFlux(long max_cells) {
    std::array<double, 5> q0 = {2.0, 0.25, 0.1, 0.1, 0.767};
    std::array<std::array<T, 5>, 4> q;
    for (int corner = 0; corner < 4; ++corner)
        for (int eqn = 0; eqn < 5; ++eqn)
            q[corner][eqn] = q0[eqn];

    std::array<HS::Point<double>, 4> node_xyz;
    node_xyz[0] = {0.0, 0.0, 0.0};
    node_xyz[1] = {1.0, 0.0, 0.0};
    node_xyz[2] = {0.0, 1.0, 0.0};
    node_xyz[3] = {0.0, 0.0, 1.0};

    auto edge_normals = HS::Element::Tet::calcDualNormals(node_xyz);

    auto one = T(1.0);
    HS::Point<T> mock_grad(one, one, one);
    HS::Point<T> ugrad = mock_grad;
    HS::Point<T> vgrad = mock_grad;
    HS::Point<T> wgrad = mock_grad;
    HS::Point<T> tgrad = mock_grad;

    std::array<HS::Point<T>, 6> edge_ugrad, edge_vgrad, edge_wgrad, edge_tgrad;
    edge_ugrad.fill(ugrad);
    edge_vgrad.fill(vgrad);
    edge_wgrad.fill(wgrad);
    edge_tgrad.fill(tgrad);

    // Hardcoded to avoid thermodynamics dependency
    T viscosity_avg(0.01);
    T thermal_conductivity_avg(0.1);
    HS::Point<T> uvw_viscosity_avg;
    uvw_viscosity_avg[0] = 1.0;
    uvw_viscosity_avg[1] = 1.1;
    uvw_viscosity_avg[2] = 1.2;

    std::array<std::array<T, 5>, 4> flux;
    int index = -1;
    int percent_done = 0;
    for (long i = 0; i < max_cells; ++i) {
        index++;
        if (index == (max_cells / 10)) {
            percent_done += 10;
            printf("[%d%%]  looping cell: %lu of %lu\n", percent_done, i, max_cells);
            index = 0;
        }
        flux = HS::ElementBasedViscousFlux<4, 5>(HS::Element::Tet::edge_to_node,
                                                 edge_normals,
                                                 edge_ugrad, edge_vgrad, edge_wgrad, edge_tgrad,
                                                 thermal_conductivity_avg,
                                                 viscosity_avg,
                                                 uvw_viscosity_avg);
    }
    printf("[%d%%] looping cell: %lu of %lu\n", 100, max_cells, max_cells);
    return flux;
}

int main(int argc, char *argv[]) {

    if (argc == 1) {
        printf("Usage: %s <number of cells to loop>\n", argv[0]);
        return 1;
    }
    std::string arg = argv[1];
    long max_cells;
    try {
        std::size_t pos;
        max_cells = long(std::stod(arg, &pos));
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << '\n';
            return 1;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
        return 1;
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
        return 1;
    }

    printf("Loop over %lu cells.\n", max_cells);

    printf("RHS:\n");
    auto real_flux = LoopCellFlux<double>(max_cells);
    printf("LHS:\n");
    auto ddt_flux = LoopCellFlux<Linearize::Ddata<double, 5>>(max_cells);
    for (int corner = 0; corner < 4; ++corner) {
        for (int eqn = 0; eqn < 5; ++eqn) {
            double difference = std::fabs(real_flux[corner][eqn] - ddt_flux[corner][eqn].value);
            if (difference > 1.e-15) return 1;
        }
    }
    return 0;
}