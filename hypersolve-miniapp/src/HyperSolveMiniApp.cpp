#include "Element.h"
#include "HyperSolveMiniApp.h"
#include "Ddata.h"
#include "TetMesh.h"

namespace HS {
    template<typename T, size_t NumEqns>
    using Solution = std::vector<std::array<T, NumEqns>>;
    template<typename T, size_t NumEqns>
    using Residual = std::vector<std::array<T, NumEqns>>;
}

template <typename ResidualType>
void zeroResidual(ResidualType& residual) {
    for (auto& res : residual)
        for (auto& R : res)
            R = 0.0;
}

void monitorCellLoop(size_t cell_id, size_t max_cells) {
    cell_id += 1;
    if (cell_id % (max_cells/10) == 0) {
        auto percent_done = cell_id * 100 / max_cells;
        printf("[%3lu%%]  looping cell: %lu of %lu\n", percent_done, cell_id, max_cells);
    }
}

template<typename T>
void LoopCellFlux(HS::Residual<T, 5>& residual, size_t cell_id,
                  const HS::Solution<double, 5> &solution, const HS::TetMesh &mesh) {
    std::array<int, 4> cell_node_ids;
    mesh.getCellNodeIds(cell_node_ids.data(), cell_id);

    std::array<std::array<T, 5>, 4> q;
    for (int corner = 0; corner < 4; ++corner)
        for (int eqn = 0; eqn < 5; ++eqn)
            q[corner][eqn] = solution[cell_node_ids[corner]][eqn];

    std::array<HS::Point<double>, 4> node_xyz;
    for (int corner = 0; corner < 4; ++corner) {
        mesh.getXYZ(node_xyz[corner].data(), cell_node_ids[corner]);
    }

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

    auto flux = HS::ElementBasedViscousFlux<4, 5>(HS::Element::Tet::edge_to_node,
                                                  edge_normals,
                                                  edge_ugrad, edge_vgrad, edge_wgrad, edge_tgrad,
                                                  thermal_conductivity_avg,
                                                  viscosity_avg,
                                                  uvw_viscosity_avg);
    for (int corner_node = 0; corner_node < 4; ++corner_node) {
        int corner_node_id = cell_node_ids[corner_node];
        for (int eqn = 0; eqn < 5; ++eqn) {
            residual[corner_node_id][eqn] += flux[corner_node][eqn];
        }
    }
}

int main(int argc, char *argv[]) {

    if (argc == 1) {
        printf("Usage: %s <number of cells to loop>\n", argv[0]);
        return 1;
    }
    std::string arg = argv[1];
    size_t max_cells;
    try {
        std::size_t pos;
        max_cells = size_t(std::stod(arg, &pos));
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

    auto tet_mesh = HS::TetMesh(max_cells);

    std::array<double, 5> q0 = {2.0, 0.25, 0.1, 0.1, 0.767};
    HS::Solution<double, 5> solution(tet_mesh.getNodeCount(), q0);

    typedef Linearize::Ddata<double, 5> ddt;
    HS::Residual<double, 5> real_residual(tet_mesh.getNodeCount());
    HS::Residual<ddt, 5> ddt_residual(tet_mesh.getNodeCount());

    printf("RHS:\n");
    zeroResidual(real_residual);
    for (size_t cell_id = 0; cell_id < tet_mesh.getCellCount(); ++cell_id) {
        LoopCellFlux<double>(real_residual, cell_id, solution, tet_mesh);
        monitorCellLoop(cell_id, max_cells);
    }
    printf("LHS:\n");
    zeroResidual(ddt_residual);
    for (size_t cell_id = 0; cell_id < tet_mesh.getCellCount(); ++cell_id) {
        LoopCellFlux<ddt>(ddt_residual, cell_id, solution, tet_mesh);
        monitorCellLoop(cell_id, max_cells);
    }

    std::array<std::array<double, 5>, 4> expected_residual_at_cell;
    expected_residual_at_cell[0] = {0.0,  6.666666666667e-03,  6.666666666667e-03,  6.666666666667e-03,  2.250000000000e+00};
    expected_residual_at_cell[1] = {0.0,  0.000000000000e+00, -3.333333333333e-03, -3.333333333333e-03, -7.833333333333e-01};
    expected_residual_at_cell[2] = {0.0, -3.333333333333e-03,  0.000000000000e+00, -3.333333333333e-03, -7.500000000000e-01};
    expected_residual_at_cell[3] = {0.0, -3.333333333333e-03, -3.333333333333e-03,  0.000000000000e+00, -7.166666666667e-01};

    double tolerance = 1.e-12;
    for (size_t cell_id = 0; cell_id < tet_mesh.getCellCount(); ++cell_id) {
        std::array<int, 4> cell_node_ids;
        tet_mesh.getCellNodeIds(cell_node_ids.data(), cell_id);
        for (int corner = 0; corner < 4; ++corner) {
            int node_id = cell_node_ids[corner];
            for (int eqn = 0; eqn < 5; ++eqn) {
                double real_difference = std::fabs(real_residual[node_id][eqn] - expected_residual_at_cell[corner][eqn]);
                if (real_difference > tolerance) {
                    printf("node: %d eqn: %d value: %e expected: %e\n", node_id, eqn, real_residual[node_id][eqn],
                           expected_residual_at_cell[corner][eqn]);
                    printf("Error: RHS difference: %e greater than tolerance: %e\n", real_difference, tolerance);
                    return 1;
                }
            }
        }
    }

    for (size_t node_id = 0; node_id < tet_mesh.getNodeCount(); ++node_id) {
        for (int eqn = 0; eqn < 5; ++eqn) {
            double difference = std::fabs(real_residual[node_id][eqn] - ddt_residual[node_id][eqn].value);
            if (difference > tolerance) {
                printf("node: %lu eqn: %d RHS: %e LHS: %e\n", node_id, eqn, real_residual[node_id][eqn],
                       ddt_residual[node_id][eqn].value);
                printf("Error: RHS and LHS don't match!\n");
                return 1;
            }
        }
    }

    return 0;
}