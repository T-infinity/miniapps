#pragma once

#include <vector>
#include <array>
#include "Point.h"

namespace HS {
    class TetMesh {
    public:
        explicit TetMesh(size_t number_of_cells_in)
                : number_of_cells(number_of_cells_in),
                  xyz_locations(buildXYZ(number_of_cells_in)),
                  cell_node_ids(buildCellNodeIds(number_of_cells_in)){}

        size_t getCellCount() const {
            return number_of_cells;
        }

        size_t getNodeCount() const {
            return number_of_cells * 4;
        }

        void getCellNodeIds(int* node_ids, int cell_id) const {
            for (int corner = 0; corner < 4; ++corner) {
                node_ids[corner] = cell_node_ids[cell_id][corner];
            }
        }

        void getXYZ(double* coordinates, int node_id) const {
            for (int i = 0; i < 3; ++i) {
                coordinates[i] = xyz_locations[node_id][i];
            }
        }

    private:
        const size_t number_of_cells;
        std::vector<Point<double>> xyz_locations;
        std::vector<std::array<int, 4>> cell_node_ids;

        std::vector<Point<double>> buildXYZ(size_t number_of_cells) const {
            std::vector<Point<double>> xyz(number_of_cells * 4);
            int i = 0;
            for (size_t cell_id = 0; cell_id < number_of_cells; ++cell_id) {
                xyz[i++] = {0.0, 0.0, 0.0};
                xyz[i++] = {1.0, 0.0, 0.0};
                xyz[i++] = {0.0, 1.0, 0.0};
                xyz[i++] = {0.0, 0.0, 1.0};
            }
            return xyz;
        }

        std::vector<std::array<int, 4>> buildCellNodeIds(size_t number_of_cells) const {
            std::vector<std::array<int, 4>> cell_node_ids(number_of_cells);
            int i = 0;
            for (size_t cell_id = 0; cell_id < number_of_cells; ++cell_id) {
                cell_node_ids[cell_id][0] = i++;
                cell_node_ids[cell_id][1] = i++;
                cell_node_ids[cell_id][2] = i++;
                cell_node_ids[cell_id][4] = i++;
            }
            return cell_node_ids;
        }
    };
}
