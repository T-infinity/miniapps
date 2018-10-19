#pragma once

#include <array>
#include <vector>
#include "Point.h"

namespace HS {
    namespace Element {
        template<size_t n_corners, size_t n_edges, size_t n_faces>
        std::array<Point<double>, n_edges>
        calcCellDualNormals(const std::array<Point<double>, n_corners> &node_coordinates,
                            const std::array<std::array<int, 2>, n_edges> &edge_to_face,
                            const Point<double> &cell_center,
                            const std::array<Point<double>, n_edges> &edge_centers,
                            const std::array<Point<double>, n_faces> &face_centers) {

            std::array<Point<double>, n_edges> edge_areas;
            for (size_t edge = 0; edge < n_edges; ++edge) {
                auto cell_center_to_edge_midpoint = cell_center - edge_centers[edge];

                int left_face = edge_to_face[edge][0];
                auto left_face_to_edge_midpoint = face_centers[left_face] - edge_centers[edge];
                auto left_area = 0.5 * Point<double>::cross(left_face_to_edge_midpoint, cell_center_to_edge_midpoint);

                int right_face = edge_to_face[edge][1];
                auto right_face_to_edge_midpoint = face_centers[right_face] - edge_centers[edge];
                auto right_area = 0.5 * Point<double>::cross(cell_center_to_edge_midpoint, right_face_to_edge_midpoint);

                edge_areas[edge] = left_area + right_area;
            }
            return edge_areas;
        };

        namespace Tet {
            static constexpr std::array<std::array<int, 2>, 6> edge_to_node = {{
                                                                                       {0, 1},
                                                                                       {1, 2},
                                                                                       {2, 0},
                                                                                       {0, 3},
                                                                                       {1, 3},
                                                                                       {2, 3}
                                                                               }};

            static constexpr std::array<std::array<int, 2>, 6> edge_to_face = {{
                                                                                       {0, 1},
                                                                                       {0, 2},
                                                                                       {0, 3},
                                                                                       {1, 3},
                                                                                       {2, 1},
                                                                                       {3, 2}
                                                                               }};
            inline std::array<Point<double>, 6>
            calcDualNormals(std::array<Point<double>, 4> node_coordinates) {
                Point<double> cell_center = 0.25 * (node_coordinates[0] + node_coordinates[1] +
                                                    node_coordinates[2] + node_coordinates[3]);
                std::array<Point<double>, 6> edge_centers{
                        0.5 * (node_coordinates[0] + node_coordinates[1]),
                        0.5 * (node_coordinates[1] + node_coordinates[2]),
                        0.5 * (node_coordinates[2] + node_coordinates[0]),
                        0.5 * (node_coordinates[0] + node_coordinates[3]),
                        0.5 * (node_coordinates[1] + node_coordinates[3]),
                        0.5 * (node_coordinates[2] + node_coordinates[3])
                };
                std::array<Point<double>, 4> face_centers{
                        (node_coordinates[0] + node_coordinates[1] + node_coordinates[2]) / 3.0,
                        (node_coordinates[0] + node_coordinates[1] + node_coordinates[3]) / 3.0,
                        (node_coordinates[3] + node_coordinates[1] + node_coordinates[2]) / 3.0,
                        (node_coordinates[0] + node_coordinates[2] + node_coordinates[3]) / 3.0
                };

                return calcCellDualNormals(node_coordinates, edge_to_face, cell_center, edge_centers, face_centers);
            }
        }
    }
}
