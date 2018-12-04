#pragma once

#include "array"

namespace HS {
    template<typename T, size_t NumEdges>
    class EdgeGradients {
    public:
        explicit EdgeGradients(const std::array<T, NumEdges> mock_gradient)
        : ux(mock_gradient), uy(mock_gradient), uz(mock_gradient),
          vx(mock_gradient), vy(mock_gradient), vz(mock_gradient),
          wx(mock_gradient), wy(mock_gradient), wz(mock_gradient),
          tx(mock_gradient), ty(mock_gradient), tz(mock_gradient) {}

        const std::array<T, NumEdges> ux;
        const std::array<T, NumEdges> uy;
        const std::array<T, NumEdges> uz;

        const std::array<T, NumEdges> vx;
        const std::array<T, NumEdges> vy;
        const std::array<T, NumEdges> vz;

        const std::array<T, NumEdges> wx;
        const std::array<T, NumEdges> wy;
        const std::array<T, NumEdges> wz;

        const std::array<T, NumEdges> tx;
        const std::array<T, NumEdges> ty;
        const std::array<T, NumEdges> tz;
    };
}


