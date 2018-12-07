
// Copyright 2018 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.
// No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
// INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE
// SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT
// SOFTWARE.  THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
// RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM
// USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING
// THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS." 
//
// RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
// AS WELL AS ANY PRIOR RECIPIENT. IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES,
// DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR
// RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE
// UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
// PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION
// OF THIS AGREEMENT.


#pragma once

#include "Point.h"
#include "EdgeGradients.h"
#include "Element.h"
#include <Eigen/Dense>

namespace HS {
    template<typename T, size_t NumEdges>
    inline Eigen::Matrix<T, 3, 3>
    calcSij(const EdgeGradients<T, NumEdges>& grad, int edge) {
        Eigen::Matrix<T, 3, 3> Sij;

        // Need for operator(*) with Ddata
        T half(0.5);

        Sij(0, 0) = grad.ux[edge];
        Sij(0, 1) = half * (grad.uy[edge] + grad.vx[edge]);
        Sij(0, 2) = half * (grad.uz[edge] + grad.wx[edge]);

        Sij(1, 0) = Sij(0, 1);
        Sij(1, 1) = grad.vy[edge];
        Sij(1, 2) = half * (grad.vz[edge] + grad.wy[edge]);

        Sij(2, 0) = Sij(0, 2);
        Sij(2, 1) = Sij(1, 2);
        Sij(2, 2) = grad.wz[edge];

        return Sij;
    }

    template<typename T>
    inline Eigen::Matrix<T, 3, 3>
    calcSij(const Eigen::Matrix<T, 3, 1> &grad_u,
            const Eigen::Matrix<T, 3, 1> &grad_v,
            const Eigen::Matrix<T, 3, 1> &grad_w) {
        Eigen::Matrix<T, 3, 3> Sij;

        // Need for operator(*) with Ddata
        T half(0.5);

        Sij(0, 0) = grad_u[0];
        Sij(0, 1) = half * (grad_u[1] + grad_v[0]);
        Sij(0, 2) = half * (grad_u[2] + grad_w[0]);

        Sij(1, 0) = Sij(0, 1);
        Sij(1, 1) = grad_v[1];
        Sij(1, 2) = half * (grad_v[2] + grad_w[1]);

        Sij(2, 0) = Sij(0, 2);
        Sij(2, 1) = Sij(1, 2);
        Sij(2, 2) = grad_w[2];

        return Sij;
    }

    template<typename T>
    Eigen::Matrix<T, 3, 3>
    calcSijBar(const Eigen::Matrix<T, 3, 3>& Sij) {
        T trace_div_3 = Sij.trace() / 3.0;
        Eigen::Matrix<T, 3, 3> SijBar = Sij;
        for (int i = 0; i < 3; ++i)
            SijBar(i,i) -= trace_div_3;
        return SijBar;
    };

    template<size_t n_corners, size_t NumEqns, size_t n_edges, typename T>
    std::array<std::array<T, NumEqns>, n_corners>
    ElementBasedViscousFlux(const std::array<std::array<int, 2>, n_edges>& edge_to_node,
                            const Element::EdgeNormals<n_edges> &edge_normals,
                            const EdgeGradients<T, n_edges> &edge_gradients,
                            const T& thermal_conductivity_avg, // pass directly because prandtl number might be different for turbulence
                            const T& viscosity_avg,
                            const Point<T>& uvw_viscosity_avg) {

        using namespace Eigen;

        Matrix<T, 3, 1> tgrad;
        Matrix<T, 3, 1> uvw_viscosity, normal;
        Matrix<T, 3, 1> heat_flux, momentum_fluxes;
        Matrix<T, 3, 1> energy_stress, energy_eqn_terms;
        Matrix<T, 3, 3> Sij, Sij_bar, tau;

        std::array<std::array<T, NumEqns>, n_corners> flux;
        for (auto &f : flux)
            for (auto &value : f)
                value = 0.0;

        uvw_viscosity << uvw_viscosity_avg[0], uvw_viscosity_avg[1], uvw_viscosity_avg[2];

        for (size_t edge = 0; edge < n_edges; ++edge) {

            tgrad << edge_gradients.tx[edge], edge_gradients.ty[edge], edge_gradients.tz[edge];

            Sij.noalias() = HS::calcSij(edge_gradients, edge);
            Sij_bar.noalias() = HS::calcSijBar(Sij);
            tau.noalias() = viscosity_avg * 2.0 * Sij_bar;

            heat_flux.noalias() = thermal_conductivity_avg * tgrad;
            energy_stress.noalias() = T(2.0) * Sij_bar * uvw_viscosity;

            energy_eqn_terms.noalias() = heat_flux + energy_stress;

            normal << T(edge_normals.nx[edge]), T(edge_normals.ny[edge]), T(edge_normals.nz[edge]);
            momentum_fluxes.noalias() = tau * normal;
            T energy_flux = (energy_eqn_terms.transpose() * normal)(0,0);
            std::array<T, 4> edge_flux = {momentum_fluxes(0,0), momentum_fluxes(1,0), momentum_fluxes(2,0), energy_flux};

            int left_node = edge_to_node[edge][0];
            int right_node = edge_to_node[edge][1];
            for (int i = 1; i < 5; ++i) {
                flux[left_node][i] += edge_flux[i-1];
                flux[right_node][i] -= edge_flux[i-1];
            }
        }
        return flux;
    };
}
