
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

#include "Tensor.h"
#include "Point.h"

namespace HS {
    template<typename T>
    inline Tensor<T, 3, 3>
    calcSij(const Point<T> &grad_u,
            const Point<T> &grad_v,
            const Point<T> &grad_w){
        Tensor<T, 3, 3> Sij;

        Sij(0, 0) = grad_u[0];
        Sij(0, 1) = 0.5 * (grad_u[1] + grad_v[0]);
        Sij(0, 2) = 0.5 * (grad_u[2] + grad_w[0]);

        Sij(1, 0) = Sij(0,1);
        Sij(1, 1) = grad_v[1];
        Sij(1, 2) = 0.5 * (grad_v[2] + grad_w[1]);

        Sij(2, 0) = Sij(0, 2);
        Sij(2, 1) = Sij(1, 2);
        Sij(2, 2) = grad_w[2];

        return Sij;
    }

    template<typename T>
    inline Tensor<T, 3, 3>
    getKroneckerDelta(){
        Tensor<T, 3, 3> delta;
        delta(0,0) = 1.0;
        delta(1,1) = 1.0;
        delta(2,2) = 1.0;
        return delta;
    };

    template<typename T>
    inline Tensor<T, 3, 3>
    calcSijBar(const Tensor<T, 3, 3>& Sij) {
        return (Sij - (Sij.trace() / 3.0) * getKroneckerDelta<T>());
    };

    template<size_t n_corners, size_t NumEqns, size_t n_edges, typename T>
    std::array<std::array<T, NumEqns>, n_corners>
    ElementBasedViscousFlux(const std::array<std::array<int, 2>, n_edges>& edge_to_node,
                            const std::array<HS::Point<double>, n_edges> &edge_normals,
                            const std::array<HS::Point<T>, n_edges> &edge_ugrad,
                            const std::array<HS::Point<T>, n_edges> &edge_vgrad,
                            const std::array<HS::Point<T>, n_edges> &edge_wgrad,
                            const std::array<HS::Point<T>, n_edges> &edge_tgrad,
                            const T& thermal_conductivity_avg, // pass directly because prandtl number might be different for turbulence
                            const T& viscosity_avg,
                            const HS::Point<T>& uvw_viscosity_avg) {

        Tensor<T, 3, 3> Sij, Sij_bar, tau;

        std::array<std::array<T, NumEqns>, n_corners> flux;
        for (auto &f : flux)
            for (auto &value : f)
                value = 0.0;

        for (size_t edge = 0; edge < n_edges; ++edge) {

            Sij = HS::calcSij(edge_ugrad[edge], edge_vgrad[edge], edge_wgrad[edge]);
            Sij_bar = HS::calcSijBar(Sij);
            tau = viscosity_avg * 2.0 * Sij_bar;

            auto heat_flux = thermal_conductivity_avg * edge_tgrad[edge];
            auto energy_stress = T(2.0) * Sij_bar * uvw_viscosity_avg;

            auto energy_eqn_terms = heat_flux + energy_stress;

            auto momentum_fluxes = tau * edge_normals[edge];
            T energy_flux = energy_eqn_terms[0] * edge_normals[edge][0]
                          + energy_eqn_terms[1] * edge_normals[edge][1]
                          + energy_eqn_terms[2] * edge_normals[edge][2];
            std::array<T, 4> edge_flux = {momentum_fluxes[0], momentum_fluxes[1], momentum_fluxes[2], energy_flux};

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
