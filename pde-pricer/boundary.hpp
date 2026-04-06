#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "grid3d.hpp"
#include <vector>
#include <algorithm>

// =============================================================================
// SECTION 7 — CONDITIONS AUX LIMITES
//
// Call FX européen : payoff = max(S - K, 0) à t=T
//
// CL en S :
//   S=0       → u = 0               (le call ne vaut rien si le sous-jacent = 0)
//   S=S_max   → u = S·e^{-rf·τ} - K·e^{-rd·τ}  (Call deep ITM → parité forward)

// =============================================================================

// Met à jour les bords S=0 et S=S_max sur toute la grille
static void applyBCinS(Grid3D& U,
                        const std::vector<double>& Svec,
                        const std::vector<double>& Rdvec,
                        const std::vector<double>& Rfvec,
                        double K, double tau)  // tau = T - t
{
    const size_t nS  = U.numS();
    const size_t nrd = U.numRd();
    const size_t nrf = U.numRf();
    const double S_max    = Svec[nS-1];

    for (size_t ird = 0; ird < nrd; ++ird) {
        const double rd = Rdvec[ird];
        for (size_t irf = 0; irf < nrf; ++irf) {
            const double rf = Rfvec[irf];

            // Bord gauche : S = 0 
            U.at(0, ird, irf) = 0.0;

            // Bord droit : S = S_max 
            U.at(nS-1, ird, irf) = std::max(
                S_max * std::exp(-rf*tau) - K * std::exp(-rd*tau),
                0.0);
        }
    }
}

// Calcule g1_bnd(ird, irf) = u(S_max) au temps t
// (à appeler avant chaque pas de temps pour mettre à jour le bord)
static Grid3D computeG1Boundary(
    const std::vector<double>& Svec,
    const std::vector<double>& Rdvec,
    const std::vector<double>& Rfvec,
    double K, double tau)
{
    const size_t nS  = Svec.size();
    const size_t nrd = Rdvec.size();
    const size_t nrf = Rfvec.size();
    const double S_max    = Svec[nS-1];

    // Dimension 1 en S (on stocke uniquement le bord S_max)
    Grid3D g1(1, nrd, nrf);
    for (size_t ird = 0; ird < nrd; ++ird) {
        const double rd = Rdvec[ird];
        for (size_t irf = 0; irf < nrf; ++irf) {
            const double rf = Rfvec[irf];
            g1.at(0, ird, irf) = std::max(
                S_max * std::exp(-rf*tau) - K * std::exp(-rd*tau),
                0.0);
        }
    }
    return g1;
}

#endif // BOUNDARY_HPP