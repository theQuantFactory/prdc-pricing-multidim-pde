#ifndef INITIAL_CONDITION_HPP
#define INITIAL_CONDITION_HPP

#include "grid3d.hpp"
#include <vector>
#include <algorithm>

// =============================================================================
// SECTION 8 — CONDITION INITIALE (payoff à t=T)
// =============================================================================
static Grid3D makeInitialCondition(
    const std::vector<double>& Svec,
    const std::vector<double>& Rdvec,
    const std::vector<double>& Rfvec,
    double K)
{
    const size_t nS  = Svec.size();
    const size_t nrd = Rdvec.size();
    const size_t nrf = Rfvec.size();

    Grid3D U(nS, nrd, nrf);
    for (size_t iS = 0; iS < nS; ++iS) {
        const double payoff = std::max(Svec[iS] - K, 0.0);
        for (size_t ird = 0; ird < nrd; ++ird)
            for (size_t irf = 0; irf < nrf; ++irf)
                U.at(iS, ird, irf) = payoff;
    }
    return U;
}

#endif // INITIAL_CONDITION_HPP