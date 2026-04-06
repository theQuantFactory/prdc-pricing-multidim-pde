#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include "grid3d.hpp"
#include <vector>
#include <algorithm>

// =============================================================================
// SECTION 10 — INTERPOLATION TRILINÉAIRE
// Interpolation dans la grille pour lire u(S0, rd0, rf0)
// =============================================================================
static double trilinearInterp(
    const Grid3D& U,
    const std::vector<double>& Svec,
    const std::vector<double>& Rdvec,
    const std::vector<double>& Rfvec,
    double S0, double rd0, double rf0,
    double dS, double drd, double drf)
{
    const size_t nS  = Svec.size();
    const size_t nrd = Rdvec.size();
    const size_t nrf = Rfvec.size();

    // Clamper dans les bornes
    S0  = std::max(Svec[0],   std::min(Svec[nS-1],   S0));
    rd0 = std::max(Rdvec[0],  std::min(Rdvec[nrd-1], rd0));
    rf0 = std::max(Rfvec[0],  std::min(Rfvec[nrf-1], rf0));

    // Indices inférieurs
    auto idx = [](double x, double x0, double h, std::size_t n) {
        return std::min(static_cast<std::size_t>((x-x0)/h), n-2);
    };
    size_t iS  = idx(S0,  Svec[0],  dS,  nS);
    size_t ird = idx(rd0, Rdvec[0], drd, nrd);
    size_t irf = idx(rf0, Rfvec[0], drf, nrf);

    // Coordonnées locales dans [0,1]
    double aS  = (S0  - Svec[iS])   / dS;
    double ard = (rd0 - Rdvec[ird]) / drd;
    double arf = (rf0 - Rfvec[irf]) / drf;

    // Interpolation trilinéaire sur les 8 sommets du cube
    return (1-aS)*(1-ard)*(1-arf)*U.at(iS,   ird,   irf  )
          +   aS *(1-ard)*(1-arf)*U.at(iS+1, ird,   irf  )
          +(1-aS)*   ard *(1-arf)*U.at(iS,   ird+1, irf  )
          +   aS *   ard *(1-arf)*U.at(iS+1, ird+1, irf  )
          +(1-aS)*(1-ard)*   arf *U.at(iS,   ird,   irf+1)
          +   aS *(1-ard)*   arf *U.at(iS+1, ird,   irf+1)
          +(1-aS)*   ard *   arf *U.at(iS,   ird+1, irf+1)
          +   aS *   ard *   arf *U.at(iS+1, ird+1, irf+1);
}

#endif // INTERPOLATION_HPP