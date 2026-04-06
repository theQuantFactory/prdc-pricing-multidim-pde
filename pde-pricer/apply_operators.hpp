#ifndef APPLY_OPERATORS_HPP
#define APPLY_OPERATORS_HPP

#include "grid3d.hpp"
#include "fd_coeffs.hpp"
#include <vector>

// =============================================================================
// SECTION 4 — APPLICATION EXPLICITE DES OPÉRATEURS
// Calcule out += A·U   (seuls les points intérieurs sont touchés)
// =============================================================================

// Applique A1·U et ajoute dans out
static void applyA1(const Grid3D& U, Grid3D& out,
                    const std::vector<double>& Svec,
                    const std::vector<double>& Rdvec,
                    const std::vector<double>& Rfvec,
                    const ModelParams& p, double dS)
{
    const size_t nS=U.numS(), nrd=U.numRd(), nrf=U.numRf();
    for (size_t ird=0; ird<nrd; ++ird)
        for (size_t irf=0; irf<nrf; ++irf)
            for (size_t iS=1; iS<nS-1; ++iS) {
                Tri c = coeffA1(Svec[iS], Rdvec[ird], Rfvec[irf], p.sigma, dS);
                out.at(iS,ird,irf) +=
                    c.lo*U.at(iS-1,ird,irf) + c.di*U.at(iS,ird,irf)
                  + c.up*U.at(iS+1,ird,irf);
            }
}

// Applique A2·U et ajoute dans out
static void applyA2(const Grid3D& U, Grid3D& out,
                    const std::vector<double>& Rdvec,
                    const ModelParams& p, double drd)
{
    const size_t nS=U.numS(), nrd=U.numRd(), nrf=U.numRf();
    for (size_t iS=0; iS<nS; ++iS)
        for (size_t irf=0; irf<nrf; ++irf)
            for (size_t ird=1; ird<nrd-1; ++ird) {
                Tri c = coeffA2(Rdvec[ird], p.theta_d, p.kappa_d, p.sigma_d, drd);
                out.at(iS,ird,irf) +=
                    c.lo*U.at(iS,ird-1,irf) + c.di*U.at(iS,ird,irf)
                  + c.up*U.at(iS,ird+1,irf);
            }
}

// Applique A3·U et ajoute dans out
static void applyA3(const Grid3D& U, Grid3D& out,
                    const std::vector<double>& Rfvec,
                    const ModelParams& p, double drf)
{
    const size_t nS=U.numS(), nrd=U.numRd(), nrf=U.numRf();
    for (size_t iS=0; iS<nS; ++iS)
        for (size_t ird=0; ird<nrd; ++ird)
            for (size_t irf=1; irf<nrf-1; ++irf) {
                Tri c = coeffA3(Rfvec[irf], p.theta_f, p.kappa_f,
                                p.sigma_f, p.sigma, p.rho_fs, drf);
                out.at(iS,ird,irf) +=
                    c.lo*U.at(iS,ird,irf-1) + c.di*U.at(iS,ird,irf)
                  + c.up*U.at(iS,ird,irf+1);
            }
}

// Applique A0·U (termes croisés) et ajoute dans out
// Discrétisation : ∂²u/∂x∂y ≈ [u(i+1,j+1)-u(i+1,j-1)-u(i-1,j+1)+u(i-1,j-1)]/(4hx·hy)
static void applyA0(const Grid3D& U, Grid3D& out,
                    const std::vector<double>& Svec,
                    const ModelParams& p,
                    double dS, double drd, double drf)
{
    const size_t nS=U.numS(), nrd=U.numRd(), nrf=U.numRf();

    for (size_t iS=1; iS<nS-1; ++iS) {
        const double s    = Svec[iS];
        // Coefficients des 3 termes croisés
        const double c_Srd  = p.rho_ds * p.sigma_d * p.sigma * s / (4.0*dS*drd);
        const double c_Srf  = p.rho_fs * p.sigma_f * p.sigma * s / (4.0*dS*drf);
        const double c_rdrf = p.rho_df * p.sigma_d * p.sigma_f   / (4.0*drd*drf);

        for (size_t ird=1; ird<nrd-1; ++ird)
            for (size_t irf=1; irf<nrf-1; ++irf) {

                // Terme 1 : ρds·σd·σ·S·∂²u/∂S∂rd
                out.at(iS,ird,irf) += c_Srd*(
                    + U.at(iS+1,ird+1,irf) - U.at(iS+1,ird-1,irf)
                    - U.at(iS-1,ird+1,irf) + U.at(iS-1,ird-1,irf));

                // Terme 2 : ρfs·σf·σ·S·∂²u/∂S∂rf
                out.at(iS,ird,irf) += c_Srf*(
                    + U.at(iS+1,ird,irf+1) - U.at(iS+1,ird,irf-1)
                    - U.at(iS-1,ird,irf+1) + U.at(iS-1,ird,irf-1));

                // Terme 3 : ρdf·σd·σf·∂²u/∂rd∂rf
                out.at(iS,ird,irf) += c_rdrf*(
                    + U.at(iS,ird+1,irf+1) - U.at(iS,ird+1,irf-1)
                    - U.at(iS,ird-1,irf+1) + U.at(iS,ird-1,irf-1));
            }
    }
}

#endif // APPLY_OPERATORS_HPP