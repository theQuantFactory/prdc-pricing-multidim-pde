#ifndef IMPLICIT_SOLVERS_HPP
#define IMPLICIT_SOLVERS_HPP

#include "grid3d.hpp"
#include "fd_coeffs.hpp"
#include "thomas_solver.hpp"
#include <vector>

// resolution implicite direction S : (I - θ·dt·A1)·U = rhs
static void solveImplicitS(Grid3D& rhs,
                            const std::vector<double>& Svec,
                            const std::vector<double>& Rdvec,
                            const std::vector<double>& Rfvec,
                            const ModelParams& p, double dS,
                            double theta, double dt,
                            const Grid3D& g1_bnd)
{
    const std::size_t nS  = rhs.numS();
    const std::size_t nrd = rhs.numRd();
    const std::size_t nrf = rhs.numRf();
    const std::size_t n   = nS - 2;  

    std::vector<double> a(n), b(n), c(n), d(n);

    for (std::size_t ird = 0; ird < nrd; ++ird) {
        for (std::size_t irf = 0; irf < nrf; ++irf) {

            // Construire la matrice (I - θ·dt·A1) pour cette ligne (ird, irf) et le RHS
            for (size_t k = 0; k < n; ++k) {
                size_t iS = k + 1;  // point intérieur
                Tri co = coeffA1(Svec[iS], Rdvec[ird], Rfvec[irf], p.sigma, dS);

                a[k] = -theta * dt * co.lo;           // sous-diag
                b[k] =  1.0   - theta * dt * co.di;   // diag    
                c[k] = -theta * dt * co.up;           // sur-diag

                d[k] = rhs.at(iS, ird, irf);          // RHS courant

                // *** GARDE : diagonale dominante ***
                // Si violée → dt trop grand → schéma instable
                if (abs(b[k]) < std::abs(a[k]) + std::abs(c[k]))
                    throw std::runtime_error(
                        "Thomas S: diag non dominante, réduire dt");
            }

            // Correction bord gauche S=0 
            d[0] -= a[0] * rhs.at(0, ird, irf);         

            // Correction bord droit S_max 
            d[n-1] -= c[n-1] * g1_bnd.at(0, ird, irf);

            // Thomas
            thomasSolve(a, b, c, d);

            // Réécrire la solution dans rhs (points intérieurs)
            for (size_t k = 0; k < n; ++k)
                rhs.at(k+1, ird, irf) = d[k];
        }
    }
}

// resolution implicite direction rd : (I - θ·dt·A2)·U = rhs

static void solveImplicitRd(Grid3D& rhs,
                             const std::vector<double>& Rdvec,
                             const ModelParams& p, double drd,
                             double theta, double dt)
{
    const std::size_t nS  = rhs.numS();
    const std::size_t nrd = rhs.numRd();
    const std::size_t nrf = rhs.numRf();
    const std::size_t n   = nrd - 2;

    std::vector<double> a(n), b(n), c(n), d(n);

    for (std::size_t iS = 0; iS < nS; ++iS) {
        for (std::size_t irf = 0; irf < nrf; ++irf) {

            for (std::size_t k = 0; k < n; ++k) {
                std::size_t ird = k + 1;
                Tri co = coeffA2(Rdvec[ird], p.theta_d, p.kappa_d, p.sigma_d, drd);

                a[k] = -theta * dt * co.lo;
                b[k] =  1.0   - theta * dt * co.di;
                c[k] = -theta * dt * co.up;
                d[k] = rhs.at(iS, ird, irf);
            }

            
            d[0]   -= a[0]   * rhs.at(iS, 0,      irf);
            d[n-1] -= c[n-1] * rhs.at(iS, nrd-1,  irf);

            thomasSolve(a, b, c, d);

            for (std::size_t k = 0; k < n; ++k)
                rhs.at(iS, k+1, irf) = d[k];

            // Propager CL Neumann aux bords
            rhs.at(iS, 0,      irf) = rhs.at(iS, 1,      irf);
            rhs.at(iS, nrd-1,  irf) = rhs.at(iS, nrd-2,  irf);
        }
    }
}

// resolution implicite direction rf : (I - θ·dt·A3)·U = rhs
static void solveImplicitRf(Grid3D& rhs,
                             const std::vector<double>& Rfvec,
                             const ModelParams& p, double drf,
                             double theta, double dt)
{
    const std::size_t nS  = rhs.numS();
    const std::size_t nrd = rhs.numRd();
    const std::size_t nrf = rhs.numRf();
    const std::size_t n   = nrf - 2;

    std::vector<double> a(n), b(n), c(n), d(n);

    for (std::size_t iS = 0; iS < nS; ++iS) {
        for (std::size_t ird = 0; ird < nrd; ++ird) {

            for (std::size_t k = 0; k < n; ++k) {
                std::size_t irf = k + 1;
                Tri co = coeffA3(Rfvec[irf], p.theta_f, p.kappa_f,
                                 p.sigma_f, p.sigma, p.rho_fs, drf);

                a[k] = -theta * dt * co.lo;
                b[k] =  1.0   - theta * dt * co.di;
                c[k] = -theta * dt * co.up;
                d[k] = rhs.at(iS, ird, irf);
            }

            d[0]   -= a[0]   * rhs.at(iS, ird, 0);
            d[n-1] -= c[n-1] * rhs.at(iS, ird, nrf-1);

            thomasSolve(a, b, c, d);

            for (std::size_t k = 0; k < n; ++k)
                rhs.at(iS, ird, k+1) = d[k];

            rhs.at(iS, ird, 0)      = rhs.at(iS, ird, 1);
            rhs.at(iS, ird, nrf-1)  = rhs.at(iS, ird, nrf-2);
        }
    }
}

#endif // IMPLICIT_SOLVERS_HPP