#ifndef DOUGLAS_STEP_HPP
#define DOUGLAS_STEP_HPP

#include "grid3d.hpp"
#include "model_params.hpp"
#include "apply_operators.hpp"
#include "implicit_solvers.hpp"
#include "boundary.hpp"
#include <vector>

// =============================================================================
// SECTION 9 — UN PAS DE TEMPS (schéma de Douglas complet)
//
// Résumé du schéma Douglas (θ = 0.5 → Crank-Nicolson généralisé) :
//
//   Y0 = U + dt·(A0 + A1 + A2 + A3)·U          ← prédicteur explicite
//
//   Y1 = Y0 - θdt·A1·U
//   Résoudre (I - θdt·A1)·Y1 = Y1               ← Thomas en S
//
//   Y2 = Y1 - θdt·A2·U
//   Résoudre (I - θdt·A2)·Y2 = Y2               ← Thomas en rd
//
//   U_new = Y2 - θdt·A3·U
//   Résoudre (I - θdt·A3)·U_new = U_new         ← Thomas en rf
// =============================================================================
static Grid3D douglasStep(
    const Grid3D& U,
    const std::vector<double>& Svec,
    const std::vector<double>& Rdvec,
    const std::vector<double>& Rfvec,
    const ModelParams& p,
    double dS, double drd, double drf,
    double theta, double dt,
    double t_current,          // temps avant ce pas
    double T, double K)        // pour mettre à jour g1_bnd
{
    const std::size_t nS  = U.numS();
    const std::size_t nrd = U.numRd();
    const std::size_t nrf = U.numRf();

    // -------------------------------------------------------------------------
    // Étape 0 : Y0 = U + dt·(A0 + A1 + A2 + A3)·U
    // -------------------------------------------------------------------------
    Grid3D Y0(nS, nrd, nrf);
    Y0.copyFrom(U);  // Y0 = U

    // tmp = A0·U  puis  Y0 += dt·tmp
    { Grid3D tmp(nS,nrd,nrf); applyA0(U,tmp,Svec,p,dS,drd,drf); Y0.addScaled(dt,tmp); }
    // tmp = A1·U  puis  Y0 += dt·tmp
    { Grid3D tmp(nS,nrd,nrf); applyA1(U,tmp,Svec,Rdvec,Rfvec,p,dS); Y0.addScaled(dt,tmp); }
    // tmp = A2·U  puis  Y0 += dt·tmp
    { Grid3D tmp(nS,nrd,nrf); applyA2(U,tmp,Rdvec,p,drd); Y0.addScaled(dt,tmp); }
    // tmp = A3·U  puis  Y0 += dt·tmp
    { Grid3D tmp(nS,nrd,nrf); applyA3(U,tmp,Rfvec,p,drf); Y0.addScaled(dt,tmp); }

    // -------------------------------------------------------------------------
    // Étape 1 : correction direction S
    //   Y1 = Y0 - θ·dt·A1·U
    //   (I - θ·dt·A1)·Y1 = Y1
    // -------------------------------------------------------------------------
    Grid3D Y1(nS,nrd,nrf); Y1.copyFrom(Y0);

    // Y1 -= θ·dt·A1·U
    { Grid3D A1U(nS,nrd,nrf); applyA1(U,A1U,Svec,Rdvec,Rfvec,p,dS);
      Y1.addScaled(-theta*dt, A1U); }

    // Mettre à jour le bord S_max pour le temps t_current - dt
    double tau_new = T - (t_current - dt);
    Grid3D g1_bnd = computeG1Boundary(Svec, Rdvec, Rfvec, K, tau_new);

    // Thomas en S (avec correction de bord g1_bnd)
    solveImplicitS(Y1, Svec, Rdvec, Rfvec, p, dS, theta, dt, g1_bnd);

    // -------------------------------------------------------------------------
    // Étape 2 : correction direction rd
    //   Y2 = Y1 - θ·dt·A2·U
    //   (I - θ·dt·A2)·Y2 = Y2
    // -------------------------------------------------------------------------
    Grid3D Y2(nS,nrd,nrf); Y2.copyFrom(Y1);

    { Grid3D A2U(nS,nrd,nrf); applyA2(U,A2U,Rdvec,p,drd);
      Y2.addScaled(-theta*dt, A2U); }

    solveImplicitRd(Y2, Rdvec, p, drd, theta, dt);

    // -------------------------------------------------------------------------
    // Étape 3 : correction direction rf → résultat final
    //   U_new = Y2 - θ·dt·A3·U
    //   (I - θ·dt·A3)·U_new = U_new
    // -------------------------------------------------------------------------
    Grid3D U_new(nS,nrd,nrf); U_new.copyFrom(Y2);

    { Grid3D A3U(nS,nrd,nrf); applyA3(U,A3U,Rfvec,p,drf);
      U_new.addScaled(-theta*dt, A3U); }

    solveImplicitRf(U_new, Rfvec, p, drf, theta, dt);

    return U_new;
}

#endif // DOUGLAS_STEP_HPP