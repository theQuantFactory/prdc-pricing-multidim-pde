#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;

#include "grid3d.hpp"
#include "model_params.hpp"
#include "fd_coeffs.hpp"
#include "apply_operators.hpp"
#include "thomas_solver.hpp"
#include "implicit_solvers.hpp"
#include "boundary.hpp"
#include "initial_condition.hpp"
#include "douglas_step.hpp"
#include "interpolation.hpp"
#include "benchmark.hpp"

int main() {
    const double S0    = 1.00;   // spot FX 
    const double K     = 1.00;   // strike ATM
    const double T     = 1.0;    // maturité 
    const double sigma = 0.15;   // vol FX

    const double rd0    = 0.03;  // taux domestique initial
    const double rf0    = 0.01;  // taux étranger initial
    const double kd     = 0;     // Mean reversion speed domestique
    const double kf     = 0;     // Mean reversion speed étranger
    const double thetad = kd * rd0;  // Long-term mean domestique
    const double thetaf = kf * rf0;  // Long-term mean étranger
    const double sig_d  = 0;     // vol rd 
    const double sig_f  = 0;     // vol rf 

    // Corrélations
    const double rho_ds = 0.20;
    const double rho_fs = -0.15;
    const double rho_df = 0.10;

    ModelParams p{ sigma, sig_d, sig_f, kd, kf, thetad, thetaf,
                   rho_ds, rho_fs, rho_df };

    // Grille spatiale (uniforme)
    GridParams g;
    g.nS  = 80;  g.S_min  = 0.0;   g.S_max  = 4.0;
    g.nrd = 30;  g.rd_min = -0.02; g.rd_max = 0.10;
    g.nrf = 30;  g.rf_min = -0.02; g.rf_max = 0.08;

    // Construire les vecteurs de nœuds
    const double dS  = (g.S_max  - g.S_min)  / (g.nS  - 1);
    const double drd = (g.rd_max - g.rd_min) / (g.nrd - 1);
    const double drf = (g.rf_max - g.rf_min) / (g.nrf - 1);

    std::vector<double> Svec(g.nS), Rdvec(g.nrd), Rfvec(g.nrf);
    for (std::size_t i=0; i<g.nS;  ++i) Svec[i]  = g.S_min  + i*dS;
    for (std::size_t i=0; i<g.nrd; ++i) Rdvec[i] = g.rd_min + i*drd;
    for (std::size_t i=0; i<g.nrf; ++i) Rfvec[i] = g.rf_min + i*drf;


    //  Paramètres temporels
    const std::size_t N_steps = 200;          
    const double dt            = T / N_steps; 
    const double theta         = 0.5;         

   
    // Condition initiale 

    Grid3D U = makeInitialCondition(Svec, Rdvec, Rfvec, K);

    cout << "Paramètres : S0=" << S0 << " K=" << K
              << " T=" << T << " σ=" << sigma
              << " rd0=" << rd0 << " rf0=" << rf0 << "\n";
    std::cout << "Grille     : nS=" << g.nS << " nrd=" << g.nrd
              << " nrf=" << g.nrf << " N_steps=" << N_steps << "\n\n";


    // Boucle temporelle de T à 0 

    double t = T;
    for (std::size_t n = 0; n < N_steps; ++n) {
        applyBCinS(U, Svec, Rdvec, Rfvec, K, T - t);
        U = douglasStep(U, Svec, Rdvec, Rfvec, p,
                        dS, drd, drf, theta, dt, t, T, K);
        t -= dt;
    }
    applyBCinS(U, Svec, Rdvec, Rfvec, K, T);


    double price_fd = trilinearInterp(U, Svec, Rdvec, Rfvec,
                                      S0, rd0, rf0,
                                      dS, drd, drf);

   
    // Benchmark modele de Garman-Kohlhagen 
   
    double price_gk = garmanKohlhagen(S0, K, rd0, rf0, sigma, T);

    cout << "--- Résultats ---\n";
    cout << "Prix FD  (ADI Douglas) : " << price_fd  << "\n";
    cout << "Prix GK  (analytique)  : " << price_gk  << "\n";
    cout << "Erreur absolue         : " << std::abs(price_fd - price_gk) << "\n";


    // Test de convergence : afficher le prix pour différents N
    
    std::cout << "--- Convergence en N_steps (grille spatiale fixe) ---\n";
    cout << setw(10) << "N_steps"
              << std::setw(14) << "Prix FD"
              << std::setw(14) << "Erreur"
                << "\n";
    cout << string(52, '-') << "\n";

    double prev_err = -1.0;
    for (size_t Ntest : {25u, 50u, 100u, 200u, 400u}) {
        const double dt_test = T / Ntest;
        Grid3D Utest = makeInitialCondition(Svec, Rdvec, Rfvec, K);
        double t_test = T;
        for (std::size_t n = 0; n < Ntest; ++n) {
            applyBCinS(Utest, Svec, Rdvec, Rfvec, K, T - t_test);
            Utest = douglasStep(Utest, Svec, Rdvec, Rfvec, p,
                                dS, drd, drf, theta, dt_test, t_test, T, K);
            t_test -= dt_test;
        }
        applyBCinS(Utest, Svec, Rdvec, Rfvec, K, T);
        double pfd = trilinearInterp(Utest, Svec, Rdvec, Rfvec,
                                     S0, rd0, rf0, dS, drd, drf);
        double err = std::abs(pfd - price_gk);
        cout << std::setw(10) << Ntest
                  << std::setw(14) << pfd
                  << std::setw(14) << err;
        cout << "\n";
        prev_err = err;
    }

    return 0;
}