#ifndef MODEL_PARAMS_HPP
#define MODEL_PARAMS_HPP

#include <cstddef>

// =============================================================================
// SECTION 2 — PARAMÈTRES DU MODÈLE ET DE LA GRILLE
// =============================================================================

// Paramètres financiers du modèle PRDc
struct ModelParams {
    double sigma;     // volatilité FX (constante)
    double sigma_d;   // vol taux domestique rd
    double sigma_f;   // vol taux étranger  rf
    double kappa_d;   // mean-reversion rd
    double kappa_f;   // mean-reversion rf
    double theta_d;   // niveau long terme rd  (θd constant ici)
    double theta_f;   // niveau long terme rf  (θf constant ici)
    double rho_ds;    // corrélation S-rd
    double rho_fs;    // corrélation S-rf
    double rho_df;    // corrélation rd-rf
};

// Paramètres de discrétisation (grille uniforme)
struct GridParams {
    // Direction S
    size_t nS;
    double S_min, S_max;

    // Direction rd
    size_t nrd;
    double rd_min, rd_max;

    // Direction rf
    size_t nrf;
    double rf_min, rf_max;
};

#endif // MODEL_PARAMS_HPP