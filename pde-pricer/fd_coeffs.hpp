#ifndef FD_COEFFS_HPP
#define FD_COEFFS_HPP

#include <cmath>

// =============================================================================
// SECTION 3 — OPÉRATEURS difference finies pour les dérivées spatiales
//
// Sur grille uniforme d'espacement h :
//   ∂u/∂x  → (u[i+1] - u[i-1]) / (2h)          central O(h²)
//   ∂²u/∂x² → (u[i-1] - 2u[i] + u[i+1]) / h²   central O(h²)
// =============================================================================

// Coefficients tridiagonaux [lower, diag, upper] pour un point intérieur
struct Tri { double lo, di, up; };

// ---- A1 : direction S -------------------------------------------------------
// A1·u = (1/2)σ²s²·∂²u/∂s² + (rd-rf)s·∂u/∂s - rd·u
static Tri coeffA1(double s, double rd, double rf,
                   double sigma, double dS)
{
    const double half_v2_s2 = 0.5 * sigma*sigma * s*s;
    const double drift       = (rd - rf) * s;
    const double inv_dS2     = 1.0 / (dS*dS);
    const double inv_2dS     = 0.5 / dS;
    return {
        half_v2_s2*inv_dS2 - drift*inv_2dS,   // lo
       -sigma*sigma*s*s*inv_dS2 - rd,           // di  (négatif !)
        half_v2_s2*inv_dS2 + drift*inv_2dS    // up
    };
}

// ---- A2 : direction rd -------------------------------------------------------
// A2·u = (1/2)σd²·∂²u/∂rd² + (θd - κd·rd)·∂u/∂rd

static Tri coeffA2(double rd, double theta_d, double kappa_d,
                   double sigma_d, double drd)
{
    const double half_v2  = 0.5 * sigma_d*sigma_d;
    const double drift    = theta_d - kappa_d*rd;
    const double inv_drd2 = 1.0 / (drd*drd);
    const double inv_2drd = 0.5 / drd;
    return {
        half_v2*inv_drd2 - drift*inv_2drd,
       -sigma_d*sigma_d*inv_drd2,
        half_v2*inv_drd2 + drift*inv_2drd
    };
}

// ---- A3 : direction rf -------------------------------------------------------
// A3·u = (1/2)σf²·∂²u/∂rf² + (θf - κf·rf - ρfs·σf·σ)·∂u/∂rf                                             
 
static Tri coeffA3(double rf, double theta_f, double kappa_f,
                   double sigma_f, double sigma, double rho_fs,
                   double drf)
{
    const double half_v2  = 0.5 * sigma_f*sigma_f;
    // Drift avec correction de Girsanov : -ρfs·σf·σ
    const double drift    = theta_f - kappa_f*rf - rho_fs*sigma_f*sigma;
    const double inv_drf2 = 1.0 / (drf*drf);
    const double inv_2drf = 0.5 / drf;
    return {
        half_v2*inv_drf2 - drift*inv_2drf,
       -sigma_f*sigma_f*inv_drf2,
        half_v2*inv_drf2 + drift*inv_2drf
    };
}

#endif // FD_COEFFS_HPP