#ifndef THOMAS_SOLVER_HPP
#define THOMAS_SOLVER_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

// =============================================================================
// SECTION 5 — ALGORITHME DE THOMAS (résolution tridiagonale) pour resoudre (I - θ·dt·A)·x = d
// =============================================================================
static void thomasSolve(std::vector<double>& a,  // sous-diag (in)
                         std::vector<double>& b,  // diag      (in/out)
                         std::vector<double>& c,  // sur-diag  (in/out)
                         std::vector<double>& d)  // RHS       (in) → solution (out)
{
    const int n = static_cast<int>(d.size());

    // --- Forward sweep : éliminer a[] ---
    for (int k = 1; k < n; ++k) {
        if (std::abs(b[k-1]) < 1e-15)
            throw std::runtime_error("Thomas: pivot nul");
        const double m = a[k] / b[k-1];  // facteur d'élimination
        b[k] -= m * c[k-1];              // mise à jour diagonale
        d[k] -= m * d[k-1];              // mise à jour RHS
    }

    // --- Backward substitution : remonter ---
    d[n-1] /= b[n-1];
    for (int k = n-2; k >= 0; --k) {
        if (std::abs(b[k]) < 1e-15)
            throw std::runtime_error("Thomas: pivot nul (backward)");
        d[k] = (d[k] - c[k]*d[k+1]) / b[k];
    }
}

#endif // THOMAS_SOLVER_HPP