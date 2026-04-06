#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include <cmath>

// =============================================================================
// SECTION 11 — BENCHMARK GARMAN-KOHLHAGEN (formule analytique)
// =============================================================================
static double normCDF(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double garmanKohlhagen(double S, double K, double rd, double rf,
                               double sigma, double T) {
    double d1 = (std::log(S/K) + (rd - rf + 0.5*sigma*sigma)*T)
                / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S*std::exp(-rf*T)*normCDF(d1) - K*std::exp(-rd*T)*normCDF(d2);
}

#endif // BENCHMARK_HPP