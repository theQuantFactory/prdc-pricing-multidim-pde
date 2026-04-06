#ifndef GRID3D_HPP
#define GRID3D_HPP

#include <vector>

// =============================================================================
// SECTION 1 — GRILLE 3D (S, rd, rf)
// Stockage : index = iS*(nrd*nrf) + ird*nrf + irf
// =============================================================================
class Grid3D {
public:
    Grid3D(size_t nS, size_t nrd, size_t nrf)
        : nS_(nS), nrd_(nrd), nrf_(nrf), data_(nS*nrd*nrf, 0.0) {}
        // construire une grille nS×nrd×nrf initialisée à 0.0

    
    double& at(size_t iS, size_t ird, size_t irf) {
        return data_[iS*nrd_*nrf_ + ird*nrf_ + irf];
    }
        // Accès à l'élément (iS, ird, irf)
    double at(size_t iS, size_t ird, size_t irf) const {
        return data_[iS*nrd_*nrf_ + ird*nrf_ + irf];
    }

    // Accès plat (pour boucles rapides)
    double& flat(size_t i)       { return data_[i]; }
    double  flat(size_t i) const { return data_[i]; }

    size_t size()  const { return data_.size(); }
    size_t numS()  const { return nS_;  }
    size_t numRd() const { return nrd_; }
    size_t numRf() const { return nrf_; }

    // U += alpha * rhs  (utilisé pour construire Y0, Y1, Y2)
    void addScaled(double alpha, const Grid3D& rhs) {
        for (size_t i = 0; i < data_.size(); ++i)
            data_[i] += alpha * rhs.data_[i];
    }

    void copyFrom(const Grid3D& src) {
        data_ = src.data_;
    }

private:
    size_t nS_, nrd_, nrf_;
    std::vector<double> data_;
};

#endif // GRID3D_HPP