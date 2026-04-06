
##################################################################################

#stockage des coefficients de corrélation pour le modèle SDEDATA 

class sdedataC():
    def __init__(self, rho_domestic_base, rho_base_fx, rho_domestic_fx):

        self.rho_domestic_base = rho_domestic_base
        self.rho_base_fx = rho_base_fx
        self.rho_domestic_fx = rho_domestic_fx
        return
##################################################################################

##################################################################################

#stocker tous les paramètres du modèle G1++(hull white) pour les deux taux d'intérêt stochastiques du modèle multi-devise.
class g1pp_dataC():
    def __init__(self, spot_FX, base_x0, domestic_x0, g1ppTA, domestic_shifttimes, g1ppTVOL, g1ppA, g1ppVOL, domestic_shiftvalues,
                g1ppdivTA, g1ppdivMEANREV, base_shifttimes, base_shiftvalues, g1ppdivTVOL, g1ppdivVOL):
        
        self.spot_FX=spot_FX
        self.base_x0 = base_x0
        self.domestic_x0 = domestic_x0
        
        self.domestic_ta = g1ppTA
        self.domestic_a = g1ppA
        self.domestic_shifttimes = domestic_shifttimes
        self.domestic_shiftvalues = domestic_shiftvalues
        self.domestic_tvol = g1ppTVOL
        self.domestic_vol = g1ppVOL

        self.base_ta = g1ppdivTA
        self.base_a = g1ppdivMEANREV
        self.base_shifttimes = base_shifttimes
        self.base_shiftvalues = base_shiftvalues
        self.base_tvol = g1ppdivTVOL
        self.base_vol = g1ppdivVOL
        return
##################################################################################

