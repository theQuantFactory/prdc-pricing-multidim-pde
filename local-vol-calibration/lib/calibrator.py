"""Local volatility calibration utilities under stochastic interest rates"""
import inspect
import math
import os
import pickle
import statistics

import numpy as np
import interpolator

def split_seq(seq, size):
    """
    Splits a given list evenly to sublists of roughly the same size

    Parameters
    ----------
    seq : list
        Original sequence
    size : int
        Number of sublists

    Returns
    -------
    newseq : list of lists

    """
    newseq = []
    splitsize = 1.0/size*len(seq)
    for i in range(size):
        newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
    return newseq

class LocalVolIRCalibration:

    def __init__(self, surface, simulator, filename=None):
        self.surface = surface
        self.simulator = simulator
        self.filename = filename

    def save_lv_surface(self, all_strikes, times, surfacevalues, surface_errs=None):
        """
        Save the local vol surface data

        Parameters
        ----------
        all_strikes : list of floats
            List of strikes
        times : list of floats
            List of times
        surfacevalues : list of lists of floats
            Vol surface data
        surface_errs : list of lists of floats
            Vol surface errors

        """
        if self.filename:
            dirname = os.path.dirname(self.filename)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            lvol_output = { 'spots': all_strikes,
                            'times': times,
                            'surfacevalues': surfacevalues,
                           }
            if surface_errs:
                lvol_output['errors'] = surface_errs
            with open(self.filename, 'wb') as f:
                pickle.dump(lvol_output, f)
            print("Written: %s" % (self.filename))

class DeterministicLocalVolDeterministicIRCalibration(LocalVolIRCalibration):


    def __init__(self, surface, simulator, filename=None):
        super(DeterministicLocalVolDeterministicIRCalibration, self).__init__(surface, simulator,
                                                                              filename=filename,
                                                                              )

    def calibrate_localvol(self, Ks, ts):
        """
        Do the calibration

        Parameters
        ----------
        Ks : list of lists of floats
            2D strike grid
        ts : list of floats
            Time grid

        Returns
        -------
        locvols : list of lists of floats
            The local volatility surface

        """
        if len(Ks) != len(ts):
            raise Exception("Strike grid length does not match times grid length")
        locvols = []
        all_strikes = Ks
        for idx, t, in enumerate(ts):
            print("Processing time: %f" % t)
            strikes = all_strikes[idx]
            locvol_slice = self.surface.evaluate_localvol_slice_det_ir(strikes, t)
            locvols.append(locvol_slice)
            self.save_lv_surface(all_strikes[:idx+1], ts[:idx+1], locvols)

        return locvols


class DeterministicLocalVolStochasticIRCalibration(LocalVolIRCalibration):


    def __init__(self, surface, simulator, filename=None):
        super(DeterministicLocalVolStochasticIRCalibration, self).__init__(surface, simulator,
                                                                           filename=filename,
                                                                           )
        self.deterministic_ir_calibration = DeterministicLocalVolDeterministicIRCalibration(surface, simulator,
                                                                                            filename=filename,
                                                                                            )

    def calibrate_localvol(self, Ks, ts):
        """
        Do the calibration

        Parameters
        ----------
        Ks : list of lists of floats
            2D strike grid
        ts : list of floats
            Time grid

        Returns
        -------
        locvols : list of lists of floats
            The local volatility surface

        """
        if len(Ks) != len(ts):
            raise excp.InvalidLengthException("Strike grid length does not match times grid length")
        locvols = []
        locvol_errs = []
        all_strikes = Ks
        for idx, t, in enumerate(ts):
            print("Processing time: %f" % t)
            strikes = all_strikes[idx]
            if idx == 0:
                locvol_slice = self.deterministic_ir_calibration.surface.evaluate_localvol_slice_det_ir(strikes, t)
                locvol_err_slice = [0.0] * len(locvol_slice)
            else:
                self.simulator.set_contracts_for_calibration(strikes, t)
                self.simulator.update_localvol(all_strikes[:idx], ts[:idx], locvols)
                self.simulator.update_domestic_bfactors(t)
                self.simulator.update_base_bfactors(t)
                self.simulator.run()
                expectations = self.simulator.get_means_from_output()
                stderrs = self.simulator.get_stderrs_from_output()
                locvol_slice, locvol_err_slice = self.surface.evaluate_localvol_slice_stoc_ir(strikes, t, expectations, stderrs)
            locvols.append(locvol_slice)
            locvol_errs.append(locvol_err_slice)
            self.save_lv_surface(all_strikes[:idx+1], ts[:idx+1], locvols, locvol_errs)

        return locvols

