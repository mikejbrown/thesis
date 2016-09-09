# -*- coding: utf-8 -*-
"""

Generate plots for thesis ch 3 on solution of 2PI HF gap equations

Modified from aip2014_poster_plots.py

@author: Michael Brown
"""

import numpy as np
import scipy as sc
import pylab as py
from matplotlib.colors import colorConverter

from common import output_fig

# Figure tweaks for LaTeX from matplotlib cookbook
# http://wiki.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
FIG_WIDTH_PT = 746.25327  # Get this from LaTeX using \showthe\columnwidth
INCHES_PER_PT = 1.0/72.27               # Convert pt to inch
GOLDEN_MEAN = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
FIG_WIDTH = FIG_WIDTH_PT*INCHES_PER_PT  # width in inches
FIG_HEIGHT = FIG_WIDTH*GOLDEN_MEAN      # height in inches
FIG_SIZE = [FIG_WIDTH, FIG_HEIGHT]
PARAMS = {'backend': 'ps',
          'lines.linewidth': 2,
          'font.size': 20,
          'axes.labelsize': 20,
          'axes.titlesize': 20,
          'text.fontsize': 20,
          'legend.fontsize': 20,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': True,
          'figure.figsize': FIG_SIZE}
py.rcParams.update(PARAMS)

def get_vev(crit_temp, temps, results):
    """ Get vev data from results in a form suitable for plotting. """
    # Combine both sets of results for the unimproved vev to cover the full
    # temps range with a single variable
    _no_si_broken_vev1 = results[:, 1]
    _no_si_broken_vev2 = results[:, 12]
    no_si_broken_vev = np.where(np.isfinite(_no_si_broken_vev1), _no_si_broken_vev1,
                                _no_si_broken_vev2)

    no_si_sym_vev = sc.zeros_like(results[:, 4])  # symmetric phase vev identically zero
    no_si_sym_vev[temps < crit_temp] = sc.nan  # no symmetric phase below the critical temperature

    return no_si_broken_vev, no_si_sym_vev

def get_higgs(crit_temp, temps, results):
    """ Get Higgs data from results in a form suitable for plotting. """
    # Combine both sets of results for the unimproved Higgs mass to cover the full
    # temps range with a single variable
    _no_si_broken_mn1 = results[:, 3]
    _no_si_broken_mn2 = results[:, 14]
    no_si_broken_mn = np.where(np.isfinite(_no_si_broken_mn1), _no_si_broken_mn1, _no_si_broken_mn2)

    _no_si_sym_mn1 = results[:, 5]
    _no_si_sym_mn2 = results[:, 16]
    no_si_sym_mn = np.where(np.isfinite(_no_si_sym_mn1), _no_si_sym_mn1, _no_si_sym_mn2)
    no_si_sym_mn[temps < crit_temp] = sc.nan  # no symmetric phase below the critical temperature

    return no_si_broken_mn, no_si_sym_mn

def get_goldstone(results):
    """ Get Goldstone data from results in a form suitable for plotting. """
    # Combine both sets of results for the unimproved Goldstone mass to cover the full
    # temps range with a single variable
    _no_si_broken_mg1 = results[:, 2]
    _no_si_broken_mg2 = results[:, 13]
    no_si_broken_mg = np.where(np.isfinite(_no_si_broken_mg1), _no_si_broken_mg1, _no_si_broken_mg2)

    return no_si_broken_mg

def gap_eq_plots(crit_temp, results, save_figs):
    """ Make various plots of gap equation solutions for thesis chapter 3. """
    temps = results[:, 0]

    # VEV
    no_si_broken_vev, no_si_sym_vev = get_vev(crit_temp, temps, results)

    # Higgs mass
    no_si_broken_mn, no_si_sym_mn = get_higgs(crit_temp, temps, results)

    # Goldstone mass
    no_si_broken_mg = get_goldstone(results)

    py.figure()
    py.plot(temps, no_si_broken_vev, 'k-.', temps, no_si_sym_vev, 'k--')
    py.plot(temps, no_si_broken_mn, 'g-')
    py.plot(temps, no_si_broken_mg, 'b-')
    py.plot(temps, no_si_sym_mn, 'r:')
    py.axis([-1, 250, -1, 550])
    py.vlines(crit_temp, 0, 550, color=colorConverter.to_rgba('0.6', alpha=0.6))
    py.legend(('$v$(broken)', '$v$(symmetric)', '$m_H$(broken)', '$m_G$(broken)',
               '$m_G=m_H$(symmetric)'), loc='best')
    py.xlabel('T (MeV)')
    py.ylabel('$v,m_G,m_H$ (MeV)')
    py.grid(False)

    py.annotate('Critical Temperature', xy=(crit_temp, 250), xycoords='data',
                xytext=(150, 0), textcoords='offset points', ha='center',
                size=PARAMS['text.fontsize'],
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                connectionstyle="angle3,angleA=90,angleB=180"),
               )

    output_fig("twopi-hf-soln.pdf", save_figs, base_path=OUTPUTDIR, dpi=144)

if __name__ == "__main__":
    SAVE_FIGS = False  # whether to show figures interactively or save to file
    OUTPUTDIR = "../images/ch3"

    VB2 = 93. ** 2 # Pi-Sigma model
    MNB2 = 500. ** 2
    LAM = 3. * MNB2 / VB2
    N = 4
    CRIT_TEMP = np.sqrt(12. * VB2 / (N + 2.))

# Numpy array file containing the results array
# see help(sym_imp_scalars_hartree.driver) for format details, or:
#    results[:, 0] = temps
#    results[:, 1] = vev computed using asymmetric branch of unimproved EOM
#    results[:, 2] = Goldstone mass computed using same
#    results[:, 3] = Higgs mass computed using same
#    results[:, 4] = Goldstone mass computed using symmetric branch of unimproved
#    results[:, 5] = Higgs mass computed using same
#    results[:, 6] = vev computed using Pilaftsis & Teresi symmetry improved EOM
#    results[:, 7] = Goldstone mass computed using same
#    results[:, 8] = Higgs mass computed using same
#    results[:, 9] = vev computed using 3PI symmetry improved (Hartree-Fock) EOM asymmetric branch
#    results[:, 10]= Goldstone mass computed using same
#    results[:, 11]= Higgs mass computed using same
#    results[:, 12]= vev computed using 2PI unimproved asymmetric branch root finder
#    results[:, 13]= Goldstone mass computed using same
#    results[:, 14]= Higgs mass computed using same
#    results[:, 15]= Goldstone mass computed using 2PI unimproved symmetric branch root finder
#    results[:, 16]= Higgs mass computed using same
#    results[:, 17]= vev computed using 3PI symmetry improved (Hartree-Fock) EOM symmetric branch
#    results[:, 18]= Goldstone mass computed using same
#    results[:, 19]= Higgs mass computed using same
    RESULTS = sc.load("hartree-results-03102014.npy")

    gap_eq_plots(CRIT_TEMP, RESULTS, SAVE_FIGS)

