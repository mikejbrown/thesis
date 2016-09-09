# -*- coding: utf-8 -*-
"""

Generate plots for thesis ch 4 on solution of SI-2PI HF gap equations

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

def gap_eq_plots(vb2, N, save_figs):
    """ Make various plots of gap equation solutions for thesis chapter 4. """
    # Numpy array file containing the results array
    # see help(sym_imp_scalars_hartree.driver) for format details, or:
#    results[:, 0] = Ts
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
#    results[:, 17] = vev computed using 3PI symmetry improved (Hartree-Fock) EOM symmetric branch
#    results[:, 18]= Goldstone mass computed using same
#    results[:, 19]= Higgs mass computed using same
    results_file = "hartree-results-03102014.npy"

    results = sc.load(results_file)


    crit_T = np.sqrt(12. * vb2 / (N + 2.))

    Ts = results[:, 0]

    # VEV

    # Combine both sets of results for the unimproved vev to cover the full
    # Ts range with a single variable
    _no_si_broken_vev1 = results[:, 1]
    _no_si_broken_vev2 = results[:, 12]
    no_si_broken_vev = np.where(np.isfinite(_no_si_broken_vev1), _no_si_broken_vev1,
                                _no_si_broken_vev2)

    no_si_sym_vev = sc.zeros_like(results[:, 4])  # symmetric phase vev identically zero
    no_si_sym_vev[Ts < crit_T] = sc.nan  # no symmetric phase below the critical temperature

    pt_vev = results[:, 6]

    #si3pi_broken_vev = results[:, 9]
    si3pi_sym_vev = results[:, 17]

    # Higgs mass

    # Combine both sets of results for the unimproved Higgs mass to cover the full
    # Ts range with a single variable
    _no_si_broken_mn1 = results[:, 3]
    _no_si_broken_mn2 = results[:, 14]
    no_si_broken_mn = np.where(np.isfinite(_no_si_broken_mn1), _no_si_broken_mn1, _no_si_broken_mn2)

#    _no_si_sym_mn1 = results[:, 5]
#    _no_si_sym_mn2 = results[:, 16]
#    no_si_sym_mn = np.where(np.isfinite(_no_si_sym_mn1), _no_si_sym_mn1, _no_si_sym_mn2)

    pt_mn = results[:, 8]

    #si3pi_broken_mn = results[:, 11]
    si3pi_sym_mn = results[:, 19]
    si3pi_sym_mn[Ts < crit_T] = sc.nan  # no symmetric phase below the critical temperature

    # Goldstone mass

    # Combine both sets of results for the unimproved Goldstone mass to cover the full
    # Ts range with a single variable
    _no_si_broken_mg1 = results[:, 2]
    _no_si_broken_mg2 = results[:, 13]
    no_si_broken_mg = np.where(np.isfinite(_no_si_broken_mg1), _no_si_broken_mg1, _no_si_broken_mg2)

#    _no_si_sym_mg1 = results[:, 4]
#    _no_si_sym_mg2 = results[:, 15]
#    no_si_sym_mn = np.where(np.isfinite(_no_si_sym_mn1), _no_si_sym_mn1, _no_si_sym_mn2)

    pt_mg = results[:, 7]

    #si3pi_broken_mg = results[:, 10]
    si3pi_sym_mg = results[:, 18]
    si3pi_sym_mg[Ts < crit_T] = sc.nan  # no symmetric phase below the critical temperature


    py.figure()
    py.plot(Ts, no_si_broken_vev, 'k-', Ts, pt_vev, 'b-.',
            Ts, si3pi_sym_vev, 'k--')
    py.vlines(crit_T, 0, 550, color=colorConverter.to_rgba('0.6', alpha=0.6))
    py.axis([-1, 250, -1, 100])
    py.legend(('2PI',
               'SI-2PI',
               'symmetric phase (all)'), loc='best')
    py.xlabel('T (MeV)')
    py.ylabel('$v$ (MeV)')
    py.grid(False)
    py.annotate('Unphysical metastable phase\n(1st order phase transition)',
                xy=(207, 60), xycoords='data',
                xytext=(-42, -75), textcoords='offset points', ha='center',
                size=PARAMS['text.fontsize'],
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                connectionstyle="angle3,angleA=90,angleB=0"),
               )
    py.annotate('2nd order phase transition', xy=(crit_T, 15), xycoords='data',
                xytext=(30, 0), textcoords='offset points',
                size=PARAMS['text.fontsize'],
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                connectionstyle="angle3,angleA=90,angleB=0"),
               )
    output_fig("si2pi-vev-comparison.pdf", save_figs, base_path=outputdir, dpi=144)

    py.figure()
    py.plot(Ts, no_si_broken_mn, 'k-', Ts, pt_mn, 'b-.',
            Ts, si3pi_sym_mn, 'k--')
    py.vlines(crit_T, 0, 550, color=colorConverter.to_rgba('0.6', alpha=0.6))
    py.axis([-1, 250, -1, 550])
    py.legend(('2PI',
               'SI-2PI',
               'symmetric phase (all)'), loc='best')
    py.xlabel('T (MeV)')
    py.ylabel('$m_H$ (MeV)')
    py.grid(False)
    py.annotate('Critical Temperature\n$T_c = \\sqrt{\\frac{12 v^2}{N + 2}}$',
                xy=(crit_T, 0), xycoords='data',
                xytext=(150, 125), textcoords='offset points', ha='center',
                size=PARAMS['text.fontsize'],
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                connectionstyle="angle3,angleA=0,angleB=45"),
               )
    output_fig("si2pi-higgs-mass-comparison.pdf", save_figs, base_path=outputdir, dpi=144)

    py.figure()
    py.plot(Ts, no_si_broken_mg, 'k-', Ts, pt_mg, 'b-.',
            Ts, si3pi_sym_mg, 'k--')
    py.vlines(crit_T, 0, 550, color=colorConverter.to_rgba('0.6', alpha=0.6))
    py.axis([-1, 250, -1, 200])
    py.legend(('2PI',
               'SI-2PI',
               'symmetric phase (all)'), loc='best')
    py.xlabel('T (MeV)')
    py.ylabel('$m_G$ (MeV)')
    py.grid(False)
    py.annotate('Goldstone theorem violated', xy=(100, 90), xycoords='data',
                xytext=(-200, 45), textcoords='offset points',
                size=PARAMS['text.fontsize'],
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                connectionstyle="angle3,angleA=0,angleB=-30"),
               )
    py.annotate('Goldstone theorem satisfied', xy=(100, 0), xycoords='data',
                xytext=(-110, 50), textcoords='offset points',
                size=PARAMS['text.fontsize'],
                #bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="fancy",
                                fc="0.6", ec="none",
                                connectionstyle="angle3,angleA=0,angleB=-30"),
               )
    output_fig("si2pi-goldstone-mass-comparison.pdf", save_figs, base_path=outputdir, dpi=144)



if __name__ == "__main__":
    save_figs = False  # whether to show figures interactively or save to file
    outputdir = "../images/ch4"

    vb2 = 93. ** 2 # Pi-Sigma model
    mnb2 = 500. ** 2
    lam = 3. * mnb2 / vb2
    N = 4
    crit_T = np.sqrt(12. * vb2 / (N + 2.))

    gap_eq_plots(vb2, N, save_figs)

