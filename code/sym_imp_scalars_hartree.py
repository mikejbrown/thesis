# -*- coding: utf-8 -*-
"""
Solve several variations of the 2PI gap equations for a O(N) scalar field
theory at the Hartree-Fock level.

Currently implemented:

- Regular, unimproved 2PI: shows a 1st order phase transition
    - have an iterative solver (no_si_solution),
    - and a scipy.optimize.root based solver (no_si_solution_root) - not as good
- Symmetry improved 2PI ala Pilaftsis & Teresi: 2nd order phase transition
    - have an iterative solver (si_2pi_hartree_solution)
- Symmetry improved 2PI imposing also the vertex Ward identity: shows a
  1st order transition
    - have an iterative solver (si_3pi_hartree_solution)

Support routines:

no_si_rhs: right hand side of gap equations

si_3pi_hartree_rhs: right hand side of gap equations

driver: main computation loop
    - The main body essentially just sets up some parameters, runs driver,
      and plots the results.

Modified for reproducibility testing purposes on Sep 10 2016
Created on Mon Aug 04 16:40:28 2014

@author: Michael Brown
"""
import time
import numpy as np
import scipy as sc

from common import thermal_tadpole

# MS-bar renormalisation point (MeV)
MU = 500.

def no_si_rhs(v2, mg2, mn2, temp, vb2, lam, N, symmetric_branch=False):
    """
    TODO: Currently breaks near the critical temperature!
    Right hand side of the equations of motion for the 2PIEA without symmetry
    improvement.

    Arguments:
    v2   = the Higgs vev squared
    mg2  = the Goldstone mass squared
    mn2  = the Higgs mass squared
    temp    = the temperature
    vb2  = the Higgs vev squared at zero temperature
    lam  = the quartic coupling constant at zero temperature
    N    = the number of field components (i.e., O(N) symmetry)
    symmetric_branch = Boolean indicating whether to compute for the symmetric
            or broken symmetry phase. Defaults to False. Note that this routine
            does not do the checking to make sure that the specified phase
            exists for the given parameter values.

    The units do not matter as long as they are all the same, e.g. GeV.

    Returns:
    A tuple (v2, mg2, mn2) of updated estimates, suitable for iteration.
    """
    assert v2 >= 0
    assert mg2 >= 0
    assert mn2 >= 0
    assert temp >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"

    if temp > 0:
        tg = thermal_tadpole(sc.sqrt(mg2), temp, MU)
        tn = thermal_tadpole(sc.sqrt(mn2), temp, MU)
    else:
        tg = 0
        tn = 0

    if symmetric_branch:
        return (0.,
                (lam / 6.) * (-vb2 + (N + 1) * tg + tn),
                (lam / 6.) * (-vb2 + (N + 1) * tg + tn))  # mn == mg in sym.
    else:
        return (vb2 - (N - 1) * tg - 3. * tn,
                lam * (tg - tn) / 3.,
                lam * v2 / 3.)


def no_si_solution(temp, vb2, lam, N, maxsteps=20, tol=1e-3,
                   symmetric_branch=False, update_weight=0.2):
    """
    Solution of the equations of motion for the 2PIEA without symmetry
    improvement.

    Arguments:
    temp    = the temperature
    vb2  = the Higgs vev squared at zero temperature
    lam  = the quartic coupling constant at zero temperature
    N    = the number of field components (i.e., O(N) symmetry)
    maxsteps = maximum number of iteration steps (defaults to 20)
    tol  = the absolute tolerance goal for the solution (component-wise)
    symmetric_branch = Boolean indicating whether to compute for the symmetric
            or broken symmetry phase. Defaults to False. Note that this routine
            does not do the checking to make sure that the specified phase
            exists for the given parameter values. If it does not exists this
            will likely fail with an AssertionError, or return with a
            RuntimeWarning complaining that maxsteps has been reached without
            convergence to a solution.
    update_weight = the weighting given to new updates in the iterative
            procedure (defaults to 0.2, must be between 0 and 1)

    The units do not matter as long as they are all the same, e.g. GeV.

    Returns:
    A tuple (v2, mg2, mn2) of
    v2  = the Higgs vev squared
    mg2 = the Goldstone mass squared
    mn2 = the Higgs mass squared
    """
    assert temp >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    assert maxsteps > 1
    assert np.floor(maxsteps) == maxsteps, "maxsteps must be an integer"
    assert 0 < update_weight <= 1.0

#    if symmetric_branch:
#        guess = (vb2, lam * vb2 / 3., lam * vb2 / 3.)
#    else:
    guess = (vb2, 0.0, lam * vb2 / 3.)

    for _ in range(maxsteps):
        new_guess = no_si_rhs(guess[0], guess[1], guess[2], temp, vb2, lam, N,
                              symmetric_branch=symmetric_branch)
        if symmetric_branch:
            new_guess = sc.array(new_guess)
            new_guess[new_guess < 0] = 0  # fix negative guesses
        if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
            return new_guess
        guess = (update_weight * sc.array(new_guess)
                 + (1. - update_weight) * sc.array(guess))

    return guess


def no_si_solution_root(temp, vb2, lam, N, maxsteps=20, tol=1e-3,
                        symmetric_branch=False, update_weight=0.2):
    """
    Solution of the equations of motion for the 2PIEA without symmetry
    improvement. Uses scipy.optimize.root instead of hand-rolled iteration.

    Arguments:
    temp    = the temperature
    vb2  = the Higgs vev squared at zero temperature
    lam  = the quartic coupling constant at zero temperature
    N    = the number of field components (i.e., O(N) symmetry)
    maxsteps = ( does nothing )
    tol  = the absolute tolerance goal for the solution (component-wise)
    symmetric_branch = Boolean indicating whether to compute for the symmetric
            or broken symmetry phase. Defaults to False. Note that this routine
            does not do the checking to make sure that the specified phase
            exists for the given parameter values. If it does not exists this
            will likely fail with an AssertionError, or return with a
            RuntimeWarning complaining that maxsteps has been reached without
            convergence to a solution.
    update_weight = the weighting given to new updates in the iterative
            procedure (defaults to 0.2, must be between 0 and 1)

    The units do not matter as long as they are all the same, e.g. GeV.

    Returns:
    A tuple (v2, mg2, mn2) of
    v2  = the Higgs vev squared
    mg2 = the Goldstone mass squared
    mn2 = the Higgs mass squared
    """
    from scipy.optimize import root
    assert temp >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    assert 0 < update_weight <= 1.0
    assert maxsteps >= 1

    guess = np.array((vb2, 20., lam * vb2 / 3.))  # mg2 = 20. just to avoid
                                                  # getting stuck at zero
    stepcount = 0

    if symmetric_branch:
        if temp < sc.sqrt(12. * vb2 / (N + 2.)):  # temp < crit_temp => no sym. phase
            return np.array([sc.nan, sc.nan, sc.nan])

        while stepcount < maxsteps:
            v2 = 0

            def residual(mn2, v2, temp, vb2, lam, N):
                """ Residual lhs - rhs for the mn2 gap equation"""
                if mn2 < 0:
                    mn2 = 0
                return mn2 - no_si_rhs(v2, mn2, mn2, temp, vb2, lam, N, symmetric_branch=True)[2]

            def jac(mn2, v2, temp, vb2, lam, N):
                """ Jacobian for the mn2 gap equation rootfinder """
                eps = 1e-12
                if mn2 <= 0:
                    mn2 = eps  # regulator for mg=0 case
                mn = sc.sqrt(mn2)
                return np.array([1.0 - (lam * (N + 2.) / (12. * mn)) *
                                 (thermal_tadpole(mn + eps, temp, MU)
                                  - thermal_tadpole(mn, temp, MU)) / eps, ])

            mn2 = root(residual, guess[2], args=(v2, temp, vb2, lam, N), jac=jac).x
            if mn2 < 0:
                mn2 = 0

            new_guess = np.array([v2, mn2, mn2])
            if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
                return new_guess
            guess = (update_weight * sc.array(new_guess)
                     + (1. - update_weight) * sc.array(guess))
    else:
        while stepcount < maxsteps:
            v2 = no_si_rhs(guess[0], guess[1], guess[2], temp, vb2, lam, N,
                           symmetric_branch=False)[0]
            mn2 = no_si_rhs(v2, guess[1], guess[2], temp, vb2, lam, N,
                            symmetric_branch=False)[2]

            def residual(mg2, v2, mn2, temp, vb2, lam, N):
                """ Residual lhs - rhs for the mg2 gap equation"""
                if mg2 < 0:
                    mg2 = 0
                return mg2 - no_si_rhs(v2, mg2, mn2, temp, vb2, lam, N,
                                       symmetric_branch=False)[1]

            def jac(mg2, v2, mn2, temp, vb2, lam, N):
                """ Jacobian for the mg2 gap equation rootfinder """
                eps = 1e-12
                if mg2 <= 0:
                    mg2 = eps  # regulator for mg=0 case
                mg = sc.sqrt(mg2)
                return np.array([1.0 - (lam / (6. * mg)) * (thermal_tadpole(mg + eps, temp, MU)
                                                            - thermal_tadpole(mg, temp, MU)) / eps, ])
            mg2 = root(residual, guess[1], args=(v2, mn2, temp, vb2, lam, N), jac=jac).x
            if mg2 < 0:
                mg2 = 0
            new_guess = np.array([v2, mg2, mn2])
            if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
                return new_guess
            guess = (update_weight * sc.array(new_guess)
                     + (1. - update_weight) * sc.array(guess))

    return guess


def si_2pi_hartree_solution(temp, vb2, lam, N,
                            maxsteps=20, tol=1e-3, update_weight=0.2):
    """
    Solution of the equations of motion for the 2PIEA with
    symmetry improvement. This is the same procedure used in:

    A. Pilaftsis and D. Teresi, Symmetry-Improved CJT Effective Action.
    Nucl. Phys. B 874, 594 (2013) http://arxiv.org/abs/1305.3221.

    Arguments:
    temp    = the temperature
    vb2  = the Higgs vev squared at zero temperature
    lam  = the quartic coupling constant at zero temperature
    N    = the number of field components (i.e., O(N) symmetry)
    maxsteps = maximum number of iteration steps (defaults to 20)
    tol  = the absolute tolerance goal for the solution (component-wise)
    symmetric_branch = whether or not to solve for the symmetric or broken
            symmetry phase solution (defaults to False)
    update_weight = the weighting given to new updates in the iterative
            procedure (defaults to 0.2, must be between 0 and 1)

    The units do not matter as long as they are all the same, e.g. GeV.

    Returns:
    A tuple (v2, mg2, mn2) of
    v2  = the Higgs vev squared
    mg2 = the Goldstone mass squared
    mn2 = the Higgs mass squared
    """
    assert temp >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    assert maxsteps > 1
    assert np.floor(maxsteps) == maxsteps, "maxsteps must be an integer"

    crit_T2 = 12. * vb2 / (N + 2.)
    mn2 = (lam * vb2 / 3.) * (1. - temp ** 2. / crit_T2)
    if temp == 0:
        return (vb2, 0., mn2)
    elif 0 < temp ** 2. <= crit_T2:
        # analytical solution exists!
        tn = thermal_tadpole(sc.sqrt(mn2), temp, MU)
        return ((3. * mn2 / lam) + ((temp ** 2.) / 12.) - tn, 0., mn2)
    else:
        # symmetric phase - same solution as the unimproved case
        return no_si_solution(temp, vb2, lam, N, maxsteps=maxsteps, tol=tol,
                              symmetric_branch=True,
                              update_weight=update_weight)


def si_3pi_hartree_rhs(v2, mg2, mn2, temp, vb2, lam, N, symmetric_branch=False):
    """
    Right hand side of the equations of motion for the 2PIEA at the
    Hartree-Fock level with 3PI symmetry improvement. More precisely,
    the vertex Ward identity is enforced in place of the Higgs equation
    of motion.

    Arguments:
    v2   = the Higgs vev squared
    mg2  = the Goldstone mass squared
    mn2  = the Higgs mass squared
    temp    = the temperature
    vb2  = the Higgs vev squared at zero temperature
    lam  = the quartic coupling constant at zero temperature
    N    = the number of field components (i.e., O(N) symmetry)

    The units do not matter as long as they are all the same, e.g. GeV.

    Returns:
    A tuple (v2, mg2, mn2) of updated estimates, suitable for iteration.
    """
    assert v2 >= 0
    assert mg2 >= 0
    assert mn2 >= 0
    assert temp >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"

    if temp > 0:
        tg = thermal_tadpole(sc.sqrt(mg2), temp, MU)
        tn = thermal_tadpole(sc.sqrt(mn2), temp, MU)
    else:
        tg = 0
        tn = 0

    if not symmetric_branch:
        # broken symmetry branch
        return (vb2 - ((N + 1.) * temp ** 2.) / 12. - tn,
                0.,
                lam * v2 / 3.)
    else:
        # symmetric branch
        return (0.,
                (lam / 6.) * (-vb2 + (N + 1) * tg + tn),
                (lam / 6.) * (-vb2 + (N + 1) * tg + tn))


def si_3pi_hartree_solution(temp, vb2, lam, N,
                            maxsteps=20, tol=1e-3, update_weight=0.2, symmetric_branch=False):
    """
    Solution of the equations of motion for the 2PIEA with 3PI symmetry
    improvement. More precisely, the vertex Ward identity is enforced in place
    of the Higgs equation of motion.

    Arguments:
    temp    = the temperature
    vb2  = the Higgs vev squared at zero temperature
    lam  = the quartic coupling constant at zero temperature
    N    = the number of field components (i.e., O(N) symmetry)
    maxsteps = maximum number of iteration steps (defaults to 20)
    tol  = the absolute tolerance goal for the solution (component-wise)
    symmetric_branch = whether or not to solve for the symmetric or broken
            symmetry phase solution (defaults to False)
    update_weight = the weighting given to new updates in the iterative
            procedure (defaults to 0.2, must be between 0 and 1)

    The units do not matter as long as they are all the same, e.g. GeV.

    Returns:
    A tuple (v2, mg2, mn2) of
    v2  = the Higgs vev squared
    mg2 = the Goldstone mass squared
    mn2 = the Higgs mass squared
    """
    assert temp >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    assert maxsteps > 1
    assert np.floor(maxsteps) == maxsteps, "maxsteps must be an integer"
    assert 0 < update_weight <= 1.0

#    if symmetric_branch:
#        return no_si_solution(temp, vb2, lam, N, maxsteps=maxsteps, tol=tol,
#                              symmetric_branch=True,
#                              update_weight=update_weight)

    crit_temp = sc.sqrt(12. * vb2 / (N + 2.))
    guess = (vb2, 0.0, lam * vb2 / 3.)

    for _ in range(maxsteps):
        new_guess = si_3pi_hartree_rhs(guess[0], guess[1], guess[2],
                                       temp, vb2, lam, N, symmetric_branch=symmetric_branch)
        new_guess = sc.array(new_guess)
        new_guess[new_guess < 0] = 0  # fix negative guesses

        if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
            # if all zeros and not at crit_temp this is a failure case
            # otherwise have converged on a solution
            if np.all(np.abs(sc.array(new_guess)) < 10 * tol) and not np.abs(temp - crit_temp) < tol:
                return (sc.nan, sc.nan, sc.nan)  # failure
            return new_guess  # success
        guess = (update_weight * sc.array(new_guess)
                 + (1. - update_weight) * sc.array(guess))

    return guess


def driver(vb2, lam, N, Ts, maxsteps=1000, weights=None):
    """
    Solve the Hartree-Fock 2PIEA equations of motion (EOM) for the given
    parameter values.

    Arguments:
    vb2     = the Higgs vev squared at zero temperature
    lam     = the quartic coupling constant at zero temperature
    N       = the number of field components (i.e., O(N) symmetry)
    Ts      = numpy array of temperatures at which to solve the EOM
    maxsteps= maximum number of steps in the iterative algorithm (default 1000)
    weights = numpy array of same shape as Ts with the update weights to use
              in the iterative algorithm. Defaults 0.2 * numpy.ones_like(Ts)

    Returns:
    results = numpy array with the shape (len(Ts), 20) with columns formatted
              as follows:
    results[:, 0] = Ts
    results[:, 1] = vev computed using asymmetric branch of unimproved EOM
    results[:, 2] = Goldstone mass computed using same
    results[:, 3] = Higgs mass computed using same
    results[:, 4] = Goldstone mass computed using symmetric branch of unimproved
    results[:, 5] = Higgs mass computed using same
    results[:, 6] = vev computed using Pilaftsis & Teresi symmetry improved EOM
    results[:, 7] = Goldstone mass computed using same
    results[:, 8] = Higgs mass computed using same
    results[:, 9] = vev computed using 3PI symmetry improved (Hartree-Fock) EOM asymmetric branch
    results[:, 10]= Goldstone mass computed using same
    results[:, 11]= Higgs mass computed using same
    results[:, 12]= vev computed using 2PI unimproved asymmetric branch root finder
    results[:, 13]= Goldstone mass computed using same
    results[:, 14]= Higgs mass computed using same
    results[:, 15]= Goldstone mass computed using 2PI unimproved symmetric branch root finder
    results[:, 16]= Higgs mass computed using same
    results[:, 17] = vev computed using 3PI symmetry improved (Hartree-Fock) EOM symmetric branch
    results[:, 18]= Goldstone mass computed using same
    results[:, 19]= Higgs mass computed using same
    """
    assert sc.all(Ts >= 0)
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    crit_temp = sc.sqrt(12. * vb2 / (N + 2.))

    if weights is None:
        weights = 0.2 * sc.ones_like(Ts)
    assert sc.all((weights > 0) & (weights <= 1.0)), "weights must all be between 0 and 1"

    _results = sc.nan * sc.zeros((len(Ts), 20))
    _results[:, 0] = Ts
    for i in range(len(Ts)):
        print("Step %d of %d. temp/Tc = %f" % (i, len(Ts), Ts[i]/crit_temp))
        temp = Ts[i]

        try:
            _v2, _mg2, _mn2 = no_si_solution(temp, vb2, lam, N,
                                             maxsteps=maxsteps,
                                             symmetric_branch=False,
                                             update_weight=weights[i])
            _results[i, 1:4] = (sc.sqrt(_v2), sc.sqrt(_mg2), sc.sqrt(_mn2))
        except AssertionError:
            pass

        try:
            _, _mg2s, _mn2s = no_si_solution(temp, vb2, lam, N,
                                             maxsteps=maxsteps,
                                             symmetric_branch=True,
                                             update_weight=weights[i])
            _results[i, 4:6] = (sc.sqrt(_mg2s), sc.sqrt(_mn2s))
        except AssertionError:
            pass

        try:
            si_v2, si_mg2, si_mn2 = si_2pi_hartree_solution(
                temp, vb2, lam, N, maxsteps=maxsteps, update_weight=weights[i])
            _results[i, 6:9] = (sc.sqrt(si_v2), sc.sqrt(si_mg2),
                                sc.sqrt(si_mn2))
        except AssertionError:
            pass

        try:
            si3_v2s, si3_mg2s, si3_mn2s = si_3pi_hartree_solution(
                temp, vb2, lam, N, maxsteps=maxsteps, update_weight=weights[i],
                symmetric_branch=False)
            _results[i, 9:12] = (sc.sqrt(si3_v2s), sc.sqrt(si3_mg2s),
                                 sc.sqrt(si3_mn2s))
        except AssertionError:
            pass

        try:
            _v2, _mg2, _mn2 = no_si_solution_root(temp, vb2, lam, N,
                                                  maxsteps=maxsteps,
                                                  symmetric_branch=False,
                                                  update_weight=weights[i])
            _results[i, 12:15] = (sc.sqrt(_v2), sc.sqrt(_mg2), sc.sqrt(_mn2))
        except AssertionError:
            pass

        try:
            _, _mg2s, _mn2s = no_si_solution_root(temp, vb2, lam, N,
                                                  maxsteps=maxsteps,
                                                  symmetric_branch=True,
                                                  update_weight=weights[i])
            _results[i, 15:17] = (sc.sqrt(_mg2s), sc.sqrt(_mn2s))
        except AssertionError:
            pass

        try:
            si3_v2s, si3_mg2s, si3_mn2s = si_3pi_hartree_solution(temp, vb2, lam, N,
                                                                  maxsteps=maxsteps,
                                                                  update_weight=weights[i],
                                                                  symmetric_branch=True)
            _results[i, 17:] = (sc.sqrt(si3_v2s), sc.sqrt(si3_mg2s),
                                sc.sqrt(si3_mn2s))
        except AssertionError:
            pass

    return _results

#%% Main (spyder cell magic)
if __name__ == "__main__":
    START = time.clock()
    # Vacuum vev and Higgs mass squared in GeV**2
    VB2 = 93. ** 2 # Pi-Sigma model
    MNB2 = 500. ** 2
    LAM = 3. * MNB2 / VB2
    _N = 4
    CRIT_TEMP = sc.sqrt(12. * VB2 / (_N + 2.))

    TEMPS = np.concatenate((sc.linspace(0, 130., 130),
                            sc.linspace(130., 140., 50),
                            sc.linspace(140., 200., 50),
                            sc.linspace(200., 225., 50),
                            sc.linspace(225., 250., 25)))

    # Stepped update_weights to focus computation near the critical temperature
    # where convergence is a bit dicey
    WEIGHTS = sc.zeros_like(TEMPS)
    WEIGHTS[TEMPS < CRIT_TEMP] = 0.9
    WEIGHTS[TEMPS >= CRIT_TEMP] = 0.1
    WEIGHTS[TEMPS > 1.5 * CRIT_TEMP] = 0.3

    # all the magic
    RESULTS = driver(VB2, LAM, _N, TEMPS, maxsteps=10000, weights=WEIGHTS)

    print("*****  Results  *****")
    print("Higgs model: vev = %.0f mH = %.0f N = %d lambda = %.3f CRIT_TEMP = %.2f" %
          (sc.sqrt(VB2), sc.sqrt(MNB2), _N, LAM, CRIT_TEMP))
    STOP = time.clock()
    print("Total time: %f" % (STOP - START))

#%% Reproducibility check (spyder cell magic)
    print("*****  Reproducibility Check  *****")

    # Numpy array file containing the RESULTS array.
    # This is the data file used to produce the thesis plots (TODO: which ones exactly?).
    # It was produced by running this code and saving manually using
    # np.save("filename.npy", RESULTS). Here it is used to check reproducibility.
    RESULTS_FILE = "hartree-results-03102014.npy"
    print("Loading results file: '%s'." % RESULTS_FILE)
    RESULTS_OLD = sc.load(RESULTS_FILE)

    # Check that RESULTS and RESULTS_OLD are the same shape
    print("Results are the same shape: %s" % (RESULTS.shape == RESULTS_OLD.shape))
    # Check that RESULTS and RESULTS_OLD are nans in all the same places
    print("Results are nans in all the same places: %s" %
          (sc.isnan(RESULTS) == sc.isnan(RESULTS_OLD)).all())
    # Check the rms residual where not nan
    RESID = sc.where(sc.isnan(RESULTS) | sc.isnan(RESULTS_OLD),
                     sc.zeros_like(RESULTS),
                     RESULTS - RESULTS_OLD)
    RMSRESID = sc.sqrt((RESID**2).sum())
    RESULTNORM = sc.sqrt((sc.where(sc.isnan(RESULTS),
                                   sc.zeros_like(RESULTS),
                                   RESULTS)**2).sum())
    print("RMS Residual = %f / %f == %f" % (RMSRESID, RESULTNORM, RMSRESID/RESULTNORM))
