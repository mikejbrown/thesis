# -*- coding: utf-8 -*-
"""
Solve several variations of the 2PI gap equations for a O(N) scalar field
theory at the Hartree-Fock level.

Currently implemented:

- Regular, unimproved 2PI: shows a 1st order phase transition
    - have an iterative solver (no_SI_solution),
    - and a scipy.optimize.root based solver (no_SI_solution_root) - not as good
- Symmetry improved 2PI ala Pilaftsis & Teresi: 2nd order phase transition
    - have an iterative solver (SI_2PI_Hartree_solution)
- Symmetry improved 2PI imposing also the vertex Ward identity: shows a
  1st order transition (!?)
    - have an iterative solver (SI_3PI_Hartree_solution)

Support routines:

no_SI_rhs: right hand side of gap equations

SI_3PI_Hartree_rhs: right hand side of gap equations

driver: main computation loop
    - The main body essentially just sets up some parameters, runs driver,
      and plots the results.

Created on Mon Aug 04 16:40:28 2014

@author: Michael
"""
from __future__ import division  # use floating point division by default on python versions < 3.0
import numpy as np
import scipy as sc
import time

from common import thermal_tadpole

mu = 500.

def no_SI_rhs(v2, mg2, mn2, T, vb2, lam, N, symmetric_branch=False):
    """
    TODO: Currently breaks near the critical temperature!
    Right hand side of the equations of motion for the 2PIEA without symmetry
    improvement.

    Arguments:
    v2   = the Higgs vev squared
    mg2  = the Goldstone mass squared
    mn2  = the Higgs mass squared
    T    = the temperature
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
    assert T >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"

    if T > 0:
        tg = thermal_tadpole(sc.sqrt(mg2), T, mu)
        tn = thermal_tadpole(sc.sqrt(mn2), T, mu)
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


def no_SI_solution(T, vb2, lam, N, maxsteps=20, tol=1e-3,
                   symmetric_branch=False, update_weight=0.2):
    """
    Solution of the equations of motion for the 2PIEA without symmetry
    improvement.

    Arguments:
    T    = the temperature
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
    import warnings
    assert T >= 0
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

    for i in range(maxsteps):
        new_guess = no_SI_rhs(guess[0], guess[1], guess[2], T, vb2, lam, N,
                              symmetric_branch=symmetric_branch)
        if symmetric_branch:
            new_guess = sc.array(new_guess)
            new_guess[new_guess < 0] = 0  # fix negative guesses
        if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
            return new_guess
        guess = (update_weight * sc.array(new_guess)
                 + (1. - update_weight) * sc.array(guess))
    warnings.warn("maxsteps reached: result may not have converged", RuntimeWarning)
    return guess


def no_SI_solution_root(T, vb2, lam, N, maxsteps=20, tol=1e-3,
                   symmetric_branch=False, update_weight=0.2):
    """
    Solution of the equations of motion for the 2PIEA without symmetry
    improvement. Uses scipy.optimize.root instead of hand-rolled iteration.

    Arguments:
    T    = the temperature
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
    import warnings
    assert T >= 0
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
        if T < sc.sqrt(12. * vb2 / (N + 2.)):  # T < crit_T => no sym. phase
            return np.array([sc.nan, sc.nan, sc.nan])

        while stepcount < maxsteps:
            v2 = 0

            def residual(mn2, v2, T, vb2, lam, N):
                if mn2 < 0:
                    warnings.warn("Found a negative mn2 in no_SI_solution: %f"%mn2, RuntimeWarning)
                    mn2 = 0
                return mn2 - no_SI_rhs(v2, mn2, mn2, T, vb2, lam, N, symmetric_branch=True)[2]
            def jac(mn2, v2, T, vb2, lam, N):
                eps = 1e-12
                if mn2 <= 0:
                    mn2 = eps  # regulator for mg=0 case
                mn = sc.sqrt(mn2)
                return np.array([1.0 - (lam * (N + 2.) / (12. * mn)) * (thermal_tadpole(mn + eps, T, mu)
                        - thermal_tadpole(mn, T, mu)) / eps, ])
            mn2 = root(residual, guess[2], args=(v2, T, vb2, lam, N), jac=jac).x
            if mn2 < 0:
                warnings.warn("Found a negative mn2 in no_SI_solution: %f"%mn2, RuntimeWarning)
                mn2 = 0

            new_guess = np.array([v2, mn2, mn2])
            if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
                return new_guess
            guess = (update_weight * sc.array(new_guess)
                     + (1. - update_weight) * sc.array(guess))
    else:
        while stepcount < maxsteps:
            v2 = no_SI_rhs(guess[0], guess[1], guess[2], T, vb2, lam, N, symmetric_branch=False)[0]
            mn2 = no_SI_rhs(v2, guess[1], guess[2], T, vb2, lam, N, symmetric_branch=False)[2]

            def residual(mg2, v2, mn2, T, vb2, lam, N):
                if mg2 < 0:
                    warnings.warn("Found a negative mg2 in no_SI_solution: %f"%mg2, RuntimeWarning)
                    mg2 = 0
                return mg2 - no_SI_rhs(v2, mg2, mn2, T, vb2, lam, N, symmetric_branch=False)[1]

            def jac(mg2, v2, mn2, T, vb2, lam, N):
                eps = 1e-12
                if mg2 <= 0:
                    mg2 = eps  # regulator for mg=0 case
                mg = sc.sqrt(mg2)
                return np.array([1.0 - (lam / (6. * mg)) * (thermal_tadpole(mg + eps, T, mu)
                        - thermal_tadpole(mg, T, mu)) / eps, ])
            mg2 = root(residual, guess[1], args=(v2, mn2, T, vb2, lam, N), jac=jac).x
            if mg2 < 0:
                    warnings.warn("Found a negative mg2 in no_SI_solution: %f"%mg2, RuntimeWarning)
                    mg2 = 0
            new_guess = np.array([v2, mg2, mn2])
            if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
                return new_guess
            guess = (update_weight * sc.array(new_guess)
                     + (1. - update_weight) * sc.array(guess))
    warnings.warn("maxsteps reached: result may not have converged", RuntimeWarning)
    return guess


def SI_2PI_Hartree_solution(T, vb2, lam, N,
                            maxsteps=20, tol=1e-3, update_weight=0.2):
    """
    Solution of the equations of motion for the 2PIEA with
    symmetry improvement. This is the same procedure used in:

    A. Pilaftsis and D. Teresi, Symmetry-Improved CJT Effective Action.
    Nucl. Phys. B 874, 594 (2013) http://arxiv.org/abs/1305.3221.

    Arguments:
    T    = the temperature
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
    assert T >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    assert maxsteps > 1
    assert np.floor(maxsteps) == maxsteps, "maxsteps must be an integer"

    crit_T2 = 12. * vb2 / (N + 2.)
    mn2 = (lam * vb2 / 3.) * (1. - T ** 2. / crit_T2)
    if T == 0:
        return (vb2, 0., mn2)
    elif 0 < T ** 2. <= crit_T2:
        # analytical solution exists!
#        print("Analytical solution")
        tn = thermal_tadpole(sc.sqrt(mn2), T, mu)
        return ((3. * mn2 / lam) + ((T ** 2.) / 12.) - tn, 0., mn2)
    else:
        # symmetric phase - same solution as the unimproved case
#        print("In symmetric phase")
        return no_SI_solution(T, vb2, lam, N, maxsteps=maxsteps, tol=tol,
                              symmetric_branch=True,
                              update_weight=update_weight)


def SI_3PI_Hartree_rhs(v2, mg2, mn2, T, vb2, lam, N, symmetric_branch=False):
    """
    Right hand side of the equations of motion for the 2PIEA at the
    Hartree-Fock level with 3PI symmetry improvement. More precisely,
    the vertex Ward identity is enforced in place of the Higgs equation
    of motion.

    Arguments:
    v2   = the Higgs vev squared
    mg2  = the Goldstone mass squared
    mn2  = the Higgs mass squared
    T    = the temperature
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
    assert T >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"

    if T > 0:
        tg = thermal_tadpole(sc.sqrt(mg2), T, mu)
        tn = thermal_tadpole(sc.sqrt(mn2), T, mu)
    else:
        tg = 0
        tn = 0

    if not symmetric_branch:
        # broken symmetry branch
        return (vb2 - ((N + 1.) * T ** 2.) / 12. - tn,
                0.,
                lam * v2 / 3.)
    else:
        # symmetric branch
        return (0.,
                (lam / 6.) * (-vb2 + (N + 1) * tg + tn),
                (lam / 6.) * (-vb2 + (N + 1) * tg + tn))


def SI_3PI_Hartree_solution(T, vb2, lam, N,
                            maxsteps=20, tol=1e-3, update_weight=0.2, symmetric_branch=False):
    """
    Solution of the equations of motion for the 2PIEA with 3PI symmetry
    improvement. More precisely, the vertex Ward identity is enforced in place
    of the Higgs equation of motion.

    Arguments:
    T    = the temperature
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
    import warnings
    assert T >= 0
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    assert maxsteps > 1
    assert np.floor(maxsteps) == maxsteps, "maxsteps must be an integer"
    assert 0 < update_weight <= 1.0

#    if symmetric_branch:
#        return no_SI_solution(T, vb2, lam, N, maxsteps=maxsteps, tol=tol,
#                              symmetric_branch=True,
#                              update_weight=update_weight)

    crit_T = sc.sqrt(12. * vb2 / (N + 2.))
    guess = (vb2, 0.0, lam * vb2 / 3.)

    for i in range(maxsteps):
        new_guess = SI_3PI_Hartree_rhs(guess[0], guess[1], guess[2],
                                       T, vb2, lam, N, symmetric_branch=symmetric_branch)
        new_guess = sc.array(new_guess)
        new_guess[new_guess < 0] = 0  # fix negative guesses
#        print(new_guess)
        if sc.all(abs(sc.array(new_guess) - sc.array(guess)) < tol):
            # if all zeros and not at crit_T this is a failure case
            # otherwise have converged on a solution
            if np.all(np.abs(sc.array(new_guess)) < 10 * tol) and not (np.abs(T - crit_T) < tol):
                return (sc.nan, sc.nan, sc.nan)  # failure
            return new_guess  # success
        guess = (update_weight * sc.array(new_guess)
                 + (1. - update_weight) * sc.array(guess))
    warnings.warn("maxsteps reached: result may not have converged", RuntimeWarning)
    return guess


def driver(vb2, mnb2, lam, N, Ts, maxsteps=1000, weights=None):
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
    import warnings
    assert sc.all(Ts >= 0)
    assert vb2 >= 0
    assert lam > 0
    assert N >= 1
    assert np.floor(N) == N, "N must be an integer"
    crit_T = sc.sqrt(12. * vb2 / (N + 2.))

    if weights == None:
        weights = 0.2 * sc.ones_like(Ts)
    assert sc.all((0 < weights) & (weights <= 1.0)), "weights must all be between 0 and 1"

    results = sc.nan * sc.zeros((len(Ts), 20))
    results[:, 0] = Ts
    for i in range(len(Ts)):
        print("Step %d of %d. T/Tc = %f" % (i, len(Ts), Ts[i]/crit_T))
        T = Ts[i]

        try:
            _v2, _mg2, _mn2 = no_SI_solution(T, vb2, lam, N,
                                             maxsteps=maxsteps,
                                             symmetric_branch=False,
                                             update_weight=weights[i])
            results[i, 1:4] = (sc.sqrt(_v2), sc.sqrt(_mg2), sc.sqrt(_mn2))
        except AssertionError:
            warnings.warn("Failed to converge on an asymmetric branch solution")

        try:
            _v2s, _mg2s, _mn2s = no_SI_solution(T, vb2, lam, N,
                                                maxsteps=maxsteps,
                                                symmetric_branch=True,
                                                update_weight=weights[i])
            results[i, 4:6] = (sc.sqrt(_mg2s), sc.sqrt(_mn2s))
        except AssertionError:
            warnings.warn("Failed to converge on a symmetric branch solution")

        try:
            si_v2, si_mg2, si_mn2 = SI_2PI_Hartree_solution(
                T, vb2, lam, N, maxsteps=maxsteps, update_weight=weights[i])
            results[i, 6:9] = (sc.sqrt(si_v2), sc.sqrt(si_mg2),
                               sc.sqrt(si_mn2))
        except AssertionError:
            warnings.warn("Failed to converge on a SI-2PI solution")

        try:
            si3_v2s, si3_mg2s, si3_mn2s = SI_3PI_Hartree_solution(
                T, vb2, lam, N, maxsteps=maxsteps, update_weight=weights[i], symmetric_branch=False)
            results[i, 9:12] = (sc.sqrt(si3_v2s), sc.sqrt(si3_mg2s),
                              sc.sqrt(si3_mn2s))
        except AssertionError:
            warnings.warn("Failed to converge on an asymmetric branch SI-3PI solution")

        try:
            _v2, _mg2, _mn2 = no_SI_solution_root(T, vb2, lam, N,
                                             maxsteps=maxsteps,
                                             symmetric_branch=False,
                                             update_weight=weights[i])
            results[i, 12:15] = (sc.sqrt(_v2), sc.sqrt(_mg2), sc.sqrt(_mn2))
        except AssertionError:
            warnings.warn("Failed to converge on an asymmetric branch solution")

        try:
            _v2s, _mg2s, _mn2s = no_SI_solution_root(T, vb2, lam, N,
                                                maxsteps=maxsteps,
                                                symmetric_branch=True,
                                                update_weight=weights[i])
            results[i, 15:17] = (sc.sqrt(_mg2s), sc.sqrt(_mn2s))
        except AssertionError:
            warnings.warn("Failed to converge on a symmetric branch solution")

        try:
            si3_v2s, si3_mg2s, si3_mn2s = SI_3PI_Hartree_solution(T, vb2, lam, N,
                maxsteps=maxsteps, update_weight=weights[i], symmetric_branch=True)
            results[i, 17:] = (sc.sqrt(si3_v2s), sc.sqrt(si3_mg2s),
                              sc.sqrt(si3_mn2s))
        except AssertionError:
            warnings.warn("Failed to converge on a symmetric branch SI-3PI solution")
#        print(results[i,:])

    return results

#%% Main (spyder cell magic)
if __name__ == "__main__":
    start = time.clock()
    # Vacuum vev and Higgs mass squared in GeV**2
    vb2  = 93. ** 2 # Pi-Sigma model
    mnb2 = 500. ** 2
    lam  = 3. * mnb2 / vb2
    N = 4
    crit_T = sc.sqrt(12. * vb2 / (N + 2.))

    Ts = np.concatenate((sc.linspace(0, 130., 130),
                         sc.linspace(130., 140., 50),
                         sc.linspace(140., 200., 50),
                         sc.linspace(200., 225., 50),
                         sc.linspace(225., 250., 25)))

    # Stepped update_weights to focus computation near the critical temperature
    # where convergence is a bit dicey
    weights = sc.zeros_like(Ts)
    weights[Ts < crit_T] = 0.9
    weights[Ts >= crit_T] = 0.1
    weights[Ts > 1.5 * crit_T] = 0.3

    # all the magic
    results = driver(vb2, mnb2, lam, N, Ts, maxsteps=10000, weights=weights)

    print("*****  Results  *****")
    print("Higgs model: vev = %.0f mH = %.0f N = %d lambda = %.3f crit_T = %.2f" % (
    sc.sqrt(vb2), sc.sqrt(mnb2), N, lam, crit_T))
    
#%% Reproducibility check (spyder cell magic)
    print("*****  Reproducibility Check  *****")
    
        # Numpy array file containing the results array
    # This is the data file used to produce the thesis plots (TODO: which ones exactly?)
    # It was produced by running this code and saving manually using np.save("filename.npy", results)
    # Here it is used to check reproducibility
    results_file = "hartree-results-03102014.npy"
    print("Loading results file: '%s'." % results_file)
    results_old = sc.load(results_file)
    
    # Check that results and results_old are the same shape
    print("Results are the same shape: %s" % (results.shape == results_old.shape))
    # Check that results and results_old are nans in all the same places
    print("Results are nans in all the same places: %s" % (sc.isnan(results) == sc.isnan(results_old)).all())
    # Check the rms residual where not nan
    resid = sc.where(sc.isnan(results) | sc.isnan(results_old), sc.zeros_like(results), (results - results_old))
    rmsresid = sc.sqrt((resid**2).sum())
    resultnorm = sc.sqrt((sc.where(sc.isnan(results), sc.zeros_like(results), results)**2).sum())
    print("RMS Residual = %f / %f == %f" % (rmsresid, resultnorm, rmsresid/resultnorm))
    
    stop = time.clock()
    print("Total time: %f" % (stop - start))
