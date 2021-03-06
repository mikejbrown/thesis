# -*- coding: utf-8 -*-
"""
Common routines

Created on Tue Aug 26 11:32:25 2014

@author: Michael Brown
"""

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import interpolate


class PrecomputeInterpolate:
    """ Interpolating decorator for hartree_integral (or any function that
    returns a list of tuples [(val1, err1), (val2, err2), ...]).

    Adding "@PrecomputeInterpolate(xmin=a, xmax=b, numpoints=c)" above a
    function definition wraps that function so that a range of c values in
    the interval [a,b] (inclusive) are precomputed and the return value
    based on an interpolation of the function over that range. Values outside
    of the precomputed range are computed by the original function.

    This speeds up functions which are repeatedly evaluated in a certain
    interval, and accuracy is preserved if the function is sufficiently smooth
    on the scale dx = (xmax - xmin)/numpoints.

    Defaults: xmin = 0, xmax = 100, numpoints = 1000"""
    def __init__(self, xmin=0, xmax=100., numpoints=1000):
        self.__xmin = xmin
        self.__xmax = xmax
        self.__numpoints = numpoints
        self.__xgrid = sc.linspace(xmin, xmax, numpoints, endpoint=True)

    def __call__(self, func):
        import time
        print("Wrapping %s with PrecomputeInterpolate" % func)
        start = time.clock()
        _tmp = np.array([func(x) for x in self.__xgrid])
        stop = time.clock()
        print("Took %f seconds." % (stop - start))
        self.__func = func
        self.__vals = _tmp[:, 0]
        self.__errs = _tmp[:, 1]
        self.__interp_vals = interpolate.InterpolatedUnivariateSpline(
            self.__xgrid, self.__vals)
        self.__interp_errs = interpolate.InterpolatedUnivariateSpline(
            self.__xgrid, self.__errs)

        def _wrapped_f(x):
            if self.__xmin <= x <= self.__xmax:
                return (self.__interp_vals(x), self.__interp_errs(x))
            else:
                return self.__func(x)

        _wrapped_f.__name__ = func.__name__
        _wrapped_f.__doc__ = func.__doc__
        _wrapped_f.__repr__ = func.__repr__
        return _wrapped_f


@PrecomputeInterpolate(numpoints=5000)
def hartree_integral(mass):
    """
    Calculate the Hartree integral
    Integrate[(4 pi x^2)/((2 pi)^3 omega[x] (exp(omega[x])-1)), {x,0,Inf}]
    where
    omega[x] = sqrt(x^2+mass^2)

    Returns a tuple (value, estimated_error).
    """
    if mass == 0:  # exact value for massless case
        return 1. / 12, 0.0
    _omega = lambda x: np.sqrt(x ** 2. + mass ** 2.)
    _pref = 1. / (2 * sc.pi ** 2)
    _integrand = lambda x: x ** 2. / (_omega(x) * (sc.exp(_omega(x)) - 1.))
    _res, _err = quad(_integrand, 0, sc.inf)
    return _pref * _res, _pref * _err


def thermal_tadpole(mass, temp, mu):
    """ Thermal tadpole graph (Hartree-Fock) with mass mass, temperature temp and
    MS-bar renormalization point mu, all in the same units. """
#    assert mass >= 0
    assert temp >= 0
    assert mu > 0

    if mass <= 0:
        return (temp ** 2) / 12.
    if temp == 0:
        return (mass ** 2) * 2 * np.log(mass / mu) / (16 * sc.pi ** 2)

    return ((temp ** 2) * hartree_integral(mass / temp)[0]
            + (mass ** 2) * 2 * np.log(mass / mu) / (16 * sc.pi ** 2))

def output_fig(name, save_figs, base_path="images", **kwargs):
    """ If save_figs == True, saves the figure as a file of name 'base_path/name'.
    If save_figs == False, display the figure interactively instead.
    kwargs are passed onto matplotlib.pyplot.savefig unmodified"""
    if save_figs:
        import os
        fname = os.path.realpath(os.path.join(base_path, name))
        dname = os.path.dirname(fname)
        if not (os.path.exists(dname) and os.path.isdir(dname)):
            os.mkdir(dname)
        try:
            plt.savefig(fname, **kwargs)
        except:
            raise RuntimeError("Could not save image {0}".format(fname))
    else:
        plt.show()
