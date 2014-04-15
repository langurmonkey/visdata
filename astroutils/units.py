"""
Small submodule that provides convenient unit conversion interface.
:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1

Requirements
------------
* `Astropy 0.3+ <http://www.astropy.org>`_
"""
from __future__ import division

import astropy.units as u


def convert(unit0, unit1, value):
    """ Converts the given value from unit0 to uni1 using astropy's unit conversion. """
    if unit0 == unit1:
        return value
    else:
        u0 = u.Unit(unit0)
        u1 = u.Unit(unit1)
        val = value * u0
        return val.to(u1).value
