# -*- coding: utf-8 -*-

try:
    __SIGGEN_SETUP__
except NameError:
    __SIGGEN_SETUP__ = False

if not __SIGGEN_SETUP__:
    __all__ = ["Detector"]

    from .Detector import Detector
    from .geometries.PPC import PPC
    from .geometries.ICPC import ICPC
    from .geometries.GEM import GEM
