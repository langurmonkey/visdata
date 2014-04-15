"""
Astrometry data module. Loads, transforms and writes astronomy (and specially astrometry) data from and to different
sources.

:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1

"""
import logging

# LOG
""" Global logger """
logger = logging.getLogger("visdata")
logger.setLevel(logging.DEBUG)
_ch = logging.StreamHandler()
_ch.setLevel(logging.DEBUG)
_fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_ch.setFormatter(_fmt)
logger.addHandler(_ch)
