visdata
=======

Python module and package to fetch and process astrometric data.
The visdata module aims to provide a flexible solution for the extraction and transformation
of astronomical data. This tool provides a fast and straightforward way of getting the right
data in the right format without having to write a single line of code.
Even though it was initially intended for the fetching of Gaia data
for visualisation purposes, it can be used for any other tasks.
It includes and uses the astroutils module, which does the 'dirty' work.

Introduction
============

This started as a straightforward module
which loaded some data from a database using the TAP client in GAVO's VOTable python 
package, transformed it and wrote it to an output file. Then, in order to make it a bit more flexible
it grew larger and larger until one could basically specify the desired output columns
in any desired units and the module would do all the conversions and transformations
for you. In the process, the astroutils package was created. It is a collection of utilities and functions
that do common astronomical calculations, conversions and transformations. 

Requirements
============

The visdata module uses the astroutils package, which was developed at the same time for this sole purpose.
The astroutils package itself uses a set of third-party libraries that have proven very reliable and useful.

-Christoph Gohlke's transformations library. A library providing matrix transformations and quaternions handling. This is already included in the astroutils package (http://www.lfd.uci.edu/gohlke/).

-Astropy's unit conversion (http://www.astropy.org).

-Numpy 1.7 (http://www.numpy.org).

-VOTable by GAVO (http://vo.ari.uni-heidelberg.de/soft/subpkgs}{http://vo.ari.uni-heidelberg.de/soft/subpkgs).


More information
================

For further information about how to use this module or the astroutils package please refer to the
documentation files under /doc.
