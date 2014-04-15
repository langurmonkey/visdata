visdata
=======

Python module and package to fetch and process astrometric data.
The VisData module aims to provide a flexible solution for the extraction and transformation
of astronomical data. This tool provides a fast and straightforward way of getting the right
data in the right format without having to write a single line of code.
Even though it was initially intended for the fetching of Gaia data
for visualisation purposes, it can be used for any other tasks.

Introduction
============

The original aim of this package was to aid in the preparation of data from ESA's Gaia
space mission for their visualisation in planetarium domes. It started as a straightforward module
which loaded some data from a database using the TAP client in GAVO's VOTable python 
package, transformed it and wrote it to an output file. Then, in order to make it a bit more flexible
it grew larger and larger until one could basically specify the desired output columns
in any desired units and the module would do all the conversions and transformations
for you. In the process, the astroutils package was created. It is a collection of utilities and functions
that do common astronomical calculations, conversions and transformations. 

More information
================

For further information about how to use this module or the astroutils package please refer to the
documentation files under /doc.
