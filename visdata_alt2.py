"""
Utility to extract data from an archive, transform it and output it in the desired format.
:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1

"""

from __future__ import division

import logging
import math
import os.path
import shutil
import time
import urlparse

import astroutils.coordinatescore as coordinateutils
import astroutils.units as units
import commonio
import visdataconfig as conf
from visdatalog import logger, _fmt


def strip_suffix(string, suffix):
    """ Strips the given suffix from the given string, if it exists, and returns it """
    return string[:len(string) - len(suffix)] if string.endswith(suffix) else string


def run(configfile):
    """ Runs the program as configured in the given configuration file.

        This is the main method to run the functionalitiy in this package. This
        loads the configuration file. Using the configuration, it loads the appropriate data,
        transforms it to the output format specified in the config file and writes it
        in a file.
    """
    conf.init_config(configfile)

    exec_timestamp = time.strftime("%d%m%Y-%H.%M.%S")

    # Configure logger to file
    hdlr = logging.FileHandler(os.path.join(conf.get_log_path(), exec_timestamp + ".log"))
    hdlr.setLevel(logging.DEBUG)
    hdlr.setFormatter(_fmt)
    logger.addHandler(hdlr)

    # Prepare data loader and writerrefsys
    data_source = conf.get_data_source()
    parts = urlparse.urlsplit(data_source)
    minus = -1
    if parts.scheme and parts.netloc:
        # accesURL is apparently an URL, a TAP service
        dataloader = commonio.TAPLoader(conf.get_data_source(), conf.get_query())
    elif parts.path:
        # data_source is apparently a file
        try:
            f = open(data_source, 'r')
            f.close()
        except OSError:
            # handle error here
            logger.error("Input file in data_source not valid: %s" % data_source)
            raise Exception("Input file in data_source not valid: %s" % data_source)

        dataloader = commonio.AsciiLoader(conf.get_data_source())
        minus = 1
    else:
        logger.error("Could not determine the type of data_source: %s" % data_source)
        raise Exception("Could not determine the type of data_source: %s" % data_source)

    # START CORE PROCESSING
    # Initialize output lists
    colnames1 = conf.get_output_columns()
    units1 = conf.get_output_units()
    fields1 = []

    # Create fields1
    for (idx, out_name) in enumerate(colnames1):
        columnfield = commonio.field(out_name, units1[idx])

        fields1.append(columnfield)

    n = len(dataloader) - minus
    logger.info("%i records found in data source: %s" % (n, data_source))

    # Prepare output writer
    try:
        outfilepath = conf.get_output_file_path()
    except:
        logger.error("Property not defined: outfile")
        return

    filename = exec_timestamp + "-%i" % n
    outfile = os.path.join(outfilepath, filename + ".txt")
    datawriter = commonio.AsciiWriter(outfile)

    try:
        # Get/update fields0
        fields0 = dataloader.fields
        if(any(fld.unit is None for fld in fields0)):
            raise Exception("Missing units in data_source: %s" % dataloader.data_source)

        # Open and write header
        datawriter.open()
        datawriter.write_header(fields1)

        granularity = max(10, 10 ** (math.floor(math.log10(n)) - 3))

        # Loop over data, which should use generators in the case of files
        for idx, line in enumerate(dataloader.get_data()):
            if (idx + 1) % (n // granularity + 1) == 0 or idx == 0:
                logger.info("Processing line %i of %i - %3.1f%% (%.3f seconds)" % (idx + 1, n, ((idx + 1) * 100 / n), time.clock()))

            # Process line
            alpha = units.convert('deg', 'rad', line[0])
            delta = units.convert('deg', 'rad', line[1])
            distance = line[2]
            mualpha = units.convert('mas/a', 'rad/s', line[3])
            mudelta = units.convert('mas/a', 'rad/s', line[4])
            radvel = line[5]

            cartesian_all = coordinateutils.spherical_to_cartesian_pm([alpha, delta, distance, mualpha, mudelta, radvel])
            cartesian_coordinates = coordinateutils.equatorial_to_galactic_v(cartesian_all[:3])
            cartesian_pm = coordinateutils.equatorial_to_galactic_v(cartesian_all[3:])

            current_row = cartesian_coordinates + cartesian_pm + line[6:]

            # Add row to matrix
            datawriter.write_line(current_row)

        # Copy configuration file to output path
        outconffile = os.path.join(outfilepath, filename + ".ini")
        logger.info("Writing metadata to file %s" % outconffile)
        shutil.copyfile(configfile, outconffile)
    except Exception, e:
        logger.exception(e)
        datawriter.delete_file()
        logger.error("An error occurred, output file removed: %s" % datawriter.data_source)
    finally:
        # Close data writer
        datawriter.close()
        logger.info("%i records written to file %s" % (datawriter.records, datawriter.data_source))

    logger.info("Program finished in %.3f seconds." % (time.clock()))
