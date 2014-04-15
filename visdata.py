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

import astroutils.coordinates as coordinates
import astroutils.units as units
import commonio
import visdataconfig as conf

from visdatalog import logger, _fmt


# Definition of common coordinates
sph_refsys_definitions = ["equatorial", "equatorial_pm", "galactic", "galactic_pm", "ecliptic", "ecliptic_pm"]
cart_refsys_definitions = ["cartesian", "cartesian_pm"]
refsys_definitions = sph_refsys_definitions + cart_refsys_definitions

# SPHERICAL
equatorial = ["alpha", "delta", "distance"]
equatorial_pm = ["mualpha", "mudelta", "radialvelocity"]

galactic = ["l", "b", "distance"]
galactic_pm = ["mul", "mub", "radialvelocity"]

ecliptic = ["lambda", "beta", "distance"]
ecliptic_pm = ["mulambda", "mubeta", "radialvelocity"]

# CARTESIAN
cartesian = ["x", "y", "z"]
cartesian_pm = ["xdot", "ydot", "zdot"]

all_spherical_pos = galactic + ecliptic + equatorial
all_spherical_vel = galactic_pm + ecliptic_pm + equatorial_pm
all_spherical = all_spherical_pos + all_spherical_vel

names_spherical = ["ecliptic", "equatorial", "galactic", "ecliptic_pm", "equatorial_pm", "galactic_pm"]
names_cartesian = ["cartesian", "cartesian_pm"]

all_cartesian = cartesian + cartesian_pm

all_coordinates = all_spherical + all_cartesian


class AutoVivification(dict):
    """Implementation of perl's auto-vivification feature for python dictionaries."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def get_values_and_fields(values, fields, field_names):
    """ Returns a list of values and fields corresponding to the given field_names

        Returns a list of values and fields corresponding to the given field_names. The fields
        list must contain all the field_names.
    """
    out_fields = []
    out_values = []

    for field_name in field_names:
        try:
            index = fields.index(field_name)
            out_fields.append(fields[index])
            out_values.append(values[index])
        except:
            logger.error("One of the output fields (%s) can not be derived from input columns: %s" % (field_name, ', '.join(map(str, fields))))
            raise Exception("One of the output fields (%s) can not be derived from input columns: %s" % (field_name, ', '.join(map(str, fields))))

    return out_fields, out_values


def get_refsys(column, inout_refsys):
    """ Works out and returns the reference system of the given column as a string.

        If the column is a spherical coordinate, then the appropriate reference
        system is returned in the following fashion:
        - alpha, delta, ra, dec = equatorial
        - l, b = galactic
        - lambda, beta = ecliptic

        If the column is a cartesian coordinate, then the value of
        inout_refsys is returned. Inout_refsys should contain the reference system of
        cartesian coordinates.

        If the column does not represent a coordinate, then the string 'na' is returned.
    """
    if column in all_cartesian:
        # CARTESIAN
        return inout_refsys
    elif column in all_spherical:
        # SPHERICAL
        for sph_refsys in sph_refsys_definitions:
            if column in globals()[sph_refsys]:
                return strip_suffix(sph_refsys, "_pm")
    else:
        return "na"


def check_refsys(column, *refsystems):
    """ Checks the column name given is in any of the lists named by the strings in refsystems. """
    if (not any(column in globals()[refsys] for refsys in refsystems)):
        raise Exception("Column name '%s' does not match with input reference system: %s" % (column, refsystems))


def get_io_coord_info(in_fields, out_fields):
    """ Returns a dictionary with information on the in/out coordinates and velocities

        This method extracts meta-information from the given input and output fields lists. It will
        get the indeces of the input positions and velocities (if any) per reference system and the
        input and output position and velocities reference systems.
    """
    coord_info = AutoVivification()
    coord_info["in"], in_pos, in_vel = _get_io_coord_info(in_fields)
    coord_info["out"], out_pos, out_vel = _get_io_coord_info(out_fields)
    return coord_info, (in_pos, in_vel), (out_pos, out_vel)


def _get_io_coord_info(fields):
    """ Gets a dictionary with meta-information on the given fields

        It parses the fields and matches them agains the spherical and cartesian coordinates
        and proper motions, constructing a dictionary with the information of which are
        the reference systems in the fields, and at which index are their values. It also
        returns information on the position and proper motion reference systems contained
        in the fields.
    """
    coord_info = AutoVivification()
    pos = None
    vel = None
    for refsys in refsys_definitions:
        if all(elem in fields for elem in globals()[refsys]):
            indices = []
            for coordinate in globals()[refsys]:
                indices.append(fields.index(coordinate))
            coord_info[refsys] = indices

            if refsys.endswith("_pm"):
                vel = refsys
            else:
                pos = refsys

    return coord_info, pos, vel


def get_coordinates(coord_info, in_tuple, out_tuple, fields0, line, cart_in_refsys):
    """ Gets a dictionary with the coordinates information and returns the coordinates object """
    if out_tuple[1] is None:
        if out_tuple[0] is not None:
            # Only positions
            refsys = in_tuple[0] if in_tuple[0] != 'cartesian' else cart_in_refsys

            indices = coord_info["in"][in_tuple[0]]
            val1 = line[indices[0]]
            val2 = line[indices[1]]
            val3 = line[indices[2]]

            unit_values = (fields0[indices[0]].unit, fields0[indices[1]].unit, fields0[indices[2]].unit)

            coord_class = getattr(coordinates, refsys.title())

            if in_tuple[0] == 'cartesian':
                return coord_class(x=val1, y=val2, z=val3, unit=unit_values), refsys
            elif in_tuple[0] == 'equatorial':
                return coord_class(ra=val1, dec=val2, distance=val3, unit=unit_values), refsys
            elif in_tuple[0] == 'ecliptic':
                return coord_class(lon=val1, lat=val2, distance=val3, unit=unit_values), refsys
            elif in_tuple[0] == 'galactic':
                return coord_class(l=val1, b=val2, distance=val3, unit=unit_values), refsys

        else:
            # No need for coordinates
            return None
    else:
        if out_tuple[0] is not None:
            # Positions and proper motions
            refsys = in_tuple[0] if in_tuple[0] != 'cartesian' else cart_in_refsys

            indices = coord_info["in"][in_tuple[0]]
            indices_pm = coord_info["in"][in_tuple[1]]
            val1 = line[indices[0]]
            val2 = line[indices[1]]
            val3 = line[indices[2]]
            val1_pm = line[indices_pm[0]]
            val2_pm = line[indices_pm[1]]
            val3_pm = line[indices_pm[2]]

            unit_values = (fields0[indices[0]].unit, fields0[indices[1]].unit, fields0[indices[2]].unit)
            unit_values_pm = (fields0[indices_pm[0]].unit, fields0[indices_pm[1]].unit, fields0[indices_pm[2]].unit)

            unit_values = unit_values + unit_values_pm

            coord_class = getattr(coordinates, refsys.title())

            if in_tuple[0] == 'cartesian':
                return coord_class(x=val1, y=val2, z=val3, xdot=val1_pm, ydot=val2_pm, zdot=val3_pm, unit=unit_values), refsys
            elif in_tuple[0] == 'equatorial':
                return coord_class(ra=val1, dec=val2, distance=val3, mualpha=val1_pm, mudelta=val2_pm, radvel=val3_pm, unit=unit_values), refsys
            elif in_tuple[0] == 'ecliptic':
                return coord_class(lon=val1, lat=val2, distance=val3, mulambda=val1_pm, mubeta=val2_pm, radvel=val3_pm, unit=unit_values), refsys
            elif in_tuple[0] == 'galactic':
                return coord_class(l=val1, b=val2, distance=val3, mul=val1_pm, mub=val2_pm, radvel=val3_pm, unit=unit_values), refsys
        else:
            # ERROR!
            raise Exception("Can not convert velocities without positions!")


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
    cart_in_refsys = conf.get_input_reference_system()
    cart_out_refsys = conf.get_output_reference_system()

    # Create fields1
    for (idx, out_name) in enumerate(colnames1):
        columnfield = commonio.field(out_name, units1[idx])
        columnfield.refsys = get_refsys(out_name, conf.get_output_reference_system())

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
        for field0 in fields0:
            field0.refsys = get_refsys(field0.name, conf.get_input_reference_system())

        # Get information on in/out coordinates and velocities
        coord_info, in_tuple, out_tuple = get_io_coord_info(fields0, fields1)

        # Open and write header
        datawriter.open()
        datawriter.write_header(fields1)

        granularity = max(10, 10 ** (math.floor(math.log10(n)) - 3))

        # Loop over data, which should use generators in the case of files
        for idx, line in enumerate(dataloader.get_data()):
            if (idx + 1) % (n // granularity + 1) == 0 or idx == 0:
                logger.info("Processing line %i of %i - %3.1f%% (%.3f seconds)" % (idx + 1, n, ((idx + 1) * 100 / n), time.clock()))

            # Prepare a value container as a dictionary: {[col_name][out_unit]: value}
            # The out_unit in the key ensures that we can have the same field in different units in the output.
            # For example, we can have a column x in pc and another column x in km
            valuesbag = AutoVivification()
            for (idx, val) in enumerate(line):
                valuesbag[fields0[idx].name][fields0[idx].refsys][fields0[idx].unit] = val

            # Parse our coord_info structure and complete valuesbag with the necessary information
            coordinates_map = {}
            original_coordinates, refsys = get_coordinates(coord_info, in_tuple, out_tuple, fields0, line, cart_in_refsys)
            coordinates_map[refsys] = original_coordinates

            for coord_key in coord_info["out"].keys():
                if coord_key in names_spherical:
                    refsys = strip_suffix(coord_key, "_pm")
                elif coord_key in names_cartesian:
                    refsys = cart_out_refsys

                if(refsys in coordinates_map):
                    # Coordinate already exists!
                    coord_transf = coordinates_map[refsys]
                else:
                    # Calculate coordinates
                    coord_transf = original_coordinates.transform_to(refsys)
                    coordinates_map[refsys] = coord_transf

                names_bag = globals()[coord_key]
                for name in names_bag:
                    out_unit = fields1[fields1.index(name)].unit
                    valuesbag[name][refsys][out_unit] = coord_transf.get_basic(name, out_unit)

            # The output list of values for this row
            current_row = []

            # For each output column check if it is in the values bag
            for (idx, field1) in enumerate(fields1):
                out_name = field1.name
                out_refsys = field1.refsys
                out_unit = field1.unit

                if(out_name in valuesbag):
                    # We have the column:
                    #    Check reference systems and add refsys transformation if necessary
                    #    Else, if units different, add unit transformation
                    if(out_refsys in valuesbag[out_name]):
                        if(out_unit in valuesbag[out_name][out_refsys]):
                            # Have it all
                            current_row.append(valuesbag[out_name][out_refsys][out_unit])
                        else:
                            # Only unit conversion
                            current_unit = valuesbag[out_name][out_refsys].keys()[0]
                            value = units.convert(current_unit, out_unit, valuesbag[out_name][out_refsys][current_unit])
                            valuesbag[out_name][out_refsys][out_unit] = value
                            current_row.append(value)
                    else:
                        # We must be in cartesian -> cartesian, add refsys transform
                        if out_refsys in coordinates_map:
                            coord_transf = coordinates_map[out_refsys]
                        else:
                            coord_transf = original_coordinates.transform_to(out_refsys)
                            coordinates_map[out_refsys] = coord_transf

                        valuesbag['x'][out_refsys][out_unit] = coord_transf.get_basic('x', out_unit)
                        valuesbag['y'][out_refsys][out_unit] = coord_transf.get_basic('y', out_unit)
                        valuesbag['z'][out_refsys][out_unit] = coord_transf.get_basic('z', out_unit)

                        current_row.append(valuesbag[out_name][out_refsys][out_unit])

                else:
                    # No way to produce desired column
                    raise Exception("Can not produce the desired column '%s' with the input data" % out_name)

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
