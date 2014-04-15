"""
Reads the configuration file into a dictionary
:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1

"""
import ConfigParser

from visdatalog import logger


def get_config_section(section, conf_file):
    dict1 = {}
    Config = ConfigParser.ConfigParser()
    Config.read(conf_file)
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                logger.debug("skip: %s" % option)
        except:
            logger.error("exception on %s!" % option)
            dict1[option] = None
    return dict1


def init_config(conf_file):
    global inprops
    global outprops
    global confprops
    inprops = get_config_section("input", conf_file)
    outprops = get_config_section("output", conf_file)
    confprops = get_config_section("conf", conf_file)


def get_data_source():
    return inprops["data_source"]


def get_query():
    return inprops["query"]


def get_output_file_path():
    return outprops["out_filepath"]


def get_output_columns():
    """ Gets the output columns specified in the config file as an array of strings. """
    return outprops["out_columns"].strip().split()


def get_output_units():
    """ Gets the output units specified in the config file as an array of strings. """
    return outprops["out_units"].strip().split()


def get_input_reference_system():
    """ Gets the reference system of the input positions """
    return inprops["in_refsys"].strip()


def get_output_reference_system():
    """ Gets the reference system of the output cartesian positions """
    return outprops["out_refsys"].strip()


def get_log_path():
    """ Gets the log file """
    return confprops["logpath"].strip()
