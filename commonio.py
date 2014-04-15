"""
Data io module. Contains classes to load data from different sources (VO using TAP, ASCII files, HDF5 files...).
Also contains a class to write ASCII files, which is the only supported format right now.
:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1

Requirements
------------
* `VOTable <http://vo.ari.uni-heidelberg.de/soft/subpkgs>`_
"""
from astropy.table import Table
from gavo import votable
import os

from visdatalog import logger


class field():

    def __init__(self, name, unit=None, refsys="na"):
        self.name = name
        self.unit = unit
        self.refsys = refsys

    # toString
    def __str__(self):
        return self.name + "[" + self.unit + "]"

    # Equals: The dictionaries are the same, all fields match.
    # When comparing to a string, we just compare self.name
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    # notEquals
    def __ne__(self, other):
        return not self.__eq__(other)


class IOData():
    """ Main abstract class """
    def __init__(self, data_source):
        self.data_source = data_source

    def get_record_count(self):
        raise NotImplementedError("Should have implemented this")

    def length(self):
        raise NotImplementedError("Should have implemented this")


class TAPLoader(IOData):

    def __init__(self, data_source, query):
        IOData.__init__(self, data_source)
        self.query = query
        # Load data
        logger.info("Connecting to TAP service:\n %s \n with query:\n %s", self.data_source, self.query)
        job = votable.ADQLTAPJob(self.data_source, self.query)
        job.run()
        dataIterator, metadata = votable.load(job.openResult())
        job.delete()

        self.datalist = list(dataIterator)

        # Prepare column names
        self.fields = []
        for fld in metadata.fields:
            self.fields.append(field(fld.name, fld.unit))

    def get_data(self):
        """ Loads the data using the TAP service specified in the configuration.

            Executes the query in the configuration file to the given TAP service
            and returns the result. It also returns the column names.
        """

        return self.datalist

    def __len__(self):
        return len(self.datalist)


class AsciiLoader(IOData):

    def __init__(self, data_source):
        IOData.__init__(self, data_source)
        # Initialize fields by reading the first commented lines of the file
        # Holds what metadata have been loaded. [0] = names, [1] = units
        descriptors = [False, False]
        self.fields = []
        f = open(self.data_source, 'r')
        for line in f:
            # Skip comments, first three comments are candidates to be names and units
            if(line.startswith("#") and not descriptors[0]):
                # Column names
                colnames = line.strip('#').split()
                for colname in colnames:
                    self.fields.append(field(colname))
                descriptors[0] = True
            elif(line.startswith("#") and descriptors[0] and not descriptors[1]):
                # Units
                units = line.strip('#').split()
                for (index, unit) in enumerate(units):
                    self.fields[index].unit = unit
                descriptors[1] = True
            elif(not line.startswith("#")):
                break
        f.close()

    def get_data(self):
        """ Returns a list of processed lines as floats
        """
        f = open(self.data_source, 'r')
        while True:
            line = f.readline()
            if not line:
                f.close()
                break
            elif not line.startswith('#'):
                yield map(float, line.split())

    def __len__(self):
        x = getattr(self, 'record_count', None)
        if x is None:
            # Calculate length
            f = open(self.data_source, 'r')
            self.record_count = sum(bl.count("\n") for bl in self._blocks(f))
            f.close()

        return self.record_count

    def _blocks(self, thefile, size=65536):
        while True:
            b = thefile.read(size)
            if not b:
                break
            yield b


class FITSLoader(IOData):
    def __init__(self, data_source):
        IOData.__init__(self, data_source)
        """ Loads the FITS file using astropy """
        self.datalist = Table.read(data_source, format='fits')

    def get_data(self):
        """ Returns the data """
        return self.datalist

    def __len__(self):
        return len(self.datalist)


class HDF5Loader(IOData):
    def __init__(self, data_source):
        """ Loads the HDF5 file using astropy """
        self.datalist = Table.read(data_source, format='hdf5')

    def get_data(self):
        """ Returns the data """
        return self.datalist

    def __len__(self):
        return len(self.datalist)


class AsciiWriter(IOData):

    def __init__(self, data_source):
        IOData.__init__(self, data_source)
        self.separator = ' '
        self.records = 0

    def open(self):
        try:
            self.f = open(self.data_source, 'w')
        except:
            raise Exception("Error opening file: %s" % self.data_source)

    def close(self):
        x = getattr(self, 'f', None)
        if not x is None and not self.f.closed:
            self.f.close()

    def delete_file(self):
        try:
            os.remove(self.data_source)
        except:
            raise Exception("Could not remove data source: %s" % self.data_source)

    def write_header(self, fields):
        """ Writes the header of the file.

            Writes the header of the output file using the currently opened file and
            the given fields. It writes three comment lines:
            *line 1: The names of the columns.
            *line 2: The units of the columns.
         """
        x = getattr(self, 'f', None)
        if x is None or self.f.closed:
            raise IOError("Attempting to write on a file that has not been opened: %s" % self.data_source)

        sep = self.separator
        try:
            # Write the column names
            self.f.write('#' + sep.join(colname.name for colname in fields) + '\n')

            # Write the units
            self.f.write('#' + sep.join((colname.unit or "NA") for colname in fields) + '\n')
        except:
            logger.error("Error writing file: %s" % self.data_source)
            raise Exception("Error writing file: %s" % self.data_source)

    def write_line(self, dataline):
        """ Saves the given data to the given output file.

            Saves the given datalist to the file adding a commented row at the top
            with the given column names.
        """
        x = getattr(self, 'f', None)
        if x is None or self.f.closed:
            raise IOError("Attempting to write on a file that has not been opened: %s" % self.data_source)

        # Write the actual data
        self.f.write(self.separator.join(map(str, dataline)) + '\n')
        self.records += 1
