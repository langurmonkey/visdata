"""
Classes representing coordinates and proper motions/velocities in different reference systems.
Conversions between the systems and between different units are included.
Distances are stored in pc and velocities are stored in km/s. All angles are in radians, and
all angular velocities are in rad/s.
:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1

Requirements
------------
* `Numpy 1.7 <http://www.numpy.org>`_
* `transformations by Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`_
* `Astropy 0.3+ <http://www.astropy.org>`_
"""
import coordinatescore
import units


def get_coordinates(refsys_name, **kwargs):
    coord_class = globals()[refsys_name.title()]
    return coord_class(**kwargs)


_refsys_types = ['equatorial', 'galactic', 'ecliptic']
_cartesian = 'cart'
_spherical = 'sph'


class CoordinatesBase():
    """ Represents base coordinates and proper motions.
    """
    _dist_names = ['distance', 'dist', 'r']
    _cart_names = ['x', 'y', 'z']
    _dist_names_pm = ['radvel', 'radialvelocity', 'rdot']
    _cart_names_pm = ['xdot', 'ydot', 'zdot']
    _pm = False

    def __init__(self, precalculate=True, **kwargs):
        # Input units
        units_all = kwargs.pop('unit')

        # Set internal units. We store angles in rad, distances in pc, velocities in km/s and angular
        # velocities in rad/s
        self._sph_unit = ['rad', 'rad', 'pc']
        self._cart_unit = ['pc', 'pc', 'pc']
        self._sph_unit_pm = ['rad/s', 'rad/s', 'km/s']
        self._cart_unit_pm = ['km/s', 'km/s', 'km/s']

        if all(cart_name in kwargs for cart_name in self._cart_names):
            # Cartesian coordinates
            self._cart = [kwargs.pop('x'), kwargs.pop('y'), kwargs.pop('z')]

            # Transform cartesian coordinates to km
            self._cart[0] = units.convert(units_all[:3][0], self._cart_unit[0], self._cart[0])
            self._cart[1] = units.convert(units_all[:3][1], self._cart_unit[1], self._cart[1])
            self._cart[2] = units.convert(units_all[:3][2], self._cart_unit[2], self._cart[2])

            if all(cart_name_pm in kwargs for cart_name_pm in self._cart_names_pm):
                # Cartesian velocities
                self._cart_pm = [kwargs.pop('xdot'), kwargs.pop('ydot'), kwargs.pop('zdot')]

                # Transform cartesian velocities to km/s
                self._cart_pm[0] = units.convert(units_all[3:][0], self._cart_unit_pm[0], self._cart_pm[0])
                self._cart_pm[1] = units.convert(units_all[3:][1], self._cart_unit_pm[1], self._cart_pm[1])
                self._cart_pm[2] = units.convert(units_all[3:][2], self._cart_unit_pm[2], self._cart_pm[2])

                self._pm = True

        elif any(latitude in kwargs for latitude in self._latitude_names) and any(longitude in kwargs for longitude in self._longitude_names) and any(dist in kwargs for dist in self._dist_names):
            # Spherical coordinates names
            long_name = next(longitude for longitude in self._longitude_names if longitude in kwargs)
            lat_name = next(latitude for latitude in self._latitude_names if latitude in kwargs)
            dist_name = next(dist for dist in self._dist_names if dist in kwargs)

            # Spherical coordinates
            self._sph = [kwargs.pop(long_name), kwargs.pop(lat_name), kwargs.pop(dist_name)]
            aux_sph_unit = units_all[:3]

            # Convert angles to rad and distance to km
            self._sph[0] = units.convert(aux_sph_unit[0], self._sph_unit[0], self._sph[0])
            self._sph[1] = units.convert(aux_sph_unit[1], self._sph_unit[1], self._sph[1])
            self._sph[2] = units.convert(aux_sph_unit[2], self._sph_unit[2], self._sph[2])

            if any(latitude_pm in kwargs for latitude_pm in self._latitude_names_pm) and any(longitude_pm in kwargs for longitude_pm in self._longitude_names_pm) and any(radvel in kwargs for radvel in self._dist_names_pm):
                # Spherical proper motions names
                long_name_pm = next(longitude_pm for longitude_pm in self._longitude_names_pm if longitude_pm in kwargs)
                lat_name_pm = next(latitude_pm for latitude_pm in self._latitude_names_pm if latitude_pm in kwargs)
                dist_name_pm = next(radvel for radvel in self._dist_names_pm if radvel in kwargs)

                # Spherical proper motions
                self._sph_pm = [kwargs.pop(long_name_pm), kwargs.pop(lat_name_pm), kwargs.pop(dist_name_pm)]
                aux_sph_unit_pm = units_all[3:]

                # Convert proper motions to rad/s and radial velocity to km/s
                self._sph_pm[0] = units.convert(aux_sph_unit_pm[0], self._sph_unit_pm[0], self._sph_pm[0])
                self._sph_pm[1] = units.convert(aux_sph_unit_pm[1], self._sph_unit_pm[1], self._sph_pm[1])
                self._sph_pm[2] = units.convert(aux_sph_unit_pm[2], self._sph_unit_pm[2], self._sph_pm[2])

                self._pm = True

        else:
            raise Exception("%s initialized with insufficient or inconsistent data: %s" % (self._name.title(), kwargs))

        if precalculate:
            # Work out cartesian or spherical right now
            self._calculate_other_coordinates()

    def transform_to(self, other_type):
        """ Transforms the current coordinates into the given other type.

            Parameters
            ----------
            other_type -- string representing the reference system to convert to: 'equatorial', 'galactic' or 'ecliptic'
        """
        if other_type not in _refsys_types:
            raise Exception("Type '%s' is not one of our coordinate system definitions: %s" % (other_type, _refsys_types))

        if other_type == self._name:
            # No conversion needed
            return self
        else:
            # Conversion needed
            converter_name = self._name + '_to_' + other_type + '_v'
            converter_method = getattr(coordinatescore, converter_name)
            out = converter_method(self._cart)
            other_class = globals()[other_type.title()]
            if self._pm:
                # Proper motions
                out_pm = converter_method(self._cart_pm)
                return other_class(x=out[0], y=out[1], z=out[2], xdot=out_pm[0], ydot=out_pm[1], zdot=out_pm[2], unit=tuple(self._cart_unit + self._cart_unit_pm))
            else:
                # No proper motions
                return other_class(x=out[0], y=out[1], z=out[2], unit=tuple(self._cart_unit))

    def __str__(self, coord_type=_spherical):
        coords = getattr(self, '_' + coord_type)
        units = getattr(self, '_' + coord_type + '_unit')
        names = self._cart_names
        if coord_type == _spherical:
            names = ['longitude', 'latitude', 'dist']
        elif coord_type == _cartesian:
            names = ['x', 'y', 'z']

        out_str = "%s (%s=%.3f %s, %s=%.3f %s, %s=%.3f %s)" % (self._name, names[0], coords[0], units[0], names[1], coords[1], units[1], names[2], coords[2], units[2])

        if self._pm:
            coords_pm = getattr(self, '_' + coord_type + '_pm')
            units_pm = getattr(self, '_' + coord_type + '_unit_pm')
            names_pm = self._cart_names_pm

            out_str += '\n'
            out_str += "  vel (%s=%.8f %s, %s=%.8f %s, %s=%.8f %s)" % (names_pm[0], coords_pm[0], units_pm[0], names_pm[1], coords_pm[1], units_pm[1], names_pm[2], coords_pm[2], units_pm[2])

        return out_str

    def to_string(self, coord_type='sph'):
        return self.__str__(coord_type)

    def get(self, **kwargs):
        """ Gets the given coordinates in the desired units

            Parameters
            ----------
            **kwargs -- Dictionary of name='unit' pairs
        """
        result = {}
        for key in kwargs.keys():
            result[key] = self.get_basic(key, kwargs[key])
        return result

    def get_basic(self, name=None, unit=None):
        """ Gets the value under the given name converting it into the given unit, if possible.

            For example:
            -- self.get_basic('x', 'km') will return the value of x in kilometers.
            -- self.get_basic('longitude', 'deg') will return the longitude angle in degrees.
        """
        if name in self._cart_names:
            # Cartesian coordinate
            i = self._cart_names.index(name)
            return units.convert(self._cart_unit[i], unit, self._cart[i])
        elif name in self._cart_names_pm:
            # Cartesian velocity
            i = self._cart_names_pm.index(name)
            return units.convert(self._cart_unit_pm[i], unit, self._cart_pm[i])

        elif name in self._longitude_names:
            return units.convert(self._sph_unit[0], unit, self._sph[0])
        elif name in self._longitude_names_pm:
            return units.convert(self._sph_unit_pm[0], unit, self._sph_pm[0])

        elif name in self._latitude_names:
            return units.convert(self._sph_unit[1], unit, self._sph[1])
        elif name in self._latitude_names_pm:
            return units.convert(self._sph_unit_pm[1], unit, self._sph_pm[1])

        elif name in self._dist_names:
            return units.convert(self._sph_unit[2], unit, self._sph[2])
        elif name in self._dist_names_pm:
            return units.convert(self._sph_unit_pm[2], unit, self._sph_pm[2])

        else:
            raise Exception("Name not recognized: %s" % name)

    def _calculate_other_coordinates(self):
        # Work out cartesian or spherical
        if hasattr(self, '_cart'):
            if self._pm:
                # With proper motions
                aux = coordinatescore.cartesian_to_spherical_pm(self._cart + self._cart_pm)
                self._sph = aux[:3]
                self._sph_pm = aux[3:]
            else:
                # Without proper motions
                self._sph = coordinatescore.cartesian_to_spherical(self._cart)

        elif hasattr(self, '_sph'):
            if self._pm:
                # With velocities
                aux = coordinatescore.spherical_to_cartesian_pm(self._sph + self._sph_pm)
                self._cart = aux[:3]
                self._cart_pm = aux[3:]
            else:
                # Without velocities
                self._cart = coordinatescore.spherical_to_cartesian(self._sph)


class Equatorial(CoordinatesBase):
    """ Represents equatorial coordinates """

    _name = 'equatorial'
    _longitude_names = ['lon', 'ra', 'alpha', 'longitude']
    _latitude_names = ['lat', 'dec', 'delta', 'latitude']
    _longitude_names_pm = ['mualpha', 'mura']
    _latitude_names_pm = ['mudelta', 'mudec']

    def __init__(self, **kwargs):
        """ Create a new Equatorial coordinates object.

            Examples
            --------
            eq_coord_1 = coordinates.Equatorial(ra=120, dec=-90, distance=10, unit=('deg', 'deg', 'pc'))
            eq_coord_2 = coordinates.Equatorial(x=10, y=10, z=10, unit=('km'))
            eq_coord_1 = coordinates.Equatorial(ra=120, dec=-90, distance=10, mualpha=2.3, mudelta=-1.43, radialvelocity=5, unit=('deg', 'deg', 'pc', 'mas/a', 'rad/s', 'km/s'))
        """
        CoordinatesBase.__init__(self, **kwargs)


class Ecliptic(CoordinatesBase):
    """ Represents ecliptic coordinates """

    _name = 'ecliptic'
    _longitude_names = ['lon', 'lambda', 'longitude']
    _latitude_names = ['lat', 'beta', 'latitude']
    _longitude_names_pm = ['mulambda']
    _latitude_names_pm = ['mubeta']

    def __init__(self, **kwargs):
        CoordinatesBase.__init__(self, **kwargs)


class Galactic(CoordinatesBase):
    """ Represents galactic coordinates """

    _name = 'galactic'
    _longitude_names = ['lon', 'l', 'longitude']
    _latitude_names = ['lat', 'b', 'latitude']
    _longitude_names_pm = ['mul', 'mulongitude']
    _latitude_names_pm = ['mub', 'mulatitude']

    def __init__(self, **kwargs):
        CoordinatesBase.__init__(self, **kwargs)
