"""
Various utility methods for converting coordinates and proper motions/velocities.
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
"""
import math
import numpy

import transformations as trans


# SOME CONSTANTS
# Obliquity of ecliptic in radians J2000 with T=0
obliquity = math.radians(23.4392916667)
# Node line galactic longitude
nodes_l = math.radians(-33)
# Inclination of galactic equator to celestial equator
incl_gal_cel = math.radians(62.9)
# Node line equatorial right ascension
nodes_ra = math.radians(282.25)


def cartesian_to_spherical(values):
    """ Transforms cartesian coordinates into spherical coordinates.

        Converts from cartesian to spherical coordinates given the x, y, and z. The
        resulting angles are in radians and the distance is in the units of x, y, z.

        Parameters
        ----------
        values -- A list with the x, y and z components in this order. They all must be in the same units.

        Returns
        -------
        A list with the alpha, delta (in rad) and distance (in the units of x, y and z).
    """
    x, y, z = values[:]
    xsq = x ** 2
    ysq = y ** 2
    zsq = z ** 2
    distance = math.sqrt(xsq + ysq + zsq)

    alpha = math.atan2(y, x)
    # Correct the value of alpha depending upon the quadrant.
    if alpha < 0:
        alpha += 2 * math.pi

    if (xsq + ysq) == 0:
        # In the case of the poles, delta is -90 or +90
        delta = math.copysign(math.pi / 2, z)
    else:
        delta = math.atan(z / math.sqrt(xsq + ysq))

    return [alpha, delta, distance]


def spherical_to_cartesian(values):
    """ Transforms spherical coordinates into cartesian coordinates.

        Converts from spehrical to cartesian coordinates given the longitude , latitude and the radius.
        Both angles must be in radians.

        Parameters
        ----------
        values -- A list with the alpha, delta, in radians, and distance.

        Returns
        -------
        A list with the x, y and z in the units of distance.
    """
    longitude, latitude, radius = values[:]
    x = radius * math.cos(latitude) * math.cos(longitude)
    y = radius * math.cos(latitude) * math.sin(longitude)
    z = radius * math.sin(latitude)

    return [x, y, z]


def spherical_to_cartesian_pm(values):
    """ Transforms spherical proper motions into cartesian velocities.

        Converts the given positions and proper motions to cartesian positions and velocities.

        Parameters
        ----------
        values -- A list containing the following elements in the given order. The distance units can be anything, here they are denoted by 'du'.
            alpha -- The right ascension in radians.
            delta -- The declination in radians.
            distance -- The distance in du.
            mualpha -- The proper motion in alpha in rad/s
            mudelta -- The proper motion in delta in rad/s
            radvel -- The radial velocity in du/s

        Returns
        -------
        A vector with x, y, z in du and xdot, ydot and zdot in du/s
    """
    alpha, delta, distance, mualpha, mudelta, radvel = values[:]
    x, y, z = spherical_to_cartesian([alpha, delta, distance])

    # In alpha, we need to multiply by cos(delta) to account for the radius from axis being smaller at higher
    # declinations
    # Transversal velocities in km/s
    vta = mualpha * math.cos(delta) * distance
    # vta = mualpha * distance
    vtd = mudelta * distance

    """ Equations of proper motion to cartesian coordinates
    vx = (vR cos(delta) cos(alpha)) - (vTA sin(alpha)) - (vTD sin(delta) cos(alpha))
    vy = (vR cos(delta) sin(alpha)) + (vTA cos(alpha)) - (vTD sin(delta) sin(alpha))
    vz = vR sin(delta) + vTD cos(delta)
    """
    xdot = (radvel * math.cos(delta) * math.cos(alpha)) - (vta * math.sin(alpha)) - (vtd * math.sin(delta) * math.cos(alpha))
    ydot = (radvel * math.cos(delta) * math.sin(alpha)) + (vta * math.cos(alpha)) - (vtd * math.sin(delta) * math.sin(alpha))
    zdot = radvel * math.sin(delta) + vtd * math.cos(delta)

    """ vx, vy, vz """
    return [x, y, z, xdot, ydot, zdot]


def cartesian_to_spherical_pm(values):
    """ Transforms cartesian velocities into spherical proper motions.

        Converts the given spherical positions, proper motions, distance and radial velocity into
        cartesian coordinates and cartesian velocities.

        Parameters
        ----------
        values -- A list containing the following elements in the given order. The distance units can be anything, here they are denoted by 'du'.
            x -- The x cartesian coordinate in du.
            y -- The y cartesian coordinate in du.
            z -- The z cartesian coordinate in du.
            xdot -- The velocity in x in du/s.
            ydot -- The velocity in y in du/s.
            zdot -- The velocity in z in du/s.

        Returns
        -------
        A vector with alpha, delta in rad, distance in du, mualpha, mudelta in rad/s and radvel in du/s
    """
    alpha, delta, distance = cartesian_to_spherical(values[:3])
    xdot, ydot, zdot = values[3:]

    """
        As yielded by Maxima:

        eqsys:[a=(z*cos(i)*cos(j))-(x*cos(i)*k*sin(j))-(y*k*sin(i)*cos(j)),b=(z*cos(i)*sin(j))+(x*cos(i)*k*cos(j))-(y*k*sin(i)*sin(j)),c=z*sin(i)+y*k*cos(i)];
        solve(eqsys,[x,y,z]);

        x=-(a*sin(j)-b*cos(j))/((cos(i)*sin(j)^2+cos(i)*cos(j)^2)*k)
        y=(c*cos(i)*sin(j)^2-b*sin(i)*sin(j)+c*cos(i)*cos(j)^2-a*sin(i)*cos(j))/(((sin(i)^2+cos(i)^2)*sin(j)^2+(sin(i)^2+cos(i)^2)*cos(j)^2)*k)
        z=(c*sin(i)*sin(j)^2+b*cos(i)*sin(j)+c*sin(i)*cos(j)^2+a*cos(i)*cos(j))/((sin(i)^2+cos(i)^2)*sin(j)^2+(sin(i)^2+cos(i)^2)*cos(j)^2)

        where:
        xdot = a
        ydot = b
        zdot = c
        alpha = j
        delta = i
        dist = k
        mualpha = x
        mudelta = y
        radvel = z
    """
    mualpha = -(xdot * math.sin(alpha) - ydot * math.cos(alpha)) / ((math.cos(delta) * math.sin(alpha) ** 2 + math.cos(delta) * math.cos(alpha) ** 2) * distance)
    mudelta = (xdot * math.cos(delta) * math.sin(alpha) ** 2 - ydot * math.sin(delta) * math.sin(alpha) + zdot * math.cos(delta) * math.cos(alpha) ** 2 - xdot * math.sin(delta) * math.cos(alpha)) / (((math.sin(delta) ** 2 + math.cos(delta) ** 2) * math.sin(alpha) ** 2 + (math.sin(delta) ** 2 + math.cos(delta) ** 2) * math.cos(alpha) ** 2) * distance)
    radvel = (zdot * math.sin(delta) * math.sin(alpha) ** 2 + ydot * math.cos(delta) * math.sin(alpha) + zdot * math.sin(delta) * math.cos(alpha) ** 2 + xdot * math.cos(delta) * math.cos(alpha)) / ((math.sin(delta) ** 2 + math.cos(delta) ** 2) * math.sin(alpha) ** 2 + (math.sin(delta) ** 2 + math.cos(delta) ** 2) * math.cos(alpha) ** 2)

    """ mualpha, mudelta [rad/s], radvel [km/s] """
    return [alpha, delta, distance, mualpha, mudelta, radvel]


def galactic_to_equatorial():
    """ Returns the transformation matrix to convert from galactic to equatorial cartesian coordinates.

        Returns the matrix to transform from galactic to equatorial cartesian coordinates.
        The inclination of the galactic equator to the celestial equator is 62.9deg.
        The intersection, or node line, of the two equators is at RA=282.25deg DEC=0deg and
        l=33deg b=0deg. So we have the Euler angles alpha=-33, beta=62.9, gamma=282.25.

        In order to use this matrix to transform the coordinates, do:
        > eq_vector = numpy.dot(gal_vector, galactic_to_equatorial().T)
    """
    return rotation_matrix_angles(nodes_l, incl_gal_cel, nodes_ra)


def galactic_to_equatorial_v(vector):
    """ Transforms the given vector from galactic to equatorial cartesian coordinates.

        Parameters
        ----------
        vector -- A list with the galactic cartesian coordinates [x, y, z].
    """
    vector4 = numpy.append(numpy.array(vector), numpy.array([1]))
    return numpy.dot(vector4, galactic_to_equatorial().T)[:3].tolist()


def galactic_to_equatorial_sph(vector):
    """ Transforms the given galactic coordinates l, b into alpha, delta.

        Parameters
        ----------
        vector -- A list with the galactic coordinates [l, b].
    """
    vector.append(1)
    gal_cart = spherical_to_cartesian(vector)
    sph_eq = cartesian_to_spherical(galactic_to_equatorial_v(gal_cart))

    return sph_eq[:2]


def equatorial_to_galactic():
    """ Returns the transformation matrix to convert from equatorial to galactic cartesian coordinates.

        In order to use this matrix to transform the coordinates, do:
        > gal_vector = numpy.dot(eq_vector, equatorial_to_galactic().T)
    """
    return rotation_matrix_angles(-nodes_ra, -incl_gal_cel, -nodes_l)


def equatorial_to_galactic_v(vector):
    """ Transforms the given vector from equatorial to galactic cartesian coordinates.

        Parameters
        ----------
        vector -- A list with the equatorial cartesian coordinates [x, y, z].
    """
    vector4 = numpy.append(numpy.array(vector), numpy.array([1]))
    return numpy.dot(vector4, equatorial_to_galactic().T)[:3].tolist()


def equatorial_to_galactic_sph(vector):
    """ Transforms the given equatorial coordinates ra, dec into l, b.

        Parameters
        ----------
        vector -- A list with the equatorial coordinates [ra, dec].
    """
    vector.append(1)
    eq_cart = spherical_to_cartesian(vector)
    sph_gal = cartesian_to_spherical(equatorial_to_galactic_v(eq_cart))

    return sph_gal[:2]


def equatorial_to_ecliptic():
    """ Returns the transformation matrix to convert from equatorial to ecliptic coordinates.

        Gets the rotation matrix to transform from equatorial to ecliptic coordinates. Since the zero point
        in both systems is the same (the vernal equinox, gamma, defined as the intersection between the equator and the ecliptic),
        alpha and gamma are zero. Beta, the angle between the up directions of both systems, is precisely the obliquity of the
        ecliptic, epsilon=23.439281 deg. So we have the Euler angles alpha=0 deg, beta=epsilon deg, gamma=0 deg.

        In order to use this matrix to transform the coordinates, do:
        > ecl_vector = numpy.dot(eq_vector, equatorial_to_ecliptic().T)
    """
    return rotation_matrix_angles(0, -obliquity, 0)


def equatorial_to_ecliptic_v(vector):
    """ Transfroms the given vector from equatorial to ecliptic coordinates.

        Parameters
        ----------
        vector -- A list with the equatorial cartesian coordinates [x, y, z].
    """
    vector4 = numpy.append(numpy.array(vector), numpy.array([1]))
    return numpy.dot(vector4, equatorial_to_ecliptic().T)[:3].tolist()


def equatorial_to_ecliptic_sph(vector):
    """ Transforms the given ra, dec (in radians) in equatorial coordinates to lambda, beta in ecliptic coordinates.

        Parameters
        ----------
        vector -- A list with the equatorial coordinates [ra, dec].
    """
    alpha = vector[0]
    delta = vector[1]

    l = math.atan2((math.sin(alpha) * math.cos(obliquity) + math.tan(delta) * math.sin(obliquity)), math.cos(alpha))
    if l < 0:
        l += math.pi * 2

    b = math.asin(math.sin(delta) * math.cos(obliquity) - math.cos(delta) * math.sin(obliquity) * math.sin(alpha))

    return {l, b}


def ecliptic_to_equatorial():
    """ Returns the transformation matrix to convert from ecliptic to equatorial coordinates

        In order to use this matrix to transform the coordinates, do:
        > eq_vector = numpy.dot(ecl_vector, ecliptic_to_equatorial().T)
    """
    return rotation_matrix_angles(0, obliquity, 0)


def ecliptic_to_equatorial_v(vector):
    """ Transfroms the given vector from ecliptic to equatorial coordinates.

        Parameters
        ----------
        vector -- A list with the ecliptic cartesian coordinates [x, y, z].
    """
    vector4 = numpy.append(numpy.array(vector), numpy.array([1]))
    return numpy.dot(vector4, ecliptic_to_equatorial().T)[0:3].tolist()


def ecliptic_to_equatorial_sph(vector):
    """ Converts lambda, beta in radians to alpha, delta.

        Parameters
        ----------
        vector -- A list with the ecliptic coordinates [lambda, beta].
    """
    l = vector[0]
    b = vector[1]

    alpha = math.atan2((math.sin(l) * math.cos(obliquity) - math.tan(b) * math.sin(obliquity)), math.cos(l))
    if alpha < 0:
        alpha += math.pi * 2

    delta = math.asin(math.sin(b) * math.cos(obliquity) + math.cos(b) * math.sin(obliquity) * math.sin(l))

    return {alpha, delta}


def ecliptic_to_galactic():
    """ Returns the transformation matrix to convert from eclpitic to galactic cartesian coordinates.

        In order to use this matrix to transform the coordinates, do:
        > gal_vector = numpy.dot(ecl_vector, equatorial_to_ecliptic().T)
    """
    # Transofrmation compositon: ECL -> GAL = ECL -> EQ + EQ -> GAL
    return numpy.dot(equatorial_to_galactic(), ecliptic_to_equatorial())


def ecliptic_to_galactic_v(vector):
    """ Transforms the given vector from ecliptic to galactic cartesian coordinates.

        Parameters
        ----------
        vector -- A list with the ecliptic cartesian coordinates [x, y, z].
    """
    vector4 = numpy.append(numpy.array(vector), numpy.array([1]))
    return numpy.dot(vector4, ecliptic_to_galactic().T)[0:3].tolist()


def ecliptic_to_galactic_sph(vector):
    """ Converts lambda, beta in radians to l, b in radians.

        Parameters
        ----------
        vector -- A list with the ecliptic coordinates [lambda, beta].
    """
    vector.append(1)
    ecl_cart = spherical_to_cartesian(vector)
    sph_gal = cartesian_to_spherical(ecliptic_to_galactic_v(ecl_cart))

    return sph_gal[:2]


def galactic_to_ecliptic():
    """ Returns the transformation matrix to convert from galactic to ecliptic cartesian coordinates.
    """
    # Transofrmation compositon: GAL -> ECL = GAL -> EQ + EQ -> ECL
    return numpy.dot(equatorial_to_ecliptic(), galactic_to_equatorial())


def galactic_to_ecliptic_v(vector):
    """ Transforms the given vector from galactic to ecliptic carteisan coordinates.

        Parameters
        ----------
        vector -- A list with the galactic cartesian coordinates [x, y, z].
    """
    vector4 = numpy.append(numpy.array(vector), numpy.array([1]))
    return numpy.dot(vector4, galactic_to_ecliptic().T)[0:3].tolist()


def galactic_to_ecliptic_sph(vector):
    """ Transforms galactic coordinates l, b in radians to ecliptic ones lambda, beta.

        Parameters
        ----------
        vector -- A list with the galactic coordinates [l, b].
    """
    vector.append(1)
    gal_cart = spherical_to_cartesian(vector)
    sph_ecl = cartesian_to_spherical(galactic_to_ecliptic_v(gal_cart))

    return sph_ecl[:2]


z = numpy.array([0, 0, 1])
x = numpy.array([1, 0, 0])


def rotation_matrix_angles(alpha, beta, gamma):
    """ Returns the rotation matrix that applies a rotation by the given Euler angles.

        It applies Rz(gamma)*Rx(beta)*Rz(alpha).
    """
    R1 = trans.rotation_matrix(gamma, z)
    R2 = trans.rotation_matrix(beta, x)
    R3 = trans.rotation_matrix(alpha, z)

    M = trans.concatenate_matrices(R1, R2, R3)

    return M
