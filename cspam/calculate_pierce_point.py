import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

def calculate_pierce_point(RA, DEC, lon, lat, height, time, screenheight = 400000):
    """
    This function etc etc etc.

    Parameters
    ----------
    RA : str
        '12h22m54.899s'
    DEC : str
        '+15d49m20.57s'
    lon : double
        31.956389 (degrees)
    lat : double
        -111.598333 (degrees)
    time : str
        '2010-01-01T12:00:00'
    screenheight : double
        Optional screen height in meters. Default at 400 km.

    Returns
    -------
    pierce_point_xyz : Astropy EarthLocation class at the moment
    """

    # Use astropy to go from equatorial system to horizontal system
    # i.e. from RA, DEC and antenna position to Az, El
    # See http://www.astropy.org/coordinates-benchmark/summary.html
    object = SkyCoord(RA, DEC)
    obstime = Time(time, scale='utc')
    antposition = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    obj_altaz = object.transform_to(AltAz(obstime=obstime, location=antposition))
    el = obj_altaz.alt.rad
    az = obj_altaz.az.rad

    # Find ECEF (Earth-Centered, Earth-Fixed) antenna positions assuming the
    # Earth is adequatly described by the WGS84 oblate spheroid.
    # The conversion is handled by astropy and relies on ERFA (Essential
    # Routines for Fundamental Astronomy) which is based on the SOFA library
    # by the International Astronomical Union.
    # See http://en.wikipedia.org/wiki/Geographic_coordinate_conversion
    antx = antposition.x.value # in meters
    anty = antposition.y.value # in meters
    antz = antposition.z.value # in meters

    # Find the unit normal vector to the ellipsoid at the antenna position.
    # Note that the WGS84 oblate spheroid is defined as:
    #
    # x^2 + y^2   z^2
    # --------- + --- = 1, with:
    #    a^2      b^2
    #
    # a = 6 378 137.0 m
    # b = 6 356 752.314 245 m
    # See http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
    a = 6378137.0
    b = 6356752.314245

    normalx = 2*antx/a**2.0
    normaly = 2*anty/a**2.0
    normalz = 2*antz/b**2.0

    unitnormalx = normalx/np.sqrt(normalx**2.0+normaly**2.0+normalz**2.0)
    unitnormaly = normaly/np.sqrt(normalx**2.0+normaly**2.0+normalz**2.0)
    unitnormalz = normalz/np.sqrt(normalx**2.0+normaly**2.0+normalz**2.0)

    # Parameterize Az, El line-of-sight towards the object in the sky.
    # Note that instead of Elevation, the inclination (polar angle) is used.
    # Also, note the minus sign because of the azimuth being definined with
    # the x-direction as north and y-direction as west.
    inclination = np.pi/2 - el

    param_x = np.sin(inclination)*np.cos(az)
    param_y = -np.sin(inclination)*np.sin(az)
    param_z = np.cos(inclination)

    # Now, rotate and translate the horizontal frame in such a way that it
    # is tangent to the Earth at the antenna position and that the
    # x-axis points towards the north in the ECEF frame.
    # This comes down to rotating and translating the parameterization of
    # the line-of-sight using theta and phi from the normal vector.
    # The angle between the north in horizontal frame and the unit normal
    # vector is theta. Phi is the angle between the x-direction in the
    # horizontal frame and the projection of the normal vector on the
    # xy-plane.
    theta = np.arccos(unitnormalz)
    phi = np.arctan(unitnormaly/unitnormalx)

    rotation_matrix_z = np.matrix( [[np.cos(phi), -np.sin(phi), 0],
                    [np.sin(phi), np.cos(phi), 0],
                    [0, 0, 1]] )
    rotation_matrix_y = np.matrix( [[np.cos(theta), 0, -np.sin(theta)],
                    [0, 1, 0],
                    [np.sin(theta), 0, np.cos(theta)]] )

    param_vector = np.array([param_x, param_y, param_z])[:,None]
    intermediate_step = np.dot(rotation_matrix_y, param_vector)
    rotated_param = np.dot(rotation_matrix_z, intermediate_step)
    rotated_param = np.squeeze(np.asarray(rotated_param))

    rotated_param_x = rotated_param[0]
    rotated_param_y = rotated_param[1]
    rotated_param_z = rotated_param[2]

    # This rotated parameterization can now be translated using the
    # antenna position. This gives a parameterization of the line-
    # of-sight in the ECEF system. The intersection with a larger
    # ellipsoid (given by the screen height) gives the pierce point.
    # The parameterization is given by:
    #     | antx |       | rotated_param_x |
    # l = | anty | + t * | rotated_param_y |
    #     | anyz |       | rotated_param_x |
    # whereas the shell representing the ionosphere is given by:
    #
    # x^2 + y^2   z^2
    # --------- + --- = 1, with:
    #   c^2      d^2
    # with c = a + height in meters
    # and d = b + height in meters

    c = a + screenheight
    d = b + screenheight

    t = -(1./(d**2.*(rotated_param_x**2. + rotated_param_y**2.) +
        c**2.*rotated_param_z**2.)) * (d**2.*(rotated_param_x*antx +
        rotated_param_y*anty) + c**2.*rotated_param_z*antz - 1/2.*np.sqrt(
        4*(d**2.*(rotated_param_x*antx + rotated_param_y*anty) +
        c**2.*rotated_param_z*antz)**2. - 4*(d**2.*(rotated_param_x**2. +
        rotated_param_y**2.) + c**2.*rotated_param_z**2.)*(d**2.*(-c**2. +
        antx**2. + anty**2.) + c**2.*antz**2.)))

    piercex = antx + t*rotated_param_x
    piercey = anty + t*rotated_param_y
    piercez = antz + t*rotated_param_z

    piercepoint_xyz = EarthLocation(x=piercex*u.m, y=piercey*u.m, z=piercez*u.m)
    piercepoint_lonlatheight = piercepoint_xyz.to_geodetic()

    return piercepoint_xyz

if __name__ == "__main__":
    RA = '12h22m54.899s'
    DEC = '+15d49m20.57s'
    lon = -111.598333
    lat = 31.956389
    height = 0
    time = '2010-01-01T12:00:00'

    print 'Found Pierce point:'
    print calculate_pierce_point(RA, DEC, lon, lat, height, time)

    print 'Comparison with existing SPAM code by Huib'
    gst = Time(time, scale='utc').sidereal_time('apparent', 'greenwich').deg
    tempobj = SkyCoord(RA, DEC)
    radeg = tempobj.ra.deg
    decdeg = tempobj.dec.deg
    antposition = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    x = antposition.x.value
    y = antposition.y.value
    z = antposition.z.value
    import sys
    sys.path.append('/net/student11/data1/kvdam/spam/python/spam')
    import sphere
    print sphere.calculate_pierce_point([x,y,z],[radeg, decdeg], gst)
