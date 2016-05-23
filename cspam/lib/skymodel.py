"""
This methods in this file come largely from the old AIPS based SPAM of 
Dr. Intema. Adapted by Kasper van Dam, Leiden Observatory 2016.
"""

###############################################################################

# import Python modules
import os
import inspect
import math
import numpy as np

###############################################################################

def is_array( a ):
  return isinstance( a, type( np.array( [ 1 ] ) ) )

###############################################################################

def amodulo( x, y ):
  if ( not is_array( x ) ):
    if ( not is_array( y ) ):
      m = x - y * math.floor( x / ( y + float( y == 0. ) ) )
    else:
      xx = x * aones( y )
      m = xx - y * math.floor( x / ( y + array( y == 0., dtype = y.dtype ) ) )
  else:
    if ( not is_array( y ) ):
      yy = y * aones( x )
      m = x - yy * math.floor( x / ( yy + array( yy == 0., dtype = yy.dtype ) ) )
    else:
      m = x - y * math.floor( x / ( y + array( y == 0., dtype = y.dtype ) ) )
  return m

###############################################################################

def calculate_angular_separation( degdeg1, degdeg2 ):
  ra1 = amodulo( degdeg1[ 0 ], 360. )
  dec1 = np.degrees( math.asin( 
      max( - 1., min( 1., math.sin( np.radians( amodulo( degdeg1[ 1 ], 360. ) ) ) ) ) ) )
  ra2 = amodulo( degdeg2[ 0 ], 360. )
  dec2 = np.degrees( math.asin( 
      max( - 1., min( 1., math.sin( np.radians( amodulo( degdeg2[ 1 ], 360. ) ) ) ) ) ) )
  dra = amodulo( ( ra2 - ra1 ) + 180., 360. ) - 180.
  ddec = dec2 - dec1
  if ( ( dra == 0. ) or ( amodulo( dec1, 180. ) == 90. ) or ( amodulo( dec2, 180. ) == 90. ) ):
    radius = abs( ddec )
    if ( dec2 >= dec1 ):
      angle = 0.
    else:
      angle = 180.
#  elif ( ddec == 0. ):
#    radius = dra * cos( radians( dec1 ) )
#    if ( radius > 0. ):
#      angle = 90.
#    else:
#      angle = 270.
#    radius = abs( radius )
  else:
    # use Haversine formula (adapted from http://www.plutoproject.com/dist_pa2.cpp)
    hav_r = ( math.sin( np.radians( ddec / 2. ) )**2 + 
        math.cos( np.radians( dec1 ) ) * math.cos( np.radians( dec2 ) ) * math.sin( np.radians( dra / 2. ) )**2 )
    radius = np.degrees( 2. * math.asin( max( - 1., min( 1., math.sqrt( hav_r ) ) ) ) )
    t = math.cos( np.radians( dec1 ) ) * math.sin( np.radians( radius ) )
    hav_a = ( math.sin( np.radians( dec1 + radius ) ) - math.sin( np.radians( dec2 ) ) ) / ( 2. * t )
    angle = np.degrees( 2. * math.asin( min( 1., math.sqrt( hav_a ) ) ) )
    if ( math.sin( np.radians( dra ) ) < 0. ):
      angle = 360. - angle
  return [ radius, angle ]

###############################################################################

def generate_source_list( radec, radius, freq, epoch = 2000.0, use_wenss = True,
    use_nvss = True, assoc_radius = 40. / 3600., spectral_index = None,
    flux_min = None, si_limits = [ -1.5, 0. ], use_vlss = False,):
# radius <= 90 degrees
# assoc_radius = 40. / 3600. # e.g. 2004MNRAS.352..909C
  
  # This is very ugly but it works
  path_of_this_file = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
  catalog_path = path_of_this_file+'/catalogs/'
  catalog_path_e = os.path.expandvars( catalog_path )
  if ( epoch == 2000.0 ):
    if ( radec[ 1 ] < -35. ):
#      nvss_file_name = catalog_path_e + 'SU00.0006'
      nvss_file_name = catalog_path_e + 'SM00.0006'
    else:
      nvss_file_name = catalog_path_e + 'NV00.0003'
    if ( not os.path.isfile( nvss_file_name ) ):
      if ( radec[ 1 ] < -35. ):
#        nvss_file_name = catalog_path_e + 'SU00.0050'
        nvss_file_name = catalog_path_e + 'SM00.0050'
      else:
        nvss_file_name = catalog_path_e + 'NV00.0030'
    wenss_file_name = catalog_path_e + 'WE00.0000'
    if ( not os.path.isfile( wenss_file_name ) ):
      wenss_file_name = catalog_path_e + 'WE00.0100'
    vlss_file_name = catalog_path_e + 'VL00.0400'
    if ( not os.path.isfile( vlss_file_name ) ):
      vlss_file_name = None
  elif ( epoch == 1950.0 ):
    nvss_file_name = catalog_path_e + 'NV50.0003'
    if ( not os.path.isfile( nvss_file_name ) ):
      nvss_file_name = catalog_path_e + 'NV50.0030'
    wenss_file_name = catalog_path_e + 'WE50.0000'
    if ( not os.path.isfile( wenss_file_name ) ):
      wenss_file_name = catalog_path_e + 'WE50.0100'
    vlss_file_name = None
  if ( ( not os.path.isfile( nvss_file_name ) ) or ( not os.path.isfile( wenss_file_name ) ) ):
    raise error( 'no catalog files available' )
  
  if ( not spectral_index is None ):
    if ( ( spectral_index < si_limits[ 0 ] ) or ( spectral_index > si_limits[ 1 ] ) ):
      raise error( 'default spectral index outside allowed spectral index range' )

  nvss_freq = 1.4e9
  if ( radec[ 1 ] < -35. ):
    nvss_freq = 843.e6
  wenss_freq = 326.e6
  if ( radec[ 1 ] < 0. ):
    wenss_freq = 352.e6
  vlss_freq = 73.8e6

  if ( not flux_min is None ):
    if ( nvss_freq > freq ):
      nvss_flux_min = flux_min * ( nvss_freq / freq )**si_limits[ 0 ]
    else:
      nvss_flux_min = flux_min * ( nvss_freq / freq )**si_limits[ 1 ]
    if ( wenss_freq > freq ):
      wenss_flux_min = flux_min * ( wenss_freq / freq )**si_limits[ 0 ]
    else:
      wenss_flux_min = flux_min * ( wenss_freq / freq )**si_limits[ 1 ]
  else:
    nvss_flux_min = 0.
    wenss_flux_min = 0.

  # determine crude selection criteria
  ra = radec[ 0 ]
  dec = radec[ 1 ]
  ra_range_contains_zero = False
  if ( dec + radius >= 90. ) or ( dec - radius <= -90. ):
    ra_min = 0.
    ra_max = 360.
    if ( dec + radius >= 90. ):
      dec_min = min( [ dec - radius, 180. - ( dec + radius ) ] )
      dec_max = 90.
    else:
      dec_min = - 90.
      dec_max = max( [ dec + radius, - 180. - ( dec - radius ) ] )
  else:
    cos_dec = math.cos( np.radians( max( [ abs( dec + radius ), abs( dec - radius ) ] ) ) )
    if ( ( radius / cos_dec ) >= 180. ):
      ra_min = 0.
      ra_max = 360.
    else:
      ra_min = amodulo( ra - ( radius / cos_dec ), 360. )
      ra_max = amodulo( ra + ( radius / cos_dec ), 360. )
      if ( ra_min > ra_max ):
        ra_range_contains_zero = True
    dec_min = dec - radius
    dec_max = dec + radius

  # get NVSS sources
  nvss_source_list = []
  if use_nvss:
    nvss_file = file( nvss_file_name, mode = 'r' )
    for line in nvss_file:
      first_index = 0
      while ( line[ first_index ] == ' ' ):
        first_index = first_index + 1
      line = line[ first_index : ]
      if ( len( line[ first_index : ] ) == 0 ):
        continue
      if ( line[ first_index ] == ';' ):
        continue
      if ( line[ first_index : first_index + 5 ] == 'Found' ):
        continue
      try:
        words = line.split()
        nvss_source = [ [ float( words[ 0 ] ), float( words[ 1 ] ) ],
            float( int( words[ 2 ] ) ) / 1.e3 ]
      except:
        try:
          # F9.5,1X,F9.5,I7,F10.4
          nvss_source = [ [ float( line[ 0 : 9 ] ), float( line[ 10 : 19 ] ) ],
              float( int( line[ 19 : 26 ] ) ) / 1.e3 ]
        except:
          continue
      if ( ( nvss_source[ 0 ][ 1 ] >= dec_min ) and 
          ( nvss_source[ 0 ][ 1 ] <= dec_max ) and
          ( nvss_source[ 1 ] >= nvss_flux_min ) ):
        if ra_range_contains_zero:
          if ( ( nvss_source[ 0 ][ 0 ] >= ra_max ) or
              ( nvss_source[ 0 ][ 0 ] <= ra_min ) ):
            [ nvss_radius, nvss_angle ] = calculate_angular_separation( 
                nvss_source[ 0 ], [ ra, dec ] )
            if ( nvss_radius <= radius ):
              nvss_source_list.append( nvss_source )
        else:
          if ( ( nvss_source[ 0 ][ 0 ] >= ra_min ) and
              ( nvss_source[ 0 ][ 0 ] <= ra_max ) ):
            [ nvss_radius, nvss_angle ] = calculate_angular_separation(
                nvss_source[ 0 ], [ ra, dec ] )
            if ( nvss_radius <= radius ):
              nvss_source_list.append( nvss_source )
    nvss_file.close()
    nvss_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )

  # get WENSS sources
  wenss_source_list = []
  if use_wenss:
    wenss_file = file( wenss_file_name, mode = 'r' )
    for line in wenss_file:
#      columns = [ column.strip() for column in line.split() ]
#      if ( len( columns ) == 3 ):
#        if ( columns[ 0 ][ 0 ] != ';' ):
#          wenss_source_data = [ float( column ) for column in columns ]
#          wenss_source = [ [ wenss_source_data[ 0 ], wenss_source_data[ 1 ] ],
#              wenss_source_data[ 2 ] / 1.e3 ]
      first_index = 0
      while ( line[ first_index ] == ' ' ):
        first_index = first_index + 1
      line = line[ first_index : ]
      if ( len( line[ first_index : ] ) == 0 ):
        continue
      if ( line[ first_index ] == ';' ):
        continue
      if ( line[ first_index : first_index + 5 ] == 'Found' ):
        continue
      try:
        words = line.split()
        wenss_source = [ [ float( words[ 0 ] ), float( words[ 1 ] ) ],
            float( int( words[ 2 ] ) ) / 1.e3 ]
      except:
        try:
          # F9.5,1X,F9.5,I7,F10.4
          wenss_source = [ [ float( line[ 0 : 9 ] ), float( line[ 10 : 19 ] ) ],
              float( int( line[ 19 : 26 ] ) ) / 1.e3 ]
        except:
          continue
      if ( ( wenss_source[ 0 ][ 1 ] >= dec_min ) and
          ( wenss_source[ 0 ][ 1 ] <= dec_max ) and
          ( wenss_source[ 1 ] >= wenss_flux_min ) ):
        if ra_range_contains_zero:
          if ( ( wenss_source[ 0 ][ 0 ] >= ra_max ) or
              ( wenss_source[ 0 ][ 0 ] <= ra_min ) ):
            [ wenss_radius, wenss_angle ] = calculate_angular_separation(
                wenss_source[ 0 ], [ ra, dec ] )
            if ( wenss_radius <= radius ):
              wenss_source_list.append( wenss_source )
        else:
          if ( ( wenss_source[ 0 ][ 0 ] >= ra_min ) and
              ( wenss_source[ 0 ][ 0 ] <= ra_max ) ):
            [ wenss_radius, wenss_angle ] = calculate_angular_separation(
                wenss_source[ 0 ], [ ra, dec ] )
            if ( wenss_radius <= radius ):
              wenss_source_list.append( wenss_source )
    wenss_file.close()
    wenss_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )

  # get VLSS sources
  vlss_source_list = []
  if ( use_vlss and ( not vlss_file_name is None ) ):
    vlss_file = file( vlss_file_name, mode = 'r' )
    for line in vlss_file:
#      columns = [ column.strip() for column in line.split() ]
#      if ( len( columns ) == 3 ):
#        if ( columns[ 0 ][ 0 ] != ';' ):
#          vlss_source_data = [ float( column ) for column in columns ]
#          vlss_source = [ [ vlss_source_data[ 0 ], vlss_source_data[ 1 ] ],
#              vlss_source_data[ 2 ] / 1.e3 ]
      first_index = 0
      while ( line[ first_index ] == ' ' ):
        first_index = first_index + 1
      line = line[ first_index : ]
      if ( len( line[ first_index : ] ) == 0 ):
        continue
      if ( line[ first_index ] == ';' ):
        continue
      if ( line[ first_index : first_index + 5 ] == 'Found' ):
        continue
      try:
        words = line.split()
        vlss_source = [ [ float( words[ 0 ] ), float( words[ 1 ] ) ],
            float( int( words[ 2 ] ) ) / 1.e3 ]
      except:
        try:
          # F9.5,1X,F9.5,I7,F10.4
          vlss_source = [ [ float( line[ 0 : 9 ] ), float( line[ 10 : 19 ] ) ],
              float( int( line[ 19 : 26 ] ) ) / 1.e3 ]
        except:
          continue
      if ( ( vlss_source[ 0 ][ 1 ] >= dec_min ) and
          ( vlss_source[ 0 ][ 1 ] <= dec_max ) ):
        if ra_range_contains_zero:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_max ) or
              ( vlss_source[ 0 ][ 0 ] <= ra_min ) ):
            [ vlss_radius, vlss_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], [ ra, dec ] )
            if ( vlss_radius <= radius ):
              vlss_source_list.append( vlss_source )
        else:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_min ) and
              ( vlss_source[ 0 ][ 0 ] <= ra_max ) ):
            [ vlss_radius, vlss_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], [ ra, dec ] )
            if ( vlss_radius <= radius ):
              vlss_source_list.append( vlss_source )
    vlss_file.close()
    vlss_source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )

  # look for WENSS sources in NVSS catalog
  # extrapolate source fluxes
  # check for additional VLSS source
  # build source list
  source_list = []
  i = 0
  for wenss_source in wenss_source_list:
    i = i + 1
    cross_source_list = []
    ra_min = wenss_source[ 0 ][ 0 ] - ( assoc_radius / math.cos( np.radians( wenss_source[ 1 ] ) ) )
    ra_max = wenss_source[ 0 ][ 0 ] + ( assoc_radius / math.cos( np.radians( wenss_source[ 1 ] ) ) )
    dec_min = wenss_source[ 0 ][ 1 ] - assoc_radius
    dec_max = wenss_source[ 0 ][ 1 ] + assoc_radius
    for nvss_source in nvss_source_list:
      if ( ( nvss_source[ 0 ][ 0 ] >= ra_min ) and ( nvss_source[ 0 ][ 0 ] <= ra_max ) and
           ( nvss_source[ 0 ][ 1 ] >= dec_min ) and ( nvss_source[ 0 ][ 1 ] <= dec_max ) ):
        [ cross_radius, cross_angle ] = calculate_angular_separation(
            nvss_source[ 0 ], wenss_source[ 0 ] )
        if ( cross_radius <= assoc_radius ):
          cross_source_list.append( nvss_source )
          nvss_source_list.remove( nvss_source )
    if ( len( cross_source_list ) > 0 ):
      nvss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
      source_sp_index = math.log( nvss_flux / wenss_source[ 1 ] ) / math.log( nvss_freq / wenss_freq )
      if ( source_sp_index < si_limits[ 0 ] ):
        source_sp_index = si_limits[ 0 ]
      elif ( source_sp_index > si_limits[ 1 ] ):
        source_sp_index = si_limits[ 1 ]
      # check for VLSS counterpart
      if ( len( vlss_source_list ) > 0 ):
        for vlss_source in vlss_source_list:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_min ) and ( vlss_source[ 0 ][ 0 ] <= ra_max ) and
             ( vlss_source[ 0 ][ 1 ] >= dec_min ) and ( vlss_source[ 0 ][ 1 ] <= dec_max ) ):
            [ cross_radius, cross_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], wenss_source[ 0 ] )
            if ( cross_radius <= assoc_radius ):
              cross_source_list.append( vlss_source )
              vlss_source_list.remove( vlss_source )
        if ( len( cross_source_list ) > 0 ):
          vlss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
          vlss_model_flux = nvss_flux * pow( vlss_freq / nvss_freq, source_sp_index )
          if ( vlss_flux < vlss_model_flux ):
            coef = 1. / ( ( 1. / vlss_flux ) - ( 1. / vlss_model_flux ) )
            source_flux_1 = nvss_flux * pow( freq / nvss_freq, source_sp_index )
            source_flux_2 = coef * pow( freq / vlss_freq, 2.5 )
            source_flux = 1. / ( ( 1. / source_flux_1 ) + ( 1. / source_flux_2 ) )
          else:
            source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
        else: # no VLSS counterfound found
          source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
      else: # no VLSS counterfound found
        source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
        source_list.append( [ wenss_source[ 0 ], source_flux ] )
    elif ( not spectral_index is None ): # no NVSS counterpart found
      if ( len( vlss_source_list ) > 0 ):
        for vlss_source in vlss_source_list:
          if ( ( vlss_source[ 0 ][ 0 ] >= ra_min ) and ( vlss_source[ 0 ][ 0 ] <= ra_max ) and
             ( vlss_source[ 0 ][ 1 ] >= dec_min ) and ( vlss_source[ 0 ][ 1 ] <= dec_max ) ):
            [ cross_radius, cross_angle ] = calculate_angular_separation(
                vlss_source[ 0 ], wenss_source[ 0 ] )
            if ( cross_radius <= assoc_radius ):
              cross_source_list.append( vlss_source )
              vlss_source_list.remove( vlss_source )
        if ( len( cross_source_list ) > 0 ):
          vlss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
          source_sp_index = math.log( vlss_flux / wenss_source[ 1 ] ) / math.log( vlss_freq / wenss_freq )
          if ( source_sp_index < si_limits[ 0 ] ):
            source_sp_index = si_limits[ 0 ]
          elif ( spectral_index > si_limits[ 1 ] ):
            source_sp_index = si_limits[ 1 ]
          source_flux = vlss_flux * pow( freq / vlss_freq, source_sp_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
        else: # no NVSS & VLSS counterpart found
          source_flux = wenss_source[ 1 ] * pow( freq / wenss_freq, spectral_index )
          source_list.append( [ wenss_source[ 0 ], source_flux ] )
      else: # no NVSS & VLSS counterpart found
        source_flux = wenss_source[ 1 ] * pow( freq / wenss_freq, spectral_index )
        source_list.append( [ wenss_source[ 0 ], source_flux ] )
  for vlss_source in vlss_source_list: # no WENSS counterpart found
    i = i + 1
    cross_source_list = []
    ra_min = vlss_source[ 0 ][ 0 ] - ( assoc_radius / math.cos( np.radians( vlss_source[ 1 ] ) ) )
    ra_max = vlss_source[ 0 ][ 0 ] + ( assoc_radius / math.cos( np.radians( vlss_source[ 1 ] ) ) )
    dec_min = vlss_source[ 0 ][ 1 ] - assoc_radius
    dec_max = vlss_source[ 0 ][ 1 ] + assoc_radius
    for nvss_source in nvss_source_list:
      if ( ( nvss_source[ 0 ][ 0 ] >= ra_min ) and ( nvss_source[ 0 ][ 0 ] <= ra_max ) and
           ( nvss_source[ 0 ][ 1 ] >= dec_min ) and ( nvss_source[ 0 ][ 1 ] <= dec_max ) ):
        [ cross_radius, cross_angle ] = calculate_angular_separation(
            nvss_source[ 0 ], vlss_source[ 0 ] )
        if ( cross_radius <= assoc_radius ):
          cross_source_list.append( nvss_source )
          nvss_source_list.remove( nvss_source )
    if ( len( cross_source_list ) > 0 ):
      nvss_flux = float( sum( [ cross_source[ 1 ] for cross_source in cross_source_list ] ) )
      source_sp_index = math.log( nvss_flux / vlss_source[ 1 ] ) / math.log( nvss_freq / vlss_freq )
      if ( source_sp_index < si_limits[ 0 ] ):
        source_sp_index = si_limits[ 0 ]
      elif ( source_sp_index > si_limits[ 1 ] ):
        source_sp_index = si_limits[ 1 ]
      source_flux = nvss_flux * pow( freq / nvss_freq, source_sp_index )
      source_list.append( [ vlss_source[ 0 ], source_flux ] )
    else: # no WENSS & NVSS counterpart
      pass
  for nvss_source in nvss_source_list: # no VLSS & WENSS counterpart found
    source_flux = nvss_source[ 1 ]
    if ( not spectral_index is None ):
      source_flux = source_flux * pow( freq / nvss_freq, spectral_index )
    source_list.append( [ nvss_source[ 0 ], source_flux ] )
  if ( len( source_list ) == 0 ):
    return source_list
  
  # apply minimum flux criteria
  source_list.sort( cmp = lambda a, b: cmp( b[ 1 ], a[ 1 ] ) )
  if ( not flux_min is None ):
    for source in source_list:
      if ( source[ 1 ] < flux_min ):
        break
    index = source_list.index( source )
    source_list = source_list[ 0 : index ]
  
  return source_list

###############################################################################

def get_primary_beam_attenuations( telescope, band, direction_radec, freq, source_radec_list, cutoff = 0.05 ):
  A_list = []
  higherparms = get_pbparms(telescope, band)
  pbparms = [ cutoff, 1. ] + higherparms
  for source_radec in source_radec_list:
    [ radius, angle ] = calculate_angular_separation( direction_radec, source_radec )
    A_list.append( calculate_pbparm_attenuation( freq, radius, pbparms ) )
  return A_list

###############################################################################

def calculate_pbparm_attenuation( freq, radius, pbparms ):
  [ cutoff, apply_pbparms, pbparm3, pbparm4, pbparm5, pbparm6, pbparm7 ] = pbparms
  if ( apply_pbparms > 0. ):
    X = ( ( freq / 1.e9 ) * ( radius * 60. ) )**2
    A = ( 1.0 + ( X * pbparm3 / 1.e3 ) + ( X**2 * pbparm4 / 1.e7 ) + 
        ( X**3 * pbparm5 / 1.e10 ) + ( X**4 * pbparm6 / 1.e13 ) + 
        ( X**5 * pbparm7 / 1.e16 ) )
    if ( A < cutoff ):
      A = 0.
    else:
      dXdR = 2 * ( ( freq / 1.e9 ) * ( radius * 60. ) ) * ( ( freq / 1.e9 ) * ( 60. ) )
      dAdX = ( ( pbparm3 / 1.e3 ) + ( 2 * X * pbparm4 / 1.e7 ) +
          ( 3 * X**2 * pbparm5 / 1.e10 ) + ( 4 * X**3 * pbparm6 / 1.e13 ) + 
          ( 5 * X**4 * pbparm7 / 1.e16 ) )
      if ( dAdX * dXdR > 0. ):
        A = 0.
  else:
    A = 1.
  return A

###############################################################################

def get_pb_attenuated_sources(radec, radius, freq, telescope, band):
  """
  RA, DEC in degrees: RA between 0 and 360, DEC between -90 and 90. Format: [ra,dec]
  Radius in degrees.
  Frequency in Hz.
  Telescope and band are both strings.
  The flux in the source list is returned in Jy.
  """
  # First generate the source list from the catalog
  source_list = generate_source_list(radec, radius, freq)
  source_radec_list = [i[0] for i in source_list]
  
  # Secondly, obtain the primary beam attenuations for this list
  attenuations = get_primary_beam_attenuations(telescope, band, radec, freq, source_radec_list)

  # Finally, multiply all fluxes in the source list with the attenuation
  corrected_source_list = []
  for i, j in zip(source_list, attenuations):
    corrected_source_list.append( [i[0], i[1]*j] )
  
  # Delete sources with a flux of zero
  corrected_source_list = [x for x in corrected_source_list if x[1] > 0]
  
  return corrected_source_list

###############################################################################

def get_pbparms(telescope, band):
  """
  This method simply returns the primary beam parameters per telescope and
  band.
  """
  if telescope == 'GMRT':
    if band == '1420':
      pbparm3 = -2.27961
      pbparm4 = 21.4611
      pbparm5 = -9.7929
      pbparm6 = 1.80153
      pbparm7 = 0
    elif band == '610':
      pbparm3 = -3.486
      pbparm4 = 47.749
      pbparm5 = -35.203
      pbparm6 = 10.399
      pbparm7 = 0
    elif band == '325':
      pbparm3 = -3.397
      pbparm4 = 47.192
      pbparm5 = -30.931
      pbparm6 = 7.803
      pbparm7 = 0
    elif band == '235':
      pbparm3 = -3.366
      pbparm4 = 46.159
      pbparm5 = -29.963
      pbparm6 = 7.529
      pbparm7 = 0
    elif band == '151':
      pbparm3 = 0
      pbparm4 = 0
      pbparm5 = 0
      pbparm6 = 0
      pbparm7 = 0
    else:
      print 'No primary beam parameters for this band.'
  elif telescope == 'EVLA':
    if band == 'P':
      pbparm3 = -0.935
      pbparm4 = 3.23
      pbparm5 = -0.378
      pbparm6 = 0
      pbparm7 = 0
    elif band == 'L':
      pbparm3 = -1.343
      pbparm4 = 6.579
      pbparm5 = -1.186
      pbparm6 = 0
      pbparm7 = 0
    else:
      print 'No primary beam parameters for this band.'
  else:
	print 'No primary beam parameters for this telescope.'
  
  return [ pbparm3, pbparm4, pbparm5, pbparm6, pbparm7 ]

###############################################################################
