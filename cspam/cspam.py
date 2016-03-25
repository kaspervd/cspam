import sys
import os
import numpy as np
import logging
import ConfigParser
import inspect

"""
# Make sure that lib can be imported (module from a relative path)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from lib import AntennaObjects
from lib import TableObjects
from lib import utils
"""

# This above only works if you import casa functionalities separately in the
# above files. Somewhat ugly but it seems to work.
# Otherwise use:
execfile('/home/kvdam/Documents/MA-project-2/cspam/cspam/lib/AntennaObjects.py')
# Contains classes: RefAntHeuristics, RefAntGeometry, RefAntFlagging

execfile('/home/kvdam/Documents/MA-project-2/cspam/cspam/lib/TableObjects.py')
# Contains classes: MSObj, STObj

def get_conf(config_file='/home/kvdam/Documents/MA-project-2/cspam/cspam/cspam.config'):
    """
    Prepare a config dictionary with all the user-defined parameters
    """

    if not os.path.isfile(config_file):
        logging.critical('Configuration file '+config_file+' not found.')
        sys.exit(1)
    
    confP = ConfigParser.ConfigParser(
		{'flag':'','cal_scans':'','tgt_scans':''})
    confP.read(config_file)

    conf = {}
    conf['steps'] = confP.get('DEFAULT','steps').replace(' ','').split(',')
    conf['debug'] = confP.getint('DEFAULT','debug')
    conf['prog_dir'] = confP.get('DEFAULT','prog_dir')
    conf['data_dir'] = confP.get('DEFAULT','data_dir')

    # creating MSs from sections
    conf['MSs'] = confP.sections()
    for MS in conf['MSs']:
        conf[MS] = {}
        conf[MS]['file_path'] = conf['data_dir']+'/'+MS
        conf[MS]['flag'] = confP.get(MS, 'flag')
        conf[MS]['cal_scans'] = confP.get(
				MS, 'cal_scans').replace(' ','').split(',')
        conf[MS]['tgt_scans'] = confP.get(
				MS, 'tgt_scans').replace(' ','').split(',')
        conf[MS]['spw'] = confP.get(MS, 'spw').replace(' ','')
        conf[MS]['central_chan_percentage'] = confP.getint(
					      MS, 'central_chan_percentage')
   
    return conf

print "Starting CSPAM version: x.x"

# Print information on MSs
conf = get_conf()
print 'Imported conf file:'
print conf
print ' '

# creating the list of MSs and setting some specific values (as flags)
MSs = []
for MS in conf['MSs']:
    MSs.append( MSObj(conf[MS], MS))


###
### Below it's only testing
###

# In this test case we only have one measurement set
mset = MSs[0]

#
# Antenna position correction
#
gencal(vis=mset.file_path, caltable=mset.dir_cal+'/'+mset.ms_name+'.pos', 
       caltype='antpos')

#
# Flagging (first chan, quack, bad ant, bad time)
#
flagdata(vis=mset.file_path, mode='manual', # from tutorial
         antenna='ea06,ea17,ea20,ea26')
flagdata(vis=mset.file_path, mode='shadow', # from tutorial
         flagbackup=False)

#spw = '0:0'
#flagdata(vis=active_ms, mode='manualflag', spw=spw, flagbackup=False)

#if badranges != {}:
#        for badant in badranges:
#            logging.debug("Flagging :"+badant+" - time: "+badranges[badant])
#            default('flagdata')
#            flagdata(vis=active_ms, mode='manualflag', antenna=badant,\
#		timerange=badranges[badant], flagbackup=False)

# quack
flagdata(vis=mset.file_path, mode='quack', quackinterval=1, quackmode='beg',
         action='apply', flagbackup=False)
    
# flag zeros
flagdata(vis=mset.file_path, mode='clip', clipzeros=True,\
         correlation='ABS_ALL', action='apply', flagbackup=False)
    
# save flag status
flagmanager(vis=mset.file_path, mode='save',
            versionname='AfterStaticFlagging')

# First RFI removal
flagdata(vis=mset.file_path, mode='tfcrop', datacolumn='data',
         timecutoff = 4., freqcutoff = 3., maxnpieces = 7,\
         action='apply', flagbackup=False)

# save flag status
flagmanager(vis=mset.file_path, mode='save',
            versionname='AfterDynamicFlagging')

#
# Just a quick test to create a cal table to test plotcal
#
gaincal(vis=mset.file_path,
        caltable=mset.dir_cal+'/'+mset.ms_name+'.initPh',
        intent='CALIBRATE_PHASE*', solint='int',
        spw='0:10~13,1;3;5~6:30~33,2:32~35,4:35~38,7:46~49',
        refant='ea24', minblperant=3,
        minsnr=3.0, calmode='p',
        gaintable=mset.dir_cal+'/'+mset.ms_name+'.pos')

# Set models

# Bandpass calibration

# Calibration

# Self Calibration

# Imaging
