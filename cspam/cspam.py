import sys
import os
import numpy as np
import logging
import ConfigParser

# Here you typically want to do something like:
#
# from lib import AntennaObjects
# from lib import TableObjects
# from lib import utils
#
# Unfortunately, you lose the casa functionalities within these
# scripts if you do this. So, instead use:

execfile('lib/AntennaObjects.py')
# Contains classes: RefAntHeuristics, RefAntGeometry, RefAntFlagging

execfile('lib/TableObjects.py')
# Contains classes: MSObj, STObj

def get_conf(config_file='cspam.config'):
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
        conf[MS]['file_name'] = conf['data_dir']+'/'+MS
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
print conf

# creating the list of MSs and setting some specific values (as flags)
MSs = []
for MS in conf['MSs']:
    MSs.append( MSObj(conf[MS]) )

