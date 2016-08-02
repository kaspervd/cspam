import sys
import os
import numpy as np
import logging
import ConfigParser
import argparse

# CSPAM Modules
from lib import TableObjects
from lib import utils
import steps

def get_conf(config_file='cspam.config'):
    """
    Prepare a config dictionary with all the user-defined parameters
    """
    if not os.path.isfile(config_file):
        logging.critical('Configuration file '+config_file+' not found.')
        sys.exit(1)
    
    confP = ConfigParser.ConfigParser()
    confP.read(config_file)

    conf = {}
    conf['steps'] = confP.get('DEFAULT','steps').replace(' ','').split(',')
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
        conf[MS]['flux_cal_field'] = confP.get(MS, 'flux_cal_field')
        conf[MS]['phase_cal_field'] = confP.get(MS, 'phase_cal_field')
        conf[MS]['leakage_cal_field'] = confP.get(MS, 'leakage_cal_field')
        conf[MS]['position_cal_field'] = confP.get(MS, 'position_cal_field')
        conf[MS]['angle_cal_field'] = confP.get(MS, 'angle_cal_field')
    
    return conf

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser(description="""CSPAM stands of CASA: 
             Source Peeling and Atmospheric Modeling""")
    
    parser.add_argument("-c", "--config", help="""Path to the config file""", 
                        type=str,required=False)
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbosity", help="""Make PNG plot""")
    
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parseCommandLineArguments()
    verbose = args['verbosity'] # Don't use it yet
    config_file = args['config']
    
    logging.basicConfig(filename='cspam.log',level=logging.DEBUG)
    
    if config_file:
        conf = get_conf(config_file=config_file)
    else: # assume default config file (cspam.config)
        conf = get_conf()

    print "Starting CSPAM version: x.x"
    print "Parsed configuration file:"
    utils.print_dict(conf)

    # Creating the list of MSs
    MSs = []
    for MSinConf in conf['MSs']:
        MSs.append(TableObjects.MSObj(conf[MSinConf]['file_path'], conf=conf[MSinConf]))
        # Note that the MSObj class might update field names in the MS

    # Carry out the wanted steps per measurement set
    for mset in MSs:
        # Set the environment
        if not os.path.isdir(mset.dir_img):
            os.makedirs(mset.dir_img)
        if not os.path.isdir(mset.dir_plot):
            os.makedirs(mset.dir_plot)
        if not os.path.isdir(mset.dir_cal):
            os.makedirs(mset.dir_cal)
        if not os.path.isdir(mset.dir_peel):
            os.makedirs(mset.dir_peel)

        # Execute the wanted steps
        if 'plots' in conf['steps']:
            print 'Step: plots'
            steps.plots(mset)

        if 'preflag' in conf['steps']:
            print 'Step: preflag'
            steps.preflag(mset)

        if 'setjy' in conf['steps']:
            print 'Step: setjy'
            steps.set_flux_density_scale(mset)

        if 'bandpass' in conf['steps']:
            print 'Step: bandpass'
            steps.bandpass_calibration(mset)

        if 'cal' in conf['steps']:
            print 'Step: cal'
            steps.calib(mset)

        if 'selfcal' in conf['steps']:
            print 'Step: selfcal'
            steps.selfcal(mset)

        if 'peeling' in conf['steps']:
            print 'Step: peeling'
            steps.peeling(mset)

        if 'createimage' in conf['steps']:
            print 'Step: createimage'
            steps.createimage(mset)


