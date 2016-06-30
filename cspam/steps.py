import os
import sys
import datetime
import logging
import numpy as np
import aplpy
import matplotlib.cm
from astropy.io import fits

# CSPAM Modules
from lib import AntennaObjects
from lib import TableObjects
from lib import utils
from lib import skymodel
from lib import peel

# The following casanova imports are somewhat ugly but importing this way
# allows for CASA-style usage of the toolkits and tasks.

# CASA Toolkits
import casac
cl = casac.casac.componentlist()

# CASA Tasks
from casat import plotants
plotants = plotants.plotants

from casat import plotms
plotms = plotms.plotms

from casat import listobs
listobs = listobs.listobs

from casat import flagdata
flagdata = flagdata.flagdata

from casat import flagmanager
flagmanager = flagmanager.flagmanager

from casat import setjy
setjy = setjy.setjy

from casat import gaincal
gaincal = gaincal.gaincal

from casat import smoothcal
smoothcal = smoothcal.smoothcal

from casat import bandpass
bandpass = bandpass.bandpass

from casat import plotcal
plotcal = plotcal.plotcal

from casat import applycal
applycal = applycal.applycal

from casat import split
split = split.split

from casat import imstat
imstat = imstat.imstat

from casat import ftw
ftw = ftw.ftw

from casat import clearcal
clearcal = clearcal.clearcal

from casat import exportfits
exportfits = exportfits.exportfits


sou_res = ['1arcsec']
sou_size = [5000]
expnoise = 1.e-6
rob=0.5


def plots(mset):
    """
    Create some plots and save them in the plot directory.
    Also create a file with the summary of the measurement set and
    save it in the plot directory.
    """
    logging.info("### CREATING PLOTS")
    plotants(vis=mset.file_path, figfile=mset.dir_plot+'/plotants.png')
    plotms(vis=mset.file_path, xaxis='time', yaxis='elevation', showgui=False,
           selectdata=True, antenna='0&1;2&3',spw='0:31', coloraxis='field', 
           plotfile=mset.dir_plot+'/el_vs_time.png', overwrite=True)
    listobs(vis=mset.file_path, verbose=True,
            listfile=mset.dir_plot+'/listobs.txt', overwrite=True)

def preflag(mset):
    """
    First flagging: quack, zeros, RFI 
    """
    logging.info("### FIRST FLAGGING")
    
    if mset.nchan == 512:
        spw='0:0'
        if mset.freq > 200e6 and mset.freq < 300e6:
            spw='0:0~130' # 235 MHz +20 border
    elif mset.nchan == 256:
        spw='0:0'
        if mset.freq > 200e6 and mset.freq < 300e6:
            spw='0:0~65' # 235 MHz +20 border
    elif mset.nchan == 128:
        spw='0:0'
        if mset.freq > 200e6 and mset.freq < 300e6:
            spw='0:0~65' # 235 MHz +20 border
    elif mset.nchan == 64:
        spw='0:0'
        if mset.freq > 200e6 and mset.freq < 300e6:
            spw='0:0~32' # 235 MHz +20 border
    
    else:
        logging.error('Cannot understand obs type.')
        sys.exit(1)

    flagdata(vis=mset.file_path, mode='manualflag', spw=spw, flagbackup=False)
    
    # quack: flag the quackinterval (in seconds) at the beginning of each scan
    flagdata(vis=mset.file_path, mode='quack', quackinterval=1,
             quackmode='beg', action='apply', flagbackup=False)
    
    # flag zeros
    flagdata(vis=mset.file_path, mode='clip', clipzeros=True,
             correlation='ABS_ALL', action='apply', flagbackup=False)
    
    # save flag status
    flagmanager(vis=mset.file_path, mode='save',
                versionname='AfterStaticFlagging', 
                comment=str(datetime.datetime.now()))

    # First RFI removal
    flagdata(vis=mset.file_path, mode='tfcrop', datacolumn='data',
             timecutoff = 4., freqcutoff = 3., maxnpieces = 7,
             action='apply', flagbackup=False)

    # save flag status
    flagmanager(vis=mset.file_path, mode='save',
                versionname='AfterDynamicFlagging',
                comment=str(datetime.datetime.now())) 

def set_flux_density_scale(mset):
    """
    Set the model visibility amp and phase of the calibrators
    """
    logging.info("### SETJY")
    
    setjy(vis=mset.file_path, scan=mset.fluxcalibrator.scans, 
          standard='Scaife-Heald 2012', usescratch=True,
          scalebychan=True)

def bandpass_calibration(mset):
    """
    This method compensates for the change of gain with frequency using a
    calibrator source.
    """

    if not os.path.isdir(mset.dir_cal+mset.fluxcalibrator.extend_dir):
        os.makedirs(mset.dir_cal+mset.fluxcalibrator.extend_dir)
    if not os.path.isdir(mset.dir_plot+mset.fluxcalibrator.extend_dir):
        os.makedirs(mset.dir_plot+mset.fluxcalibrator.extend_dir)
        
    # In order for plotcal to work a symbolic link is needed.
    # Plotcal assumes the measurement set is in the same directory
    # as the cal table.
    syscommand = 'ln -s '+mset.file_path+' '+\
                 mset.dir_cal+mset.fluxcalibrator.extend_dir+'/'+mset.ms_name+\
                 '.ms'
    os.system(syscommand)

    # Select a spectral window: 0:10~20 mean window 0, channels 10 to 20
    if mset.nchan == 512: initspw = '0:240~260'
    elif mset.nchan == 256: initspw = '0:120~130'
    elif mset.nchan == 128: initspw = '0:70~80'
    elif mset.nchan == 64: initspw = '0:30~50'

    for step in ['cycle1','final']:
        # Do two bandpass calibration steps per calibrator.
        logging.info("Start bandpass step: "+step)

        gaintables=[]
        interp=[]

        refAntObj = AntennaObjects.RefAntHeuristics(vis=mset.file_path, 
                                   field=mset.fluxcalibrator.field_name, geometry=True, 
                                   flagging=True)
        refAnt = refAntObj.calculate()[0]
        
        # Select minimal signal to noise ratio
        if mset.freq < 500e6:
            minsnr=3.0
        else:
            minsnr=5.0
        
        caltablepath = mset.dir_cal+mset.fluxcalibrator.extend_dir+'/'
        
        # Start with phase calibration
        # The 'B' solutions are limited by the signal-to-noise ratio 
        # available per channel, which may be quite small. It is therefore 
        # important that the data be coherent over the time-range of the 
        # 'B' solutions. As a result, 'B' solutions are almost always
        # preceded by an initial 'G' or 'T' solve using gaincal.
        # 'G' is polarization dependent gain, 'T' is polarization
        # independent gain. So boot stands for startup.
        logging.info("Phase calibration")
        gaincal(vis=mset.file_path, caltable=caltablepath+step+'-boot.Gp', 
                field=mset.fluxcalibrator.field_name, selectdata=True, uvrange='>50m',
                scan=mset.fluxcalibrator.scans, spw=initspw, solint='int', combine='', 
                refant=refAnt, minblperant=mset.minBL_for_cal, minsnr=0, 
                calmode='p')
        bootGp = TableObjects.STObj(caltablepath+step+'-boot.Gp')
        mset.fluxcalibrator.derived_cal_tables.append(bootGp)
        
        # Smoothing solutions
        smoothcal(vis=mset.file_path, tablein=bootGp.file_path,
                  caltable=caltablepath+step+'-boot.Gp-smooth')
        bootGpsmooth = TableObjects.STObj(caltablepath+step+'-boot.Gp-smooth')
        mset.fluxcalibrator.derived_cal_tables.append(bootGpsmooth)
        
        # Initial bandpass correction.
        # Calibration type 'B' differs from 'G' only in that it is 
        # determined for each channel in each spectral window.
        logging.info("Bandpass calibration 1")
        bandpass(vis=mset.file_path, caltable=caltablepath+step+'-boot.B', 
                 field=mset.fluxcalibrator.field_name, selectdata=True, uvrange='>100m', 
                 scan=mset.fluxcalibrator.scans, solint='inf', combine='scan,field', 
                 refant=refAnt, minblperant=mset.minBL_for_cal, 
                 minsnr=minsnr, solnorm=True, bandtype='B', 
                 gaintable=[bootGpsmooth.file_path], 
                 interp=['linear'])
        bootB = TableObjects.STObj(caltablepath+step+'-boot.B')
        mset.fluxcalibrator.derived_cal_tables.append(bootB)

        # Find leftover time-dependent delays, Type 'K' solves for simple
        # antenna-based delays via Fourier transforms of the spectra on
        # baselines to the reference antenna. 
        logging.info("BP: Delay calibration")
        gaincal(vis=mset.file_path, caltable=caltablepath+step+'.K', 
                field=mset.fluxcalibrator.field_name, selectdata=True, uvrange='>100m',
                scan=mset.fluxcalibrator.scans, solint='int',combine='', refant=refAnt, 
                interp=interp+['nearest,nearestflag'], gaintype='K',
                minblperant=mset.minBL_for_cal, minsnr=minsnr,
                gaintable=gaintables+[bootB.file_path])
        K = TableObjects.STObj(caltablepath+step+'.K')
        mset.fluxcalibrator.derived_cal_tables.append(K)
        
        # flag outliers
        utils.FlagCal(K.file_path, sigma = 5, cycles = 3)
        
        # plot
        K.plot(mset.dir_plot+mset.fluxcalibrator.extend_dir)
        
        gaintables.append(K.file_path)
        interp.append('linear')
        
        # Find time-dependant gains, now taking both amplitude and phase
        # into account.
        logging.info("BP: Gain calibration")
        gaincal(vis=mset.file_path, caltable=caltablepath+step+'.Gap', 
                field=mset.fluxcalibrator.field_name, selectdata=True, uvrange='>100m', 
                scan=mset.fluxcalibrator.scans, solint='int',combine='', refant=refAnt, 
                interp=interp+['nearest,nearestflag'],
                minblperant=mset.minBL_for_cal, minsnr=minsnr,  
                gaintype='G', calmode='ap',
                gaintable=gaintables+[bootB.file_path])
        Gap = TableObjects.STObj(caltablepath+step+'.Gap')
        mset.fluxcalibrator.derived_cal_tables.append(K)
        
        # flag outliers
        utils.FlagCal(Gap.file_path, sigma = 3, cycles = 3)
        
        # plot
        Gap.plot(mset.dir_plot+mset.fluxcalibrator.extend_dir)
        
        gaintables.append(caltablepath+step+'.Gap')
        interp.append('linear')

        # Find cross-K, type 'KCROSS' solves for global cross-hand delays,
        # so this is again to find delays.
        logging.info("BP: Kcross calibration")
        gaincal(vis=mset.file_path, caltable=caltablepath+step+'.Kcross',
                field=mset.fluxcalibrator.field_name, selectdata=True, uvrange='>100m',
                scan=mset.fluxcalibrator.scans, solint='inf',combine='scan,field', 
                refant=refAnt, interp=interp+['nearest,nearestflag'],
                minblperant=mset.minBL_for_cal, minsnr=minsnr,
                gaintype='KCROSS',
                gaintable=gaintables+[bootB.file_path])
        Kcross = TableObjects.STObj(caltablepath+step+'.Kcross')
        mset.fluxcalibrator.derived_cal_tables.append(Kcross)
                
        # plot (custom plot, not in STObj class)
        plotcal(caltable = Kcross.file_path, xaxis = 'antenna', 
                yaxis = 'delay', showgui=False,
                figfile= mset.dir_plot+mset.fluxcalibrator.extend_dir+'/'+ \
                step+'.Kcross.png' )
        
        gaintables.append(caltablepath+step+'.Kcross')
        interp.append('nearest')

        # Recalculate bandpass taking delays into account.
        logging.info("Bandpass calibration 2")
        bandpass(vis=mset.file_path, caltable=caltablepath+step+'.B',
                 field=mset.fluxcalibrator.field_name, selectdata=True, uvrange='>100m', 
                 scan=mset.fluxcalibrator.scans, solint='inf', combine='scan,field',
                 refant=refAnt, interp=interp, 
                 minblperant=mset.minBL_for_cal, minsnr=minsnr, 
                 solnorm=True, bandtype='B', gaintable=gaintables)
        B = TableObjects.STObj(caltablepath+step+'.B')
        mset.fluxcalibrator.derived_cal_tables.append(B)
                 
        # plot
        B.plot(mset.dir_plot+mset.fluxcalibrator.extend_dir)
                        
        gaintables.append(B.file_path)
        interp.append('nearest,nearestflag')

        # Apply the gaintables from this bandpass calibration cycle.
        logging.info("Apply bandpass")
        applycal(vis=mset.file_path, selectdata=True, field=mset.fluxcalibrator.field_name,
                 scan=mset.fluxcalibrator.scans, gaintable=gaintables, calwt=False,
                 flagbackup=False, interp=interp)
        
        if step != 'final':
            # clip on residuals
            utils.clipresidual(mset.file_path, f=mset.fluxcalibrator.field_name, s=mset.fluxcalibrator.scans)


    # remove K, amp from gaintables, we keep B and Kcross which are global 
    # and T-indep
    gaintables=[caltablepath+'final.B', caltablepath+'final.Kcross']
    interp=['nearest,nearestflag','nearest,nearestflag']
    
    utils.statsFlag(mset.file_path, note='Before apply bandpass')

    # apply final bandpass to target(s)
    for target in mset.targetsources:
        applycal(vis=mset.file_path, selectdata=True, field=target.field_name, 
                 scan=target.scans, gaintable=gaintables, calwt=False, 
                 flagbackup=False, interp=interp)
        # calibrator is already corrected (also with G and K, not a big deal)

    utils.statsFlag(mset.file_path, note='After apply bandpass, before rflag')

    # run the final flagger
    flagdata(vis=mset.file_path, mode='rflag', ntime='scan',
             combinescans=False, datacolumn='corrected', winsize=3,
             timedevscale=5, freqdevscale=5, action='apply', flagbackup=False)

    # flag statistics after flagging
    utils.statsFlag(mset.file_path, note='After rflag')

def calib(mset):
    """
    This method compensates for the change of gain with time using a
    calibrator source.
    """
    
    logging.info("### CALIB")
    
    caltablepathg = mset.dir_cal+mset.phasecalibrator.extend_dir+'/'
    caltablepathf = mset.dir_cal+mset.fluxcalibrator.extend_dir+'/'

    if not os.path.isdir(mset.dir_cal+mset.phasecalibrator.extend_dir):
        os.makedirs(mset.dir_cal+mset.phasecalibrator.extend_dir)
    if not os.path.isdir(mset.dir_plot+mset.phasecalibrator.extend_dir):
        os.makedirs(mset.dir_plot+mset.phasecalibrator.extend_dir)
    if not os.path.isdir(mset.dir_img+mset.phasecalibrator.extend_dir):
        os.makedirs(mset.dir_img+mset.phasecalibrator.extend_dir)
        
    # In order for plotcal to work a symbolic link is needed.
    # Plotcal assumes the measurement set is in the same directory
    # as the cal table.
    syscommand = 'ln -s '+mset.file_path+' '+mset.dir_cal+ \
                 mset.phasecalibrator.extend_dir+'/'+mset.ms_name+'.ms'
    os.system(syscommand)

    gaintables=[]
    interp=[]
    n_cycles = 3
    for cycle in xrange(n_cycles):
        # Three cycles of phase calibration, time delay calibration,
        # (another) phase calibration and amplitude calibration.
        logging.info("Start CALIB cycle: "+str(cycle))

        refAntObj = AntennaObjects.RefAntHeuristics(vis=mset.file_path, 
                                   field=mset.phasecalibrator.field_name, geometry=True, 
                                   flagging=True)
        refAnt = refAntObj.calculate()[0]
        
        # Reset these for this cycle
        gaintables=[caltablepathf+'final.B', caltablepathf+'final.Kcross']
        interp=['nearest,nearestflag','nearest,nearestflag']
        
        # Gain cal phase
        if mset.freq < 500e6:
            minsnr=2.0
        else:
            minsnr=4.0
        gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                str(cycle)+'-boot.Gp', field=mset.phasecalibrator.field_name, selectdata=True,
                uvrange='>100m', scan=mset.phasecalibrator.scans, solint='int', refant=refAnt, 
                interp=interp, minblperant=mset.minBL_for_cal, 
                minsnr=minsnr, calmode='p', gaintable=gaintables)
        bootGp = TableObjects.STObj(caltablepathg+'gain'+str(cycle)+'-boot.Gp')
        mset.phasecalibrator.derived_cal_tables.append(bootGp)

        # Find leftover time-dependent delays, Type 'K' solves for simple
        # antenna-based delays via Fourier transforms of the spectra on
        # baselines to the reference antenna. 
        gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                str(cycle)+'.K', field=mset.phasecalibrator.field_name, selectdata=True,
                uvrange='>100m', scan=mset.phasecalibrator.scans, solint='int', refant=refAnt,
                minblperant=mset.minBL_for_cal, minsnr=minsnr, gaintype='K',
                interp=interp+['linear'],gaintable=gaintables+ \
                [bootGp.file_path])
        K = TableObjects.STObj(caltablepathg+'gain'+str(cycle)+'.K')
        mset.phasecalibrator.derived_cal_tables.append(K)
        
        # Flag
        utils.FlagCal(K.file_path, sigma = 5, cycles = 3)
        
        # Plot
        K.plot(mset.dir_plot+mset.phasecalibrator.extend_dir)

        # Redo gain cal phase (now with K)
        gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                str(cycle)+'.Gp', field=mset.phasecalibrator.field_name, selectdata=True,
                uvrange='>100m', scan=mset.phasecalibrator.scans, solint='int', refant=refAnt,
                minblperant=mset.minBL_for_cal, minsnr=minsnr, calmode='p', 
                gaintable=gaintables+[K.file_path],
                interp=interp+['linear'])
        Gp = TableObjects.STObj(caltablepathg+'gain'+str(cycle)+'.Gp')
        mset.phasecalibrator.derived_cal_tables.append(Gp)
                
        # Smoothing solutions
        smoothcal(vis=mset.file_path, tablein=caltablepathg+'gain'+ \
                  str(cycle)+'.Gp', caltable=caltablepathg+'gain'+ \
                  str(cycle)+'.Gp-smooth')
        Gpsmooth = TableObjects.STObj(caltablepathg+'gain'+str(cycle)+'.Gp-smooth')
        mset.phasecalibrator.derived_cal_tables.append(Gpsmooth)
                  
        #Plot
        Gpsmooth.plot(mset.dir_plot+mset.phasecalibrator.extend_dir, phase_only=True)
                          
        gaintables.append(Gpsmooth.file_path)
        interp.append('linear')

        # Gain cal amp
        if mset.freq < 500e6:
            minsnr=3.0
        else:
            minsnr=5.0
        gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                str(cycle)+'.Ga', field=mset.phasecalibrator.field_name, selectdata=True,
                uvrange='>100m', scan=mset.phasecalibrator.scans, solint='60s', minsnr=minsnr,
                refant=refAnt, minblperant=mset.minBL_for_cal, calmode='a', 
                gaintable=gaintables)
        Ga = TableObjects.STObj(caltablepathg+'gain'+str(cycle)+'.Ga')
        mset.phasecalibrator.derived_cal_tables.append(Ga)
        
        # Flag        
        utils.FlagCal(Ga.file_path, sigma = 3, cycles = 3)
        
        # Plot
        Ga.plot(mset.dir_plot+mset.phasecalibrator.extend_dir, amp_only=True)
        
        gaintables.append(Ga.file_path)
        interp.append('linear')

        applycal(vis=mset.file_path, field=mset.phasecalibrator.field_name, scan=mset.phasecalibrator.scans,
                 gaintable=gaintables, interp=interp, calwt=False, 
                 flagbackup=False)
        
        # clip of residuals not on the last cycle (useless and prevent
        # imaging of calibrator)
        if cycle != n_cycles-1:
            utils.clipresidual(mset.file_path, f=mset.phasecalibrator.field_name, s=mset.phasecalibrator.scans)

    # make a test img of the gain cal to check that everything is fine
    parms = {'vis':mset.file_path, 'field':mset.phasecalibrator.field_name,
             'imagename':mset.dir_img+mset.phasecalibrator.extend_dir+'/'+ \
             mset.phasecalibrator.field_name+'_gcal',
             'gridmode':'widefield', 'wprojplanes':128,
             'mode':'mfs', 'nterms':2, 'niter':1000, 'gain':0.1,
             'psfmode':'clark', 'imagermode':'csclean', 'imsize':512,
             'cell':sou_res, 'weighting':'briggs', 'robust':0,
             'usescratch':False}
    utils.cleanmaskclean(parms, makemask=False)
    
    # Apply to all targets
    # Use a different cycle to compensate for messing up with uvsub during the 
    # calibration of other sources in this way the CRRECTED_DATA are OK for 
    # all fields
    for target in mset.targetsources:
        applycal(vis=mset.file_path, field=target.field_name, scan=target.scans,
                 gaintable=gaintables, interp=interp, calwt=False,
                 flagbackup=False)

def selfcal(mset):
    """
    This method compensates for the change of gain with time using the data
    itself (and possibly a skymodel).
    """
    logging.info("### SELFCAL")

    # First, set the width to be used for splitting the target from the
    # measurement set. This width is the number of channels to average to form
    # one output channel.
    if mset.freq > 1000e6: width = 16
    if mset.freq > 550e6 and mset.freq < 650e6: width = 16
    if mset.freq > 300e6 and mset.freq < 350e6: width = 8
    if mset.freq > 200e6 and mset.freq < 300e6: width = 8
    # renormalize if chans were not 512, force int to prevent bug in split() 
    # if width is a numpy.int64
    width = int(width / (512/mset.nchan))
    logging.info("Average with width="+str(width))
    
    # Perform self calibration on each target field individually
    for target in mset.targetsources:
    
        # Create the directories
        if not os.path.isdir(mset.dir_cal+target.extend_dir):
            os.makedirs(mset.dir_cal+target.extend_dir)
        if not os.path.isdir(mset.dir_plot+target.extend_dir):
            os.makedirs(mset.dir_plot+target.extend_dir)
        if not os.path.isdir(mset.dir_img+target.extend_dir):
            os.makedirs(mset.dir_img+target.extend_dir)
    
        # Split the target field from the old measurement set
        target_file_path = mset.file_path[:-3]+'_target_'+target.field_name+'.ms'
        split(vis=mset.file_path, outputvis=target_file_path, scan=target.scans,
              field=target.field_name, width=width, datacolumn='corrected',
              keepflags=False)
        
        # Create newinstance of MSObj for this target and tell the original 
        # mset that this target has been split off
        target_mset = TableObjects.MSObj(target_file_path)
        mset.targetmsets.append(target_mset)

        # In order for plotcal to work a symbolic link is needed.
        # Plotcal assumes the measurement set is in the same directory
        # as the cal table.
        syscommand = 'ln -s '+mset.file_path[:-3]+'_target_'+target.field_name+\
                     '.ms '+mset.dir_cal+target.extend_dir+'/'+\
                     mset.ms_name+'_target_'+target.field_name+'.ms'
        os.system(syscommand)
    
        # Perform several self calibration cycles.
        for cycle in xrange(6):
     
            logging.info("Start SELFCAL cycle: "+str(cycle))
            
            # save flag for recovering
            flagmanager(vis=target_file_path, mode='save', 
                        versionname='selfcal-c'+str(cycle))

            ts = str(expnoise*10*(5-cycle))+' Jy' # expected noise this cycle
            
            # Store different iterations in different directories
            if not os.path.isdir(mset.dir_img+target.extend_dir+'/'\
                                 +target.field_name+str(cycle)):
                os.makedirs(mset.dir_img+target.extend_dir+'/'\
                            +target.field_name+str(cycle))
            
            # The first step is to create an image, create a mask for the
            # background and re-image with the background masked. This is done
            # because clean tends to perform better, and is less likely to
            # diverge, if the clean component placement is limited by a mask
            # where real emission is expected to be.
            parms = {'vis':target_file_path, 'imagename':mset.dir_img+ \
                     target.extend_dir+'/'+target.field_name+str(cycle)+\
                     '/im',
                     'gridmode':'widefield', 'wprojplanes':512, 'mode':'mfs', 
                     'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 
                     'imagermode':'csclean', 'imsize':sou_size, 'cell':sou_res, 
                     'weighting':'briggs', 'robust':rob, 'usescratch':True, 
                     'mask':'', 'threshold':ts, 'multiscale':[]}
            utils.cleanmaskclean(parms)

            # Get the rms from this image, this is needed to compare different 
            # cycles. "< 1" is to invert the mask, so that the rms is caculated
            # only for the background. The syntax is very weird, see:
            # https://casa.nrao.edu/aips2_docs/notes/223/index.shtml
            mask = '"'+mset.dir_img+target.extend_dir+'/'\
                   +target.field_name+str(cycle)+'/im.newmask'+'"'+' < 1'
            rms = imstat(imagename=mset.dir_img+target.extend_dir+'/'\
                         +target.field_name+str(cycle)+'/im-masked.image.tt0',
                         mask=mask)['rms'][0]
                         
            # If there was a previous cycle, check if there is improvement,
            # otherwise quit the self calibration.
            if cycle != 0 and old_rms * 1.1 < rms:
                logging.warning('Image rms noise ('+str(rms)+' Jy/b) is higher \
                than previous cycle ('+str(old_rms)+' Jy/b). Apply old cal \
                tables and quitting selfcal.')

                # get previous flags
                flagmanager(vis=target_file_path, mode='restore',
                            versionname='selfcal-c'+str(cycle-1))
                
                # for the first cycle just remove all calibration i.e. no selfcal
                # for the others get the previous cycle tables.
                if cycle == 1:
                    clearcal(vis=target_file_path)
                # For the first few cycles only phase calibration is applied.
                elif cycle < 4:
                    Gp_prev = target.self_cal_gp_tables[-2]
                    Gp_prev.plot(mset.dir_plot+target.extend_dir+'/Gp_cycle'+str(cycle-2), phase_only=True)
                    
                    applycal(vis=target_file_path, gaintable=gaintable, 
                             interp=['linear','linear'], calwt=False,
                             flagbackup=False)
                # For the later cycles also the amplitude is used.    
                elif cycle >= 4: 
                    Gp_prev = target.self_cal_gp_tables[-2]
                    Gp_prev.plot(mset.dir_plot+target.extend_dir+'/Gp_cycle'+str(cycle-2))
                    
                    Ga_prev = target.self_cal_ga_tables[-2]
                    Ga_prev.plot(mset.dir_plot+target.extend_dir+'/Gp_cycle'+str(cycle-2), amp_only=True)
                    
                    applycal(vis=target_file_path, gaintable=gaintable, 
                             interp=['linear','linear'], calwt=False, 
                             flagbackup=False)
                # Since there is no improvement, quit the self cal loop.
                break
            # Log the new rms if it's better (or within a factor 1.1)
            elif cycle != 0:
                logging.info('Rms noise change: '+str(old_rms)+' Jy/b -> '\
                             +str(rms)+' Jy/b.')

            # Don't do one more useless calibration if we're in the last cycle
            if cycle == 5: break
            
            # If this is the first cycle: start with a source model created 
            # from a catalog (e.g. NVSS, WENSS)
            if cycle == 0:
                # Obtain the catalog
                direction = mset.get_direction_from_tgt_field_id(target.field_id)
                m0 = direction['m0']['value'] # radians
                m1 = direction['m1']['value'] # radians
                ra = np.degrees(m0)
                dec = np.degrees(m1)
                catalog_list = skymodel.get_pb_attenuated_sources([ra,dec], 3,
                                        mset.freq, mset.telescope, mset.band)

                # Create an image showing the earlier image with
                # the source from the catalog (just for inspection)
                fluxlist = [i[1] for i in catalog_list]
                normalizedfluxlist = fluxlist / max(fluxlist)
                colors = matplotlib.cm.cool(normalizedfluxlist)
                fitsfig = aplpy.FITSFigure(mset.dir_img+target.extend_dir+'/'+target.field_name+\
                                           str(cycle)+'/im.image.tt0.fits')
                fitsfig.show_grayscale()
                fitsfig.show_circles([i[0][0] for i in catalog_list], 
                                     [i[0][1] for i in catalog_list], 0.02, 
                                     color=colors, lw=2.5)
                fitsfig.save(mset.dir_img+target.extend_dir+'/'\
                             +target.field_name+str(cycle)+\
                             '/im.image.tt0-skymodel.png')
                
                # Add each source in the catalog to the component list
                for i in catalog_list:
                    
                    # Obtain direction in radians
                    m0source = np.radians(i[0][0])
                    m1source = np.radians(i[0][1])

                    # convert direction to casa format
                    direc = {'m0': {'value': m0source, 'unit': 'rad'},  
                             'm1': {'value': m1source, 'unit': 'rad'},  
                             'refer': 'J2000',
                             'type': 'direction'}
                    flux = i[1]
                    
                    # Add this source to the model component list
                    # dir is a python keyword: use a dictionary to
                    # parse it as a variable name
                    parms = {'flux':flux, 'fluxunit':'Jy', 'shape':'point',
                             'dir':direc, 'freq':mset.freq}
                    cl.addcomponent(**parms)
                
                # Delete catalog if it already exists
                syscommand='rm -rf '+mset.dir_img+target.extend_dir+\
                           '/'+target.field_name+str(cycle)+'/catalog_component.cl'
                os.system(syscommand)
                
                # Save the catalog
                cl.rename(mset.dir_img+target.extend_dir+'/'\
                          +target.field_name+str(cycle)+'/catalog_component.cl')
                           
                # Close the catalog
                cl.close()
                
                # Create the source model
                ftw(vis=target_file_path, complist=mset.dir_img+target.extend_dir+'/'+target.field_name+str(cycle)+\
                    '/catalog_component.cl', wprojplanes=512, usescratch=True)
            
            # If this is not the first cycle and there is improvement (we
            # didn't quit earlier), use the earlier masked image as a source 
            # model for the self calibration
            else:
                ftw(vis=target_file_path,
                    model=[mset.dir_img+target.extend_dir+'/'\
                           +target.field_name+str(cycle)+'/im-masked.model.tt0',
                           mset.dir_img+target.extend_dir+'/'\
                           +target.field_name+str(cycle)+'/im-masked.model.tt1'],
                    nterms = 2, wprojplanes=512, usescratch=True)
            
            ## Perform the actual self calbration
            # recalibrating
            refAntObj = AntennaObjects.RefAntHeuristics(vis=target_file_path,
                                       field='0', geometry=True, flagging=True)
            refAnt = refAntObj.calculate()[0][0]

            # Gaincal - phases
            if cycle==0:
                solint='30s'
                minsnr=4
            if cycle==1:
                solint='15s'
                minsnr=3
            if cycle==2:
                solint='int'
                minsnr=3
            if cycle==3:
                solint='int'
                minsnr=2
            if cycle==4:
                solint='int'
                minsnr=2

            if mset.freq < 400e6:
                minsnr -= 1.

            gaincal(vis=target_file_path, 
                    caltable=mset.dir_cal+target.extend_dir+'/'\
                             +target.field_name+str(cycle)+'.Gp',
                    solint=solint, minsnr=minsnr, selectdata=True, 
                    uvrange='>50m', refant=refAnt, 
                    minblperant=mset.minBL_for_cal,
                    gaintable=[], calmode='p')
            Gp = TableObjects.STObj(mset.dir_cal+target.extend_dir+'/'+target.field_name+str(cycle)+'.Gp')
            target.self_cal_gp_tables.append(Gp)
           
            # Gaincal - amp (only for later cycles)
            if cycle >= 3:        
                if cycle==3: 
                    solint='30s'
                    minsnr = 4.
                if cycle==4: 
                    solint='15s'
                    minsnr = 3.
                if mset.freq < 400e6:
                    minsnr -= 1.
                gaincal(vis=target_file_path,
                        caltable=mset.dir_cal+target.extend_dir+'/'\
                                 +target.field_name+str(cycle)+'.Ga',
                        selectdata=True, uvrange='>50m', solint=solint,
                        minsnr=minsnr, refant=refAnt, gaintable=[],
                        minblperant=mset.minBL_for_cal, calmode='a')
                Ga = TableObjects.STObj(mset.dir_cal+target.extend_dir+'/'+target.field_name+str(cycle)+'.Ga')
                target.self_cal_ga_tables.append(Ga)
                        
                utils.FlagCal(Ga.file_path, sigma = 3, cycles = 3)
     
            # plot gains (amps only in later cycles)
            if cycle >= 3:
                Gp.plot(mset.dir_plot+target.extend_dir+'/Gp_cycle'+str(cycle))
            else:
                Gp.plot(mset.dir_plot+target.extend_dir+'/Gp_cycle'+str(cycle), phase_only=True)

            # add to gaintable (add amps only in later cycles)
            if cycle >= 3:
                gaintable=[Gp.file_path,Ga.file_path]
            else:
                gaintable=[Gp.file_path,]

            # Apply the newly found gains
            applycal(vis=target_file_path, field = '', gaintable=gaintable, 
                     interp=['linear','linear'], calwt=False, flagbackup=False)           
            utils.statsFlag(target_file_path, 
                            note='After apply selfcal (cycle: '+str(cycle)+')')
                            
            # For the next cycle, the current rms becomes the old_rms
            old_rms = rms

def peeling(mset):

    # TEMP
    mset.targetmsets.append(TableObjects.MSObj('/data2/kvdam/test_data/GMRT_testdata/XSshaped_610_target_XSSHAPED.ms'))
    # TEMP

    # Perform peeling on each target mset individually
    for target_mset in mset.targetmsets:
        # This measurement set only has one target
        tgt_fieldname = target_mset.targetsources[0].field_name
    
        # Determine the total number of integrations in this dataset
        total_integration_time = target_mset.summary['IntegrationTime']
        integration_times = []
        for cal_scan_id in target_mset.scansummary.keys():
            integration_times.append(target_mset.summary['scan_%i' % int(cal_scan_id)]
                                     ['0']['IntegrationTime'])
        # The integration time per scan may vary, the median is a good measure
        integration_time = np.median(integration_times)
        number_of_integrations = total_integration_time / integration_time
    
        extend_dir = '/'+tgt_fieldname
        # Create the directories
        if not os.path.isdir(mset.dir_peel+extend_dir):
            os.makedirs(mset.dir_peel+extend_dir)
        if not os.path.isdir(mset.dir_peel+extend_dir+'/pre'):
            os.makedirs(mset.dir_peel+extend_dir+'/pre')

        # Create the pre-peeling image before peeling
        path_to_currentimage = mset.dir_peel+extend_dir+'/pre/pre-peeling'
        parms = {'vis':target_mset.file_path, 'imagename':path_to_currentimage, 
                 'gridmode':'widefield', 'wprojplanes':512, 'mode':'mfs', 
                 'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 
                 'imagermode':'csclean', 'imsize':sou_size, 'cell':sou_res, 
                 'weighting':'briggs', 'robust':rob, 'usescratch':True,
                 'uvtaper':True, 'threshold':str(expnoise)+' Jy'}
        #utils.cleanmaskclean(parms)
        
        # Obtain a list of source in this image
        # execute this in another python session since importing casac in 
        # casanova messes up pybdsm
        syscommand = 'lib/find_sources.py '+path_to_currentimage+'-masked.image.tt0'+\
                     ' -c '+mset.dir_peel+extend_dir+\
                     '/list_of_sources.fits --atrous_do'
        #os.system(syscommand)
        
        peelcatalog = fits.open(mset.dir_peel+extend_dir+\
                                '/list_of_sources.fits')
        peel_cat_data = peelcatalog[1].data
        peel_cat_dec = peel_cat_data['DEC']
        peel_cat_ra = peel_cat_data['RA']
        
        # Use these parameters to select sources that are suitable for peeling
        peel_cat_tflux = peel_cat_data['Total_flux']
        peel_cat_pflux = peel_cat_data['Peak_flux']
        peel_cat_isl_rms = peel_cat_data['Isl_rms']

        # Use total_flux, peak_flux, source_code, rms to determine if the
        # source is suitable.
        peelsources = []
        radeclist = []
        
        peel_radius_in_arcmin = 7
        
        for ra, dec, tflux, pflux, rms in zip(peel_cat_ra, peel_cat_dec, 
                                              peel_cat_tflux, peel_cat_pflux,
                                              peel_cat_isl_rms):
            pnr = pflux/rms
            snr = tflux/rms
            
            qualityfactor = 0.5
            qualitymeasure = pnr**(qualityfactor)*snr**(1-qualityfactor)

            if qualitymeasure > 50:
                # Now, also check if the peel patches don't overlap
                diff_large_enough = True
                for othersource in radeclist:
                    other_ra = othersource[0]
                    other_dec = othersource[1]
                    difference = ((ra-other_ra)**2.0 + (dec-other_dec)**2.0)**.5
                    
                    # I demand more than 2/3 separation in diameter
                    if difference < (2./3.) * (2*peel_radius_in_arcmin) * (1./60):
                        diff_large_enough = False
                
                if diff_large_enough:
                    peelsources.append([ra, dec, snr, pnr, rms, pflux, tflux])
                
                radeclist.append([ra, dec])

        # Sort the peel sources in order of decreasing peak flux to noise ratio
        peelsources.sort(key=lambda x: x[2], reverse=True)        
        # Create an image showing the earlier image with the potential sources to peel
        #fitsfig = aplpy.FITSFigure(path_to_currentimage+'-masked.image.tt0.fits')
        #fitsfig.show_grayscale()
        #fitsfig.show_circles([i[0] for i in peelsources],
                             #[i[1] for i in peelsources], peel_radius_in_arcmin*(1./60), lw=2.5) # (1./60) degrees is an arcmin
        #fitsfig.save(path_to_currentimage+'-masked.image.tt0.fits.potential_sources.png')

        # Now that we have a list of sources to peel, calculate for each source
        # the solution interval. Determine the solution interval needed to obtain
        # a certain SNR but also reject the source if the solution interval is
        # too large.
        min_SNR_per_solution_interval = 15
        fudge_factor = 2.0**0.5
        peelsourcesupdated = []
        for source in peelsources:
            snr = source[2]
            pnr = source[3]
            rms = source[4]
            pflux = source[5]
            tflux = source[6]
            
            noise_per_interval = rms * (number_of_integrations)**0.5
            snr_per_interval = tflux/noise_per_interval
            time_steps_needed = ((min_SNR_per_solution_interval/fudge_factor)/snr_per_interval)**2.0 # SNR grows with sqrt(time)
            time_needed_in_sec = integration_time*time_steps_needed
            
            if time_needed_in_sec < 120:
                solint = '%.1fs' % time_needed_in_sec
                source.append(solint)
                peelsourcesupdated.append(source)
        
        # Create an image showing the earlier image with the final sources to peel
        #fitsfig = aplpy.FITSFigure(path_to_currentimage+'-masked.image.tt0.fits')
        #fitsfig.show_grayscale()
        #fitsfig.show_circles([i[0] for i in peelsourcesupdated],
                             #[i[1] for i in peelsourcesupdated], peel_radius_in_arcmin*0.01667, lw=2.5) # 0.01667 degrees is an arcmin
        #fitsfig.save(path_to_currentimage+'-masked.image.tt0.fits.final_sources.png')
        
        # Create a residual image needed for peeling (i.e.: add peel source to
        # residual image, self calibrate, remove peel source and use updated
        # residual image for further peeling.
        prepeelingmodels = [path_to_currentimage+'-masked.model.tt0',
                            path_to_currentimage+'-masked.model.tt1']
        peel.subtract(target_mset.file_path, prepeelingmodels, wprojplanes=512)
        residualMSpath = mset.dir_peel+extend_dir+'/pre/'+target_mset.ms_name+'_residual.ms'
        split(vis=target_mset.file_path, outputvis=residualMSpath)
        residualMS = TableObjects.MSObj(residualMSpath)
        
        # Peel all sources
        for i, source in enumerate(peelsourcesupdated):
            ra = source[0]
            dec = source[1]
            rms = source[4]
            solint = source[7]
            solint_data = {'solint':solint, 'numint':number_of_integrations, 
                           'minsnr':min_SNR_per_solution_interval,
                           'ff':fudge_factor, 'inttime':integration_time}
            
            if not os.path.isdir(mset.dir_peel+extend_dir+'/region'+str(i+1)):
                os.makedirs(mset.dir_peel+extend_dir+'/region'+str(i+1))

            # For each source create a directory named region# with a region 
            # file in it describing the source                
            with open(mset.dir_peel+extend_dir+'/region'+str(i+1)+'/region.crtf', 'w') as regionfile:
                string_format = utils.deg2HMS(ra=ra, dec=dec)
                ra_string = string_format[0]
                dec_string = string_format[1]
                regionfile.write('#CRTFv0')
                regionfile.write('\n')
                regionfile.write('circle[['+ra_string+', '+dec_string+'], '+str(peel_radius_in_arcmin)+'arcmin]')
        
            residualMS = peel.peel(residualMS, target_mset, prepeelingmodels, mset.dir_peel+extend_dir+'/region'+str(i+1), solint_data, rms)

        

def createimage(mset):
    logging.info("### CREATE FINAL IMAGE")

    # Create image for each target field individually
    for target_mset in mset.targetmsets:
        # This measurement set only has one target
        tgt_fieldname = target_mset.targetsources[0].field_name

        parms = {'vis':target_mset.file_path, 'imagename':target_mset.file_path[:-3]+'_final', 'gridmode':'widefield', 'wprojplanes':512,
                 'mode':'mfs', 'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',
                 'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True,
                 'uvtaper':True, 'threshold':str(expnoise)+' Jy'}
        utils.cleanmaskclean(parms, makemask=False)
    
        # Create fits file
        exportfits(imagename=target_mset.file_path[:-3]+'_final.image.tt0',
                   fitsimage=target_mset.file_path[:-3]+'_final.fits',
                   history=False, overwrite=True)
