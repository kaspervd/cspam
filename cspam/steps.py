import datetime
import logging
import os
import sys
from lib import AntennaObjects
from lib import utils

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
    plotms(vis=mset.file_path, xaxis='time', yaxis='elevation',
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
    #
    # This comes from Fransesco's GMRT pipeline 
    # (two times a number of channels?)
    #
    #elif len(mset.nchan) == 2 and mset.nchan[0] == 128 and mset.nchan[1] == 128:
    #    spw='0:0,1:0' # TODO: is also the chan 0 of the second spw to be flagged?
    #    if freq > 200e6 and freq < 300e6: spw='0:0~65,1:0' # 235 MHz +20 low border
    #
    
    # CHECK IF THIS IS CORRECT
    elif mset.nchan == 64:
        spw='0:0'
        if mset.freq > 200e6 and mset.freq < 300e6:
            spw='0:0~32' # 235 MHz +20 border
    # CHECK IF ABOVE IS CORRECT
    
    else:
        logging.error('Cannot understand obs type.')
        sys.exit(1)

    flagdata(vis=mset.file_path, mode='manualflag', spw=spw, flagbackup=False)
    
    # TO DO:
    # DEAL WITH BAD ANTENNAS IN CONFIG FILE
    
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
    
    scans = ','.join(mset.cal_scan_ids)
    setjy(vis=mset.file_path, scan=scans, 
          standard='Perley-Butler 2010', usescratch=True,
          scalebychan=True)

def bandpass_calibration(mset):
    """
    This method compensates for the change of gain with frequency using a
    calibrator source.
    """
    done = []
    for cal_field_id in mset.cal_field_ids:
        # Use each calibrator to do bandpass calibration.
        cal_fieldname = mset.get_field_name_from_field_id(cal_field_id)
        
        if not os.path.isdir(mset.dir_cal+'/flux_cal_'+cal_fieldname):
            os.makedirs(mset.dir_cal+'/flux_cal_'+cal_fieldname)
        if not os.path.isdir(mset.dir_plot+'/flux_cal_'+cal_fieldname):
            os.makedirs(mset.dir_plot+'/flux_cal_'+cal_fieldname)
            
        # In order for plotcal to work a symbolic link is needed.
        # Plotcal assumes the measurement set is in the same directory
        # as the cal table.
        syscommand = 'ln -s '+mset.file_path+' '+mset.dir_cal+'/flux_cal_'+ \
                     cal_fieldname+'/'+mset.ms_name+'.ms'
        os.system(syscommand)

        # Select a spectral window: 0:10~20 mean window 0, channels 10 to 20
        if mset.nchan == 512: initspw = '0:240~260'
        elif mset.nchan == 256: initspw = '0:120~130'
        elif mset.nchan == 128: initspw = '0:70~80'
        elif mset.nchan == 64: initspw = '0:40~50' # CHECK IF THIS IS CORRECT
        #
        # This comes from Fransesco's GMRT pipeline 
        # (two times a number of channels?)
        #
        #elif len(mset.nchan) == 2 and mset.nchan[0] == 128 and mset.nchan[1] == 128: initspw = '0:70~80, 1:70~80'

        for step in ['cycle1','final']:
            # Do two bandpass calibration steps per calibrator.
            logging.info("Start bandpass step: "+step)

            gaintables=[]
            interp=[]
    
            refAntObj = AntennaObjects.RefAntHeuristics(vis=mset.file_path, 
                                       field=cal_fieldname, geometry=True, 
                                       flagging=True)
            refAnt = refAntObj.calculate()[0]
            
            # Select minimal signal to noise ratio
            if mset.freq < 500e6:
                minsnr=3.0
            else:
                minsnr=5.0
            
            caltablepath = mset.dir_cal+'/flux_cal_'+cal_fieldname+'/'
            scans = ','.join(mset.cal_scan_ids)
            
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
                    field=cal_fieldname, selectdata=True, uvrange='>50m',
                    scan=scans, spw=initspw, solint='int', combine='', 
                    refant=refAnt, minblperant=mset.minBL_for_cal, minsnr=0, 
                    calmode='p')
            # Smoothing solutions
            smoothcal(vis=mset.file_path, tablein=caltablepath+step+'-boot.Gp',
                      caltable=caltablepath+step+'-boot.Gp-smooth')
            
            # Initial bandpass correction.
            # Calibration type 'B' differs from 'G' only in that it is 
            # determined for each channel in each spectral window.
            logging.info("Bandpass calibration 1")
            bandpass(vis=mset.file_path, caltable=caltablepath+step+'-boot.B', 
                     field=cal_fieldname, selectdata=True, uvrange='>100m', 
                     scan=scans, solint='inf', combine='scan,field', 
                     refant=refAnt, minblperant=mset.minBL_for_cal, 
                     minsnr=minsnr, solnorm=True, bandtype='B', 
                     gaintable=[caltablepath+step+'-boot.Gp-smooth'], 
                     interp=['linear'])

            # Find leftover time-dependent delays, Type 'K' solves for simple
            # antenna-based delays via Fourier transforms of the spectra on
            # baselines to the reference antenna. 
            logging.info("BP: Delay calibration")
            gaincal(vis=mset.file_path, caltable=caltablepath+step+'.K', 
                    field=cal_fieldname, selectdata=True, uvrange='>100m',
                    scan=scans, solint='int',combine='', refant=refAnt, 
                    interp=interp+['nearest,nearestflag'], gaintype='K',
                    minblperant=mset.minBL_for_cal, minsnr=minsnr,
                    gaintable=gaintables+[caltablepath+step+'-boot.B'])
            # flag outliers
            utils.FlagCal(caltablepath+step+'.K', sigma = 5, cycles = 3)
            # plot
            utils.plotGainCal(caltablepath+step+'.K', mset.dir_plot+ \
                              '/flux_cal_'+cal_fieldname, delay=True)
            gaintables.append(caltablepath+step+'.K')
            interp.append('linear')
            
            # Find time-dependant gains, now taking both amplitude and phase
            # into account.
            logging.info("BP: Gain calibration")
            gaincal(vis=mset.file_path, caltable=caltablepath+step+'.Gap', 
                    field=cal_fieldname, selectdata=True, uvrange='>100m', 
                    scan=scans, solint='int',combine='', refant=refAnt, 
                    interp=interp+['nearest,nearestflag'],
                    minblperant=mset.minBL_for_cal, minsnr=minsnr,  
                    gaintype='G', calmode='ap',
                    gaintable=gaintables+[caltablepath+step+'-boot.B'])
            # flag outliers
            utils.FlagCal(caltablepath+step+'.Gap', sigma = 3, cycles = 3)
            # plot
            utils.plotGainCal(caltablepath+step+'.Gap', mset.dir_plot+ \
                              '/flux_cal_'+cal_fieldname, amp=True, phase=True)
            gaintables.append(caltablepath+step+'.Gap')
            interp.append('linear')

            # Find cross-K, type 'KCROSS' solves for global cross-hand delays,
            # so this is again to find delays.
            logging.info("BP: Kcross calibration")
            gaincal(vis=mset.file_path, caltable=caltablepath+step+'.Kcross',
                    field=cal_fieldname, selectdata=True, uvrange='>100m',
                    scan=scans, solint='inf',combine='scan,field', 
                    refant=refAnt, interp=interp+['nearest,nearestflag'],
                    minblperant=mset.minBL_for_cal, minsnr=minsnr,
                    gaintype='KCROSS',
                    gaintable=gaintables+[caltablepath+step+'-boot.B'])
            # plot
            plotcal(caltable = caltablepath+step+'.Kcross', xaxis = 'antenna', 
                    yaxis = 'delay', showgui=False,
                    figfile= mset.dir_plot+'/flux_cal_'+cal_fieldname+'/'+ \
                    step+'.Kcross.png' )
            gaintables.append(caltablepath+step+'.Kcross')
            interp.append('nearest')

            # Recalculate bandpass taking delays into account.
            logging.info("Bandpass calibration 2")
            bandpass(vis=mset.file_path, caltable=caltablepath+step+'.B',
                     field=cal_fieldname, selectdata=True, uvrange='>100m', 
                     scan=scans, solint='inf', combine='scan,field',
                     refant=refAnt, interp=interp, 
                     minblperant=mset.minBL_for_cal, minsnr=minsnr, 
                     solnorm=True, bandtype='B', gaintable=gaintables)
            # plot
            utils.plotBPCal(caltablepath+step+'.B', mset.dir_plot+ \
                            '/flux_cal_'+cal_fieldname, amp=True, phase=True)
            gaintables.append(caltablepath+step+'.B')
            interp.append('nearest,nearestflag')

            # Apply the gaintables from this bandpass calibration cycle.
            logging.info("Apply bandpass")
            applycal(vis=mset.file_path, selectdata=True, field=cal_fieldname,
                     scan=scans, gaintable=gaintables, calwt=False,
                     flagbackup=False, interp=interp)
            
            if step != 'final':
                # clip on residuals
                utils.clipresidual(mset.file_path, f=cal_fieldname, s=scans)

        # end of 2 bandpass cycles
        done.append(cal_fieldname)
    # end of loop through all calibrator sources

    # remove K, amp from gaintables, we keep B and Kcross which are global 
    # and T-indep
    gaintables=[caltablepath+'final.B', caltablepath+'final.Kcross']
    interp=['nearest,nearestflag','nearest,nearestflag']
    
    utils.statsFlag(mset.file_path, note='Before apply bandpass')

    # apply final bandpass to target
    for tgt_field_id in mset.tgt_field_ids:
        tgt_fieldname = mset.get_field_name_from_field_id(tgt_field_id)
        scans = ','.join(mset.tgt_scan_ids)
        
        applycal(vis=mset.file_path, selectdata=True, field=tgt_fieldname, 
                 scan=scans, gaintable=gaintables, calwt=False, 
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
    This method ...
    """
    
    logging.info("### CALIB")
    
    # Define these variables early since they are needed in all loops
    gaintables=[]
    interp=[]
    
    for cal_field_id in mset.cal_field_ids:
        cal_fieldname = mset.get_field_name_from_field_id(cal_field_id)
        caltablepathg = mset.dir_cal+'/gain_cal_'+cal_fieldname+'/'
        caltablepathf = mset.dir_cal+'/flux_cal_'+cal_fieldname+'/'

        if not os.path.isdir(mset.dir_cal+'/gain_cal_'+cal_fieldname):
            os.makedirs(mset.dir_cal+'/gain_cal_'+cal_fieldname)
        if not os.path.isdir(mset.dir_plot+'/gain_cal_'+cal_fieldname):
            os.makedirs(mset.dir_plot+'/gain_cal_'+cal_fieldname)
        if not os.path.isdir(mset.dir_img+'/gain_cal_'+cal_fieldname):
            os.makedirs(mset.dir_img+'/gain_cal_'+cal_fieldname)
            
        # In order for plotcal to work a symbolic link is needed.
        # Plotcal assumes the measurement set is in the same directory
        # as the cal table.
        syscommand = 'ln -s '+mset.file_path+' '+mset.dir_cal+'/gain_cal_'+ \
                     cal_fieldname+'/'+mset.ms_name+'.ms'
        os.system(syscommand)

        n_cycles = 3
        for cycle in xrange(n_cycles):
    
            logging.info("Start CALIB cycle: "+str(cycle))
    
            refAntObj = AntennaObjects.RefAntHeuristics(vis=mset.file_path, 
                                       field=cal_fieldname, geometry=True, 
                                       flagging=True)
            refAnt = refAntObj.calculate()[0]
            
            # Reset these for this cycle
            gaintables=[caltablepathf+'final.B', caltablepathf+'final.Kcross']
            interp=['nearest,nearestflag','nearest,nearestflag']
            
            scans = ','.join(mset.cal_scan_ids)
            
            # Gain cal phase
            if mset.freq < 500e6:
                minsnr=2.0
            else:
                minsnr=4.0
            gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                    str(cycle)+'-boot.Gp', field=cal_fieldname, selectdata=True,
            	    uvrange='>100m', scan=scans, solint='int', refant=refAnt, 
            	    interp=interp, minblperant=mset.minBL_for_cal, 
            	    minsnr=minsnr, calmode='p', gaintable=gaintables)

            # Find leftover time-dependent delays, Type 'K' solves for simple
            # antenna-based delays via Fourier transforms of the spectra on
            # baselines to the reference antenna. 
            gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                    str(cycle)+'.K', field=cal_fieldname, selectdata=True,
                    uvrange='>100m', scan=scans, solint='int', refant=refAnt,
                    minblperant=mset.minBL_for_cal, minsnr=minsnr, gaintype='K',
                    interp=interp+['linear'],gaintable=gaintables+ \
                    [caltablepathg+'gain'+str(cycle)+'-boot.Gp'])
            utils.FlagCal(caltablepathg+'gain'+str(cycle)+'.K', sigma = 5, 
                          cycles = 3)
            utils.plotGainCal(caltablepathg+'gain'+str(cycle)+'.K', 
                              mset.dir_plot+'/gain_cal_'+cal_fieldname,
                              delay=True)

            # Redo gain cal phase (now with K)
            gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                    str(cycle)+'.Gp', field=cal_fieldname, selectdata=True,
            	    uvrange='>100m', scan=scans, solint='int', refant=refAnt,
                    minblperant=mset.minBL_for_cal, minsnr=minsnr, calmode='p', 
                    gaintable=gaintables+[caltablepathg+'gain'+str(cycle)+'.K'],
                    interp=interp+['linear'])
            # Smoothing solutions
            smoothcal(vis=mset.file_path, tablein=caltablepathg+'gain'+ \
                      str(cycle)+'.Gp', caltable=caltablepathg+'gain'+ \
                      str(cycle)+'.Gp-smooth')
            utils.plotGainCal(caltablepathg+'gain'+str(cycle)+'.Gp-smooth',
                              mset.dir_plot+'/gain_cal_'+cal_fieldname,
                              phase=True)
            gaintables.append(caltablepathg+'gain'+str(cycle)+'.Gp-smooth')
            interp.append('linear')
    
            # Gain cal amp
            if mset.freq < 500e6:
                minsnr=3.0
            else:
                minsnr=5.0
            gaincal(vis=mset.file_path, caltable=caltablepathg+'gain'+ \
                    str(cycle)+'.Ga', field=cal_fieldname, selectdata=True,
                    uvrange='>100m', scan=scans, solint='60s', minsnr=minsnr,
                    refant=refAnt, minblperant=mset.minBL_for_cal, calmode='a', 
                    gaintable=gaintables)
            utils.FlagCal(caltablepathg+'gain'+str(cycle)+'.Ga', sigma = 3, 
                          cycles = 3)
    
            """
            # if gain and flux cal are the same the fluxscale cannot work
            # do it only in the last cycle, so the next clip can work, otherwise the uvsub subtract
            # a wrong model for (amp==1) for the gain_cal if it had been rescaled
            if s.g != s.f and cycle == n_cycles-1:
                # fluxscale
                logging.debug("Rescale gaincal sol with fluxcal sol.")
                #default('fluxscale')
                fluxscale(vis=active_ms, caltable=caltablepath+'gain'+str(cycle)+'.Ga',\
                	fluxtable=caltablepath+'gain'+str(cycle)+'.Ga_fluxscale', reference=s.f, transfer=s.g)

                plotGainCal(caltablepath+'gain'+str(cycle)+'.Ga_fluxscale', amp=True)
                gaintables.append(caltablepath+'gain'+str(cycle)+'.Ga_fluxscale')
                interp.append('linear')
            else:
                plotGainCal(caltablepath+'gain'+str(cycle)+'.Ga', amp=True)
                gaintables.append(caltablepath+'gain'+str(cycle)+'.Ga')
                interp.append('linear')
            """
            
            # In this case gain and flux calibrators are one and the same
            utils.plotGainCal(caltablepathg+'gain'+str(cycle)+'.Ga',
                              mset.dir_plot+'/gain_cal_'+cal_fieldname,
                              amp=True)
            gaintables.append(caltablepathg+'gain'+str(cycle)+'.Ga')
            interp.append('linear')
            
            # BLcal TODO: do BLcal on the fluxcal?
            #if s.f in s.fmodel:
            #    print "WARNING: flux_cal has a model and its being used for BLCAL, model must be superprecise!" 
            #blcal(vis=active_ms, caltable='cal/'+s.name+'/gain'+str(cycle)+'.BLap',  field=s.f,\
            #    scan=s.fscan, combine='', solint='inf', calmode='ap', gaintable=gaintables, solnorm=True)
            #FlagBLcal('cal/'+s.name+'/gain'+str(cycle)+'.BLap', sigma = 3)
            #plotGainCal('cal/'+s.name+'/gain'+str(cycle)+'.BLap', amp=True, phase=True, BL=True)
            #gaintables.append('cal/'+s.name+'/gain'+str(cycle)+'.BLap')
            #interp.append('nearest')

            applycal(vis=mset.file_path, field=cal_fieldname, scan=scans,
                     gaintable=gaintables, interp=interp, calwt=False, 
                     flagbackup=False)
            
            # clip of residuals not on the last cycle (useless and prevent imaging of calibrator)
            if cycle != n_cycles-1:
                utils.clipresidual(mset.file_path, f=cal_fieldname, s=scans)

        # make a test img of the gain cal to check that everything is fine
        parms = {'vis':mset.file_path, 'field':cal_fieldname,
                 'imagename':mset.dir_img+'/gain_cal_'+cal_fieldname+'/'+ \
                 cal_fieldname+'_gcal',
                 'gridmode':'widefield', 'wprojplanes':128,
              	 'mode':'mfs', 'nterms':2, 'niter':1000, 'gain':0.1,
              	 'psfmode':'clark', 'imagermode':'csclean', 'imsize':512,
              	 'cell':sou_res, 'weighting':'briggs', 'robust':0,
              	 'usescratch':False}
        utils.cleanmaskclean(parms, makemask=False)

    # use a different cycle to compensate for messing up with uvsub during the calibration of other sources
    # in this way the CRRECTED_DATA are OK for all fields
    for field_id in mset.cal_field_ids+mset.tgt_field_ids:
        fieldname = mset.get_field_name_from_field_id(field_id)
        # apply B, Gp, Ga
        if field_id in mset.cal_field_ids:
            # Calibrator
            scans = ','.join(mset.cal_scan_ids)    
        elif field_id in mset.tgt_field_ids:
			# target
			scans = ','.join(mset.tgt_scan_ids)
        
        applycal(vis=mset.file_path, field=fieldname, scan=scans,
                 gaintable=gaintables, interp=interp, calwt=False,
                 flagbackup=False)

def selfcal(mset):
    logging.info("### SELFCAL")

    if mset.freq > 1000e6: width = 16
    if mset.freq > 550e6 and mset.freq < 650e6: width = 16
    if mset.freq > 300e6 and mset.freq < 350e6: width = 8
    if mset.freq > 200e6 and mset.freq < 300e6: width = 8
    # renormalize if chans were not 512, force int to prevent bug in split() 
    # if width is a numpy.int64
    width = int(width / (512/mset.nchan))
    logging.info("Average with width="+str(width))
    
    for tgt_field_id in mset.tgt_field_ids:
        tgt_fieldname = mset.get_field_name_from_field_id(tgt_field_id)
    
        if not os.path.isdir(mset.dir_cal+'/self_cal_'+tgt_fieldname):
            os.makedirs(mset.dir_cal+'/self_cal_'+tgt_fieldname)
        if not os.path.isdir(mset.dir_plot+'/self_cal_'+tgt_fieldname):
            os.makedirs(mset.dir_plot+'/self_cal_'+tgt_fieldname)
        if not os.path.isdir(mset.dir_img+'/self_cal_'+tgt_fieldname):
            os.makedirs(mset.dir_img+'/self_cal_'+tgt_fieldname)
    
        # From now on, only use the target field -> create new MS with target
        target_file_path = mset.file_path[:-3]+'_target.ms'
        split(vis=mset.file_path, outputvis=target_file_path,
              field=tgt_fieldname, width=width, datacolumn='corrected',
              keepflags=False)

        # In order for plotcal to work a symbolic link is needed.
        # Plotcal assumes the measurement set is in the same directory
        # as the cal table.
        syscommand = 'ln -s '+mset.file_path[:-3]+'_target.ms '+mset.dir_cal+\
                     '/self_cal_'+ tgt_fieldname+'/'+mset.ms_name+'_target.ms'
        os.system(syscommand)
    
        for cycle in xrange(6):
     
            logging.info("Start SELFCAL cycle: "+str(cycle))
            
            # save flag for recovering
            flagmanager(vis=target_file_path, mode='save', 
                        versionname='selfcal-c'+str(cycle))

            ts = str(expnoise*10*(5-cycle))+' Jy' # expected noise this cycle
            
            # Store different iterations in different directories
            if not os.path.isdir(mset.dir_img+'/self_cal_'+tgt_fieldname+'/'\
                                 +tgt_fieldname+str(cycle)):
                os.makedirs(mset.dir_img+'/self_cal_'+tgt_fieldname+'/'\
                            +tgt_fieldname+str(cycle))
            parms = {'vis':target_file_path, 'imagename':mset.dir_img+ \
                     '/self_cal_'+tgt_fieldname+'/'+tgt_fieldname+str(cycle)+\
                     '/im',
                     'gridmode':'widefield', 'wprojplanes':512, 'mode':'mfs', 
                     'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 
                     'imagermode':'csclean', 'imsize':sou_size, 'cell':sou_res, 
                     'weighting':'briggs', 'robust':rob, 'usescratch':True, 
                     'mask':'', 'threshold':ts, 'multiscale':[]}
            
            utils.cleanmaskclean(parms)

            # Get img rms and if it higher apply old gaintables/flags and quit
            # "<1" is to invert the mask, this syntax is very weird, see
            # https://casa.nrao.edu/aips2_docs/notes/223/index.shtml
            mask = '"'+mset.dir_img+'/self_cal_'+tgt_fieldname+'/'\
                   +tgt_fieldname+str(cycle)+'/im.newmask'+'"'+' < 1'
            rms = imstat(imagename=mset.dir_img+'/self_cal_'+tgt_fieldname+'/'\
                         +tgt_fieldname+str(cycle)+'/im-masked.image.tt0',
                         mask=mask)['rms'][0]
                         
            if cycle != 0 and old_rms * 1.1 < rms:
                logging.warning('Image rms noise ('+str(rms)+' Jy/b) is higher \
                than previous cycle ('+str(old_rms)+' Jy/b). Apply old cal \
                tables and quitting selfcal.')

                # Not sure how I will setup the peeling part
                # rename last image so peeling doesn't use it
                #os.system('cd img/'+s.name+' && rename s/self'+str(cycle)+'/badimage/ *')

                # get previous flags
                flagmanager(vis=target_file_path, mode='restore',
                            versionname='selfcal-c'+str(cycle-1))
                
                # for the first cycle just remove all calibration i.e. no selfcal
                # for the others get the previous (i.e. cycle-2) cycle tables
                if cycle == 1:
                    clearcal(vis=target_file_path)
                elif cycle < 4:
                    utils.plotGainCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+\
                                      '/'+tgt_fieldname+str(cycle-2)+'.Gp', 
                                      mset.dir_plot+'/self_cal_'+tgt_fieldname+\
                                      '/Gp_cycle'+str(cycle-2), phase=True)
                    applycal(vis=target_file_path, gaintable=gaintable, 
                             interp=['linear','linear'], calwt=False,
                             flagbackup=False)           
                elif cycle >= 4: 
                    utils.plotGainCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+\
                                      '/'+tgt_fieldname+str(cycle-2)+'.Gp', 
                                      mset.dir_plot+'/self_cal_'+tgt_fieldname+\
                                      '/Gp_cycle'+str(cycle-2), phase=True)
                    utils.plotGainCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+\
                                      '/'+tgt_fieldname+str(cycle-2)+'.Ga', 
                                      mset.dir_plot+'/self_cal_'+tgt_fieldname+\
                                      '/Ga_cycle'+str(cycle-2), amp=True)
                    applycal(vis=target_file_path, gaintable=gaintable, 
                             interp=['linear','linear'], calwt=False, 
                             flagbackup=False)           

                break

            elif cycle != 0:
                logging.info('Rms noise change: '+str(old_rms)+' Jy/b -> '\
                             +str(rms)+' Jy/b.')

            if cycle == 5: break # don't do one more useless calibration

            old_rms = rms

            # ft() model back - if clean doesn't converge clean() fails to 
            # put the model, better do it by hand
            ftw(vis=target_file_path, 
                model=[mset.dir_img+'/self_cal_'+tgt_fieldname+'/'\
                       +tgt_fieldname+str(cycle)+'/im-masked.image.tt0',
                       mset.dir_img+'/self_cal_'+tgt_fieldname+'/'\
                       +tgt_fieldname+str(cycle)+'/im-masked.image.tt1'],
                nterms=2, wprojplanes=512, usescratch=True)
            
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
                solint='5s'
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
                    caltable=mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                             +tgt_fieldname+str(cycle)+'.Gp',
                    solint=solint, minsnr=minsnr, selectdata=True, 
                    uvrange='>50m', refant=refAnt, 
                    minblperant=mset.minBL_for_cal,
                    gaintable=[], calmode='p')
           
            # Gaincal - amp
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
                        caltable=mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                                 +tgt_fieldname+str(cycle)+'.Ga',
                	    selectdata=True, uvrange='>50m', solint=solint,
                	    minsnr=minsnr, refant=refAnt, gaintable=[],
                        minblperant=mset.minBL_for_cal, calmode='a')
                utils.FlagCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                              +tgt_fieldname+str(cycle)+'.Ga',
                              sigma = 3, cycles = 3)
     
            # plot gains
            if cycle >= 3: 
                utils.plotGainCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                                  +tgt_fieldname+str(cycle)+'.Gp', 
                                  mset.dir_plot+'/self_cal_'+tgt_fieldname+\
                                  '/Gp_cycle'+str(cycle),
                                  phase=True)
                utils.plotGainCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                                  +tgt_fieldname+str(cycle)+'.Ga', 
                                  mset.dir_plot+'/self_cal_'+tgt_fieldname+\
                                  '/Ga_cycle'+str(cycle),
                                  amp=True)
            else:
                utils.plotGainCal(mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                                  +tgt_fieldname+str(cycle)+'.Gp', 
                                  mset.dir_plot+'/self_cal_'+tgt_fieldname+\
                                  '/Gp_cycle'+str(cycle),
                                  phase=True)

            # add to gaintable
            if cycle >= 3: 
                gaintable=[mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                           +tgt_fieldname+str(cycle)+'.Gp',
                	       mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                           +tgt_fieldname+str(cycle)+'.Ga']
            else:
                gaintable=[mset.dir_cal+'/self_cal_'+tgt_fieldname+'/'\
                           +tgt_fieldname+str(cycle)+'.Gp',]

            applycal(vis=target_file_path, field = '', gaintable=gaintable, 
                     interp=['linear','linear'], calwt=False, flagbackup=False)           
            utils.statsFlag(target_file_path, 
                            note='After apply selfcal (cycle: '+str(cycle)+')') 

        # end of selfcal loop
    
    # end of cycle on sources


def createimage(mset):
    logging.info("### CREATE FINAL IMAGE")

    target_file_path = mset.file_path[:-3]+'_target.ms'

    parms = {'vis':target_file_path, 'imagename':target_file_path[:-3]+'_final', 'gridmode':'widefield', 'wprojplanes':512,
           	 'mode':'mfs', 'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',
             'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True,
             'uvtaper':True, 'threshold':str(expnoise)+' Jy'}
    utils.cleanmaskclean(parms, makemask=False)
    
    exportfits(imagename=target_file_path[:-3]+'_final.image.tt0',fitsimage=target_file_path[:-3]+'_final.fits',history=False)

    #for s in sources:
        #check_rm('img/'+s.name+'/lowres*')

        #parms = {'vis':s.ms, 'imagename':'img/'+s.name+'/lowres', 'gridmode':'widefield', 'wprojplanes':512,\
           	#'mode':'mfs', 'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 'imagermode':'csclean',\
            #'imsize':sou_size, 'cell':sou_res, 'weighting':'briggs', 'robust':rob, 'usescratch':True, 'mask':s.mask, \
            #'uvtaper':True, 'outertaper':[taper], 'threshold':str(s.expnoise)+' Jy', 'multiscale':s.multiscale}
        #cleanmaskclean(parms, s)
        
        ## pbcorr
        #correctPB('img/'+s.name+'/lowres-masked.image.tt0', freq, phaseCentre=None)
