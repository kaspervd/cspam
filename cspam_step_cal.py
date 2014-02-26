#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2014 - Francesco de Gasperin - Huib Intema
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import logging

def cspam_step_cal(MSs, conf):
    """
    Calibration step
    """
   
    logging.info('### STARTING DATA CALIBRATION STEP')

    for MS in MSs:
        logging.info('## WORKING ON '+MS.file_name)

        # SetJy
        for cal_scan in MS.cal_scans:
            # check if there's a specific model
            #if flux_cal in MS.models:
            #    logging.info("Using model "+models[flux_cal]+" for fux_cal "+str(flux_cal))
            #    default('ft')
            #    ft(vis=MS.file_name, field=cal_scan, complist=models[flux_cal], usescratch=True)
            #else:
            logging.info("Using default model for fux_cal "+str(flux_cal))
            default('setjy')
            setjy(vis=MS.file_name, field=cal_scan, standard='Scaife-Heald 2012', usescratch=True, scalebychan=True)

        # check for dead antennas/swapped pol


        # Bandpass calibration
        logging.info("# STARTING BANDPASS CALIBRATION")
        for cal_scan in MS.cal_scans:
            for step in ['preflag','postflag','final']:

                gaintables=[]
    
	            refAntObj = RefAntHeuristics(vis=active_ms, field=flux_cal, geometry=True, flagging=True)
	            refAnt = refAntObj.calculate()[0]
	            logging.debug("Refant: " + refAnt)
	        
	            # Delay on a narrow set of chan for BP and flagging
	            default('gaincal')
                gaincal(vis=MS.file_, caltable=MS.dir_cal'/cal'+cal_scan+'-init_'+step+'.K', gaintype = 'K',\
                    selectdata=True, uvrange=MS.uvrange, scan=cal_scan, spw=MS.central_chans,\
                    solint='int', combine='', refant=refAnt, minblperant=MS.minBL_for_cal, minsnr=4, calmode='p')

	            gaintables.append(MS.dir_cal'/cal'+cal_scan+'-init_'+step+'.K')
	            plot_cal_table(MS.dir_cal'/cal'+cal_scan+'-init_'+step+'.K', MS=MS)

	            # Phase gaincal on a narrow set of chan for BP and flagging
	            default('gaincal')
	            gaincal(vis=MS.file_, caltable=MS.dir_cal'/cal'+cal_scan+'-init_'+step+'.Gp', gaintype = 'G',\
	            	selectdata=True, uvrange=MS.uvrange, scan=cal_scan, spw=MS.central_chans,\
	                solint='int', combine='', refant=refAnt, minblperant=MS.minBL_for_cal, minsnr=4, calmode='p')
	    
	            gaintables.append(MS.dir_cal'/cal'+cal_scan+'-init_'+step+'.Gp')
	            plot_cal_table(MS.dir_cal'/cal'+cal_scan+'-init_'+step+'.Gp', MS=MS)
	    
	    
	            # Gain cal phase TODO: amp and ph cal -> CLCAL -> clipping (narrower 3 times)
	
	            # init bandpass correction
	            if freq < 500e6:
	                minsnr=2.0
	            else:
	                minsnr=5.0
	            default('bandpass')
	            bandpass(vis=active_ms, caltable='cal/'+str(i)+step+'.B', field=flux_cal, selectdata=True,\
	            	uvrange='>100m', scan=flux_cal_scan, solint='inf', combine='scan,field', refant=refAnt,\
	            	minblperant=minBL_for_cal, minsnr=minsnr, solnorm=True, bandtype='B', gaintable=gaintables)

                flag_bp()
	    
	            # Plot bandpass
	            plotBPCal('cal/'+str(i)+step+'.B', amp=True, phase=True)
	            
	            # Apply cal
	            if step == 'preflag':
	                field=flux_cal
	                scan=flux_cal_scan
	            if step == 'postflag' or step == 'final':
	                field=flux_cal+','+gain_cal+','+sou
	                scan=",".join(filter(None, [flux_cal_scan,gain_cal_scan,sou_scan]))
	            default('applycal')
	            applycal(vis=active_ms, selectdata=True, field=field, scan=scan,\
	            	gaintable=['cal/'+str(i)+step+'.B'], calwt=False, flagbackup=False)
	         
	            # Run an rflag after the first cycle
	            # to remove most obvious RFI
	            if step == 'preflag':
	                default('flagdata')
	                flagdata(vis=active_ms, mode='rflag', field=flux_cal, scan=flux_cal_scan,\
	                	ntime='scan', combinescans=False, datacolumn='corrected', winsize=3,\
	                	timedevscale=4.0, freqdevscale=4.0, action='apply', flagbackup=False)
	                default('flagdata')
	                flagdata(vis=active_ms, mode='extend', field=flux_cal, scan=flux_cal_scan, flagbackup=False)
	
	                   
	            # Flag with aoflagger at the second round,
	            # then redo the bandpass for the third and last time
	            if step == 'postflag' and i == len(obs)-1:
	                    
	                # reload only static initial flags
	                default('flagmanager')
	                flagmanager(vis=active_ms, mode='restore', versionname='AfterInitialFlagging')
	                    
	                # run aoflagger
	                syscommand = '~/opt/src/aoflagger/build/src/aoflagger -column CORRECTED_DATA -strategy ~/phd/obs/GMRT/rfi_GMRT610.rfis -indirect-read ' + active_ms
	                os.system(syscommand)
	                    
	                # flag statistics after flagging
	                stats_flag(active_ms)
	                    
	                default('flagmanager')
	                flagmanager(vis=active_ms, mode='save', versionname='AfterDeepFlagging', comment=str(datetime.datetime.now()))
	                logging.info("Saved flags in AfterDeepFlagging")
	    
	        # end of 3 bandpass cycles
	 
