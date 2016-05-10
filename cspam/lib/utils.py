import casac
import json
import logging
import numpy as np
import math
import os

af = casac.casac.agentflagger()
tb = casac.casac.table()
ms = casac.casac.ms()

from casat import plotcal
plotcal = plotcal.plotcal

from casat import flagdata
flagdata = flagdata.flagdata

from casat import uvsub
uvsub = uvsub.uvsub

from casat import clean
clean = clean.clean

def print_dict(dictionary):
    """
    Print a dictionary in a more readable format
    """
    print json.dumps(dictionary, indent=4)

def angularSeparationOfDirectionsArcsec(dir1,dir2,returnComponents=False):
    """
    Accepts two direction dictionaries and returns the separation in arcsec.
    It computes great circle angle using the Vincenty formula.
    Todd Hunter
    """
    retval = angularSeparationOfDirections(dir1, dir2, returnComponents)
    if (returnComponents):
        retval = np.array(retval) * 180*3600 / np.pi
    else:
        retval *= 180*3600 / np.pi
    return(retval)

def angularSeparationOfDirections(dir1,dir2,returnComponents=False):
    """
    Accepts two direction dictionaries and returns the separation in radians.
    It computes great circle angle using the Vincenty formula.
    --Todd Hunter
    """
    rad = angularSeparationRadians(dir1['m0']['value'], dir1['m1']['value'], dir2['m0']['value'], dir2['m1']['value'],returnComponents)
    return(rad)

def angularSeparationRadians(ra0,dec0,ra1,dec1,returnComponents=False):
    """
    Computes the great circle angle between two celestial coordinates.
    using the Vincenty formula (from wikipedia) which is correct for all
    angles, as long as you use atan2() to handle a zero denominator.
    See  http://en.wikipedia.org/wiki/Great_circle_distance
    Input and output are in radians.  It also works for the az,el coordinate system.
    returnComponents=True will return: [separation, raSeparation, decSeparation, raSeparationCosDec]
    See also angularSeparation()
    -- Todd Hunter
    """
    result = angularSeparation(ra0*180/math.pi, dec0*180/math.pi, ra1*180/math.pi, dec1*180/math.pi,returnComponents)
    if (returnComponents):
        return(np.array(result)*math.pi/180.)
    else:
        return(result*math.pi/180.)

def angularSeparation(ra0,dec0,ra1,dec1, returnComponents=False):
    """
    Computes the great circle angle between two celestial coordinates.
    using the Vincenty formula (from wikipedia) which is correct for all
    angles, as long as you use atan2() to handle a zero denominator.
    See  http://en.wikipedia.org/wiki/Great_circle_distance
    ra,dec must be given in degrees, as is the output.
    It also works for the az,el coordinate system.
    Component separations are field_0 minus field_1.
    See also angularSeparationRadians()
    returnComponents: if True, then also compute angular separation in both
         coordinates and the position angle of the separation vector on the sky
    -- Todd Hunter
    """
    ra0 *= math.pi/180.
    dec0 *= math.pi/180.
    ra1 *= math.pi/180.
    dec1 *= math.pi/180.
    deltaLong = ra0-ra1
    argument1 = (((math.cos(dec1)*math.sin(deltaLong))**2) +
                 ((math.cos(dec0)*math.sin(dec1)-math.sin(dec0)*math.cos(dec1)*math.cos(deltaLong))**2))**0.5
    argument2 = math.sin(dec0)*math.sin(dec1) + math.cos(dec0)*math.cos(dec1)*math.cos(deltaLong)
    angle = math.atan2(argument1, argument2) / (math.pi/180.)
    if (angle > 360):
        angle -= 360
    if (returnComponents):
        cosdec = math.cos((dec1+dec0)*0.5)
        radegreesCosDec = np.degrees(ra0-ra1)*cosdec
        radegrees = np.degrees(ra0-ra1)
        decdegrees = np.degrees(dec0-dec1)
        if (radegrees > 360):
            radegrees -= 360
        if (decdegrees > 360):
            decdegrees -= 360
#       positionAngle = -math.atan2(decdegrees*math.pi/180., radegreesCosDec*math.pi/180.)*180/math.pi
        retval = angle,radegrees,decdegrees, radegreesCosDec
    else:
        retval = angle
    return(retval)

def FlagCal(caltable, sigma = 5, cycles = 3):
    """
    Flag sol outside n sigmas
    Better high number of cycles (3) at high sigma (5)
    """
    tb.open(caltable, nomodify=False)
    if 'CPARAM' in tb.colnames():
        pars=tb.getcol('CPARAM')
    elif 'FPARAM' in tb.colnames():
        pars=tb.getcol('FPARAM')
    else:
        logging.error("Cannot flag "+caltable+". Unknown type.")
        return
    flags=tb.getcol('FLAG')
    ants=tb.getcol('ANTENNA1')
    totflag_before = sum(flags.flatten())
    for c in xrange(cycles):
        for ant in set(ants):
            parant = pars[:,:, np.where( ants == ant ) ]
            flagant = flags[:,:, np.where( ants == ant ) ]
            good = np.logical_not(flagant)
            if sum(good.flatten()) == 0: continue # all flagged antenna, continue
            flagant[ np.abs( parant - np.mean(parant[good]) ) > sigma * np.std( parant[good] ) ] = True
            flags[:,:, np.where( ants == ant ) ] = flagant
    tb.putcol('FLAG', flags)
    totflag_after = sum(flags.flatten())
    logging.debug(caltable+": Flagged "+str(totflag_after-totflag_before)+" points out of "+str(len(flags.flatten()))+".")
    tb.close()

def getMaxAmp(caltable):
    """
    Return maximum unflagged amp for plotting purposes
    """
    tb.open(caltable)
    cpar=tb.getcol('CPARAM')
    flgs=tb.getcol('FLAG')
    tb.close()
    amps=np.abs(cpar)
    good=np.logical_not(flgs)
    maxamp=np.max(amps[good])
    return maxamp

def plotGainCal(calt, plotdirectory, amp=False, phase=False, BL=False, delay=False):
    """
    Do the standard plot of gain solutions
    """
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)
    
    tb.open( '%s/ANTENNA' % calt)
    nameAntenna = tb.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tb.close()
    nplots=int(numAntenna/3)
    if amp == True:
        plotmax = getMaxAmp(calt)
        for ii in range(nplots):
            filename=plotdirectory+'/'+'a_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            if BL: xaxis = 'antenna2'
            else: xaxis = 'time'
            if BL: plotsymbol = 'o'
            else: plotsymbol = 'o-'
            plotcal(caltable=calt,xaxis=xaxis,yaxis='amp',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,0,plotmax],plotsymbol=plotsymbol,plotcolor='red',\
                markersize=5.0,fontsize=8.0,showgui=False,figfile=filename)
    if phase == True:
        for ii in range(nplots):
            filename=plotdirectory+'/'+'p_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            if BL: xaxis = 'antenna2'
            else: xaxis = 'time'
            plotcal(caltable=calt,xaxis=xaxis,yaxis='phase',antenna=antPlot,subplot=311,\
                overplot=False,clearpanel='Auto',iteration='antenna',plotrange=[0,0,-180,180],\
                plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=8.0,showgui=False,\
                figfile=filename)
    if delay == True:
        for ii in range(nplots):
            filename=plotdirectory+'/'+'d_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            plotcal(caltable=calt,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                overplot=False,clearpanel='Auto',iteration='antenna',plotrange=[],\
                plotsymbol='o-',markersize=5.0,fontsize=8.0,showgui=False,\
                figfile=filename)

def plotBPCal(calt, plotdirectory, amp=False, phase=False):
    """
    Do the standard plot of bandpass solutions
    """
    tb.open(calt)
    dataVarCol = tb.getvarcol('CPARAM')
    flagVarCol = tb.getvarcol('FLAG')
    tb.close()
    rowlist = dataVarCol.keys()
    nrows = len(rowlist)
    maxmaxamp = 0.0
    maxmaxphase = 0.0
    for rrow in rowlist:
        dataArr = dataVarCol[rrow]
        flagArr = flagVarCol[rrow]
        amps=np.abs(dataArr)
        phases=np.arctan2(np.imag(dataArr),np.real(dataArr))
        good=np.logical_not(flagArr)
        tmparr=amps[good]
        if (len(tmparr)>0):
            maxamp=np.max(amps[good])
            if (maxamp>maxmaxamp):
                maxmaxamp=maxamp
        tmparr=np.abs(phases[good])
        if (len(tmparr)>0):
            maxphase=np.max(np.abs(phases[good]))*180./np.pi
            if (maxphase>maxmaxphase):
                maxmaxphase=maxphase
    ampplotmax=maxmaxamp
    phaseplotmax=maxmaxphase

    tb.open( '%s/ANTENNA' % calt)
    nameAntenna = tb.getcol( 'NAME' )
    numAntenna = len(nameAntenna)
    tb.close()
    nplots=int(numAntenna/3)

    if amp == True:
        for ii in range(nplots):
            filename=plotdirectory+'/'+'BP_a_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            #default('plotcal')
            plotcal(caltable=calt,xaxis='freq',yaxis='amp',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,0,ampplotmax],showflags=False,plotsymbol='o',\
                plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

    if phase == True:
        for ii in range(nplots):
            filename=plotdirectory+'/'+'BP_p_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            #default('plotcal')
            plotcal(caltable=calt,xaxis='freq',yaxis='phase',antenna=antPlot,subplot=311,\
                iteration='antenna',plotrange=[0,0,-phaseplotmax,phaseplotmax],showflags=False,\
                plotsymbol='o',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)


def clipresidual(active_ms, f='', s=''):
    """
    Create residuals in the CORRECTED_DATA (then unusable!)
    and clip at 5 times the total flux of the model
    NOTE: the ms CORRECTED_DATA will be corrupted!
    """

    #default('uvsub')
    uvsub(vis=active_ms)

    # flag statistics before flagging
    statsFlag(active_ms, note='Before BL flag')

    logging.debug("Removing baselines with high residuals:")
    import itertools
    ms.open(active_ms, nomodify=False)
    metadata = ms.metadata()
    flag = {}
    # datadesc ids are usually one per spw, but ms can also be splitted in corr
    for datadescid in metadata.datadescids():
        flag[datadescid] = {}
        logging.debug("Working on datadesc: "+str(datadescid))
        ms.msselect({'field':f, 'scan':s})
        d = ms.getdata(['corrected_amplitude','flag','antenna1','antenna2','axis_info'], ifraxis=True)
        # cycle on corr
        for corr in xrange(len(d['corrected_amplitude'])):
            flag[datadescid][corr] = {}
            logging.debug("Working on corr: "+d['axis_info']['corr_axis'][corr])
            # cycle on channels
            for chan in xrange(len(d['corrected_amplitude'][corr])):
                flag[datadescid][corr][chan] = {}
                meds = []
                # cycle on bl
                for bl in xrange(len(d['corrected_amplitude'][corr][chan])):
                    amp = d['corrected_amplitude'][corr][chan][bl][~d['flag'][corr][chan][bl]] # get unflagged data
                    if amp != []: meds.append(np.median( amp[(amp == amp)] ))
                if meds != []:
                    med = np.mean(meds)
                    rms = np.std(meds)
                for bl in xrange(len(d['corrected_amplitude'][corr][chan])):
                    flag[datadescid][corr][chan][bl] = False
                    amp = d['corrected_amplitude'][corr][chan][bl][~d['flag'][corr][chan][bl]] # get unflagged data
                    if amp != []: 
                        bl_med = np.median( amp[(amp == amp)] )
                        # if BL residuals are 3 times out of med rms, flag
                        if abs(bl_med - med) > 3*rms:
                            logging.debug("Flagging corr: "+d['axis_info']['corr_axis'][corr]+" - chan:"+str(chan)+" - BL: "+d['axis_info']['ifr_axis']['ifr_name'][bl])
                            flag[datadescid][corr][chan][bl] = True
    ms.close()

    # extend flags to all scans
    ms.open(active_ms, nomodify=False)
    metadata = ms.metadata()
    # datadesc ids are usually one per spw, but ms can also be splitted in corr
    for datadescid in metadata.datadescids():
        ms.msselect()
        w = ms.getdata(['flag'], ifraxis=True)
        for corr in xrange(len(w['flag'])):
            # TODO: extend flags on all chan for BLs which appear often
            for chan in xrange(len(w['flag'][corr])):
                for bl in xrange(len(w['flag'][corr][chan])):
					
                    #print '--'
                    #print len(d['corrected_amplitude'][corr][chan])
                    #print len(w['flag'][corr][chan])
                    #print len(flag[datadescid][corr][chan])
                    #print '--'
                    
                    # So, sometimes len(w['flag'][corr][chan]) > len(flag[datadescid][corr][chan])
                    # CHECK IF THIS IS A GOOD SOLUTION
                    bl = min(bl, len(flag[datadescid][corr][chan])-1)
					
                    if flag[datadescid][corr][chan][bl] == True:
                        w['flag'][corr][chan][bl] = True
                        
        ms.putdata({'flag':w['flag']})
    ms.close()

    statsFlag(active_ms, note='After clipping')

def statsFlag(active_ms, field='', scan='', note=''):
    t = flagdata(vis=active_ms, mode='summary', field=field, scan=scan, action='calculate')
    log = 'Flag statistics ('+note+'):'
    log += '\nAntenna, '
    for k in sorted(t['antenna']):
        log += k +': %.2f%% - ' % (100.*t['antenna'][k]['flagged']/t['antenna'][k]['total'])
    log += '\nCorrelation, '
    for k, v in t['correlation'].items():
        log += k +': %.2f%% - ' % (100.*v['flagged']/v['total'])
    log += '\nSpw, '
    for k, v in t['spw'].items():
        log += k +': %.2f%% - ' % (100.*v['flagged']/v['total'])
    log += '\nTotal: %.2f%%' % (100.*t['flagged']/t['total'])
    correctedlog = log.replace(' - \n','\n')
    logging.debug(correctedlog)
    return correctedlog

def cleanmaskclean(parms, makemask=True):
    """
    Clean then make a mask and clean again
    parms: dict of parameters for the clean task
    makemask: if false quit after first clean
    """
    clean(**parms)
    if not makemask: return

    # make mask and re-do image
    if parms['nterms'] > 1:
        img = parms['imagename']+'.image.tt0' # tt = Taylor terms
    else: img = parms['imagename']+'.image'

    """
    if s.extended:
        os.system(pipdir+'/setpp.sh make_mask.py '+img+' -m'+parms['imagename']+\
                  '.newmask --threshpix=6 --threshisl=3 --atrous_do')
    else:
        os.system(pipdir+'/setpp.sh make_mask.py '+img+' -m'+parms['imagename']+\
                  '.newmask --threshpix=6 --threshisl=3')
    """
    # Let's make the mask with a trous wavelet decomposition in any case
    # see: ftp://ftp.hs.uni-hamburg.de/pub/outgoing/rafferty/PyBDSM/PyBDSM_1.8.pdf
    # section 3.2.2
    
    # fetch to directory, make_mask.py is in the same directory as utils.py
    libdirectory = os.path.dirname(os.path.realpath(__file__))
    # execute this in another python session since importing casac in casanova
    # messes up pybdsm
    os.system(libdirectory+'/make_mask.py '+img+' -m'+parms['imagename']+'.newmask --threshpix=6 --threshisl=3 --atrous_do')
    
    """
    if s.mask_faint != '':
        parms['mask']=[parms['imagename']+'.newmask',s.mask_faint]
    else:
        parms['mask']=parms['imagename']+'.newmask'
    """
    # Let's assume no faint mask
    parms['mask']=parms['imagename']+'.newmask'
    parms['imagename']=parms['imagename']+'-masked'
    parms['niter']=parms['niter']/3 # reduce number if clean iterations in masked mode
    
    clean(**parms)
