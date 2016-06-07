# Description:
# ------------
# This file contains two classes that deal with CASA tables:
# MSObj (Measurement Set Object) and STObj (Solution Table Object).
# 
# These classes come originally from cspam version 0.1 by Francesco
# de Gasperin and Huib Intema. Edited by Kasper van Dam, 2016.
#
# These objects rely heavily on the CASA environment interfaced using casanova.

import os
import sys
import logging
import numpy as np

# CSPAM Modules
import SourceObjects
import utils

# CASA Toolkits
import casac
ms = casac.casac.ms()
tb = casac.casac.table()

# CASA Tasks
from casat import plotcal
plotcal = plotcal.plotcal

class MSObj:
    """
    The MSObj class provides information on CASA Measurement Sets

    Available instance attributes:
    ------------------------------
    name                    type    info

    file_path             - str   - absolute path to the file
    ms_name               - str   - name of the measurement set (without .ms)
    summary               - dict  - summary of the measurement set
    scansummary           - dict  - summary of the main table
    dir_img               - str   - absolute path to the img directory
    dir_cal               - str   - absolute path to the cal directory
    dir_plot              - str   - absolute path to the plot directory
    dir_peel              - str   - absolute path to the peel directory
    minBL_for_cal         - int   - minimum number of baselines needed
    nchan                 - int   - number of channels
    freq                  - float - reference frequency
    telescope             - str   - name of the telescope
    band                  - str   - name of the band
    fluxcalibrator        - user  - user defined class with flux cal info
    phasecalibrator       - user  - user defined class with phase cal info
    leakcalibrator        - user  - user defined class with leakeage cal info
    poscalibrator         - user  - user defined class with position cal info
    anglecalibrator       - user  - user defined class with angle cal info
    targetsources         - list  - list of target sources (user defined class)
    *spw                   - str   - spectral window
    *central_chans         - str   - central channels
    *flag                  - dict  - flag from config in casa style

    * = not used at the moment

    Available public instance methods:
    ----------------------------------
    get_scan_ids_from_field_id(field, user_scan_ids = [])
    
    get_field_name_from_field_id(field):
    
    get_field_name_from_scan_id(scan):
    
    get_field_id_from_scan_id(scan):

    get_direction_from_tgt_field_id(tgt)

    These are pretty much self explanatory. 
    """

    def __init__(self, path, conf = None):
		
        # Names and paths
        self.file_path = path
        try:
            occurence_last_slash = path.rfind('/')
            self.ms_name = path[occurence_last_slash+1:]
        except:
            self.ms_name = path
        self.ms_name = self.ms_name.replace('.ms','')
        self.ms_name = self.ms_name.replace('.MS','')

        # Measurement set summaries
        ms.open(self.file_path)
        self.summary = ms.summary()
        self.scansummary = ms.getscansummary()
        ms.close()

        # MS reduction related directories
        self.dir_img = self.file_path.replace('.ms','')
        self.dir_img = self.dir_img.replace('.MS','')
        self.dir_img = self.dir_img+'-img'
        
        self.dir_cal = self.file_path.replace('.ms','')
        self.dir_cal = self.dir_cal.replace('.MS','')
        self.dir_cal = self.dir_cal+'-cal'
        
        self.dir_plot = self.file_path.replace('.ms','')
        self.dir_plot = self.dir_plot.replace('.MS','')
        self.dir_plot = self.dir_plot+'-plot'

        self.dir_peel = self.file_path.replace('.ms','')
        self.dir_peel = self.dir_peel.replace('.MS','')
        self.dir_peel = self.dir_peel+'-peel'
                
        # Minimum number of baselines
        self.minBL_for_cal = self._get_minBL_for_cal()
        
        # Number of channels
        self.nchan = self._get_nchan()
        
        # Frequency
        self.freq = self._get_freq()
        
        # Telescope
        self.telescope = self._get_telescope()
        assert self.telescope == 'GMRT' or self.telescope == 'EVLA'
        
        # Band
        self.band = self._get_band()
        
        # Separate the scans for the calibrator and targets
        # Note that cal_scan_ids is an empty list if there are only targets
        cal_scan_ids = self._determine_cal_scan_ids()
        tgt_scan_ids = list(set(self.scansummary.keys())-set(cal_scan_ids))
        
        cal_field_id = self.get_field_id_from_scan_id(cal_scan_ids)

        # Start with the assumption that there is one calibrator for all types
        # Override these parameters later if different settings exist in the
        # configuration file.
        
        # Check if there are calibrators
        if len(cal_field_id) > 0:
            # Check if there are multiple calibrators
            if len(cal_field_id) > 1:
                # If there are more calibrators simply take the first and
                # remove the scans associated with the other calibrators
                flux_cal_field_id = cal_field_id[0]
                flux_cal_field_name = self.get_field_name_from_field_id(
                                           flux_cal_field_id)
                flux_cal_scan_ids = self.get_scan_ids_from_field_id(
                                    flux_cal_field_id, user_scan_ids=cal_scan_ids)
            else:
                # there is only one calibrator
                flux_cal_field_id = cal_field_id[0]
                flux_cal_field_name = self.get_field_name_from_field_id(
                                           flux_cal_field_id)
                flux_cal_scan_ids = cal_scan_ids

            # One calibrator assumption so equalize everything
            phase_cal_field_id = flux_cal_field_id
            phase_cal_field_name = self.get_field_name_from_field_id(
                                        phase_cal_field_id)
            phase_cal_scan_ids = flux_cal_scan_ids
            
            leak_cal_field_id = flux_cal_field_id
            leak_cal_field_name = self.get_field_name_from_field_id(
                                       leak_cal_field_id)
            leak_cal_scan_ids = flux_cal_scan_ids
        
            pos_cal_field_id = flux_cal_field_id
            pos_cal_field_name = self.get_field_name_from_field_id(pos_cal_field_id)
            pos_cal_scan_ids = flux_cal_scan_ids
        
            angle_cal_field_id = flux_cal_field_id
            angle_cal_field_name = self.get_field_name_from_field_id(
                                        angle_cal_field_id)
            angle_cal_scan_ids = flux_cal_scan_ids
        
        # Check if there is configuration file and possibly override
        # If this mset only contains targets there is no configuration file.
        if conf is not None:
			
            # See if tgt_scans/cal_scans/cal_types are set and override			
            if conf['cal_scans']:
                # Remove the calibrator scans that are unwanted
                cal_scan_ids = list(set(cal_scan_ids).intersection(
                                                      conf['cal_scans']))
            if conf['tgt_scans']:
                # Remove the target scans that are unwanted
                tgt_scan_ids = list(set(tgt_scan_ids).intersection(
                                                      conf['tgt_scans']))
			
            if conf['flux_cal_field']:
                flux_cal_field_id = conf['flux_cal_field']
                flux_cal_scan_ids = self.get_scan_ids_from_field_id(
                                         flux_cal_field_id, 
                                         user_scan_ids=cal_scan_ids)

            if conf['phase_cal_field']:
                phase_cal_field_id = conf['phase_cal_field']
                phase_cal_scan_ids = self.get_scan_ids_from_field_id(
                                         phase_cal_field_id, 
                                         user_scan_ids=cal_scan_ids)

            if conf['leakage_cal_field']:
                leak_cal_field_id = conf['leakage_cal_field']
                leak_cal_scan_ids = self.get_scan_ids_from_field_id(
                                         leak_cal_field_id, 
                                         user_scan_ids=cal_scan_ids)
			
            if conf['position_cal_field']:
                pos_cal_field_id = conf['position_cal_field']
                pos_cal_scan_ids = self.get_scan_ids_from_field_id(
                                         pos_cal_field_id, 
                                         user_scan_ids=cal_scan_ids)
			
            if conf['angle_cal_field']:
                angle_cal_field_id = conf['angle_cal_field']
                angle_cal_scan_ids = self.get_scan_ids_from_field_id(
                                         angle_cal_field_id, 
                                         user_scan_ids=cal_scan_ids)
			
			## The following variables (flag, spw, central_chans) aren't used yet
            # Flags (manual flag to be used with flagdata command)
            self.flag = self._convert_flag(conf['flag']) 
		
            # Spectral windows
            self.spw = conf['spw']
            
            # Central channel
            central_nchan = self.nchan*conf['central_chan_percentage']/100.
            self.central_chans = str(int(round(self.nchan/2.-central_nchan/2.))) \
                                 + '~' + \
                                 str(int(round(self.nchan/2.+central_nchan/2.)))

        # Check if there are calibrators
        if len(cal_field_id) > 0:
            # Calibrators
            self.fluxcalibrator = SourceObjects.CalibratorSource('flux', 
                                  flux_cal_field_name, flux_cal_field_id,
                                  flux_cal_scan_ids)
            self.phasecalibrator = SourceObjects.CalibratorSource('phase',
                                   phase_cal_field_name, phase_cal_field_id,
                                   phase_cal_scan_ids)
            self.leakcalibrator = SourceObjects.CalibratorSource('leakage',
                                  leak_cal_field_name, leak_cal_field_id, 
                                  leak_cal_scan_ids)
            self.poscalibrator = SourceObjects.CalibratorSource('position',
                                 pos_cal_field_name, pos_cal_field_id, 
                                 pos_cal_scan_ids)
            self.anglecalibrator = SourceObjects.CalibratorSource('angle',
                                   angle_cal_field_name, angle_cal_field_id, 
                                   angle_cal_scan_ids)

        # Targets (list of an arbitrary number of target sources)
        self.targetsources = []
        tgt_field_ids = self.get_field_id_from_scan_id(tgt_scan_ids)
        for tgt_field_id in tgt_field_ids:
            tgt_field_name = self.get_field_name_from_field_id(tgt_field_id)
            tgt_scan_ids = self.get_scan_ids_from_field_id(tgt_field_id)
            target = SourceObjects.TargetSource(tgt_field_name, tgt_field_id, 
                                                tgt_scan_ids)
            self.targetsources.append(target)
        
        # List of splitt off target measurment sets (which are also instances
        # of this class). Note that for a target mset this list is empty.
        self.targetmsets = []


    # Private Methods
    # (actually private methods don't exist in Python, the _ is a convention)

    def _get_nchan(self):
        """
        Return: the number of channels
        """
        tb.open(self.file_path+'/SPECTRAL_WINDOW')
        nchan = tb.getcol('NUM_CHAN')[0]
        tb.close()
        return nchan

    def _get_telescope(self):
        """
        Return: the telscope name
        """
        tb.open(self.file_path+'/OBSERVATION')
        telescope = tb.getcol('TELESCOPE_NAME')[0]
        tb.close()
        return telescope

    def _get_band(self):
        """
        Return telescope band
        Note that if you add more bands later on, you also need to add more
        primary beam attenuations in skymodel.py
        """
        if self.telescope == 'GMRT':
            if self.freq > 650e6: return '1420'
            if self.freq > 550e6 and self.freq < 650e6: return '610'
            if self.freq > 250e6 and self.freq < 350e6: return '325'
            if self.freq > 200e6 and self.freq < 300e6: return '235'
            if self.freq < 200e6: return '151'
        elif self.telescope == 'EVLA':
            if self.freq < 1e9: return 'P'
            if self.freq >= 1e9: return 'L'

    def _get_freq(self):
        """
        Return: the reference frequency
        """
        tb.open(self.file_path+'/SPECTRAL_WINDOW')
        freq = tb.getcol('REF_FREQUENCY')[0]
        tb.close()
        return freq

    def _get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        tb.open( '%s/ANTENNA' % self.file_path)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def _get_minBL_for_cal(self):
        """
        Return: estimate the minimum BL for calibration steps
        """
        num_antenna = len(self._get_antenna_names())
        return max(3,int(num_antenna/4.0))

    def _convert_flag(self, flag_string=''):
        """
        Convert flag command from the config file to a casa command
        e.g.: C14,E03,E04,S01,W01=; =22:30:00~22:43:00; C03=22:52:30~22:55:30
        """
        flag = {}
        if flag_string != '':
            for flag_group in flag_string.replace(' ','').split(';'):
                ant = flag_group.split('=')[0]
                time = flag_group.split('=')[1]
                flag[ant] = time
            
        return flag

    def _determine_cal_scan_ids(self):
        """
        Save the calibrator scans in a list
        If empy list given, then check for coords
        """
        known_cals = {
        '3C147':{'m0': {'unit': 'rad', 'value': 1.49488177653836},
                'm1': {'unit': 'rad', 'value': 0.87008056907685105},
                'refer': 'J2000',
                'type': 'direction'},
        '3C196':{'m0': {'unit': 'rad', 'value': 2.1537362969610028},
                'm1': {'unit': 'rad', 'value': 0.841554132080366},
                'refer': 'J2000',
                'type': 'direction'},
        '3C286':{'m0': {'unit': 'rad', 'value': 3.5392586514514845},
                'm1': {'unit': 'rad', 'value': 0.53248541037303654},
                'refer': 'J2000',
                'type': 'direction'},
        '3C295':{'m0': {'unit': 'rad', 'value': 3.7146787856873482},
                'm1': {'unit': 'rad', 'value': 0.91111035090915105},
                'refer': 'J2000',
                'type': 'direction'},
        '3C380':{'m0': {'unit': 'rad', 'value': 4.8412379124131713},
                'm1': {'unit': 'rad', 'value': 0.85078013643188044},
                'refer': 'J2000',
                'type': 'direction'},
        '3C48':{'m0': {'unit': 'rad', 'value': 0.42624576436309852},
                'm1': {'unit': 'rad', 'value': 0.57874633182450852},
                'refer': 'J2000',
                'type': 'direction'}
        }
                
        ms.open(self.file_path)
        cal_scan_ids = []
        for cal_scan_id in self.scansummary.keys():
            cal_field_id = self.get_field_id_from_scan_id(cal_scan_id)
            direc = ms.getfielddirmeas(fieldid=int(cal_field_id))
            # if distance with known cal is < than 60" then add it
            for known_cal, known_cal_dir in known_cals.iteritems():
                if utils.angularSeparationOfDirectionsArcsec(
                         direc, known_cal_dir) <= 60:
                    logging.info('Found '+known_cal+' in scan: '
                                 +cal_scan_id)
                    cal_scan_ids.append(cal_scan_id)

                    # update field name for SetJy
                    tb.open('%s/FIELD' % self.file_path, nomodify=False)
                    tb.putcell('NAME', int(cal_field_id), known_cal)
                    source_id = tb.getcell('SOURCE_ID', int(cal_field_id))
                    tb.close()
                    tb.open('%s/SOURCE' % self.file_path, nomodify=False)
                    tb.putcell('NAME', source_id, known_cal)
                    tb.close()

                    break # cal found, useless to keep on iterating
        ms.close()

        if cal_scan_ids == []:
            logging.info('No calibrators found, this ms only contains targets')

        return cal_scan_ids

    # Public Methods

    def get_scan_ids_from_field_id(self, field_id, user_scan_ids = []):
        """
        This method returns a list of scan ids for a given field_id while
        possibly also excluding scans that are unwanted.
        Note that user_scan_ids is a list of wanted scans.
        """
        scans = []
        for scan_id in self.scansummary.keys():
            test_field_id = self.get_field_id_from_scan_id(scan_id)
            if int(test_field_id) == int(field_id):
                scans.append(scan_id)
        if len(user_scan_ids) > 0:
            scans = list(set(scans).intersection(user_scan_ids))
        return scans

    def get_field_name_from_field_id(self, field):
        field_name = self.summary['field_'+str(field)]['name']
        return field_name
 
    def get_field_name_from_scan_id(self, scan):
        """
        this method returns the field name of a scan id 
        or returns a list of field names of a list of scan ids
        """
        if isinstance(scan, list):
            fieldnames = []
            for i in scan:
                field_name = self.summary['scan_'+str(i)]['0']['FieldName']
                if field_name not in fieldnames:
                    fieldnames.append(field_name)
            return fieldnames
        else:
            field_name = self.summary['scan_'+str(scan)]['0']['FieldName']
            return field_name

    def get_field_id_from_scan_id(self, scan):
        """
        This method returns the field id of a scan id 
        or returns a list of field ids of a list of scan ids
        """
        if isinstance(scan, list):
            fieldids = []
            for i in scan:
                field_id = self.summary['scan_'+str(i)]['0']['FieldId']
                if field_id not in fieldids:
                    fieldids.append(field_id)
            return fieldids
        else:
            field_id = self.summary['scan_'+str(scan)]['0']['FieldId']
            return str(field_id)
        
    def get_direction_from_tgt_field_id(self, tgt):
        direction = self.summary['field_'+str(tgt)]['direction']
        return direction



class STObj:
    """
    Class used to provide information on Solution Tables
    """
    def __init__(self, file_path):
        self.file_path = file_path
        self.st_type = self._get_type()

    # Private Methods
    # (actually private methods don't exist in Python, the _ is a convention)

    def _get_type(self):
        """
        Return the Cal Table type
        """
        tb.open(self.file_path)
        st_type = tb.getkeyword('VisCal')
        tb.close()
        return st_type

    # Public Methods

    def plot(self, plotdirectory, phase_only = False, amp_only = False):
        """
        For G Jones tables plot both amp and phase, if phase_only is True, plot
        only phase.
        """
        if not os.path.isdir(plotdirectory):
            os.makedirs(plotdirectory)

        if self.st_type == 'K Jones':
            # Plot delay
            tb.open( '%s/ANTENNA' % self.file_path)
            nameAntenna = tb.getcol( 'NAME' )
            numAntenna = len(nameAntenna)
            tb.close()
            nplots=int(numAntenna/3)
            for ii in range(nplots):
                filename=plotdirectory+'/'+'d_'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                plotcal(caltable=self.file_path,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                        overplot=False,clearpanel='All',iteration='antenna',plotrange=[],\
                        plotsymbol='o-',markersize=5.0,fontsize=8.0,showgui=False,\
                        figfile=filename)
                    
        if self.st_type == 'G Jones':
            tb.open( '%s/ANTENNA' % self.file_path)
            nameAntenna = tb.getcol( 'NAME' )
            numAntenna = len(nameAntenna)
            tb.close()
            nplots=int(numAntenna/3)
            # Plot amp
            if not phase_only:
                tb.open(self.file_path)
                cpar=tb.getcol('CPARAM')
                flgs=tb.getcol('FLAG')
                tb.close()
                amps=np.abs(cpar)
                good=np.logical_not(flgs)
                plotmax=np.max(amps[good])
                BL = False # Not used in this pipeline
                for ii in range(nplots):
                    filename=plotdirectory+'/'+'a_'+str(ii)+'.png'
                    syscommand='rm -rf '+filename
                    os.system(syscommand)
                    antPlot=str(ii*3)+'~'+str(ii*3+2)
                    if BL: xaxis = 'antenna2'
                    else: xaxis = 'time'
                    if BL: plotsymbol = 'o'
                    else: plotsymbol = 'o-'
                    plotcal(caltable=self.file_path,xaxis=xaxis,yaxis='amp',antenna=antPlot,subplot=311,\
                            iteration='antenna',plotrange=[0,0,0,plotmax],plotsymbol=plotsymbol,plotcolor='red',\
                            markersize=5.0,fontsize=8.0,showgui=False,figfile=filename,clearpanel='All')
            # Plot phase
            if not amp_only:
                for ii in range(nplots):
                    filename=plotdirectory+'/'+'p_'+str(ii)+'.png'
                    syscommand='rm -rf '+filename
                    os.system(syscommand)
                    antPlot=str(ii*3)+'~'+str(ii*3+2)
                    BL = False # Not used in this pipeline
                    if BL: xaxis = 'antenna2'
                    else: xaxis = 'time'
                    plotcal(caltable=self.file_path,xaxis=xaxis,yaxis='phase',antenna=antPlot,subplot=311,\
                            overplot=False,clearpanel='All',iteration='antenna',plotrange=[0,0,-180,180],\
                            plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=8.0,showgui=False,\
                            figfile=filename)

        if self.st_type == 'B Jones':
            # Plot bandpass
            tb.open(self.file_path)
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

            tb.open( '%s/ANTENNA' % self.file_path)
            nameAntenna = tb.getcol( 'NAME' )
            numAntenna = len(nameAntenna)
            tb.close()
            nplots=int(numAntenna/3)

            for ii in range(nplots):
                filename=plotdirectory+'/'+'BP_a_'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                plotcal(caltable=self.file_path,xaxis='freq',yaxis='amp',antenna=antPlot,subplot=311,\
                        iteration='antenna',plotrange=[0,0,0,ampplotmax],showflags=False,plotsymbol='o',\
                        plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

            for ii in range(nplots):
                filename=plotdirectory+'/'+'BP_p_'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                plotcal(caltable=self.file_path,xaxis='freq',yaxis='phase',antenna=antPlot,subplot=311,\
                        iteration='antenna',plotrange=[0,0,-phaseplotmax,phaseplotmax],showflags=False,\
                        plotsymbol='o',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)
