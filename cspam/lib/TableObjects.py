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
from scipy import interpolate

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
        
    def get_all_available_timestamps(self):
        tb.open(self.file_path)
        timestamps = tb.getcol('TIME')
        tb.close()
        return timestamps



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
            #tb.open( '%s/ANTENNA' % self.file_path)
            #nameAntenna = tb.getcol( 'NAME' )
            #numAntenna = len(nameAntenna)
            #tb.close()
            
            tb.open(self.file_path)
            ants1 = tb.getcol('ANTENNA1') # Antenna
            nameAntenna = np.sort(np.unique(ants1))
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
                        plotsymbol='o-',markersize=5.0,fontsize=10.0,showgui=False,\
                        figfile=filename)
                    
        if self.st_type == 'G Jones':
            #tb.open( '%s/ANTENNA' % self.file_path)
            #nameAntenna = tb.getcol( 'NAME' )
            #numAntenna = len(nameAntenna)
            #tb.close()
            
            tb.open(self.file_path)
            ants1 = tb.getcol('ANTENNA1') # Antenna
            nameAntenna = np.sort(np.unique(ants1))
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
                            markersize=5.0,fontsize=10.0,showgui=False,figfile=filename,clearpanel='All')
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
                            iteration='antenna',plotrange=[0,0,-180,180],showflags=False,\
                            plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,\
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

            #tb.open( '%s/ANTENNA' % self.file_path)
            #nameAntenna = tb.getcol( 'NAME' )
            #numAntenna = len(nameAntenna)
            #tb.close()
            
            tb.open(self.file_path)
            ants1 = tb.getcol('ANTENNA1') # Antenna
            nameAntenna = np.sort(np.unique(ants1))
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

    def invert_table(self):
        """
        Invert a calibration table and return the new calibration table path
        """
        syscommand = "cp -r "+self.file_path+" "+self.file_path+"_inv"
        os.system(syscommand)
        caltab = self.file_path+"_inv"
        tb.open(caltab, nomodify=False) # open the caltable
        gVals = tb.getcol('CPARAM')#, startrow=start, nrow=incr) # get the values from the GAIN column
        mask = abs(gVals) > 0.0 # only consider non-zero values
        gVals[mask] = 1.0 / gVals[mask] # do the inversion
        tb.putcol('CPARAM', gVals)#, startrow=start, nrow=incr) # replace the GAIN values with the inverted values
        tb.close() # close the table
        return caltab

    def re_reference_table(self, refant=1):
        """
        CASA's gaincal changes the reference antenna if the previous reference
        antenna is flagged. This method sets the varying reference antennas to 
        one antenna.
        """

        # Create a copy of the calibration table
        syscommand = 'rm -rf '+self.file_path+"_reref"
        os.system(syscommand)
        syscommand = "cp -rf "+self.file_path+" "+self.file_path+"_reref"
        os.system(syscommand)
        caltab = self.file_path+"_reref"

        # Open the calibration table and fetch needed parameters
        tb.open(caltab, nomodify=False)
        times = tb.getcol('TIME')
        gains = tb.getcol('CPARAM')
        gainsRR = gains[0,0,:] # RR polarization
        gainsLL = gains[1,0,:] # LL polarization
        ants1 = tb.getcol('ANTENNA1') # Antenna
        ants2 = tb.getcol('ANTENNA2') # Current reference antenna
        flags = tb.getcol('FLAG')
        flagsRR = flags[0,0,:] # RR polarization
        flagsLL = flags[1,0,:] # LL polarization

        # Create empty lists to be filled
        updatedRefant = []
        updatedgainRR = []
        updatedgainLL = []
        updatedflagRR = []
        updatedflagLL = []
        
        # Create numpy array for easy access (drawback is that all values in
        # this 2D array are complex numbers now)
        table = np.column_stack((times, ants1, ants2, gainsRR, gainsLL, flagsRR, flagsLL))
        for time, ant1, ant2, gRR, gLL, fRR, fLL in table:
            if int(np.real(ant2)) is not refant:
                # This gain needs to be re-referenced. In order to do this we
                # need the gain corresponding with the new reference antenna.
                # I.e.:  g    = g    . g*   / |g   |
                #         i,n    i,o    n,o     n,o
                # where n is the new reference antenna, o the old one, | |
                # denotes the absolute value and * is the complex conjugate.

                # Find the gain relating the old reference antenna to the new
                # one.
                rerefgainRR = None
                rerefgainLL = None
                rerefflagRR = None
                rerefflagLL = None
                for time_2, ant1_2, ant2_2, gRR_2, gLL_2, fRR_2, fLL_2 in table:
                    if time_2 == time and int(np.real(ant1_2)) == refant and ant2_2 == ant2:
                        # This is the gain we need for re-referencing.
                        rerefgainRR = gRR_2
                        rerefgainLL = gLL_2
                        rerefflagRR = fRR_2
                        rerefflagLL = fLL_2
                        break # No need to search any further
                
                # Quit if no re-reference gain was found.
                if rerefgainRR == None:
                    print 'Unable to find correct antenna for re-referencing'
                    sys.exit()
                
                # Calculate new gains with respect to new reference antenna
                newgainRR = gRR * np.conj(rerefgainRR) / np.absolute(rerefgainRR)
                newgainLL = gLL * np.conj(rerefgainLL) / np.absolute(rerefgainLL)
                
                # Set flags to false if one of antennas is flagged
                newflagRR = True
                if bool(np.real(rerefflagRR)) is False and bool(np.real(fRR)) is False:
                    newflagRR = False
                
                newflagLL = True
                if bool(np.real(rerefflagLL)) is False and bool(np.real(fLL)) is False:
                    newflagLL = False
                
                # Correct for the fact that np.column_stack casts all values tp
                # complex numbers.
                updatedRefant.append( refant )
                updatedgainRR.append( newgainRR )
                updatedgainLL.append( newgainLL )
                updatedflagRR.append( newflagRR )
                updatedflagLL.append( newflagLL )
            else:
                # This gain is OK
                # Correct for the fact that np.column_stack casts all values tp
                # complex numbers.
                updatedRefant.append( int(np.real(ant2)) )
                updatedgainRR.append( gRR )
                updatedgainLL.append( gLL )
                updatedflagRR.append( bool(np.real(fRR)) )
                updatedflagLL.append( bool(np.real(fLL)) )
                        
        # Add new data to existing calibration table
        gains[0,0,:] = updatedgainRR
        gains[1,0,:] = updatedgainLL
        flags[0,0,:] = updatedflagRR
        flags[1,0,:] = updatedflagLL

        tb.putcol('CPARAM', gains)
        tb.putcol('ANTENNA2', updatedRefant)
        tb.putcol('FLAG', flags)

        tb.close()

        return caltab

    def normalize_reference_antenna(self):
        """
        Somehow, the calibration tables created by CASA's gaincal have a
        changing constant offset with respect to the reference antenna. For
        example, one would expect the reference antenna to have a phase of zero
        everywhere (after all, your measure the phase with respect to the phase
        of the reference antenna), but this is not the case. The phase of the
        reference antenna changes with constant values in time. This method
        corrects for this.
        """

        # Create a copy of the calibration table
        syscommand = 'rm -rf '+self.file_path+"_normref"
        os.system(syscommand)
        syscommand = "cp -rf "+self.file_path+" "+self.file_path+"_normref"
        os.system(syscommand)
        caltab = self.file_path+"_normref"

        # Open the calibration table and fetch needed parameters
        tb.open(caltab, nomodify=False)
        times = tb.getcol('TIME')
        gains = tb.getcol('CPARAM')
        gainsRR = gains[0,0,:] # RR polarization
        gainsLL = gains[1,0,:] # LL polarization
        ants1 = tb.getcol('ANTENNA1') # Antenna
        ants2 = tb.getcol('ANTENNA2') # Current reference antenna

        # Create empty lists to be filled
        updatedgainRR = gainsRR
        updatedgainLL = gainsLL

        # Create numpy array for easy access (drawback is that all values in
        # this 2D array are complex numbers now)
        table = np.column_stack((times, ants1, ants2, gainsRR, gainsLL))
        updated_table = table
        for time, ant1, ant2, gRR, gLL in table:
            if ant1 == ant2:
                if gRR != (1+0j) or gLL != (1+0j):
                    # Gains with this timestamp need to be fixed. Find the
                    # the gains with corresponding times.
                    for i, [time_2, ant1_2, ant2_2, gRR_2, gLL_2] in enumerate(table):
                        if time == time_2:
                            # Check which polarizations need to be fixed
                            if gRR != (1+0j):
                                # Fix RR polarization
                                updated_gRR_2 = gRR_2 * np.conjugate(gRR) / np.absolute(gRR)
                                updatedgainRR[i] = updated_gRR_2
                        
                            if gLL != (1+0j):
                                # Fix LL polarization
                                updated_gLL_2 = gLL_2 * np.conjugate(gLL) / np.absolute(gLL)
                                updatedgainLL[i] = updated_gLL_2

        # Add new data to existing calibration table
        gains[0,0,:] = updatedgainRR
        gains[1,0,:] = updatedgainLL
        tb.putcol('CPARAM', gains)
        
        tb.close()

        return caltab

    def resample_solutions(self, interpolationtimes, interp_type = 'spline'):
        """
        Calibration tables created by CASA's gaincal can have a smaller number
        of timestamps than the original measurement set (e.g. if the solution
        interval was larger than 'int'). This method resamples the calibration
        table to the original number of timestamps available in the measurement
        set by interpolation.
        
        Note that the calibration table created by this functions does not have
        correct data in the following columns: FIELD_ID, INTERVAL, 
        OBSERVATION_ID, PARAMERR, SCAN_NUMBER, SNR, SPECTRAL_WINDOW_ID (these
        data are simply not interpolated). So that means that using these
        columns in applycal is not possible.
        """
        
        # Only use the unique timestamps
        newtimes = np.sort(np.unique(interpolationtimes))

        # Create a copy of the calibration table
        syscommand = 'rm -rf '+self.file_path+'_resamp'
        os.system(syscommand)
        syscommand = "cp -rf "+self.file_path+" "+self.file_path+"_resamp"
        os.system(syscommand)
        caltab = self.file_path+'_resamp'

        # Open the calibration table and fetch needed parameters
        tb.open(self.file_path, nomodify=False)
        
        tbkeyw = tb.getkeywords()
        tbcoltime = tb.getcolkeywords('TIME')
        tbcolant1 = tb.getcolkeywords('ANTENNA1')
        tbcolant2 = tb.getcolkeywords('ANTENNA2')
        tbcolcpar = tb.getcolkeywords('CPARAM')
        tbcolflag = tb.getcolkeywords('FLAG')
        
        tbcolfieldid = tb.getcolkeywords('FIELD_ID')
        tbcolint = tb.getcolkeywords('INTERVAL')
        tbcolobsid = tb.getcolkeywords('OBSERVATION_ID')
        tbcolpar = tb.getcolkeywords('PARAMERR')
        tbcolscann = tb.getcolkeywords('SCAN_NUMBER')
        tbcolsnr = tb.getcolkeywords('SNR')
        tbcolspwid = tb.getcolkeywords('SPECTRAL_WINDOW_ID')
        tbcolwei = tb.getcolkeywords('WEIGHT')
        
        tbinfo = tb.info()
        times = tb.getcol('TIME')
        gains = tb.getcol('CPARAM')
        gainsRR = gains[0,0,:] # RR polarization
        gainsLL = gains[1,0,:] # LL polarization
        ants1 = tb.getcol('ANTENNA1') # Antenna
        ants2 = tb.getcol('ANTENNA2') # Current reference antenna
        flags = tb.getcol('FLAG')
        flagsRR = flags[0,0,:] # RR polarization
        flagsLL = flags[1,0,:] # LL polarization
        tabdesc = tb.getdesc()  
        dminfo  = tb.getdminfo()
                
        tb.close()
        
        # Check if all reference antennas are the same
        uniqueants = np.unique(ants2)
        if len(uniqueants) > 1:
            # Multiple reference antennas in calibration table
            print 'Multiple reference antennas in calibration table'
            sys.exit()
        else:
            refant = uniqueants[0]

        # Separate different antennas: info_per_ant is a list containing
        # for each antenna a list with the time, gains and flags for both
        # polarizations
        uniqueants = np.unique(ants1)
        info_per_ant = []
        for unique_ant in uniqueants:
            this_ants_time = []
            this_ants_gRR = []
            this_ants_gLL = []
            this_ants_fRR = []
            this_ants_fLL = []
            for time, gRR, gLL, ant, flagRR, flagLL in zip(times, gainsRR, gainsLL, ants1, flagsRR, flagsLL):
                if ant == unique_ant:
                    this_ants_time.append(time)
                    this_ants_gRR.append(gRR)
                    this_ants_fRR.append(flagRR)
                    this_ants_gLL.append(gLL)
                    this_ants_fLL.append(flagLL)
            
            info_per_ant.append([this_ants_time, this_ants_gRR, this_ants_gLL, this_ants_fRR, this_ants_fLL])
        
        # For each antenna interpolate the data
        updated_ant_info = []
        for ant_info in info_per_ant:
            ant_time = ant_info[0]
            ant_gRR = ant_info[1]
            ant_gLL = ant_info[2]
            ant_fRR = ant_info[3]
            ant_fLL = ant_info[4]
        
            ant_phase_RR = np.angle(ant_gRR)
            ant_mag_RR = np.abs(ant_gRR)
            ant_phase_LL = np.angle(ant_gLL)
            ant_mag_LL = np.abs(ant_gLL)

            
            ### REMOVE FLAGGED VALUES ###
                    
            # Although if data is flagged, gains are still stored in the table.
            # This messes up the interpolation so remove these values first.
            updated_ant_timeRR = []
            updated_ant_phaseRR = []
            updated_ant_magRR = []
            for time, phaseRR, magRR, fRR in zip(ant_time, ant_phase_RR, ant_mag_RR, ant_fRR):
                if not fRR:
                    updated_ant_timeRR.append(time)
                    updated_ant_phaseRR.append(phaseRR)
                    updated_ant_magRR.append(magRR)

            updated_ant_timeLL = []
            updated_ant_phaseLL = []
            updated_ant_magLL = []
            for time, phaseLL, magLL, fLL in zip(ant_time, ant_phase_LL, ant_mag_LL, ant_fLL):
                if not fLL:
                    updated_ant_timeLL.append(time)
                    updated_ant_phaseLL.append(phaseLL)
                    updated_ant_magLL.append(magLL)

            ### UNWRAP PHASES ###

            updated_ant_phaseRR = np.unwrap(updated_ant_phaseRR)
            updated_ant_phaseLL = np.unwrap(updated_ant_phaseLL)


            ### INTERPOLATE VALUES ###
            
            if interp_type == 'spline':
                # Create the interpolation functions and new values
                if len(updated_ant_timeRR) > 4:
                    # Total number of points must be greater than the degree of the
                    # spline
                    splineRepPhaseRR = interpolate.splrep(updated_ant_timeRR, updated_ant_phaseRR, s=0)
                    splineRepMagRR = interpolate.splrep(updated_ant_timeRR, updated_ant_magRR, s=0)
                    newphasesRR = interpolate.splev(newtimes, splineRepPhaseRR, der=0)
                    newmagRR = interpolate.splev(newtimes, splineRepMagRR, der=0)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesRR = np.linspace(0,0,len(newtimes))
                    newmagRR = np.linspace(1,1,len(newtimes))
                
                if len(updated_ant_timeLL) > 4:
                    # Total number of points must be greater than the degree of the
                    # spline
                    splineRepPhaseLL = interpolate.splrep(updated_ant_timeLL, updated_ant_phaseLL, s=0)
                    splineRepMagLL = interpolate.splrep(updated_ant_timeLL, updated_ant_magLL, s=0)
                    newphasesLL = interpolate.splev(newtimes, splineRepPhaseLL, der=0)
                    newmagLL = interpolate.splev(newtimes, splineRepMagLL, der=0)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesLL = np.linspace(0,0,len(newtimes))
                    newmagLL = np.linspace(1,1,len(newtimes))

            elif interp_type == 'linear':
                # Create the interpolation functions and new values
                if len(updated_ant_timeRR) > 2:
                    newphasesRR = np.interp(newtimes, updated_ant_timeRR, updated_ant_phaseRR)
                    newmagRR = np.interp(newtimes, updated_ant_timeRR, updated_ant_magRR)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesRR = np.linspace(0,0,len(newtimes))
                    newmagRR = np.linspace(1,1,len(newtimes))
                
                if len(updated_ant_timeLL) > 2:
                    newphasesLL = np.interp(newtimes, updated_ant_timeLL, updated_ant_phaseLL)
                    newmagLL = np.interp(newtimes, updated_ant_timeLL, updated_ant_magLL)
                else:
                    # (almost) everything is flagged, so just fill it with dummy
                    # values
                    newphasesLL = np.linspace(0,0,len(newtimes))
                    newmagLL = np.linspace(1,1,len(newtimes))

            new_ant_gRR = newmagRR * np.exp(1j * newphasesRR)
            new_ant_gLL = newmagLL * np.exp(1j * newphasesLL)

            
            ### FIX INTERPOLATION BETWEEN FLAGS ###
                    
            # Fix interpolation between flagged values. This comes basically
            # down to linear interpolation between flag values (zero and one).
            # Values < 1 will next be regarded as flagged.
            ant_fRR_num = map(int, ant_fRR)
            new_ant_fRR = np.interp(newtimes, ant_time, ant_fRR_num)
            new_ant_fRR = np.ceil(new_ant_fRR)
            new_ant_fRR = map(bool, new_ant_fRR)
            
            ant_fLL_num = map(int, ant_fLL)
            new_ant_fLL = np.interp(newtimes, ant_time, ant_fLL_num)
            new_ant_fLL = np.ceil(new_ant_fLL)
            new_ant_fLL = map(bool, new_ant_fLL)


            ### FLAG LARGE GAPS ###

            # Sometimes, for larger gaps, timestamps exist, flags are not
            # true but data is not available. This means that for the entire
            # gap the interpolater can do what it wants. So these large gaps
            # need to be flagged as well.
            # First RR
            for i in xrange(len(updated_ant_timeRR)-1):
                # Find out how many interpolated steps are between two data
                # points.
                time = updated_ant_timeRR[i]
                nexttime = updated_ant_timeRR[i+1]
                
                # Since we're dealing with a calibration table with less
                # timestamps, we need to find the timestamps in the measurement
                # set which correspond best with those in the cal table.
                time_matches = []
                nexttime_matches = []
                for match in newtimes:
                    time_matches.append([match-time, match])
                    nexttime_matches.append([match-nexttime, match])
                time_matches.sort(key=lambda x: x[0])
                nexttime_matches.sort(key=lambda x: x[0])
                best_match_with_time = time_matches[0][1]
                best_match_with_nexttime = nexttime_matches[0][1]
                
                # Calculate how many intermediate steps there are
                counter = 0
                start_counter = False
                for interp_time in newtimes:
                    if interp_time == best_match_with_time:
                        start_counter = True
                    if start_counter:
                        counter += 1
                    if interp_time == best_match_with_nexttime:
                        start_counter = False
            
                # If the numbers of steps is too large, flag the gap.
                expected_steps = int(len(newtimes)/np.float32(len(updated_ant_timeRR)) + 2)
                if counter > expected_steps:
                    # Flag this gap
                    do_flag = False
                    for i, interp_time in enumerate(newtimes):
						# Only flag values between the two data steps.
                        if interp_time == best_match_with_nexttime:
                            do_flag = False
                        if do_flag:
                            new_ant_fRR[i] = True
                        if interp_time == best_match_with_time:
                            do_flag = True
            # Also do LL
            for i in xrange(len(updated_ant_timeLL)-1):
                # Find out how many interpolated steps are between two data
                # points.
                time = updated_ant_timeLL[i]
                nexttime = updated_ant_timeLL[i+1]
                
                # Since we're dealing with a calibration table with less
                # timestamps, we need to find the timestamps in the measurement
                # set which correspond best with those in the cal table.
                time_matches = []
                nexttime_matches = []
                for match in newtimes:
                    time_matches.append([match-time, match])
                    nexttime_matches.append([match-nexttime, match])
                time_matches.sort(key=lambda x: x[0])
                nexttime_matches.sort(key=lambda x: x[0])
                best_match_with_time = time_matches[0][1]
                best_match_with_nexttime = nexttime_matches[0][1]
                
                # Calculate how many intermediate steps there are
                counter = 0
                start_counter = False
                for interp_time in newtimes:
                    if interp_time == best_match_with_time:
                        start_counter = True
                    if start_counter:
                        counter += 1
                    if interp_time == best_match_with_nexttime:
                        start_counter = False
            
                # If the numbers of steps is too large, flag the gap.
                expected_steps = int(len(newtimes)/np.float32(len(updated_ant_timeLL)) + 2)
                if counter > expected_steps:
                    # Flag this gap
                    do_flag = False
                    for i, interp_time in enumerate(newtimes):
						# Only flag values between the two data steps.
                        if interp_time == best_match_with_nexttime:
                            do_flag = False
                        if do_flag:
                            new_ant_fLL[i] = True
                        if interp_time == best_match_with_time:
                            do_flag = True


            ### SET FLAGGED VALUES TO (1+0j) ###

            # This is just the casa convention. Doesn't really matter since
            # the data is flagged.
            for i, flag in enumerate(new_ant_fRR):
                if flag:
                    new_ant_gRR[i] = (1+0j)
            for i, flag in enumerate(new_ant_fLL):
                if flag:
                    new_ant_gLL[i] = (1+0j)

            # Store this new interpolation data
            updated_ant_info.append([newtimes, new_ant_gRR, new_ant_gLL, new_ant_fRR, new_ant_fLL])

        # Put the updated info back in the correct places in a new cal table
        correct_format_times = []
        correct_format_antids = []
        correct_format_antrefids = []
        correct_format_gRR = []
        correct_format_gLL = []
        correct_format_fRR = []
        correct_format_fLL = []
        for time in newtimes:
            for ant_id, updated_info in enumerate(updated_ant_info):
                ant_times = updated_info[0]
                ant_gRR = updated_info[1]
                ant_gLL = updated_info[2]
                ant_fRR = updated_info[3]
                ant_fLL = updated_info[4]
                for i, ant_time in enumerate(ant_times):
                    if ant_time == time:
                        correct_format_times.append(ant_time)
                        correct_format_antids.append(ant_id)
                        correct_format_antrefids.append(refant)
                        correct_format_gRR.append(ant_gRR[i])
                        correct_format_gLL.append(ant_gLL[i])
                        correct_format_fRR.append(ant_fRR[i])
                        correct_format_fLL.append(ant_fLL[i])

        updatedgains = np.zeros([2,1,len(correct_format_gRR)], dtype=np.complex)
        updatedflags = np.zeros([2,1,len(correct_format_gRR)], dtype=bool)
        updatedgains[0,0,:] = correct_format_gRR
        updatedgains[1,0,:] = correct_format_gLL
        updatedflags[0,0,:] = correct_format_fRR
        updatedflags[1,0,:] = correct_format_fLL
        
        # Copy existing table information
        tb.create(caltab, tabdesc, dminfo=dminfo)
        tb.addrows(len(correct_format_times))
        tb.putinfo(tbinfo)
        tb.putkeywords(tbkeyw)
        tb.putcolkeywords('TIME', tbcoltime)
        tb.putcolkeywords('ANTENNA1', tbcolant1)
        tb.putcolkeywords('ANTENNA2', tbcolant2)
        tb.putcolkeywords('CPARAM', tbcolcpar)
        tb.putcolkeywords('FLAG', tbcolflag)
        tb.putcolkeywords('FIELD_ID', tbcolfieldid)
        tb.putcolkeywords('INTERVAL', tbcolint)
        tb.putcolkeywords('OBSERVATION_ID', tbcolobsid)
        tb.putcolkeywords('PARAMERR', tbcolpar)
        tb.putcolkeywords('SCAN_NUMBER', tbcolscann)
        tb.putcolkeywords('SNR', tbcolsnr)
        tb.putcolkeywords('SPECTRAL_WINDOW_ID', tbcolspwid)
        tb.putcolkeywords('WEIGHT', tbcolwei)
        
        # Useful stuff
        tb.putcol('TIME', correct_format_times)
        tb.putcol('CPARAM', updatedgains)
        tb.putcol('ANTENNA1', correct_format_antids)
        tb.putcol('ANTENNA2', correct_format_antrefids)
        tb.putcol('FLAG', updatedflags)
        
        # Needed but not useful
        emptylist = np.linspace(0,0,len(correct_format_times))
        secondemptylist = np.zeros_like(updatedgains, dtype=np.float)
        tb.putcol('FIELD_ID', emptylist)
        tb.putcol('INTERVAL', emptylist)
        tb.putcol('OBSERVATION_ID', emptylist)
        tb.putcol('PARAMERR', secondemptylist) # different format
        tb.putcol('SCAN_NUMBER', emptylist)
        tb.putcol('SNR', secondemptylist) # different format
        tb.putcol('SPECTRAL_WINDOW_ID', emptylist)
        # The WEIGHT COLUMN WAS EMPTY SO LEAVE IT EMPTY
        tb.close()

        return caltab
