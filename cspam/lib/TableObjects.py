# Description:
# ------------
# This file contains two classes that deal with CASA tables:
# MSObj (Measurement Set Object) and STObj (Solution Table Object).
# 
# These classes come originally from cspam version 0.1 by Francesco
# de Gasperin and Huib Intema. Edited by Kasper van Dam.
#
# These objects rely heavily on the CASA environment.

import logging
import numpy as np
#import casa, casac

"""
# The following is needed to run this outside the CASA
# environment, apparently CASA doesn't let you import
# modules while maintaining CASA functionalitites within these modules
# and I want this to be a submodule of the main code. 
casapath = '/software/casa/casa-release-4.5.2-el6/'
import os
import sys
sys.path.insert(0,casapath+'lib/python2.7/') # contains all tasks
# sys.path.insert(0,casapath+'lib/python2.7/site-packages')
# the above path contains all casa approved 3rd party modules
# use system python installation instead.
sys.path.insert(0,casapath+'lib/python2.7/__casac__') # contains all toolkits
sys.path.insert(0,casapath+'xml/')
from __casac__ import ms
ms = ms.ms() # We need an instance of ms
#
# End of this hacky work-around
"""

class MSObj:
    """
    Class used to provide information on MSs
    """
    def __init__(self, conf, name):
        self.file_path = conf['file_path']
        self.ms_name = name.replace('.ms','')
        self.ms_name = self.ms_name.replace('.MS','')

        # save summary info which are often used
        ms.open(self.file_path)
        self.summary = ms.summary()
        self.scansummary = ms.getscansummary()
        ms.close()

        self.dir_img = self.file_path.replace('.ms','')
        self.dir_img = self.dir_img.replace('.MS','')
        self.dir_img = self.dir_img+'-img'
        
        # annoying workaround to keep caltables in the same dir of MSs, so plotcal works
        self.dir_cal = os.path.dirname(self.file_path) 
        #self.dir_cal = self.file_path.replace('.MS','')+'-cal/'
        
        self.dir_plot = self.file_path.replace('.ms','')
        self.dir_plot = self.dir_plot.replace('.MS','')
        self.dir_plot = self.dir_plot+'-plot'
        self.minBL_for_cal = self.get_minBL_for_cal()
        self.nchan = self.get_nchan()
        self.freq = self.get_freq()
        self.telescope = self.get_telescope()
        assert self.telescope == 'GMRT' or self.telescope == 'EVLA'
        self.band = self.get_band()
        self.spw = conf['spw']
        central_nchan = self.nchan*conf['central_chan_percentage']/100.
        self.central_chans = str(int(round(self.nchan/2.-central_nchan/2.))) + '~' + str(int(round(self.nchan/2.+central_nchan/2.)))

        self.flag = self.set_flag(conf['flag']) # manual flag to be used with flagdata command
        self.set_cal_scan_ids(conf['cal_scans']) # sets cal_scan_ids cal_scan_names and cal_scan_unwanted
        self.tgt_scan_dict = self.set_tgt_scan_dict(conf['tgt_scans']) # dict of cal scans contining groups of relative cal+tgt scans
        
        if self.telescope == 'GMRT': self.uvrange = '>1000m'
        if self.telescope == 'EVLA': self.uvrange = ''

    def get_nchan(self):
        """
        Return: the number of channels
        """
        #LocTb = casa.table
        tb.open(self.file_path+'/SPECTRAL_WINDOW')
        nchan = tb.getcol('NUM_CHAN')[0]
        tb.close()
        return nchan

    def get_telescope(self):
        """
        Return: the telscope name
        """
        #LocTb = casa.table
        tb.open(self.file_path+'/OBSERVATION')
        telescope = tb.getcol('TELESCOPE_NAME')[0]
        tb.close()
        return telescope

    def get_band(self):
        """
        Return telescope band
        """
        if self.telescope == 'GMRT':
            if self.freq > 650e6: return '1420'
            if self.freq > 600e6 and self.freq < 650e6: return '610'
            if self.freq > 300e6 and self.freq < 350e6: return '325'
            if self.freq > 200e6 and self.freq < 300e6: return '235'
            if self.freq < 200e6: return '151'
        elif self.telescope == 'EVLA':
            if self.freq < 1e9: return 'P'
            if self.freq >= 1e9: return 'L'

    def get_freq(self):
        """
        Return: the reference frequency
        """
        #LocTb = casa.table
        tb.open(self.file_path+'/SPECTRAL_WINDOW')
        freq = tb.getcol('REF_FREQUENCY')[0]
        tb.close()
        return freq

    def get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        #LocTb = casa.table
        tb.open( '%s/ANTENNA' % self.file_path)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def get_minBL_for_cal(self):
        """
        Return: estimate the minimum BL for calibration steps
        """
        num_antenna = len(self.get_antenna_names())
        return max(3,int(num_antenna/4.0))

    def set_flag(self, flag_string=''):
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

    def set_cal_scan_ids(self, usr_cal_scan_ids=['']):
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
        cal_field_ids = []
        cal_scan_ids_unwanted = []
        for cal_scan_id in self.scansummary.keys():
            cal_field_id = self.get_field_id_from_scan_id(cal_scan_id)
            dir = ms.getfielddirmeas(fieldid=int(cal_field_id))
            # if distance with known cal is < than 60" then add it
            for known_cal, known_cal_dir in known_cals.iteritems():
                if  True: #au.angularSeparationOfDirectionsArcsec(dir, known_cal_dir) <= 60:
                    if cal_scan_id not in usr_cal_scan_ids and usr_cal_scan_ids != ['']:
                         # user do not want this scan
                         logging.info('Found '+known_cal+' in scan: '+cal_scan_id+' *** Ignored')
                         cal_scan_ids_unwanted.append(cal_scan_id)
                    else:
                        logging.info('Found '+known_cal+' in scan: '+cal_scan_id)
                        cal_scan_ids.append(cal_scan_id)
                        cal_field_ids.append(str(cal_field_id))

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
            logging.critical('No calibrators found')
            sys.exit(1)

        self.cal_scan_ids = cal_scan_ids
        self.cal_field_ids = list(set(cal_field_ids))
        self.cal_scan_ids_unwanted = cal_scan_ids_unwanted


    def set_tgt_scan_dict(self, usr_tgt_scan_ids=['']):
        """
        Set the target scans in a dict associating each calibrator to the closest in time {cal1:[tgt1,tgt2],cal2:[tgt3]}
        if tgt_scans is given, then restrict to those scans
        """

        # getting calibrators central time
        cal_times = {}
        for cal_scan_id in self.cal_scan_ids:
            begin = self.scansummary[cal_scan_id]['0']['BeginTime']
            end = self.scansummary[cal_scan_id]['0']['EndTime']
            cal_times[cal_scan_id] = begin+(end-begin)/2.

        # finding the closest cal in time for each target
        tgt_scans_dict = {}
        for tgt_scan_id, tgt_scan_data in self.scansummary.iteritems():

            if tgt_scan_id not in usr_tgt_scan_ids and usr_tgt_scan_ids != ['']: continue # user do not want this scan
            if tgt_scan_id in self.cal_scan_ids_unwanted: continue # it's a discarded calibrator

            if tgt_scan_id in self.cal_scan_ids: continue # cal scan, not a tgt scan
                
            begin = tgt_scan_data['0']['BeginTime']
            end = tgt_scan_data['0']['EndTime']
            tgt_time = begin+(end-begin)/2.

            min_time = np.inf
            for cal_scan_id, cal_time in cal_times.iteritems():
                if abs(cal_time - tgt_time) < abs(min_time - tgt_time):
                    min_time = cal_time
                    min_cal = cal_scan_id

            # add the calibrator to the dict tgt_scans_dict
            if not min_cal in tgt_scans_dict: tgt_scans_dict[min_cal] = [min_cal]
            tgt_scans_dict[min_cal].append(tgt_scan_id)

        # printing the results
        for cal_scan_id in tgt_scans_dict:
            logging.info('Calibrator '+cal_scan_id+' ('+self.get_field_name_from_scan_id(cal_scan_id)+'):')
            for tgt_scan_id in tgt_scans_dict[cal_scan_id]:
                logging.info('\t * '+tgt_scan_id+' ('+self.get_field_name_from_scan_id(tgt_scan_id)+') ')

        return tgt_scans_dict


    def get_field_name_from_field_id(self, field):
        """
        Return: the source name associated with a given field id
        """
        field_name = self.summary['field_'+str(field)]['name']
        return field_name
 
    def get_field_name_from_scan_id(self, scan):
        """
        Return: the source name associated with a given scan id
        """
        field_name = self.summary['scan_'+str(scan)]['0']['FieldName']
        return field_name

    def get_field_id_from_scan_id(self, scan):
        """
        Return: the field id associated with a given scan id
        """
        field_id = self.summary['scan_'+str(scan)]['0']['FieldId']
        return str(field_id)



class STObj:
    """
    Class used to provide information on Solution Tables
    """
    def __init__(self, file_path):
        self.file_path = file_path
        self.st_type = self.get_type()

    def get_type(self):
        """
        Return the Cal Table type
        """
        tb.open(self.file_path)
        st_type = tb.getkeyword('VisCal')
        tb.close()
        return st_type

    def get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        tb.open( '%s/ANTENNA' % self.file_path)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def get_antenna_coords(self):
        """
        Return: list of antenna coords
        """
        tb.open( '%s/ANTENNA' % self.file_path)
        antenna_pos = tb.getcol( 'POSITION' )
        tb.close()
        return antenna_pos

    def get_source_dir(self):
        """
        Return: source direction
        """

    def get_col(self, colname=''):
        """
        Return: list of colname column content
        non-default colnames = VAL, ERR, ANT, SPW, SCAN
        """
        tb.open(self.file_path)
        if colname == 'VAL':
            if self.st_type == 'K Jones':
                val = tb.getcol('FPARAM')
            else:
                val = tb.getcol('CPARAM')
        elif colname == 'ERR':
            if self.st_type == 'K Jones':
                val = tb.getcol('FPARMERR')
            else:
                val = tb.getcol('CPARMERR')
        elif colname == 'TIME':
            val = tb.getcol('TIME')
            u = tb.getcolkeywords('TIME')['QuantumUnits'][0]
            if u != 's':
                val = qa.getvalue(qa.convert(qa.quantity(val, u), 's'))
        elif colname == 'SPW':
            val = tb.getcol('SPECTRAL_WINDOW_ID')
        elif colname == 'SCAN':
            val = tb.getcol('SCAN_NUMBER')
        elif colname == 'ANT':
            val = tb.getcol('ANTENNA1')
        else:
            val = tb.getcol(colname)
        
        tb.close()
        return val

    def put_col(self, colname, val):
        """
        Set a column value
        non-default colnames = VAL, ERR, ANT, SPW, SCAN
        """
        tb.open(self.file_path, nomodify=False)
        if colname == 'VAL':
            if self.st_type == 'K Jones':
                val = tb.putcol('FPARAM', val)
            else:
                val = tb.putcol('CPARAM', val)
        elif colname == 'ERR':
            if self.st_type == 'K Jones':
                val = tb.putcol('FPARMERR', val)
            else:
                val = tb.putcol('CPARMERR', val)
        elif colname == 'SPW':
            val = tb.putcol('SPECTRAL_WINDOW_ID', val)
        elif colname == 'SCAN':
            val = tb.putcol('SCAN_NUMBER', val)
        elif colname == 'ANT':
            val = tb.putcol('ANTENNA1', val)
        else:
            val = tb.putcol(colname, val)
        
        tb.close()


#    def get_val_arr(self, antenna=[], spw=[], scan=[], pol=[]):
#        """
#        Return: multidim-array with unflagged values per antenna and spw
#        """
#        ant = self.get_col('ANT')
#        ant_u = np.unique(ant).tolist()
#        vf_ant = np.vectorize(lambda x: ant_u.index(x), otypes='i')
#        spw = self.get_col('SPW')
#        spw_u = np.unique(spw).tolist()
#        vf_spw = np.vectorize(lambda x: spw_u.index(x), otypes='i')
#        scan = self.get_col('SCAN')
#        scan_u = np.unique(scan).tolist()
#        vf_scan = np.vectorize(lambda x: scan_u.index(x), otypes='i')
#        time = self.get_col('TIME')
#        time_u = np.unique(time).tolist()
#        vf_time = np.vectorize(lambda x: time_u.index(x), otypes='i')
#
#        # extract pol/chan
#        
#
#        # create a 4-d matrix of ant/spw/scan/time
#        val = self.get_col('VAL')
#        val_matrix = np.zeros((len(np.unique(ant)), len(np.unique(spw)), len(np.unique(scan)), len(np.unique(time))), np.complex)
#        val_matrix[vf_ant(ant), vf_spw(spw), vf_scan(scan), vf_time(time)] = val
#        return val_matrix


    def get_val_iter(self, return_axes=['VAL'],
    iter_axes=['SPECTRAL_WINDOW_ID','SCAN_NUMBER','ANTENNA1']):
        """
        Return an iterator which iterates along SPW, SCAN, ANT
        return_axis = 'VAL' : the returned axis
        iter_axes = specified another set of iteration axes

        Return: dict of coords and np.array of vals
        """
        import itertools
        # get return axis values
        return_axes_val = {}
        for return_axis in return_axes:
            assert return_axis not in iter_axes
            return_axes_val[return_axis] = self.get_col(return_axis)

        # get iter_axes unique and complete values
        iter_axes_val = []
        iter_axes_val_u = []
        for iter_axis in iter_axes:
            iter_axes_val.append(self.get_col(iter_axis))
            iter_axes_val_u.append(np.unique(iter_axes_val[-1]))

        for this_coords in list(itertools.product(*iter_axes_val_u)):
            coords = {}

            # create mask
            mask = np.ones(len(return_axes_val[return_axes[0]]))
            for i, iter_axis in enumerate(iter_axes):
                mask = np.logical_and(mask, iter_axes_val[i] == this_coords[i])
                coords[iter_axis] = this_coords[i]

            # apply masks to all return_axes
            return_axes_val_m = {}
            for return_axis in return_axes:
                return_axes_val_m[return_axis] = return_axes_val[return_axis][:,:,mask]

            yield coords, return_axes_val_m


    def put_val(self, coord, val, axis=''):
        """
        Set the subset of the "axis" column identified with coords (a dict as returned by get_val_iter)
        coord: a dict as returned by get_val_iter
        val: the values
        axis: the col name (e.g. FLAG, VAL)
        """
        # fetch the whole column
        val_all = self.get_col(axis)
        # create a mask to update only selected coords
        mask = np.zeros(shape = val_all.shape )
        for c, cval in coord.iteritems():
            mask[ np.where( self.get_col(c) == cval ) ] = 1
            
        # update the new values    
        val_all[mask] = val
        self.put_col(axis, val_all)


    def get_min_max(self, aptype = 'a'):
        """
        Return minumum and maximum unflagged val
        aptype can be: 'a' or 'p' (return val in rad)
        """
        assert aptype == 'a' or aptype == 'p'
        if 'a' == aptype: val = np.abs( self.get_col('VAL') )
        elif 'p' == aptype: val = np.arctan2(np.imag( self.get_col('VAL') ),np.real( self.get_col('VAL') ))
        good = np.logical_not(self.get_col('FLAG'))
        maxval = np.max(val[good])
        minval = np.min(val[good])
        return minval, maxval
