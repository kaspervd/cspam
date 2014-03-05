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
import numpy as np

class STobj():
    """
    Class used to provide information on Solution Tables
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self.st_type = self.get_type()

    def get_type(self):
        """
        Return the Cal Table type
        """
        tb.open(self.file_name)
        st_type = tb.getkeyword('VisCal')
        tb.close()
        return st_type

    def get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        tb.open( '%s/ANTENNA' % self.file_name)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def get_antenna_coords(self):
        """
        Return: list of antenna coords
        """
        tb.open( '%s/ANTENNA' % self.file_name)
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
        colnames = VAL, ERR, TIME, FLAG, ANT, SPW, SCAN
        """
        tb.open(self.file_name)
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
        elif colname == 'FLAG':
            val = tb.getcol('FLAG')
        elif colname == 'TIME':
            val = tb.getcol('TIME')
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


    def get_val_iter(self, flag=False, return_axis='VAL', iter_axes=['SPW','SCAN','ANT'], **parm):
        """
        Return an iterator which iterates along SPW, SCAN, ANT
        flag = True|False : return only unflagged data
        return_axis = 'VAL' : the returned axis
        iter_axes = specified another set of iteration axes
        **parm = for each axes can specify restrictions

        Return: dict of coords and np.array of vals
        """
        import itertools
        assert return_axis not in iter_axes

        # get return axis values
        return_axis_val = self.get_col(return_axis)

        # get iter_axes unique and complete values
        iter_axes_val = []
        iter_axes_val_u = []
        for iter_axis in iter_axes:
            iter_axes_val.append(self.get_col(iter_axis))
            iter_axes_val_u.append(np.unique(iter_axes_val[-1]))

        # create mask
        for this_coords in list(itertools.product(*iter_axes_val_u)):
            coords = {}
            mask = np.ones(len(return_axis_val[-1]))
            for i, iter_axis in enumerate(iter_axes):
                mask = np.logical_and(mask, iter_axes_val[i] == this_coords[i])
                coords[iter_axis] = this_coords[i]
            yield coords, return_axis_val[:,:,mask]


    def get_min_max(self, aptype = 'a'):
        """
        Return minumum and maximum unflagged val
        aptype can be: 'a' or 'p' (return val in rad)
        """
        if 'a' in ctype: val = np.abs( self.get_col('VAL') )
        if 'p' in ctype: val = np.arctan2(np.imag( self.get_col('VAL') ),np.real( self.get_col('VAL') ))
        good = np.logical_not(self.get_col('FLAG'))
        maxval = np.max(val[good])
        minval = np.min(val[good])
        return minval, maxval


def unwrap_phase( x, window = 10, alpha = 0.01, iterations = 3,
    clip_range = [ 170., 180. ] ): 
    """
    Unwrap the x array, if it is shorter than 2*window, use np.unwrap()
    """

    if len(x) < 2*window: return np.unwrap(x)

    xx = array( x, dtype = float64 )
#   a = zeros( ( window ), dtype = float64 )
#   a[ -1 ] = 1.
    if ( len( clip_range ) == 2 ):
        o = clip_range[ 0 ]
        s = ( clip_range[ 1 ] - clip_range[ 0 ] ) / 90.
    a = ones( ( window ), dtype = float64 ) / float( window )
    xs = xx[ 0 ]
    for j in range( 2 * iterations ):
        for k in range( window, len( x ) ): 
            xi = xx[ k - window : k ]
            xp = dot( xi, a )
            e = xx[ k ] - xp
            e = amodulo( e + 180., 360. ) - 180.
            if ( len( clip_range ) == 2 ):
                if ( abs( e ) > o ):
                    e = sign( e ) * ( s * degrees( atan( radians( ( abs( e ) - o ) / s ) ) ) + o )  
            xx[ k ] = xp + e
#           a = a + xx[ k - window : k ] * alpha * e / 360.
            a = a + xi * alpha * e / ( dot( xi, xi ) + 1.e-6 )
        xx = xx[ : : -1 ].copy()
    xx = xx - xx[ 0 ] + xs 
    return xx


def sol_filter(calt, ctype=''):
    """
    Do complex filtering on solution tables
    """
    ST = STobj(calt) 
    
    # find the best reference antenna
    # looking at flagged data
    

    # interpolate for missing points


    # de-trend
    scipy.signal.detrend

    # filtering



def plot_cal_table(calt, MS):
    """
    Do the standard plot of gain solutions
    """
    ST = STobj(calt, ctype) 
    nplots = len(ST.get_antenna_names())/3

    if 'G' in ctype:
        if 'a' in ctype:
            plotmin, plotmax = ST.get_min_max('a')
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_a'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='time',yaxis='amp',antenna=antPlot,subplot=311,\
                    iteration='antenna',plotrange=[0,0,plotmin,plotmax],plotsymbol='o-',plotcolor='red',\
                    markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

        if 'p' in ctype:
            for ii in range(nplots):
                plotmin, plotmax = ST.get_min_max('p')*180/np.pi
                filename=MS.dir_plot+calt.split('/')[-1]+'_p'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='time',yaxis='phase',antenna=antPlot,subplot=311,\
                    overplot=False,clearpanel='Auto',iteration='antenna',plotrange=[0,0,plotmin,plotmax],\
                    plotsymbol='o-',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,\
                    figfile=filename)

    if ctype == 'K':
        for ii in range(nplots):
            filename=MS.dir_plot+calt.split('/')[-1]+'_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            default('plotcal')
            plotcal(caltable=calt,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                iteration='antenna',plotsymbol='o-',plotcolor='green',\
                markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

    if 'B' in ctype:
        if 'a' in ctype:
            plotmin, plotmax = ST.get_min_max('a')
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_a'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='freq',yaxis='amp',antenna=antPlot,subplot=311,\
                    iteration='antenna',plotrange=[0,0,plotmin,plotmax],showflags=False,plotsymbol='o',\
                    plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

        if 'p' in ctype:
            plotmin, plotmax = ST.get_min_max('p')*180/np.pi
            for ii in range(nplots):
                filename=MS.dir_plot+calt.split('/')[-1]+'_p'+str(ii)+'.png'
                syscommand='rm -rf '+filename
                os.system(syscommand)
                antPlot=str(ii*3)+'~'+str(ii*3+2)
                default('plotcal')
                plotcal(caltable=calt,xaxis='freq',yaxis='phase',antenna=antPlot,subplot=311,\
                    iteration='antenna',plotrange=[0,0,plotmin,plotmax],showflags=False,\
                    plotsymbol='o',plotcolor='blue',markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)
 

def flag_low_Ba(calt, p = 5):
    """
    Flag from the given Ba table the channels below "p" percentage from the mean
    """
    tb.open(calt, nomodify=False)
    cpar = np.abs(tb.getcol('CPARAM'))
    flags = tb.getcol('FLAG')
    flags_sum = np.sum(flags)
    antennas = tb.getcol('ANTENNA1')
    # cpar has shape pol:chan:ant
    npol, nchan, nant = cpar.shape
    for pol in xrange(npol):
        for ant in xrange(nant):
            bad_chans = np.where(cpar[pol,:,ant] < p/100.*max(cpar[pol,:,ant]))
            flags[pol,bad_chans,ant] = 1
    tb.putcol('FLAG', flags)
    tb.close()
    logging.info("Flagged channel below "+str(p)+"%: "+str((np.sum(flags)-flags_sum)/float(flags_sum)*100)+"%")


