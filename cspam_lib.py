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

class MSobj():
    """
    Class used to provide information on MSs
    """
    def __init__(self, file_name):
        self.file_name = file_name
        #self.telescope = self.get_telescope()
        #self.ncahn = self.get_nchan()
        self.flag = {}

    def get_nchan(self):
        """
        Return: the number of channels
        """
        tb.open(self.file_name+'/SPECTRAL_WINDOW')
        nchan = tb.getcol('NUM_CHAN')[0]
        tb.close()
        return nchan

    def get_telescope(self):
        """
        Return: the telscope name
        """
        tb.open(self.file_name+'/OBSERVATION')
        telescope = tb.getcol('TELESCOPE_NAME')[0]
        tb.close()
        return telescope

    def get_freq(self):
        """
        Return: the reference frequency
        """
        tb.open(self.file_name+'/SPECTRAL_WINDOW')
        freq = tb.getcol('REF_FREQUENCY')[0]
        tb.close()
        return freq

    def get_antenna_names(self):
        """
        Retrun: list of antenna names
        """
        tb.open( '%s/ANTENNA' % active_ms)
        antenna_names = tb.getcol( 'NAME' )
        tb.close()
        return antenna_names

    def get_minBL_for_cal(self):
        """
        Return: estimate the minimum BL for calibration steps
        """
        num_antenna = len(self.get_antenna_names())
        return max(3,int(numAntenna/4.0))

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
            
        self.flag = flag



