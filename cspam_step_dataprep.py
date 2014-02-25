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

def cspam_step_dataprep(MSs, conf):
    """
    Data preparation step
    """

    logging.info('### STARTING DATA PREPARATION STEP')

    for i, MS in enumerate(MSs):

        # setting environment
        if os.path.exists('img/ms'+str(i)):
            logging.warning('Cleaning ./img/ms'+str(i)+'/ directory')
            os.system('rm -r img/ms'+str(i))
        os.makedirs('img/ms'+str(i))
        if os.path.exists('cal/ms'+str(i)):
            logging.warning('Cleaning ./cal/ms'+str(i)+'/ directory')
            os.system('rm -r cal/ms'+str(i))
        os.makedirs('cal/ms'+str(i))
        if os.path.exists('plots/ms'+str(i)):
            logging.warning('Cleaning ./plots/ms'+str(i)+'/ directory')
            os.system('rm -r plots/ms'+str(i))
        os.makedirs('plots/ms'+str(i))

        # create listobs.txt for references
        default('listobs')
        if not os.path.isfile('listobs-ms'+str(i)+'.txt'):
            listobs(vis=MS.file_name, verbose=True, listfile='listobs-ms'+str(i)+'.txt')

        # plot ants
        default('plotants')
        plotants(vis=MS.file_name, figfile='plots/ms'+str(i)+'/plotants.png')

        # plot elevation
        default('plotms')
        plotms(vis=MS.file_name, xaxis='time', yaxis='elevation', selectdata=True, antenna='0&1;2&3',\
            spw='0:31', coloraxis='field', plotfile='plots/ms'+str(i)+'/el_vs_time.png', overwrite=True)

        # report initial statistics
        statsflags = getStatsflag(MS.file_name)
        logging.info("Initial flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%")

        # Flag bad channel for GMRT
        if MS.telescope == 'GMRT':
            if MS.nchan == 512:
                if freq > 600e6 and freq < 650e6: spw='0:0~10,0:502~511' # 610 MHz
                if freq > 300e6 and freq < 350e6: spw='0:0~10,0:502~511' # 325 MHz
                if freq > 200e6 and freq < 300e6: spw='0:0~130,0:450~511' # 235 MHz +20 border
            elif MS.nchan == 256:
                if freq > 600e6 and freq < 650e6: spw='0:0~5,0:251~255' # 610 MHz
                if freq > 300e6 and freq < 350e6: spw='0:0~5,0:251~255' # 325 MHz
                if freq > 200e6 and freq < 300e6: spw='0:0~65,0:225~255' # 235 MHz +20 border

            default('flagdata')
            flagdata(vis=MS.file_name, mode='manualflag', spw=spw, flagbackup=False)

        if MS.flag != {}:
            for badant in MS.flag:
                print "* Flagging :", badant, " - time: ", MS.flag[badant]
                default('flagdata')
                flagdata(vis=MS.file_name, mode='manualflag', antenna=badant,\
                    timerange=MS.flag[badant], flagbackup=False)

        # quack
        default('flagdata')
        flagdata(vis=MS.file_name, mode='quack', quackinterval=1, quackmode='beg', action='apply', flagbackup=False)

        # flag zeros
        default('flagdata')
        flagdata(vis=MS.file_name, mode='clip', clipzeros=True,\
            correlation='ABS_ALL', action='apply', flagbackup=False)

        # flag statistics after pre-flag
        statsflags = getStatsflag(MS.file_name)
        logging.info("After initial-flagging flag percentage: " + str(statsflags['flagged']/statsflags['total']*100.) + "%")

        # save flag status
        default('flagmanager')
        flagmanager(vis=MS.file_name, mode='save', versionname='AfterInitialFlagging', comment=str(datetime.datetime.now()))
        logging.info("Saved flags in AfterInitialFlagging")

    return 0
