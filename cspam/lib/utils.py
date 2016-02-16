def stats_flag(ms, spw='', field=''):
    """
    Print (and return) the falg statistics
    """
    from casa import agentflagger as af
    af.open(ms)
    af.selectdata(field=field, spw=spw)
    agentSummary={'mode':'summary'}
    af.parseagentparameters(agentSummary)
    af.init()
    summary = af.run()
    af.done()
    del af

    array_flag = summary['report0']['array']['0']
    logging.info("Flag percentage: " + str(array_flag['flagged']/array_flag['total']*100.) + "%")
    return summary
    
def unwrap_phase( x, window = 10, alpha = 0.01, iterations = 3,
    clip_range = [ 170., 180. ] ): 
    """
    Unwrap the x array, if it is shorter than 2*window, use np.unwrap()
    """

    if len(x) < 2*window: return np.unwrap(x)

    xx = np.array( x, dtype = np.float64 )
#   a = zeros( ( window ), dtype = np.float64 )
#   a[ -1 ] = 1.
    if ( len( clip_range ) == 2 ):
        o = clip_range[ 0 ]
        s = ( clip_range[ 1 ] - clip_range[ 0 ] ) / 90.
    a = np.ones( ( window ), dtype = np.float64 ) / float( window )
    xs = xx[ 0 ]
    for j in range( 2 * iterations ):
        for k in range( window, len( x ) ): 
            xi = xx[ k - window : k ]
            xp = np.dot( xi, a )
            e = xx[ k ] - xp
            e = np.mod( e + 180., 360. ) - 180.
            if ( len( clip_range ) == 2 ):
                if ( abs( e ) > o ):
                    e = sign( e ) * ( s * degrees( atan( radians( ( abs( e ) - o ) / s ) ) ) + o )  
            xx[ k ] = xp + e
#           a = a + xx[ k - window : k ] * alpha * e / 360.
            a = a + xi * alpha * e / ( np.dot( xi, xi ) + 1.e-6 )
        xx = xx[ : : -1 ].copy()
    xx = xx - xx[ 0 ] + xs 
    return xx


# TODO: what if reference antenna changes?
def sol_filter_G(calt, window_ph = 60., window_amp = 60., order = 1, max_gap = 5.,
    max_ncycles = 10, max_rms_amp = 3., max_rms_ph=3.):
    """
    Do complex filtering on G Jones solution tables
    return:
    dict with rms values. Dict is inexed by a tuple of (scan, ant, spw, pol, chan)
    """
    ST = STobj(calt) 
    rms_dict = {}

    for coord, val in ST.get_val_iter(return_axes=['TIME','VAL','FLAG','ANTENNA2'],
        iter_axes=['SCAN_NUMBER', 'SPECTRAL_WINDOW_ID', 'ANTENNA1']):

        npol, nchan, __null = val['VAL'].shape
        # check if the reference antenna is always the same
        if not ( val['ANTENNA2'][0] == val['ANTENNA2'] ).all():
            raise "Reference antenna changing!"

        allamp = np.log10(abs(val['VAL'])) # amp in log so to have gaussian noise
        allph = np.arctan2(np.imag( val['VAL'] ),np.real( val['VAL'] ))
        new_flag = np.copy(val['FLAG'])

        # cycle across chan and polarizations
        for pol in xrange(npol):
            for chan in xrange(nchan):

                print "working on pol "+str(pol)+" - chan "+str(chan)
                print coord

                # restrict to a pol/chan and remove flagged data
                flag = val['FLAG'][pol][chan]
                amp = allamp[pol][chan][~flag]
                ph = allph[pol][chan][~flag]
                time = val['TIME'][~flag]

                # unwrap phases
                ph = unwrap_phase( ph, alpha = 0.01 )
                if ( np.nan == ph ).any:
                    ph = unwrap_phase( ph, alpha = 0.001 )
                if ( np.nan == ph ).any or max( np.fabs( ph ) ) > 1.e4:
                    ph = np.unwrap( ph )
        
                # filtering
                new_flag[pol][chan][~flag], rms_amp = outlier_rej(amp, time, max_ncycles, max_rms_amp, window_amp, order, max_gap)
                new_flag[pol][chan][~flag], rms_ph  = outlier_rej(ph, time, max_ncycles, max_rms_ph, window_ph, order, max_gap)

                # creating return dict (scan, ant, spw, pol, chan)
                rms_dict[tuple(['amp'] + [c for col, c in coord.iteritems()] + [pol] + [chan])] = np.sum(rms_amp)
                rms_dict[tuple(['ph'] + [c for col, c in coord.iteritems()] + [pol] + [chan])] = np.sum(rms_ph)

        # writing back
        ST.put_val(coord, new_flag, axis='FLAG')

    return rms_dict


def smooth(data, times, window = 60., order = 1, max_gap = 5. ):
    """
    Remove a trend from the data
    window: in sec, sliding phase window dimension
    order: 0: remove avg, 1: remove linear, 2: remove cubic
    max_gap: maximum allawed gap in minutes

    return: detrendized data array
    """
    
    final_data = np.copy(data)

    # loop over solution times
    for time in times:

        # get data to smooth (values inside the time window)

        data_array = data[ np.where( abs(times - time) <= window / 2. ) ]
        data_offsets = times[ np.where( abs(times - time) <= window / 2. ) ] - time

        # check and remove big gaps in data
        if ( len( data_offsets ) > 1 ):
          ddata_offsets = data_offsets[ 1 : ] - data_offsets[ : -1 ]
          sel = np.where( ddata_offsets > max_gap * 60. )[0]
          if ( len( sel ) > 0 ):
            min_data_index = 0
            max_data_index = len( data_offsets )
            this_time_index = np.where( abs( data_offsets ) == abs( data_offsets ).min() )[0]
            # find min and max good indexes for this window
            for s in sel:
              if ( s < this_time_index ):
                min_data_index = s + 1
              if ( s >= this_time_index ):
                max_data_index = s + 1
                break
            # redefine data arrays
            data_array = data_array[ min_data_index : max_data_index ]
            data_offsets = data_offsets[ min_data_index : max_data_index ]

        # smooth
        if len( data_array ) > 1:
          dim = min( len( data_array ) - 2, order )
          if ( dim == 0 ):
            smooth_data = np.median( data_array )
          else:
            P = np.zeros( ( len( data_offsets ), dim + 1 ), dtype = data_offsets.dtype )
            P[ : , 0 ] = 1.
            if ( dim >= 1 ):
                P[ : , 1 ] = data_offsets
            if ( dim >= 2 ):
                P[ : , 2 ] = data_offsets**2
            Pt = np.transpose( P )
            smooth_data = np.dot( np.linalg.inv( np.dot( Pt, P ) ), np.dot( Pt, data_array ) )[ 0 ]
          n = np.where( times == time )[0]
          final_data[n] = smooth_data
        
    return final_data


def outlier_rej(val, time, max_ncycles = 10, max_rms = 3., window = 60., order = 1, max_gap = 5.):
    """
    Reject outliers using a running median
    val = the array (avg must be 0)
    max_ncycles = maximum number of cycles
    max_rms = number of rms times for outlier flagging
    return: flags array and final rms
    """
    flags = np.zeros(shape=val.shape, dtype=np.bool)
    val_detrend = np.zeros(shape=val.shape)

    # DEBUG plotting
    #import matplotlib.pyplot as plt
    #import matplotlib.cm as cm
    #colors = iter(cm.rainbow(np.linspace(0, 1, max_ncycles)))
    #plt.plot( time[~flags], val[~flags], 'o-', color='k' )

    for i in xrange(max_ncycles):
        
        # smoothing (input with no flags!)
        val_smoothed = smooth(val[~flags], time[~flags], window, order, max_gap)
        val_detrend[~flags] = val[~flags] - val_smoothed

        # DEBUG plotting
        #color = next(colors)
        #plt.plot( time[~flags], val_smoothed, 'x:', label=str(i), color=color )

        # median calc
        rms =  1.4826 * np.median(abs(val_detrend[~flags]))

        # rejection  
        flags[ ~flags ] = np.logical_or(flags[~flags], abs(val_detrend[~flags]) > max_rms * rms )

        # DEBUG plotting
        #plt.plot( time[ flags ], val[ flags ], '*', markersize=20, color=color )

        if (flags == True).all():
            rms == 0.
            break
            
        # median calc
        this_rms =  1.4826 * np.median(abs(val_detrend[~flags]))
        if rms - this_rms == 0.:
            break

    # DEBUG plotting
    #plt.legend()
    #plt.show()
    return flags, rms

def plot_cal_table(calt, MS):
    """
    Do the standard plot of gain solutions
    """
    ST = STobj(calt) 
    nplots = len(ST.get_antenna_names())/3

    if 'G' == ST.st_type()[0]:
        if 'a' in calt[-3,-1]:
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

        if 'p' in calt[-3,-1]:
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

    if 'K' == ST.st_type()[0]:
        for ii in range(nplots):
            filename=MS.dir_plot+calt.split('/')[-1]+'_'+str(ii)+'.png'
            syscommand='rm -rf '+filename
            os.system(syscommand)
            antPlot=str(ii*3)+'~'+str(ii*3+2)
            default('plotcal')
            plotcal(caltable=calt,xaxis='time',yaxis='delay',antenna=antPlot,subplot=311,\
                iteration='antenna',plotsymbol='o-',plotcolor='green',\
                markersize=5.0,fontsize=10.0,showgui=False,figfile=filename)

    if 'B' == ST.st_type()[0]:
        if 'a' in calt[-3,-1]:
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

        if 'p' in calt[-3,-1]:
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
