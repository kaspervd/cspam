#!/usr/bin/env python

# use bdsm to write a catalog of sources (used for peeling)

def find_sources_to_peel(image_name, threshpix=5, threshisl=3, atrous_do=False, catalog_name=None, rmsbox=(55,12)):

    import lofar.bdsm as bdsm

    # wavelets are required to fit gaussians
    if atrous_do: stop_at = None
    else: stop_at = 'isl'

    # DO THE SOURCE DETECTION
    img = bdsm.process_image(image_name, mean_map='zero', rms_box=rmsbox, \
        thresh_pix=int(threshpix), thresh_isl=int(threshisl), atrous_do=atrous_do, atrous_jmax=3, \
        adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(80,20), \
        stop_at=stop_at, blank_limit=1e-5, quiet=True)

    # SAVE THE CATALOG
    img.write_catalog(outfile=catalog_name, format='fits', clobber=True)


if __name__=='__main__':
    import optparse
    opt = optparse.OptionParser(usage='%prog [-v|-V] imagename \n Kasper van Dam', version='1.0')
    opt.add_option('-p', '--threshpix', help='Threshold pixel (default=5)', type='int', default=5)
    opt.add_option('-i', '--threshisl', help='Threshold island (default=3)', type='int', default=3)
    opt.add_option('-t', '--atrous_do', help='BDSM extended source detection (default=False)', action='store_true', default=False)
    opt.add_option('-c', '--catalog', help='Name of the catalog (default=None -> automatically named)', default=None)
    opt.add_option('-r', '--rmsbox', help='rms box size (default=55,12)', default='55,12')
    (options, args) = opt.parse_args()
    
    rmsbox = (int(options.rmsbox.split(',')[0]),int(options.rmsbox.split(',')[1]))
    find_sources_to_peel(args[0].rstrip('/'), options.threshpix, options.threshisl, options.atrous_do, options.catalog, rmsbox)

