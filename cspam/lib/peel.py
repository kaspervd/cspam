import os
import numpy as np
import sys

# CSPAM Modules
import utils
import AntennaObjects
import TableObjects

# CASA Toolkits
import casac
ia = casac.casac.image()
rg = casac.casac.regionmanager()
tb = casac.casac.table()
ms = casac.casac.ms()

# CASA Tasks
from casat import immath
immath = immath.immath

from casat import ftw
ftw = ftw.ftw

from casat import ft
ft = ft.ft

from casat import uvsub
uvsub = uvsub.uvsub

from casat import clean
clean = clean.clean

from casat import fixvis
fixvis = fixvis.fixvis

from casat import gaincal
gaincal = gaincal.gaincal

from casat import applycal
applycal = applycal.applycal

from casat import split
split = split.split

from casat import exportfits
exportfits = exportfits.exportfits

from casat import delmod
delmod = delmod.delmod

from casat import imfit
imfit = imfit.imfit

from casat import imstat
imstat = imstat.imstat

sou_res = ['1arcsec']
sou_size = [5000]
expnoise = 1.e-6
rob=0.5

def extrModel(modelimg, region, cycle, compl=False, extend=None):
    """
    Extract only the part described by the region file
    from one or more (nterms>1) model img
    """
    blankedmodelimg = []

    # create a new region with a large ellipses i.e. "expand" the region
    if extend != None:
        os.system('cp '+region+' '+region.replace('.crtf','-ext.crtf'))
        region = region.replace('.crtf','-ext.crtf')
        with open(region, 'w') as f:
            f.write('#CRTFv0')
            f.write('\n')
            f.write('ellipse [['+extend[0]+', '+extend[1]+'], [900arcsec, 900arcsec], 90.00000000deg] coord=J2000, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=magenta, font=Ubuntu, fontsize=11, fontstyle=normal, usetex=false')

    for i, modelimgtt in enumerate(modelimg):
        if compl:
            # copy model
            os.system("cp -r "+modelimgtt+" "+region.replace('.crtf','')+"_compl.model.tt"+str(i))
            ia.open(region.replace('.crtf','')+"_compl.model.tt"+str(i))
            reg = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
            # set to 0 all the pixels in the region,
            # so the rest of the field is untouched
            ia.set(pixels='0', region=reg)
            ia.close()

            blankedmodelimg.append(region.replace('.crtf','')+"_compl.model.tt"+str(i))

        else:
            os.system("rm -rf "+region.replace('.crtf','')+'_model'+str(cycle)+'.tt'+str(i))            
            immath(imagename = modelimgtt, mode = 'evalexpr', expr = 'IM0',
                   region = region, outfile = region.replace('.crtf','')+'_model'+str(cycle)+'.tt'+str(i))

            blankedmodelimg.append(region.replace('.crtf','')+"_model"+str(cycle)+".tt"+str(i))

    return blankedmodelimg

def subtract(active_ms, modelimg, wprojplanes=0):
    """General function to call the necessary steps to subtract point sources
    the modelimg must have only point sources one wants to keep into the field.
    second_ms_path: MS with calibrated data in DATA
    modelimg: model of the whole sky (array of tt)
    wprojplanes: number of w-projection planes, if 0 a direct ft() will be used (best for small field)
    """
    if wprojplanes == 0:
        ft(vis=active_ms, model=modelimg, nterms=len(modelimg), usescratch=True)
    else:
        ftw(vis=active_ms, model=modelimg, nterms=len(modelimg), wprojplanes=wprojplanes, usescratch=True)
    uvsub(vis=active_ms)

def findShape(img):
    """Find a minimal shape for the source to peel
    """
    ia.open(img)
    csys = ia.coordsys()
    shape1 = ia.shape()[csys.findcoordinate('direction')['pixel'][0]]
    shape2 = ia.shape()[csys.findcoordinate('direction')['pixel'][1]]
    cell = str(int(abs(csys.increment()['numeric'][csys.findcoordinate('direction')['pixel'][0]]*180./np.pi*3600.)))+'arcsec'
    ia.close()
    shape = max(shape1, shape2)*1.0 # add extra space (TODO: DEBUG put at 1.5)
    # good image shapes
    goodvalues = np.array([6400,6144,5600,5400,5184,5000,4800,4608,4320,4096,3840,3600,3200,3072,2880,2560,2304,2048, 1600, 1536, 1200, 1024, 800, 512, 256, 128, 64])
    shape = min(goodvalues[np.where(goodvalues>=shape)])
    return shape, cell

def peel(residualMS, target_mset, modelimg, current_peel_directory, solint_data, prev_rms):

    if not os.path.isdir(current_peel_directory+'/cal'):
        os.makedirs(current_peel_directory+'/cal')
    if not os.path.isdir(current_peel_directory+'/img'):
        os.makedirs(current_peel_directory+'/img')
    if not os.path.isdir(current_peel_directory+'/plot'):
        os.makedirs(current_peel_directory+'/plot')

    # Declare some needed parameters
    solint = solint_data['solint']
    number_of_integrations = solint_data['numint']
    min_SNR_per_solution_interval = solint_data['minsnr']
    fudge_factor = solint_data['ff']
    integration_time = solint_data['inttime']
    
    # Set the peel region
    region = current_peel_directory+'/region.crtf'

    # Obtain directions for later
    orig_center = target_mset.get_direction_from_tgt_field_id(0)
    ra = orig_center['m0']['value'] # NOTE I assume radians
    dec = orig_center['m1']['value']
    orig_center_string = 'J2000 '+str(ra)+'rad '+str(dec)+'rad'

    ia.open(modelimg[0])
    reg = rg.fromtextfile(filename=region, shape=ia.shape(), csys=ia.coordsys().torecord())
    peel_source_center = reg['center']
    ra_p = peel_source_center['*1']['value']
    dec_p = peel_source_center['*2']['value']
    peel_source_center_string = 'J2000 '+str(ra_p)+'rad '+str(dec_p)+'rad'

    # Cutout everything except the peel source from the model image
    peel_source_model = extrModel(modelimg, region, 0)

    # Remove the on-the-fly model in the header (we only use the MODEL_DATA
    # column).
    delmod(vis=residualMS.file_path)
    # Add this source to the residual image
    ftw(vis=residualMS.file_path, model=peel_source_model,
        nterms=len(peel_source_model), usescratch=True, wprojplanes=512)
    uvsub(vis=residualMS.file_path, reverse = True)

    # Split off the data set in order to make sure that the DATA column
    # holds the residual map + recently added peel source
    peelMSfilepath = current_peel_directory+'/img/'+residualMS.ms_name+'_and_peelsource.ms'
    split(residualMS.file_path, peelMSfilepath)
    peelMS = TableObjects.MSObj(peelMSfilepath)

    # The cutout is smaller than the original image and placed outside the
    # center. This means that a phase shift has to be applied.
    # Shift the phase center to the peel source position
    fixvis(vis=peelMS.file_path, outputvis=peelMS.file_path,
           phasecenter = peel_source_center_string)

    # Get some values for clean
    shape, cell = findShape(peel_source_model[0])
 
    # Determine the best reference antenna
    refAntObj = AntennaObjects.RefAntHeuristics(vis=peelMS.file_path,        
                               geometry=True, flagging=True)                              
    refAnt = refAntObj.calculate()[0] 

    # In order for plotcal to work a symbolic link is needed.
    # Plotcal assumes the measurement set is in the same directory
    # as the cal table.
    syscommand = 'ln -s '+peelMS.file_path+' '+current_peel_directory+\
            '/cal/'+peelMS.ms_name+'.ms'
    os.system(syscommand)

    # Use the model from this clean session for the first selfcal cycle
    delmod(vis=peelMS.file_path)
    ft(vis=peelMS.file_path, model=peel_source_model,
       nterms=len(peel_source_model), usescratch=True)
    
    # Start the selfcalibration
    numberofcycles = 10
    usedcycles = -1
    for i in xrange(numberofcycles):
        # Create new directories
        if not os.path.isdir(current_peel_directory+'/img/cycle'+str(i)):
            os.makedirs(current_peel_directory+'/img/cycle'+str(i))
        if not os.path.isdir(current_peel_directory+'/plot/cycle'+str(i)):
            os.makedirs(current_peel_directory+'/plot/cycle'+str(i))

        # Do selfcal
        gaincal(vis=peelMS.file_path, caltable=current_peel_directory+\
                '/cal/cycle'+str(i)+'.Gp', solint=solint, minsnr=0.01, selectdata = True,
                uvrange='>50m', refant=refAnt, calmode='p')
        Gp = TableObjects.STObj(current_peel_directory+'/cal/cycle'+str(i)+'.Gp')
        Gp.plot(current_peel_directory+'/plot/cycle'+str(i), phase_only=True)

        applycal(vis=peelMS.file_path, gaintable=Gp.file_path,
                 calwt=False, flagbackup=False)
        
        # Image the result
        parms = {'vis':peelMS.file_path, 'imagename':current_peel_directory+\
                '/img/cycle'+str(i)+'/im',
                 'gridmode':'widefield', 'wprojplanes':512, 'mode':'mfs', 
                 'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 
                 'imagermode':'csclean', 'imsize':[shape], 'cell':cell, 
                 'weighting':'briggs', 'robust':rob, 'usescratch':True, 
                 'mask':''}
        utils.cleanmaskclean(parms)

        # Get the flux and rms from this image, this is needed to compare different 
        # cycles. "< 1" is to invert the mask, so that the rms is caculated
        # only for the background. The syntax is very weird, see:
        # https://casa.nrao.edu/aips2_docs/notes/223/index.shtml
        mask = '"'+current_peel_directory+'/img/cycle'+str(i)+'/im.newmask'+'"'+' < 1'
        rms = imstat(imagename=current_peel_directory+'/img/cycle'+str(i)+'/im-masked.image.tt0',
                     mask=mask)['rms'][0]
        fit = imfit(imagename=current_peel_directory+'/img/cycle'+str(i)+'/im-masked.image.tt0')
        print 'I ASSUMED THAT COMPONENT0 IS OUR PEEL SOURCE (BRIGHTEST)'
        print fit
        tflux = fit['results']['component0']['flux']['value'][0] # Jy, assume component0 is the brightest

        # Do a check and update solint
        if rms < 1.1*prev_rms:# or i < 2:
            print '---'
            print '---'
        
            noise_per_interval = rms * (number_of_integrations)**0.5
            snr_per_interval = tflux/noise_per_interval
            time_steps_needed = ((min_SNR_per_solution_interval/fudge_factor)/snr_per_interval)**2.0 # SNR grows with sqrt(time)
            time_needed_in_sec = integration_time*time_steps_needed
            solint = '%.1fs' % time_needed_in_sec
        
            print solint
            print '---'
            print '---'
        
            # Set the new model into the MODEL_DATA column (only select region)
            new_model = [current_peel_directory+'/img/cycle'+str(i)+'/im-masked.model.tt0', 
                             current_peel_directory+'/img/cycle'+str(i)+'/im-masked.model.tt1']
            updated_peel_source_model = extrModel(new_model, region, i+1)
            ft(vis=peelMS.file_path, model=updated_peel_source_model, nterms=len(updated_peel_source_model),
               usescratch=True)

            usedcycles = i
            prev_rms = rms
        else:
            # No improvement and no need to proceed with further cycles
            break
           

    # Now, after enough self-cal cycles, we have a solution table and model for
    # this peel source. Subtract this new model in order to get an updated
    # residual measurement set.
    if usedcycles == -1:
        # No improvent at all
        best_updated_model = peel_source_model
    else:
        best_updated_model = [current_peel_directory+'/img/cycle'+str(usedcycles)+'/im-masked.model.tt0', 
                              current_peel_directory+'/img/cycle'+str(usedcycles)+'/im-masked.model.tt1']
    delmod(vis=peelMS.file_path)
    subtract(peelMS.file_path, best_updated_model)
    
    # Split off the data set in order to make sure that the DATA column
    # only holds the residual map
    updatedResMSfilepath = current_peel_directory+'/updated_residual.ms'
    split(peelMS.file_path, updatedResMSfilepath)
    updatedResidualMS = TableObjects.MSObj(updatedResMSfilepath)
    
    if not usedcycles == -1:
        # Apply the inverse phase solutions for this particular peel source
        allTimestamps = updatedResidualMS.get_all_available_timestamps()
        Gp = TableObjects.STObj(current_peel_directory+'/cal/cycle'+str(usedcycles)+'.Gp')
        print '-- refant --'
        print refAnt
        Gp_reref_path = Gp.re_reference_table(refant=1)
        Gp_reref = TableObjects.STObj(Gp_reref_path)
        Gp_normref_path = Gp_reref.normalize_reference_antenna()
        Gp_normref = TableObjects.STObj(Gp_normref_path)
        Gp_resamp_path = Gp_normref.resample_solutions(allTimestamps, interp_type = 'linear')
        Gp_resamp = TableObjects.STObj(Gp_resamp_path)
        Gp_resamp.plot(current_peel_directory+'/plot/cycle'+str(usedcycles)+'/resamp', phase_only=True)
        Gp_inv_path = Gp_resamp.invert_table()
        applycal(vis=updatedResidualMS.file_path, gaintable=Gp_inv_path, calwt=False, 
                 flagbackup=False)

    # Set the phase center again to the initial position
    fixvis(vis=updatedResidualMS.file_path, outputvis=updatedResidualMS.file_path,
           phasecenter = orig_center_string)
    
    # CHECK IF EVERYTHING WORKS
    parms = {'vis':updatedResidualMS.file_path, 'imagename':current_peel_directory+'/updated_res',
             'gridmode':'widefield', 'wprojplanes':512, 'mode':'mfs', 
             'nterms':2, 'niter':10000, 'gain':0.1, 'psfmode':'clark', 
             'imagermode':'csclean', 'imsize':[5000], 'cell':cell, 
             'weighting':'briggs', 'robust':rob, 'usescratch':True, 
             'mask':''}
    utils.cleanmaskclean(parms, makemask = False)
    exportfits(imagename=current_peel_directory+'/updated_res.image.tt0',
               fitsimage=current_peel_directory+'/updated_res.image.tt0.fits',
               history=False, overwrite=True)
               
    residualrms = imstat(imagename=current_peel_directory+'/updated_res.image.tt0')['rms'][0]
    print '--'
    print '--'
    print '--'
    print 'residual rms: ', residualrms
    print '--'
    print '--'
    print '--'
    
    return updatedResidualMS
