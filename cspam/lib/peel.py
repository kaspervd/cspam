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

sou_res = ['1arcsec']
sou_size = [5000]
expnoise = 1.e-6
rob=0.5

def extrModel(modelimg, region, compl=False, extend=None):
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
            immath(imagename = modelimgtt, mode = 'evalexpr', expr = 'IM0',
                   region = region, outfile = region.replace('.crtf','')+'_model.tt'+str(i))

            blankedmodelimg.append(region.replace('.crtf','')+"_model.tt"+str(i))

    return blankedmodelimg

def subtract(active_ms, modelimg, wprojplanes=0):
    """General function to call the necessary steps to subtract point sources
    the modelimg must have only point sources one wants to keep into the field.
    second_ms_path: MS with calibrated data in DATA
    modelimg: model of the whole sky (array of tt)
    wprojplanes: number of w-projection planes, if 0 a direct ft() will be used (best for small field)
    """
    ftw(vis=active_ms, model=modelimg, nterms=len(modelimg), wprojplanes=wprojplanes, usescratch=True)
    uvsub(vis=active_ms)

def invertTable(caltab):
    """Invert a calibration table
    """
    syscommand = "cp -r "+caltab+" "+caltab+"_inv"
    os.system(syscommand)
    caltab = caltab+"_inv"
    tb.open(caltab, nomodify=False) # open the caltable
    gVals = tb.getcol('CPARAM')#, startrow=start, nrow=incr) # get the values from the GAIN column
    mask = abs(gVals) > 0.0 # only consider non-zero values
    gVals[mask] = 1.0 / gVals[mask] # do the inversion
    tb.putcol('CPARAM', gVals)#, startrow=start, nrow=incr) # replace the GAIN values with the inverted values
    tb.close() # close the table
    return caltab

def peel(residualMS, modelimg, current_peel_directory, solint):

    if not os.path.isdir(current_peel_directory+'/cal'):
        os.makedirs(current_peel_directory+'/cal')
    if not os.path.isdir(current_peel_directory+'/img'):
        os.makedirs(current_peel_directory+'/img')
    if not os.path.isdir(current_peel_directory+'/plot'):
        os.makedirs(current_peel_directory+'/plot')

    region = current_peel_directory+'/region.crtf'

    # Obtain directions for later
    orig_center = residualMS.get_direction_from_tgt_field_id(0)
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
    peel_source_model = extrModel(modelimg, region, compl=False)

    # This cutout is smaller than the original image and placed outside the
    # center. This means that a phase shift has to be applied.
    # Shift the phase center to the peel source position
    fixvis(vis=residualMS.file_path, outputvis=residualMS.file_path,
           phasecenter = peel_source_center_string)

    # Add this source to the residual image
    ftw(vis=residualMS.file_path, model=peel_source_model, 
        nterms=len(peel_source_model), wprojplanes=512, usescratch=True)
    uvsub(vis=residualMS.file_path, reverse=True)

    # Create a new measurement set with now only the residual image + peel source
    peel_source_ms_path = current_peel_directory+'/'+residualMS.ms_name+'_peel.ms'
    split(vis=residualMS.file_path, outputvis=peel_source_ms_path)
    peelsourceMS = TableObjects.MSObj(peel_source_ms_path)



    clean(vis=peelsourceMS.file_path, imagename=peelsourceMS.file_path+'_test1', 
          gridmode='widefield', wprojplanes=512, mode='mfs',
          niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', 
          interactive=False, imsize=sou_size, cell=sou_res,
          stokes='I', nterms=2, weighting='briggs', robust=rob,
          usescratch=True)
    exportfits(imagename=peelsourceMS.file_path+'_test1.image.tt0',
               fitsimage=peelsourceMS.file_path+'_test1.image.tt0.fits',
               history=False, overwrite=True)
               


    # Use the new image (with only peel source) as a model for self calibration
    ft(vis=peelsourceMS.file_path, model=peel_source_model,
       nterms=len(peel_source_model), usescratch=True)
  
    # Determine the best reference antenna
    refAntObj = AntennaObjects.RefAntHeuristics(vis=peelsourceMS.file_path,        
                               geometry=True, flagging=True)                              
    refAnt = refAntObj.calculate()[0] 

    # In order for plotcal to work a symbolic link is needed.
    # Plotcal assumes the measurement set is in the same directory
    # as the cal table.
    syscommand = 'ln -s '+peelsourceMS.file_path+' '+current_peel_directory+\
            '/cal/'+peelsourceMS.ms_name+'_peel.ms'
    os.system(syscommand)

    # selfcal cycle 1                                                          
    gaincal(vis=peelsourceMS.file_path, caltable=current_peel_directory+\
            '/cal/peel1.Gp', solint=solint, refant=refAnt, minsnr=1,
            minblperant=peelsourceMS.minBL_for_cal, calmode='p', uvrange='>50m')
    Gp = TableObjects.STObj(current_peel_directory+'/cal/peel1.Gp')
    Gp.plot(current_peel_directory+'/plot/peel1', phase_only=True)

    applycal(vis=peelsourceMS.file_path, gaintable=Gp.file_path, 
             calwt=False, flagbackup=False)





    clean(vis=peelsourceMS.file_path, imagename=peelsourceMS.file_path+'_test2', 
          gridmode='widefield', wprojplanes=512, mode='mfs',
          niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', 
          interactive=False, imsize=sou_size, cell=sou_res,
          stokes='I', nterms=2, weighting='briggs', robust=rob,
          usescratch=True)
    exportfits(imagename=peelsourceMS.file_path+'_test2.image.tt0',
               fitsimage=peelsourceMS.file_path+'_test2.image.tt0.fits',
               history=False, overwrite=True)






    sys.exit()
    # Cutout the source to peel
    rest_of_field = extrModel(modelimg, region, compl=True, extend=[str(ra_p)+'rad',str(dec_p)+'rad'])

    # Convert this new model (without peel source) into visibilities and
    # and instert into the MODEL_DATA column.
    # Next, subtract these new model visibility data (i.e. MODEL_DATA) from
    # the corrected visibility data (i.e. CORRECTED_DATA or simply the
    # antenna data with gains applied).
    subtract(peelmset.file_path, rest_of_field, wprojplanes=512)
    
    # Make a measurement set with now only the peel source
    peelmsetfilepath2 = peelmset.file_path+'_peelsource'
    split(vis=peelmset.file_path, outputvis=peelmsetfilepath2)
    
    # Cutout everything except the peel source from the model image
    peel_source = extrModel(modelimg, region, compl=False)

    # This cutout is smaller than the original image and placed outside the
    # center. This means that a phase shift has to be applied.
    # Shift the phase center to the peel source position
    fixvis(vis=peelmsetfilepath2, outputvis=peelmsetfilepath2,
           phasecenter = peel_source_center_string)

    # Use the new image (with only peel source) as a model for self calibration
    ft(vis=peelmsetfilepath2, model=peel_source, nterms=len(peel_source), usescratch=True)
  
    # Determine the best reference antenna
    refAntObj = AntennaObjects.RefAntHeuristics(vis=peelmsetfilepath2,        
                               geometry=True, flagging=True)                              
    refAnt = refAntObj.calculate()[0] 

    # In order for plotcal to work a symbolic link is needed.
    # Plotcal assumes the measurement set is in the same directory
    # as the cal table.
    syscommand = 'ln -s '+peelmsetfilepath2+' '+current_peel_directory+\
            '/cal/'+peelmset.ms_name+'.ms_peelsource'
    os.system(syscommand)

    # selfcal cycle 1                                                          
    gaincal(vis=peelmsetfilepath2, caltable=current_peel_directory+\
            '/cal/peel1.Gp',
            solint='int', refant=refAnt, minsnr=1, minblperant=peelmset.minBL_for_cal,
            calmode='p', uvrange='>50m')
    Gp = TableObjects.STObj(current_peel_directory+'/cal/peel1.Gp')
    Gp.plot(current_peel_directory+'/plot/peel1', phase_only=True)
    gaincal(vis=peelmsetfilepath2, caltable=current_peel_directory+\
            '/cal/peel1.Ga',
            solint='int', refant=refAnt, minsnr=1, minblperant=peelmset.minBL_for_cal,
            calmode='a', uvrange='>50m')
    Ga = TableObjects.STObj(current_peel_directory+'/cal/peel1.Ga')
    Ga.plot(current_peel_directory+'/plot/peel1', amp_only=True)

    applycal(vis=peelmsetfilepath2, gaintable=[Gp.file_path, Ga.file_path], 
             calwt=False, flagbackup=False)

    # selfcal cycle 2
    gaincal(vis=peelmsetfilepath2, caltable=current_peel_directory+\
            '/cal/peel2.Gp',
            solint='int', refant=refAnt, minsnr=1, minblperant=peelmset.minBL_for_cal,
            calmode='p', uvrange='>50m')
    Gp2 = TableObjects.STObj(current_peel_directory+'/cal/peel2.Gp')
    Gp2.plot(current_peel_directory+'/plot/peel2', phase_only=True)
    gaincal(vis=peelmsetfilepath2, caltable=current_peel_directory+\
            '/cal/peel2.Ga',
            solint='int', refant=refAnt, minsnr=1, minblperant=peelmset.minBL_for_cal,
            calmode='a', uvrange='>50m')
    Ga2 = TableObjects.STObj(current_peel_directory+'/cal/peel2.Ga')
    Ga2.plot(current_peel_directory+'/plot/peel2', amp_only=True)

    applycal(vis=peelmsetfilepath2, gaintable=[Gp2.file_path, Ga2.file_path], 
             calwt=False, flagbackup=False)

    clean(vis=peelmsetfilepath2, imagename=current_peel_directory+\
          '/img/peelsource', 
          gridmode='widefield', wprojplanes=512, mode='mfs',
          niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', 
          interactive=False, imsize=sou_size, cell=sou_res,
          stokes='I', nterms=2, weighting='briggs', robust=rob,
          usescratch=True, mask=region)
    exportfits(imagename=current_peel_directory+'/img/peelsource.image.tt0',
               fitsimage=current_peel_directory+'/img/peelsource.image.tt0.fits',
               history=False, overwrite=True)

    # IDEA: APPLY SOLUTIONS, DELETE PEEL SOURCE IN FOURIER SPACE, APPLY 
    # INVERSE SOLUTIONS

    # Now that we have the solutions, apply these solutions on the original
    # image. For this we first need to subtract the rest_of_field visibilities.
    subtract(peelmset.file_path+'_orig', rest_of_field, wprojplanes=512)

    # Split into a new measurement set
    split(vis=peelmset.file_path+'_orig', outputvis=peelmset.file_path+'_peelsource_calibrated')

    # Apply the solutions obtained earlier
    applycal(vis=peelmset.file_path+'_peelsource_calibrated', gaintable=[Gp2.file_path, Ga2.file_path], calwt=False, flagbackup=False)

    # Create an image for inspection
    clean(vis=peelmset.file_path+'_peelsource_calibrated', imagename=current_peel_directory+\
          '/img/peelsource_calibrated', 
          gridmode='widefield', wprojplanes=512, mode='mfs',
          niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', 
          interactive=False, imsize=sou_size, cell=sou_res,
          stokes='I', nterms=2, weighting='briggs', robust=rob,
          usescratch=True)
    exportfits(imagename=current_peel_directory+'/img/peelsource_calibrated.image.tt0',
               fitsimage=current_peel_directory+'/img/peelsource_calibrated.image.tt0.fits',
               history=False, overwrite=True)

    # SHOULD I DO THIS?
    #invcaltaba = invertTable(Ga2.file_path)
    #Ga2inv = TableObjects.STObj(invcaltaba)
    #Ga2inv.plot(current_peel_directory+'/plot/inv',amp_only=True)
    
    #invcaltabp = invertTable(Gp2.file_path)
    #Gp2inv = TableObjects.STObj(invcaltabp)
    #Gp2inv.plot(current_peel_directory+'/plot/inv',phase_only=True)

    # Split off a new measurement set in which the other sources will be added
    # back in.
    split(vis=peelmset.file_path+'_peelsource_calibrated', outputvis=peelmset.file_path+'_source_added_back')
    
    # SHOULD I DO THIS?
    #applycal(vis=peelmset.file_path+'_source_added_back', gaintable=[Ga2inv.file_path,Gp2inv.file_path], calwt=False, flagbackup=False)

    # Add in the other sources
    ftw(vis=peelmset.file_path+'_source_added_back', model=rest_of_field, nterms=len(rest_of_field), wprojplanes=512, usescratch=True)
    uvsub(vis=peelmset.file_path+'_source_added_back', reverse=True)

    # Create an image for inspection
    clean(vis=peelmset.file_path+'_source_added_back', imagename=current_peel_directory+\
          '/img/source_added_back', 
          gridmode='widefield', wprojplanes=512, mode='mfs',
          niter=5000, gain=0.1, psfmode='clark', imagermode='csclean', 
          interactive=False, imsize=sou_size, cell=sou_res,
          stokes='I', nterms=2, weighting='briggs', robust=rob,
          usescratch=True)
    exportfits(imagename=current_peel_directory+'/img/source_added_back.image.tt0',
               fitsimage=current_peel_directory+'/img/source_added_back.image.tt0.fits',
               history=False, overwrite=True)
