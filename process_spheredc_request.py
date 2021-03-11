#!/usr/bin/env python3.6

import argparse
# from glob import glob
import os
import numpy as np
import subprocess
import time
from pathlib import Path
from astropy.io import fits
import matplotlib.pyplot as plt

"""
This script needs to be updated following some changes done at the Data Center in the structure of the 
download script.
A new version of a similar script exists at
/Users/millij/Documents/DataCenter/Beltramo_project/process_spheredc_request.py

"""

def create_new_script(files,\
                      text_in_lines_to_comment=['.json','reduced_image.png',\
                                                'unsharp13.fits',\
                                                'TDB.tar.gz',\
                                                '/reduced_image_median.fits',\
                                                'MISC']): #'allcomb.tar.gz'
    """
    This function takes in input a download script file (sphere_dl_scipt.sh) or 
    a list of files (sphere_dl*.sh) and creates a copy of this script after commenting
    some irrelevant wget lines (to avoid downloading useless data). 

    It also returns the interesting image files and tar.gz files and returns them 
    in the form of 2 lists    
    """

    # the line below specifies which type of IRDIS images are returned by the function
    # (to be later opened in DS9)
    images_to_return_IRDIS = ['TLOCI','PCA','cADI','no_ADI']
    # the line below specifies which type of IFS images are returned by the function
    # (to be later opened in DS9). We only put here TLOCI and PCA as IFS is not
    # really relevant for CADI and no-ADI
    images_to_return_IFS = ['TLOCI','PCA']

    list_images_filenames = []
    list_text_filenames = []
    list_modified_scripts = []
    list_specal_charac_tar_files = []
    list_contrast_curves_tar_files = []
    list_allcomb_tar_files = []
    
    for f in files:
        h = open(f)
        f2 = f.replace('.sh','_commented.sh')
        list_modified_scripts.append(f2)
        if os.path.exists(f2):
            os.remove(f2)    
        h2 = open(f2, 'w')
        for l in h.readlines():
            if 'wget' in l:

                # 1st check if this is actually not a file you don't want to download
                text_is_there = np.asarray([t in l for t in text_in_lines_to_comment],dtype=bool)
                comment_line = np.any(text_is_there)

                # 2nd: check that the file actually does not already exist
                start_index = l.index('SPHERE_DC_DATA')
                end_index = l.index('" "http://')         
                filename = l[start_index:end_index]

                if comment_line:
                    l = l.replace('wget', '#wget')
                elif os.path.isfile(filename):
                    comment_line = True
                    print(filename,'already downloaded. Nothing to be done')
                    l = l.replace('wget', '#wget')

                algo_type_list_IRDIS = np.array([algo in l for algo in images_to_return_IRDIS],dtype=bool)
                algo_type_list_IFS   = np.array([algo in l for algo in images_to_return_IFS],dtype=bool)
                
                # for IRDIS level 1 images:                    
                conditions_IRDIS_level1 = ('IRD_SCIENCE_REDUCED_MASTER_CUBE-center_im.fits' in l or \
                                           'IRD_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits' in l)
                # for IRDIS level 2 images:
                condition_IRDIS_level2 = ('ird_specal_dc-IRD_SPECAL_MEDIAN_RESIDUAL_CUBE-cube_reduced_image_median.fits' in l and \
                    np.any(algo_type_list_IRDIS))
                    
                # for IFS level 1 images:
                conditions_IFS_level1 = ('IFS_SCIENCE_REDUCED_SPECTRAL_MASTER_CUBE-center_im.fits' in l or \
                    'IFS_SCIENCE_PSF_MASTER_CUBE-median_unsat.fits' in l)                                
                    
                # for IFS level 2 images:                    
                conditions_IFS_level2 = ('ifs_specal_dc-IFS_SPECAL_MEDIAN_RESIDUAL_STACK-reduced_image_median.fits' in l and \
                    np.any(algo_type_list_IFS))
                    
                # for text files:
                condition_text_files = ('convert_log.txt' in l or 'ATMO_CONDITIONS-seeing.txt' in l)
                    
                if condition_IRDIS_level2 or conditions_IFS_level1 or conditions_IRDIS_level1 or conditions_IFS_level2:
                    list_images_filenames.append(filename)      
                if condition_text_files:
                    list_text_filenames.append(filename)
                if 'specalcarac_outputs.tar.gz' in l:
                    list_specal_charac_tar_files.append(filename)      
                if 'contrast_curves.tar.gz' in l:
                    list_contrast_curves_tar_files.append(filename)                                              
                if 'allcomb.tar.gz' in l:
                    list_allcomb_tar_files.append(filename)
            h2.write(l)
        h.close()
        h2.close()
        print('{0:s} -> {1:s}'.format(f, f2))    
    return list_images_filenames,list_specal_charac_tar_files,list_contrast_curves_tar_files,list_allcomb_tar_files,list_text_filenames,list_modified_scripts

if __name__ == "__main__":
    # We parse the input
    parser = argparse.ArgumentParser(description='Type the file or files you want to process')

    parser.add_argument('files', type=str, help='file(s) you want to process', nargs='*')
    parser.add_argument('-d', help='download data', action='store_true') #if used True (we download) else False: no download
    parser.add_argument('-untar', help='untar archives and open', action='store_true')
    parser.add_argument('-ds9', help='visualize data in DS9', action='store_true')

    args = parser.parse_args()
    list_images_filenames,list_specal_charac_tar_files,list_contrast_curves_tar_files,list_allcomb_tar_files,list_text_filenames,list_modified_scripts = create_new_script(args.files)
    download = args.d
    untar = args.untar
    ds9=args.ds9
    
    if download:
        # we execute the file
        for script in list_modified_scripts:
            cmd = 'sh {0:s}'.format(script) 
            print(cmd)
            process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()

    if ds9:
        # we load the images in DS9
        if len(list_images_filenames)>0:
            print('Opening the files:')
            print('\n'.join(list_images_filenames))
            # #Method 1 with VIP
            # from astropy.io import fits
            # array_images = []
            # array_images = []
            # for img in list_images_filenames:
            #     array_images.append(fits.getdata(img))
            # import vip_hci as vip
            # ds9 = vip.ds9=vip.Ds9Window()
            # ds9.display(*array_images)

            # Method 2
            for f in list_images_filenames:
                process = subprocess.Popen(['open',f], stdout=subprocess.PIPE)
                output, error = process.communicate()
                time.sleep(2)
                
    if untar:
        if len(list_specal_charac_tar_files)>0:
            for f in list_specal_charac_tar_files:
                print('Untar',f)
                dirname = os.path.dirname(f)
                filename = os.path.basename(f)                
                process = subprocess.Popen(['tar','xvfz',f,'-C',dirname], stdout=subprocess.PIPE)
                output, error = process.communicate()

                file_to_open = os.path.join(dirname,'Reduction_0000/PointSource_astro_photo/comp_astro_photo_tdb.csv')
                print('Open',file_to_open)
                process = subprocess.Popen(['open',file_to_open], stdout=subprocess.PIPE)
                output, error = process.communicate()

        if len(list_contrast_curves_tar_files)>0:
            for icontrastfile,f in enumerate(list_contrast_curves_tar_files):
                print('Untar',f)
                dirname = os.path.dirname(f)
                filename = os.path.basename(f)                
                process = subprocess.Popen(['tar','xvfz',f,'-C',dirname], stdout=subprocess.PIPE)
                output, error = process.communicate()
                if 'ifs' in dirname:
                    ins='IFS'
                elif 'ird' in dirname:
                    ins='IRD'
                else:
                    print('Instrument not recognized')
                path_results = Path(dirname)
                list_contrast_curves = list(path_results.glob('*/*/*contrast_curve_tab.fits'))
                print('Displaying\n',['   - '+str(c) for c in list_contrast_curves])
                if len(list_contrast_curves)==1:
                    for icontrast_curve,contrast_curve in enumerate(list_contrast_curves):
                        hdu_list = fits.open(contrast_curve)
                        table = hdu_list[1].data
                        ncurves,nsep = table['NSIGMA_CONTRAST'].shape
                        plt.figure(10*icontrastfile+icontrast_curve)
                        if len(np.unique(table['LAM']))>1:
                            if ins == 'IFS':
                                generator_lambda = range(3,ncurves,3)
                            else:
                                generator_lambda = range(0,ncurves)                                
                            for k in generator_lambda:
                                plt.semilogy(table['SEPARATION'][k,:],table['NSIGMA_CONTRAST'][k,:],label=table['LAM'][k])
                        else:
                            plt.semilogy(table['SEPARATION'][0,:],table['NSIGMA_CONTRAST'][0,:],label=table['LAM'][0])                            
                        plt.legend(frameon=False,loc='best')
                        plt.xlabel('Separation in "')
                        plt.ylabel('$5 \\sigma$ contrast')
                        plt.title(contrast_curve.parent.parent)
                        plt.ylim(1e-7,1e-3)
                        if ins == 'IFS':
                            plt.xlim(0,2)
                        else:
                            plt.xlim(0,8)
                        plt.grid()
                        plt.savefig(contrast_curve.parent.joinpath('contrast_curve_tab.pdf'))
                        hdu_list.close()
                elif len(list_contrast_curves)>1:
                    nb_contrast_curves = len(list_contrast_curves)
                    sep = []
                    con = []
                    lab = []
                    for icontrast_curve,contrast_curve in enumerate(list_contrast_curves):
                        hdu_list = fits.open(contrast_curve)
                        table = hdu_list[1].data
                        ncurves,nsep = table['NSIGMA_CONTRAST'].shape
                        if len(np.unique(table['LAM']))>1:
                            print('{0:d} contrast curves available, we only plot one'.format(ncurves))
                        sep.append(table['SEPARATION'][-1,:])
                        con.append(table['NSIGMA_CONTRAST'][-1,:])
                        lab.append(contrast_curve.parent.parent.name)
                    plt.figure(0)
                    for k in range(nb_contrast_curves):
                        plt.semilogy(sep[k],con[k],label=lab[k])                        
                    plt.legend(frameon=False,loc='best')
                    plt.xlabel('Separation in "')
                    plt.ylabel('$5 \\sigma$ contrast')
                    plt.title(path_results.name)
                    plt.ylim(1e-7,1e-3)
                    if ins == 'IFS':
                        plt.xlim(0,2)
                    else:
                        plt.xlim(0,8)
                    plt.grid()
                    plt.savefig(path_results.joinpath('contrast_curve_tab.pdf'))
                    hdu_list.close()
                else:
                    print('No contrast curve to plot')

        plt.show()

        if len(list_allcomb_tar_files)>0:
            for iallcomb,f in enumerate(list_allcomb_tar_files):
                print('Untar',f)
                dirname = os.path.dirname(f)
                filename = os.path.basename(f)                
                process = subprocess.Popen(['tar','xvfz',f,'-C',dirname], stdout=subprocess.PIPE)
                output, error = process.communicate()


    if len(list_text_filenames)>0:
        for f in list_text_filenames:
            print('Opening {0:s}'.format(f))
            h3 = open(f)
            # print(h3.readlines())
            for l in h3.readlines():
                print(l[:-2])
            # print('\n')
