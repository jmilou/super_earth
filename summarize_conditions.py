#!/usr/bin/env python3.6

# import argparse
from glob import glob
import os
import numpy as np
# import subprocess
# import time
from astropy.io import fits
import pandas as pd

def return_info_level1(file):
    """
    """
    h = open(file)
    for l in h.readlines():
        if 'mkdir' in l:
            directory  = l[l.index('SPHERE_DC_DATA/'):-2]
            process = int(directory[-7:-1])
            tmp = directory[directory.index('/')+1:]
            target_name_DC = tmp[0:tmp.index('_')]
        if "MASTER_CUBE-center_im.fits" in l:
            filename = get_filename_from_line(l)
            header = fits.getheader(filename)
            date_start = header['HIERARCH ESO OBS START']
            date_end = header['DATE-OBS']
            ob_id = int(header['HIERARCH ESO OBS ID'])
            dit = np.round(float(header['HIERARCH ESO DET SEQ1 DIT']),2)
            irdis_filter = header['HIERARCH ESO INS1 OPTI2 NAME']
            ifs_filter = header['HIERARCH ESO INS2 OPTI2 NAME']
            dimm_seeing_start = np.round(float(header['HIERARCH ESO TEL AMBI FWHM START']),2)
            dimm_seeing_end = np.round(float(header['HIERARCH ESO TEL AMBI FWHM END']),2)
            nframes = int(header['NAXIS3'])
            prog_id = header['HIERARCH ESO OBS PROG ID']
        if "PARA_ROTATION_CUBE-rotnth.fits" in l:
            filename = get_filename_from_line(l)
            parang_list = fits.getdata(filename)
            delta_parang = np.round(np.max(parang_list) - np.min(parang_list),2)
        if 'ATMO_CONDITIONS-seeing.txt' in l:
            filename = get_filename_from_line(l)
            atmf = open(filename)
            line = atmf.readline()
            while not line.startswith('-'):
                line = atmf.readline()
            line = atmf.readline()
            seeing,tau0,airmass=[float(a) for a in line.split()] 
            seeing = np.round(seeing,2)
            airmass = np.round(airmass,2)
            tau0 = np.round(tau0*1000,2)
            atmf.close()
    h.close()
    dico = {'process':process,'target name DC':target_name_DC,'directory':directory,\
            'date start':date_start,'date end':date_end,'ob id':ob_id,'DIT [s]':dit,\
            'IRDIS filter':irdis_filter,'IFS disperser':ifs_filter,'MASS DIMM seeing start ["]':dimm_seeing_start,\
            'MASS DIMM seeing end ["]':dimm_seeing_end,'nb frames':nframes,'delta parang':delta_parang,\
            'MASS DIMM seeing avg ["]':seeing,'MASS DIMM tau0 avg [ms]':tau0,'airmass avg':airmass,\
            'program ID':prog_id}
    return dico

def return_info_irdis_level1(file):
    """
    """
    h = open(file)
    for l in h.readlines():
        if 'PARAMETER_LOG-convert_log.txt' in l:
            filename = get_filename_from_line(l)
            log = open(filename)
            for line in log.readlines() :
                if 'plate scale calibration' in line:
                    psc_date = line[line.index('is :')+4:].split()[0]
                if 'The pixscale' in line:
                    px = float(line[line.index('is :')+4:])
                if 'Date of true North calibration' in line:
                    TN_calibration_date = line[line.index('is :')+4:].split()[0]
                if 'TN correction' in line:
                    TN_correction = float(line[line.index('is :')+4:].split()[0])
                if 'Error on TN' in line:
                    TN_error = float(line[line.index('is :')+4:].split()[0])
                    if TN_error==0:
                        TN_error = np.nan
                if 'photometric variation' in line:
                    phot_variation_percentage = np.round(float(line[line.index('[%]')+4:]),2)
                if 'DIT (s)' in line:                    
                    DIT_log = float(line[line.index(': ')+1:])
                if 'MEAN Empirical SR / RTC SR' in line:
                    strehl_measured,strehl_RTC = [np.round(float(sr),2) for sr in line[line.index(':')+4:].split()]
                if 'lambda:' in line:
                    wavelengths = [float(v) for v in line[line.index('lambda:')+7:].split()]
            log.close()                                     
    h.close()
    dico = {'platescale':px,'platescale calibration date':psc_date,\
            'true north calibration date':TN_calibration_date,'true north':TN_correction,\
            'true north error':TN_error,'wavelength channels [microns]':wavelengths,\
            'photometric variations [%]':phot_variation_percentage,'DIT PSF [s]':DIT_log,\
            'Strehl (measured)':strehl_measured,'Strehl (RTC)':strehl_RTC,'camera':'IRDIS'}
    return dico

def return_info_ifs_level1(file):
    """
    """
    h = open(file)
    for l in h.readlines():
        if 'PARAMETER_LOG-convert_log.txt' in l:
            filename = get_filename_from_line(l)
            log = open(filename)
            for line in log.readlines() :
                if 'TN correction' in line:
                    TN_correction = float(line[line.index('is :')+4:].split()[0])
                if 'pixscale' in line:
                    px = float(line[line.index('is :')+4:])
                if 'filter corono:' in line:
                    wavelengths = [float(v) for v in line[line.index('lambda:')+7:].split()]
                if 'photometric variation' in line:
                    phot_variation_percentage = np.round(float(line[line.index('[%]')+4:]),2)
            log.close()                         

            log = open(filename)
            line = log.readline()
            while not line.startswith('filter corono:'):
                line = log.readline()
            wavelengths = [float(v) for v in line[line.index('lambda:')+7:].split()]
            line = log.readline()
            while len(wavelengths)<39:
                wavelengths = np.concatenate((wavelengths,[float(v) for v in line.split()]))
                line = log.readline()
            log.close()                         

    h.close()
    dico = {'platescale':px,'platescale calibration date':'NA',\
            'true north calibration date':'NA','true north':TN_correction,\
            'true north error':np.nan,'wavelength channels [microns]':wavelengths,\
            'photometric variations [%]':phot_variation_percentage,'DIT PSF [s]':np.nan,\
            'Strehl (measured)':np.nan,'Strehl (RTC)':np.nan,'camera':'IFS'}
    # dico = {'platescale':px,\
    #         'true north':TN_correction,\
    #         'wavelength channels [microns]':wavelengths,\
    #         'photometric variations [%]':phot_variation_percentage,'camera':'IFS'}
    return dico

def get_filename_from_line(line):
    """
    Return the filename (relative path) corresponding to a given lin of the download script
    """
    start_index = line.index('SPHERE_DC_DATA')
    end_index = line.index('" "http://')         
    filename = line[start_index:end_index]
    return filename

if __name__ == "__main__":

    files_IRDIS = glob('*IRDIS_level1.sh')    
    nb_obs_irdis = len(files_IRDIS)

    files_IFS = [file_IRDIS.replace('IRDIS','IFS') for file_IRDIS in files_IRDIS]      
    files_IRDIS_level2 = [file_IRDIS.replace('level1','level2') for file_IRDIS in files_IRDIS]      
    files_IFS_level2 = [file_IFS.replace('level1','level2') for file_IFS in files_IFS]      

    for i,file_IFS in enumerate(files_IFS):
        if not os.path.exists(file_IFS):
            print(file_IFS,'does not exist')

    for i,file_IRDIS_level2 in enumerate(files_IRDIS_level2):
        if not os.path.exists(file_IRDIS_level2):
            print(file_IRDIS_level2,'does not exist')
                
    for i,file_IFS_level2 in enumerate(files_IFS_level2):
        if not os.path.exists(file_IFS_level2):
            print(file_IFS_level2,'does not exist')


    dico = {'target name':[],\
            'target name DC':[],\
            'camera':[],\
            'program ID':[],\
            'date start':[],\
            'date end':[],\
            'ob id':np.ndarray((nb_obs_irdis*2),dtype=int),\
            'process':np.ndarray((nb_obs_irdis*2),dtype=int),\
            'directory':[],\
            'DIT [s]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'IRDIS filter':[],\
            'IFS disperser':[],\
            'MASS DIMM seeing start ["]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'MASS DIMM seeing end ["]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'nb frames':np.ndarray((nb_obs_irdis*2),dtype=int),\
            'delta parang':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'MASS DIMM seeing avg ["]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'MASS DIMM tau0 avg [ms]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'airmass avg':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'platescale':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'platescale calibration date':[],\
            'true north':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'true north calibration date':[],\
            'true north error':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'wavelength channels [microns]':[],\
            'photometric variations [%]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'DIT PSF [s]':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'Strehl (measured)':np.ndarray((nb_obs_irdis*2),dtype=float),\
            'Strehl (RTC)':np.ndarray((nb_obs_irdis*2),dtype=float)}

    for i,file_IRDIS_level1 in enumerate(files_IRDIS):
        
        # file_IRDIS_level2 = files_IRDIS_level2[i]
        # file_IFS_level2 = files_IFS_level2[i]

        dirname_IRDIS_level1 = os.path.dirname(file_IRDIS_level1)
        filename_IRDIS_level1 = os.path.basename(file_IRDIS_level1)                

        name_tmp = file_IRDIS_level1[file_IRDIS_level1.index('script_')+7:]
        target_name = name_tmp[0:name_tmp.index('_')]

        print('processing',file_IRDIS_level1)
        dico_common = return_info_level1(file_IRDIS_level1)
        dico_common['target name'] = target_name
        for k, v in dico_common.items():
            if type(v) not in (float,int,np.float,np.float32,np.float64):
                dico[k].append(v)
            else:
                dico[k][i] = v

        dico_irdis = return_info_irdis_level1(file_IRDIS_level1)
        for k, v in dico_irdis.items():
            if type(v) not in (float,int,np.float,np.float32,np.float64):
                dico[k].append(v)
            else:
                dico[k][i] = v

    for i,file_IFS_level1 in enumerate(files_IFS):
        
        # file_IRDIS_level2 = files_IRDIS_level2[i]
        # file_IFS_level2 = files_IFS_level2[i]

        dirname_IFS_level1 = os.path.dirname(file_IFS_level1)
        filename_IFS_level1 = os.path.basename(file_IFS_level1)                

        name_tmp = file_IFS_level1[file_IFS_level1.index('script_')+7:]
        target_name = name_tmp[0:name_tmp.index('_')]

        print('processing',file_IFS_level1)
        dico_common = return_info_level1(file_IFS_level1)
        dico_common['target name'] = target_name
        for k, v in dico_common.items():
            if type(v) not in (float,int,np.float,np.float32,np.float64):
                dico[k].append(v)
            else:
                dico[k][nb_obs_irdis+i] = v

        dico_IFS = return_info_ifs_level1(file_IFS_level1)
        for k, v in dico_IFS.items():
            if type(v) not in (float,int,np.float,np.float32,np.float64):
                dico[k].append(v)
            else:
                dico[k][nb_obs_irdis+i] = v
    
    # for i,file_IRDIS_level1 in enumerate(files_IRDIS):
        
    #     file_IFS_level1 = files_IFS[i]
    #     file_IRDIS_level2 = files_IRDIS_level2[i]
    #     file_IFS_level2 = files_IFS_level2[i]

    #     dirname_IRDIS_level1 = os.path.dirname(file_IRDIS_level1)
    #     filename_IRDIS_level1 = os.path.basename(file_IRDIS_level1)                

    #     name_tmp = file_IRDIS_level1[file_IRDIS_level1.index('script_')+7:]
    #     target_name = name_tmp[0:name_tmp.index('_')]

    #     dico_common = return_info_level1(file_IRDIS_level1)
    #     dico_common['target name'] = target_name
    #     for k, v in dico_common.items():
    #         if type(v) not in (float,int,np.float,np.float32,np.float64):
    #             dico[k].append(v)
    #             dico[k].append(v)
    #         else:
    #             dico[k][2*i] = v
    #             dico[k][2*i+1] = v

    #     print('processing',file_IRDIS_level1)
    #     dico_irdis = return_info_irdis_level1(file_IRDIS_level1)
    #     for k, v in dico_irdis.items():
    #         if type(v) not in (float,int,np.float,np.float32,np.float64):
    #             dico[k].append(v)
    #         else:
    #             dico[k][2*i] = v

    #     print('processing',file_IFS_level1)
    #     dico_ifs = return_info_ifs_level1(file_IFS_level1)
    #     for k, v in dico_ifs.items():
    #         if type(v) not in (float,int,np.float,np.float32,np.float64):
    #             dico[k].append(v)
    #         else:
    #             dico[k][2*i+1] = v
        
    pd_statistics = pd.DataFrame(dico)    
    pd_statistics.to_csv('/Volumes/Ast/super_earths/csv_tables/statistics.csv',index=None)
    
    pd_for_wiki  = pd.DataFrame(pd_statistics,columns=['target name',\
            'target name DC','camera','program ID','date start','date end',\
            'ob id','process','IRDIS filter','IFS disperser','nb frames','DIT [s]','MASS DIMM seeing avg ["]',\
            'MASS DIMM tau0 avg [ms]','delta parang','airmass avg','Strehl (measured)', 'Strehl (RTC)'])   
    pd_for_wiki.to_csv('/Volumes/Ast/super_earths/csv_tables/statistics_for_wiki.csv',index=None)
                                                      
    unique_names = np.unique(dico['target name'])
    pd_name_only = pd.DataFrame({'target name':unique_names}) 
    pd_name_only.to_csv('/Volumes/Ast/super_earths/csv_tables/names.csv',index=None)