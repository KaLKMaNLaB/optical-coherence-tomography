# -*- coding: utf-8 -*-
"""
Author = JdeWit
Date: 2019-08-15
email: J.deWit-1@tudelft.nl

This is a series of functions to read a .oct files from Thorlabs OCT software
and extract the raw data as well as imaging parameters from these files.
"""

import numpy as np
import zipfile
from bs4 import BeautifulSoup

def OCTgetDataCombined(filepath,spectrumindex=0):
    '''
    This function combines the functions below and gives the header, raw data, 
    spectrum and the field of view from the filepath.
    input parameters:
        filepath: path where the data file can be found including the file name
        spectrumindex:determine which Bscan to load the data from. This is 
            especially relevant for recordings with multiple Bscans or 3D datasets.
    output parameters:
        header: header containing the metadata of the measurement
        rawdata: file with the raw data (2D) from the OCT spectrometer 
            (compensated for the camera offset)
        spectrum: spectrum of the reference arm
        FOV: field of view, list with order [FOVz,FOVx,FOVy] in the 1D and 2D 
            mode it is a list with 2 items
    '''
    directory=OCTfileOpen(filepath)
    header=OCTreadHeader(directory)
    rawdata,spectrum=OCTgetRawData(directory,header,spectrumindex)
    FOV=OCTgetFOV(header)
    return [header,rawdata,spectrum,FOV]

def OCTfileOpen(filepath):
    '''
    this function creates pointer into the zip file where the header and the 
    data can be extracted from
    '''
    zip_ref = zipfile.ZipFile(filepath,'r')
    return zip_ref

def OCTreadHeader(directory):
    '''
    this function reads the header from the xml file in the directory that 
    is given as input
    '''
    file=directory.open('Header.xml')
    header=BeautifulSoup(file,'lxml-xml')
    return header

def OCTgetRawData(directory,header,spectrumindex=0):
    ''' 
    This function obtains the raw interference data with the indicated 
    spectrum index as well as the apodization spectrum for the measurement.
    Both the raw data and the apodization spectrum are corrected with the 
    offset
    input parameters:
        directory: pointer into the zip file 
        header: header corresponding to the measurements
        spectrumindex: determine which Bscan to load the data from. This is 
            especially relevant for recordings with multiple Bscans or 3D datasets.
    output parameters:
        raw2: the interference spectrum compensated for the camera offset.
        ApodizationSpectrum: the spectrum from the reference arm compensated 
            for the camera offset
    '''    
    # get offset
    offset_obj=directory.open('data/OffsetErrors.data')
    offset=offset_obj.read()
    offset=np.frombuffer(offset,dtype=np.float32)
    
    # get raw data
    bbPixel = int(header.Ocity.Instrument.BytesPerPixel.string)
    isSigned = header.Ocity.Instrument.RawDataIsSigned.string;
    if bbPixel==2:
        if isSigned=='False':
            dtype = np.uint16
        else:
            dtype = np.int16
    BinaryToElectronCountScaling=np.double(header.Ocity.Instrument.BinaryToElectronCountScaling.string)
    
    Raw_Data_File=header.find("DataFile",Type='Raw',string="data\Spectral"+str(spectrumindex)+".data")
    if not Raw_Data_File:
        print('Error: the desired spectrum is not found')
    else: 
        size=[int(Raw_Data_File['SizeX']),int(Raw_Data_File['SizeZ'])]
        ScanRegionStart0=int(Raw_Data_File['ScanRegionStart0'])
        try: 
            NumApos=int(Raw_Data_File['ApoRegionEnd0'])
        except:
            NumApos=0
             
        raw_data_obj=directory.open('data/Spectral'+str(spectrumindex)+'.data')
        raw_data=raw_data_obj.read()
        raw=np.frombuffer(raw_data,dtype)
        raw2=np.reshape(raw,size)*BinaryToElectronCountScaling
        for i in range(size[0]):
            raw2[i,:]=raw2[i,:]-offset
                    
    # get reference spectrum
    if NumApos==0:
        ApodizationSpectrum_obj=directory.open('data/ApodizationSpectrum.data')
        ApodizationSpectrum=ApodizationSpectrum_obj.read()
        ApodizationSpectrum=np.frombuffer(ApodizationSpectrum,dtype=np.float32)-offset
    else:
        ApodizationSpectrum=np.mean(raw2[0:NumApos,:],axis=0)
            
    # select raw data that makes up the scan
    if ScanRegionStart0>0:
        raw2=raw2[ScanRegionStart0:,:] #-1
    
    return [raw2, ApodizationSpectrum]

def OCTgetFOV(header):
    '''
    this file extract the Field of View from the header and converts it into
    meters
    '''
    FOV=np.zeros(3)
    FOV[0]=np.double(header.Ocity.Image.SizeReal.SizeZ.string)
    if header.Ocity.MetaInfo.AcquisitionMode.string=="Mode2D":
        FOV[1]=np.double(header.Ocity.Image.SizeReal.SizeX.string)
        FOV=FOV[0:2]
    elif header.Ocity.MetaInfo.AcquisitionMode.string=="Mode3D":
        FOV[1]=np.double(header.Ocity.Image.SizeReal.SizeX.string)
        FOV[2]=np.double(header.Ocity.Image.SizeReal.SizeY.string)
    else:
        FOV=FOV[0:2]
    FOV=FOV*1e-3
    return FOV


