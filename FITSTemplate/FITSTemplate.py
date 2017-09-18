#!/opt/local/bin/python

import sys
import os
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
#import scipy.ndimage.interpolation as spr
import skimage.transform as skt
from astropy import wcs
import ConfigParser

# Global dictionary of fits header types
headers = {}
headers["TAN"] = ["RA---TAN","DEC--TAN"]
headers["AIT"] = ["RA---AIT","DEC--AIT"]
headers["GTAN"] = ["GLAT-TAN","GLON-TAN"]
headers["GAIT"] = ["GLAT-AIT","GLON-AIT"]

class FITSCore (object):

    def __init__ (self,\
                 paramfile="",\
                 xsize=None,\
                 ysize=None,\
                 RA=None,\
                 dec=None,\
                 px_size=None,\
                 projection="TAN"):
        """
        Class to hold basic data on creating a FITS template. This class must
        be instantiated before any templates are generated.

        Constructor Arguments:
        ---------------------
        paramfile: Name of the parameter file holding information on the FITS file.

        If paramfile is not specified, other arguments may be given explicitly.
        
        xsize: Size of the FITS file, in pixels.
        ysize: Size of the FITS file in pixels.
        RA: Decimal representation of the target right ascension.
        dec: Decimal representation of the target declination.
        px_size: Angular size of the FITS pixels in degrees.
        projection: Projection on to the sky plane to be used:
           "TAN": Tangential projection, used for small areas
                "RA---TAN", "DEC--TAN"
           "AIT": Aitoff projection, used by Fermi-LAT
                "RA---AIT", "DEC--AIT"
           "GTAN": Tangential projection with galactic coordinates
                "GLAT---TAN", "GLON--TAN"
           "GAIT": Aitoff projection with galactic coordinates
                "GLAT---AIT", "GLON--AIT"

        Note that if both paramfile *and* the explicit arguments are given,
        the values in paramfile take priority.

        Methods:
        --------
        DiskTemplate: Generate a circular source with constant intensity.
        GaussianTemplate: Generate a circular source with a Gaussian intensity.
        EllipticalGaussianTemplate: Generate an elliptical source with a Gaussian
             intensity.
        """
        if (os.path.isfile(paramfile)):
            config=ConfigParser.RawConfigParser()
            config.read(paramfile)
            self.xsize = config.getint("FOV","xsize")
            self.ysize = config.getint("FOV","ysize")
            self.RA_val = config.getfloat("FOV","RA_val")
            self.dec_val = config.getfloat("FOV","dec_val")
            self.px_size = config.getfloat("FOV","px_size")
            self.projection = ""
        else:
            self.xsize = int(xsize)
            self.ysize = int(ysize)
            self.RA_val = float(RA)
            self.dec_val = float(dec)
            self.px_size = float(px_size)
            self.projection = str(projection)            

        try:
            self.header = headers[self.projection]
        except KeyError:
            print "Header type {0} not recognized. Defaulting to TAN projection.".format(self.header)
            self.header = headers['TAN']
            
        # This should be a dictionary, shouldn't it?
        self.field=np.zeros((self.xsize,self.ysize))
        self.x_c = self.xsize/2.0-0.5
        self.y_c = self.ysize/2.0-0.5

        
    def _SaveToFITS (self, outname, overwrite):
        # FIXME: Galactic coordinate inputs not handled correctly
        # FIXME: Cos dec size distortion issue? Might be upstream
        w = wcs.WCS(naxis=2)
        xvals = np.arange(0,self.xsize)
        yvals = np.arange(0,self.ysize)
        w.wcs.crval = [self.RA_val,self.dec_val]
        w.wcs.cdelt = np.array([self.px_size*-1.0,self.px_size])
        w.wcs.crpix = [0.5+self.xsize/2.0,0.5+self.ysize/2.0]
        w.wcs.ctype = self.header
        header = w.to_header()
        hdu = fits.PrimaryHDU(self.field,header=header)
        try:
            hdu.writeto(outname,clobber=overwrite)
        except IOError:
            print "Could not create file! {0} already exists. Try again with overwrite set to True.".format(outname)

            
    def DiskTemplate (self, radius, intensity, outname, overwrite=False):
        """
        Generate a circular source with constant intensity.

        Function Call:
        -------------
        FITSCore.DiskTemplate(radius,intensity,outname)

        Arguments:
        ---------
        radius: Radius of the disk, in pixels.
        intensity: Intensity of the source, in arbitrary units.
        outname: Name of the output file.
        overwrite(optional): If the output file exists, should it be overwritten?
             Defaults to False.
        """
        for x in range(0,self.xsize):
            for y in range(0,self.ysize):
                if ((x-self.x_c)**2 + (y-self.y_c) **2 <= radius**2):
                    self.field[x][y] = intensity
        if outname:
            self._SaveToFITS(outname,overwrite)

        
    def SplitDiskTemplate(self, radius, intensity_a, intensity_b, split_angle, outname, overwrite=False):
        """
        Generate a disk split between two different intensities.

        Function Call:
        -------------
        FITSCore.DiskTemplate(radius,intensity_a,intensity_b,split_angle,outname)

        Arguments:
        ---------
        radius: Radius of the disk, in pixels.
        intensity_a: Intensity of one half of the disk, in arbitrary units.
        intensity_b: Intensity of the other half of the disk, in arbitrary units.
        split_angle: Angle at which the dividing line between the two halves of the
         disk will be, in degrees clockwise from straight north-south.
        outname: File output name.
        overwite: If the output file exists, should it be overwritten? Defaults to False.
        """
        temp_field = np.zeros
        for x in range(0,self.xsize):
            for y in range(0,self.ysize):
                if (x-self.x_c)**2 + (y-self.y_c)**2 <= radius**2:
                    if y >= self.y_c:
                        self.field[x][y] = intensity_a
                    elif y < self.y_c:
                        self.field[x][y] = intensity_b
        # Now rotate
        # spr.rotate(input = self.field,\
        #                                 angle = split_angle,\
        #                                 reshape = False,\
        #                                 output = self.field)
        # FIXME: Rotate function results in values that are not the two intended values; maybe just use a rotation matrix here
        self.field=skt.rotate(image=self.field, angle=split_angle,center=(self.x_c,self.y_c))
        # Kill any junk outside the intended radius
        for x in range(0,self.xsize):
            for y in range(0,self.ysize):
                if (x-self.x_c)**2 + (y-self.y_c)**2 > radius**2:
                    self.field[x][y] = 0
        if outname:
            self._SaveToFITS(outname,overwrite)

            
    def EllipticalGaussianTemplate(self, FWHM, peak, epsilon, theta, outname, overwrite=False):
        """
        Generate an elliptical source with a Gaussian intensity.

        Function Call:
        -------------
        FITSCore.EllipticalGaussianTemplate(FWHM,peak,epsilon,theta,outname)

        Arguments:
        ---------
        FWHM: Full-width half-maximum of the source, in pixels.
        peak: Peak value of the source brightness.
        epsilon: Eccentricity of the source. Throws ValueError if out of bounds
             (>1 or <0).
        theta: Source rotation angle, in degrees.
        outname: Name of the output file.
        overwrite(optional): If the output file exists, should it be overwritten?
             Defaults to False.        
        """
        # Throw ValueError if we're out of bounds for the eccentricity.
        if (epsilon > 1 or epsilon < 0):
            raise ValueError('Epsilon must be a value between 0 and 1.')
            
        # As defined on the CXC Sherpa website.
        FWHM=float(FWHM)
        peak=float(peak)
        theta=np.radians(theta)
        gauss_c = 2.7725887 # 4 log 2; relationship between FWHM and sigma
        ellip = float(1-epsilon)
        for x in range(0,self.xsize):
            for y in range(0,self.ysize):
                x_o = (x-self.x_c)*np.cos(theta) + (y-self.y_c)*np.sin(theta)
                y_o = (y-self.y_c)*np.cos(theta) - (x-self.x_c)*np.sin(theta)
                r = np.sqrt(((x_o**2 * ellip**2) + y_o**2)/ellip)
                self.field[x][y]=peak*np.e**(-gauss_c*(r/FWHM)**2)
        if outname:           
            self._SaveToFITS(outname, overwrite)
        
    def GaussianTemplate (self, FWHM, peak, outname,overwrite=False):
        """
        Generate a circular source with a Gaussian intensity.

        Function Call:
        -------------
        FITSCore.GaussianTemplate(FWHM,peak,outname)

        Arguments:
        ---------
        FWHM: Full-width half-maximum of the source, in pixels.
        peak: Peak value of the source brightness.
        outname: Name of the output file.
        overwrite(optional): If the output file exists, should it be overwritten?
             Defaults to False.
        """
        self.EllipticalGaussianTemplate(FWHM,peak,0,0,outname,overwrite)

    def AnnularTemplate(self,r_in,r_out,intensity,outname,overwrite=False):
        """
        Generate an annular source.

        Function Call:
        -------------
        FITSCore.AnnularTemplate(r_in,r_out,intensity,outname,[outname=False])
        

        Arguments:
        ---------
        r_in: annulus inner radius
        r_out: annulus outer radius
        intensity: value to set template pixels to
        outname: Name of the outputfile
        overwite(optional): If the output file exists, should it be overwritten?
        
        """
        for x in range(0,self.xsize):
            for y in range(0,self.ysize):
                r = np.sqrt(((x-self.x_c)**2) + ((y-self.y_c)**2))
                if (r <= r_out and r >= r_in):                    
                    self.field[x][y] = intensity
        if outname:
            self._SaveToFITS(outname,overwrite)

    def DiffusionProfileTemplate(self, r_s, r_d, norm, outname, overwrite=False):
        """
        Generate a symmetric template with a diffusion radial brightness profile.
        The brightness is constant until r_s, then falls off in a roughly Gaussian
        manner with a characteristic scale length of r_d.
        """        
        for x in range(0, self.xsize):
            for y in range(0, self.ysize):
                r = np.sqrt(((x-self.x_c)**2) + ((y-self.y_c)**2))
                if (r <= r_s):
                    self.field[x][y] = norm
                else:
                    r = r - r_s
                    self.field[x][y] = norm/(r_d*(r+0.06*r_d))*np.exp(-(r**2)/(r_d**2))
        if outname:
            self._SaveToFITS(outname,overwrite)

# FITSTemplate.py ends here
