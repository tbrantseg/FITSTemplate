#!/opt/local/bin/python2.7

from FITSTemplate import *
import argparse

# Expected parameters
DefaultFits=FITSCore()

# The following parameters are not defined in FITSCore
intensity = 1.0
output='test.fits'
epsilon=1.0
theta=0.0

# Establish what parameters we are expecting
parser = argparse.ArgumentParser(description="FITS Template Generator")
parser.add_argument('-R','--RA', type=float, help='right ascension (default=%6.3f)'%DefaultFits.RA_val,required=False)
parser.add_argument('-D','--Dec', type=float, help='declination (default=%+6.3f)'%DefaultFits.dec_val,required=False)
parser.add_argument('-I','--intensity', type=float, help='peak intensity (default=%3.2f)'%intensity,required=False)
parser.add_argument('-O','--output', help='output file name (default=\'%s\')'%output, required=False)
parser.add_argument('-P','--pixsize', type=float, help='pixel size (degrees, default=%5.4f)'%DefaultFits.px_size,required=False)
parser.add_argument('-T','--type', help='template type (\'disk\',\'gaus\',\'asymgaus\')',required=True)
parser.add_argument('-X','--xsize', type=int, help='number of pixels in x (default=%d)'%DefaultFits.xsize,required=False)
parser.add_argument('-Y','--ysize', type=int, help='number of pixels in y (default=%d)'%DefaultFits.ysize,required=False)
parser.add_argument('-o','--overwrite', type=bool, help='overwrite(\'True\',\'False\')', default=False, required=False)
parser.add_argument('-r','--radius', type=float, help='radius of disk model in degrees', required=False)
parser.add_argument('-s','--sigma', type=float, help='sigma of Gaussian model in degrees', required=False)
parser.add_argument('-e','--epsilon', type=float, help='eccentricity of asymmetric Gaussian model(0->1, default=%f)'%epsilon, required=False)
parser.add_argument('-t','--theta', type=float, help='rotation angle of asymmetric Gaussian model (degrees, default=%d)'%theta, required=False)
parser.add_argument('-p','--projection', type=string, help='Projection type', required=False)

# Parse the input arguments
args = parser.parse_args()

# Setup the passed values
if (args.RA is not None):
    DefaultFits.RA_val=args.RA
if (args.Dec is not None):
    DefaultFits.dec_val=args.Dec
if (args.intensity is not None):
    intensity=args.intensity
if (args.output is not None):
    output=args.output
if (args.pixsize is not None):
    DefaultFits.px_size=args.pixsize
if (args.xsize is not None):
    DefaultFits.xsize=args.xsize
if (args.ysize is not None):
    DefaultFits.ysize=args.ysize
if (args.epsilon is not None):
    epsilon=args.epsilon
if (args.theta is not None):
    theta=args.theta
if (args.projection is not None):
    projection=args.projection

ThisFits=FITSCore("",DefaultFits.xsize, DefaultFits.ysize, DefaultFits.RA_val, DefaultFits.dec_val, DefaultFits.px_size,projection)

# Print some information about the model that is about to be produced
print ''
print 'Input model parameters:'
print '   Output    : %s' % output
print '   overwrite : %s' % args.overwrite
print '   xsize     : %d pixels' % ThisFits.xsize
print '   ysize     : %d pixels' % ThisFits.ysize
print '   pixel size: %f degrees' % ThisFits.px_size
print '   intensity : %f' % intensity
print '   RA        : %f deg' % ThisFits.RA_val
print '   Dec       : %+f deg' % ThisFits.dec_val
print '   Type      : %s' % args.type
print '   Projection: {0}'.format(projection)
if (args.type=='disk'):
    print '   radius    : %f degrees' % args.radius
elif(args.type=='gaus'):
    print '   sigma     : %f degrees' % args.sigma
elif(args.type=='asymgaus'):
    print '   sigma     : %f degrees' % args.sigma
    print '   epsilon   : %f' % epsilon
    print '   theta     : %f' % theta
print ''
################################
# Handle the disk model
################################
if (args.type == 'disk'):
    # make sure a radius was given
    if ((args.radius=='None')or(args.radius<0.0)):
        print "[ERROR] Radius must be positive!"
    else:
        # Convert the radius from degrees to pixels
        radius = args.radius/ThisFits.px_size
        # Create the disk template
        print 'Creating disk model with radius=%f pixels' % radius
        FITSCore.DiskTemplate(ThisFits, radius, intensity, output, args.overwrite)

################################
# Handle Gaussian model
################################
elif (args.type == 'gaus'):
    # Make sure that a sigma was supplied
    if ((args.sigma is None) or (args.sigma<=0.0)):
        print "[ERROR] Gaussian sigma must be greater than 0!"
    else:
        # Convert the sigma to FWHM in units of pixels
        FWHM = (args.sigma*2.7725887)/ThisFits.px_size
        
        print 'Creating Gaussian model with sigma=%f (FWHM=%f pixels)' % (args.sigma, FWHM)
        # Create the Gaussian template
        FITSCore.GaussianTemplate(ThisFits, FWHM, intensity, output, args.overwrite)

#################################
# Handle Asymmetric Gaussian
#################################
elif (args.type == 'asymgaus'):
    # make sure an appropriate value for sigma was supplied
    if ((args.sigma is None)or(args.sigma<=0.0)):
        print "[ERROR] Gaussian sigma must be greater than 0!"
    elif (args.epsilon is None):
        print "[ERROR] Epsilon must be defined!"
    else:
        # Convert he sigma to FWHM in units of pixels
        FWHM = (args.sigma*2.7725887)/ThisFits.px_size
        # Create the asymmetric Gaussian template
        print 'Creating Gaussian model with sigma=%f (FWHM=%f pixels), epsilon=%f, theta=%f' % (args.sigma, FWHM, epsilon, theta)
        FITSCore.EllipticalGaussianTemplate(ThisFits, FWHM, intensity, epsilon, theta, output, args.overwrite)


# Finish GenTemplate.py
