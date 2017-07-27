import FITSTemplates
import numpy as np
tmpl = FITSTemplates.FITSCore(xsize=20, ysize=20, RA=83.63, dec=22.01, px_size=0.1)
seed_array = np.array([[1,1],[1,1]])
tmpl.AddArray(tmpl.field,seed_array,len(tmpl.field[1])/2,len(tmpl.field)/2)
'''tmpl.AddArray(tmpl.field,...,...,...)'''
tmpl.Around(tmpl.field)
tmpl._SaveToFITS("template.fits",True)
