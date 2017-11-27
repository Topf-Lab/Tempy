from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.class_arg import TempyParser
import sys,os, glob, numpy as np

tp = TempyParser()
tp.generate_args()
#map input
m = tp.args.inp_map
if m is None: 
    m = tp.args.inp_map1
    if m is None: sys.exit('Input a map, its resolution(required) and contour(optional)')
    elif tp.args.res1 is None: sys.exit('Input a map, is resolution(required) and contour(optional)')
    else: r = tp.args.res1
elif tp.args.res is None: sys.exit('Input a map, its resolution(required) and contour(optional)')
else: r = tp.args.res

#model input
pdir = tp.args.pdir
if pdir is None:
    plist = tp.args.plist
    if plist is None: sys.exit('Input path to the directory with fitted models (-pdir) or list of pdbs (-plist)')
    elif not os.path.isfile(plist): sys.exit('file with list of models not found')
    else:
        list_to_check = []
        with open(plist,'r') as pf:
            for ln in pf:
                if os.path.isfile(ln[:-1]): list_to_check.append(os.path.abspath(ln[:-1]))
                else: print 'File not found', l[:-1]
            
else:
    DATADIR = pdir
    list_models = glob.glob1(DATADIR,'*.pdb') #TODO check for pdb and cif extensions?
    list_to_check = [os.path.abspath(os.path.join(pdir,mo)) for mo in list_models]

#select score
score_sel = 'ccc' #TODO use argument input
list_styles = [':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-']
print '#', sys.argv

#calculate map contour
def map_contour(m,t=-1.):
  mName = os.path.basename(m).split('.')[0]
  #print 'reading map'
  emmap=MapParser.readMRC(m)
  c1 = None
  if t != -1.0:
    print 'calculating map contour'
    zeropeak,ave,sigma1 = emmap._peak_density()
    if not zeropeak is None: c1 = zeropeak+(t*sigma1)
    else:
      c1 = 0.0
  return mName,emmap,c1
#calculate model contour
def model_contour(p,res=4.0,emmap=False,t=-1.):
  pName,modelmap = blur_model(p,res,emmap)
  c1 = None
  if t != -1.0:
    print 'calculating model derived map contour'
    c1 = t*emmap.std()#0.0
  return pName,modelmap,c1
def blur_model(p,res=4.0,emmap=False):
  pName = os.path.basename(p).split('.')[0]
  print 'reading the model'
  structure_instance=PDBParser.read_PDB_file(pName,p,hetatm=False,water=False)
  print 'filtering the model'
  blurrer = StructureBlurrer()
  if res is None: sys.exit('Map resolution required..')
  #emmap = blurrer.gaussian_blur(structure_instance, res,densMap=emmap_1,normalise=True)
  modelmap = blurrer.gaussian_blur_real_space(structure_instance, res,sigma_coeff=0.187,densMap=emmap,normalise=True) 
  return pName,modelmap

#optional args
c1 = tp.args.thr
if c1 is None: c1 = tp.args.thr1
if c1 is None:
    Name1,emmap1,c1 = map_contour(m,t=1.5)
else:
    Name1 = os.path.basename(m).split('.')[0]
    emmap1=MapParser.readMRC(m)


dict_scores_hits = {}
list_models_calc = []
for pfile in list_to_check:
    Name2,emmap2,c2 = model_contour(pfile,res=r,emmap=emmap1,t=0.5)
    if None in [Name1,Name2]: sys.exit('Calculation failed, check input map and model files')
    print '#Scoring...', Name2
    sc = ScoringFunctions()
    #OVR
    try:
        ccc_mask,ovr = sc.CCC_map(emmap1,emmap2,c1,c2,3,meanDist=True)
        print 'Percent overlap:', ovr
        if ovr < 0.0: ovr = 0.0
    except:
        print 'Exception for lccc and overlap score'
        print_exc()
        ovr = 0.0
    if ovr < 0.02:
        print "Maps do not overlap: ", Name2
        continue
    #SCCC
    print 'Local correlation score: ', ccc_mask
    if ccc_mask < -1.0 or ccc_mask > 1.0:
        ccc_mask = 0.0
    #LMI
    try:
        mi_mask = sc.MI(emmap1,emmap2,c1,c2,3)
        print 'Local Mutual information score: ', mi_mask
        if mi_mask < 0.0: mi_mask = 0.0
    except:
        print 'Exception for MI score'
        print_exc()
        mi_mask = 0.0
      #NMI
    try:
        nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
        print 'Normalized Mutual information score:', nmi
        if nmi < 0.0: nmi = 0.0
    except:
        print 'Exception for NMI score'
        print_exc()
        nmi = 0.0
        
    list_models_calc.append(Name2)
    try: dict_scores_hits['local_correlation'].append(ccc_mask)
    except KeyError: dict_scores_hits['local_correlation'] = [ccc_mask]
    try: dict_scores_hits['local_mi'].append(mi_mask)
    except KeyError: dict_scores_hits['local_mi'] = [mi_mask] 
    try: dict_scores_hits['overlap'].append(ovr)
    except KeyError: dict_scores_hits['overlap'] = [ovr]
    try: dict_scores_hits['nmi'].append(nmi)
    except KeyError: dict_scores_hits['nmi'] = [nmi]

dict_scores_hits['mi_ov'] = sc.scale_median(dict_scores_hits['overlap'],dict_scores_hits['local_mi'])
dict_scores_hits['ccc_ov'] = sc.scale_median(dict_scores_hits['overlap'],dict_scores_hits['local_correlation'])

'''
for i in range(len(list_models_calc)):    
    med_dev_ov = np.median(np.absolute(dict_scores_hits['overlap'] - np.median(dict_scores_hits['overlap'])))
    #MI OVR    
    med_dev = np.median(np.absolute(dict_scores_hits['local_mi'] - np.median(dict_scores_hits['local_mi'])))
    scale_ov = med_dev/med_dev_ov
    shift_ov = np.median(dict_scores_hits['local_mi'])-(scale_ov*np.median(dict_scores_hits['overlap']))
    if (max(dict_scores_hits['overlap']) - min(dict_scores_hits['overlap'])) > 0.1: 
        dict_scores_hits['mi_ov'].append(((scale_ov*dict_scores_hits['overlap'][i]+shift_ov) + dict_scores_hits['local_mi'][i])/2.)
    else: dict_scores_hits['mi_ov'].append(dict_scores_hits['local_mi'][i])
    
    #CCC OVR
    med_dev = np.median(np.absolute(dict_scores_hits['local_correlation'] - np.median(dict_scores_hits['local_correlation'])))
    scale_ov = med_dev/med_dev_ov
    shift_ov = np.median(dict_scores_hits['local_correlation'])-(scale_ov*np.median(dict_scores_hits['overlap']))
    if (max(dict_scores_hits['overlap']) - min(dict_scores_hits['overlap'])) > 0.1: 
        dict_scores_hits['ccc_ov'].append(((scale_ov*dict_scores_hits['overlap'][i]+shift_ov) + dict_scores_hits['local_correlation'][i])/2.)
    else: dict_scores_hits['ccc_ov'].append(dict_scores_hits['local_correlation'][i])
'''

list_miov = dict_scores_hits['mi_ov'][:]
sorted_ind = sorted(range(len(list_miov)), key=list_miov.__getitem__)
#list_miov.sort()
print '**Models ranked based on LMI_OV score'
print '{:15} {:6} {:6} {:6}'.format('Model', 'LMI_OV', 'LMI', 'OVR')
for i in reversed(sorted_ind):
    print '{:15} {:6.4f} {:6.4f} {:6.4f}'.format(list_models_calc[i], list_miov[i], dict_scores_hits['local_mi'][i] \
                                                 , dict_scores_hits['overlap'][i])
list_cccov = dict_scores_hits['ccc_ov'][:]
sorted_ind = sorted(range(len(list_cccov)), key=list_cccov.__getitem__)
print '**Models ranked based on SCCC_OV score'
print '{:15} {:6} {:6} {:6}'.format('Model', 'SCCC_OV', 'SCCC', 'OVR')
for i in reversed(sorted_ind):
    print '{:15} {:6.4f} {:6.4f} {:6.4f}'.format(list_models_calc[i], list_cccov[i], dict_scores_hits['local_correlation'][i] \
    , dict_scores_hits['overlap'][i])

list_nmi = dict_scores_hits['nmi'][:]
sorted_ind = sorted(range(len(list_nmi)), key=list_nmi.__getitem__)
print '**Models ranked based on NMI score'
print '{:15} {:6}'.format('Model', 'NMI')
for i in reversed(sorted_ind):
    print '{:15} {:6.4f}'.format(list_models_calc[i], dict_scores_hits['nmi'][i])
    
try: import matplotlib.pyplot as plt
except ImportError:flagread = 0 
try: plt.style.use('ggplot')
except AttributeError: pass
from matplotlib import pylab
#axes = plt.gca()
#fig, axes = plt.subplots()
#fig.tight_layout()
plt.xlabel = 'Model'
plt.ylabel = 'Scores'
x = range(1,len(list_models_calc)+1)
plt.plot(x,list_miov,linewidth=3.0,label='LMI_OV')
plt.plot(x,dict_scores_hits['local_mi'],linewidth=3.0,label='LMI')
plt.plot(x,list_cccov,linewidth=3.0,label='SCCC_OV')
plt.plot(x,dict_scores_hits['local_correlation'],linewidth=3.0,label='SCCC')
plt.plot(x,dict_scores_hits['overlap'],linewidth=3.0,label='OVR')
#axes.set_xticklabels(list_models_calc)
plt.xticks(x,list_models_calc, rotation='vertical')
leg = plt.legend(loc='upper right')
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)
pylab.savefig("scoreplot.png",bbox_inches='tight')
plt.close()
        
    
    