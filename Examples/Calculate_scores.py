from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
import os,sys
from TEMPy.class_arg import TempyParser
from traceback import print_exc

EXAMPLEDIR = 'Test_Files'
print 'use --help for help'
print '-m/-m1 [map] for input map, -m1,-m2 for two input maps'
print '-p/-p1 [pdb] for input pdb'
print '-r [resolution]; -r1,r2 for the two map resolutions'
print '-t [density threshold]; -t1,t2 for the two map thresholds'
print '\n\n###################################'

sc = ScoringFunctions()
tp = TempyParser()
tp.generate_args()
#COMMAND LINE OPTIONS
m1 = tp.args.inp_map1
m2 = tp.args.inp_map2
m = tp.args.inp_map
r1 = tp.args.res1
r2 = tp.args.res2
r = tp.args.res
c1 = tp.args.thr1
c2 = tp.args.thr2
c = tp.args.thr
p = tp.args.pdb
p1 = tp.args.pdb1
p2 = tp.args.pdb2
list_to_check = []
#multiple models
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

    
#EXAMPLE RUN
flag_example = False
if len(sys.argv) == 1:
  path_example=os.path.join(os.getcwd(),EXAMPLEDIR)
  if os.path.exists(path_example)==True:
    print "%s exists" %path_example
  #else: sys.exit('No input')
  #os.chdir(path_out)
  flag_example = True

#calculate map contour
def map_contour(m,t=-1.):
  mName = os.path.basename(m).split('.')[0]
  #print 'reading map'
  emmap=MapParser.readMRC(m)
  c1 = None
  if t != -1.0:
    #print 'calculating contour', mName
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
    #print 'calculating contour'
    c1 = t*emmap.std()#0.0
  return pName,modelmap,c1
def blur_model(p,res=4.0,emmap=False):
  pName = os.path.basename(p).split('.')[0]
  #print 'reading the model', pName
  structure_instance=PDBParser.read_PDB_file(pName,p,hetatm=False,water=False)
  #print 'filtering the model', pName
  blurrer = StructureBlurrer()
  if res is None: sys.exit('Map resolution required..')
  #emmap = blurrer.gaussian_blur(structure_instance, res,densMap=emmap_1,normalise=True)
  modelmap = blurrer.gaussian_blur_real_space(structure_instance, res,sigma_coeff=0.187,densMap=emmap,normalise=True) 
  return pName,modelmap
  
def lpfilter(emmap,r):
    cutoff = emmap.apix/float(r)
    #lowpass filter
    #print 'lowpass'
    mapfilt = emmap._tanh_lowpass(cutoff)
    return mapfilt
  
def match_grid(emmap1,emmap2,c1,c2):
  # DETERMINE A COMMON ALIGNMENT BOX : fill minvalue for extra voxel pads
    spacing = emmap2.apix
    if emmap2.apix < emmap1.apix: spacing = emmap1.apix
    grid_shape, new_ori = emmap1._alignment_box(emmap2,spacing)
    # INTERPOLATE TO NEW GRID
    try: emmap_1 = emmap1._interpolate_to_grid1(grid_shape,spacing,new_ori)
    except: emmap_1 = emmap1._interpolate_to_grid(grid_shape,spacing,new_ori)
    try: c1 = emmap_1._find_level(np.sum(emmap1.fullMap>c1)*(emmap1.apix**3))
    except: pass
    del emmap1.fullMap
    del emmap1
    try: emmap_2 = emmap2._interpolate_to_grid1(grid_shape,spacing,new_ori)
    except: emmap_2 = emmap2._interpolate_to_grid(grid_shape,spacing,new_ori)
    try: c2 = emmap_2._find_level(np.sum(emmap2.fullMap>c2)*(emmap2.apix**3))
    except: pass
    del emmap2.fullMap
    del emmap2
    return emmap_1, emmap_2


#GET INPUT DATA
if flag_example:
  
  p = os.path.join(path_example,'1J6Z.pdb')
  m = os.path.join(path_example,'emd_5168_monomer.mrc')
  res = 6.6
  Name1 = os.path.basename(m).split('.')[0]
  Name2 = os.path.basename(p).split('.')[0]
  emmap1=MapParser.readMRC(m)
  structure_instance=PDBParser.read_PDB_file(Name2,p,hetatm=False,water=False)
  blurrer = StructureBlurrer()
  emmap2 = blurrer.gaussian_blur(structure_instance, res,densMap=emmap1)
  c1 = 9.7
  c2 = 1.0
elif all(x is None for x in [m,m1,m2]):
    # for 2 models
    if None in [p1,p2]:
        sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
    Name1,emmap1,c1 = model_contour(p1,res=4.0,emmap=False,t=0.5)
    r1 = r2 = r = 4.0
    if c2 is None: Name2,emmap2,c2 = model_contour(p2,res=r,emmap=False,t=0.5)
    else: Name2,emmap2 = blur_model(p2,res=r,emmap=False)
    flag_filt = False
    flag_scale = False
    list_compare = [emmap2]
    list_names = [Name2]
    list_contours = [c2]
elif None in [m1,m2]:
    # for one map and model
    if m1 is None: m = tp.args.inp_map
    else: m = m1
    m = tp.args.inp_map
    #print 'reading map'
    if c is None: 
        if c1 is None: Name1,emmap1,c1 = map_contour(m,t=1.5)
        else: 
            Name1 = os.path.basename(m).split('.')[0]
            emmap1=MapParser.readMRC(m)
    else:
        Name1 = os.path.basename(m).split('.')[0]
        emmap1=MapParser.readMRC(m)
    if r1 is None and r is None: sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
    elif r1 is None: r1 = r
    if len(list_to_check) == 0:
        if all(x is None for x in [p,p1,p2]): sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
        elif p is None: p = p2
        if p is None: sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
        Name2,emmap2,c2 = model_contour(p,res=r1,emmap=emmap1,t=0.5)
        list_compare = [emmap2]
        list_names = [Name2]
        list_contours = [c2]
    #multiple models
    else: 
        list_compare = []
        list_names = []
        list_contours = []
        for p in list_to_check:
            l_details = model_contour(p,res=r1,emmap=emmap1,t=0.5)
            list_compare.append(l_details[1])
            list_names.append(l_details[0])
            list_contours.append(l_details[2])
        
    r2 = 3.0
    #print 'reading model'
    
else: 
    # For 2 input maps
    if None in [r1,r2]: sys.exit('Input two maps, their resolutions(required) and contours(optional)')
    print 'reading map1'
    if c1 is None:
        Name1,emmap1,c1 = map_contour(m1,t=1.5)
    else:
        Name1 = os.path.basename(m1).split('.')[0]
        emmap1=MapParser.readMRC(m1)
    print 'reading map2' 
    if c2 is None:
        Name2,emmap2,c2 = map_contour(m2,t=1.5)
    else:
        Name2 = os.path.basename(m2).split('.')[0]
        emmap2=MapParser.readMRC(m2)
        
    if not sc.mapComparison(emmap1,emmap2):
        if None in [r1,r2]: sys.exit('Input two maps, their resolutions(required) and contours(optional)')
        emmap1._crop_box(c1,0.5)
        emmap2._crop_box(c2,0.5)
        
        if r1 > 1.25*r2: 
            emmap_2 = lpfilter(emmap2,r1)
            emmap1, emmap2 = match_grid(emmap1,emmap_2,c1,c2)
        elif r2 > 1.25*r1:
            emmap_1 = lpfilter(emmap1,r2)
            emmap1, emmap2 = match_grid(emmap_1,emmap2,c1,c2)
        else:
            emmap1, emmap2 = match_grid(emmap1,emmap2,c1,c2)
    list_compare = [emmap2]
    list_names = [Name2]
    list_contours = [c2]    

scores = {}
sc = ScoringFunctions()

for i in range(len(list_compare)):
  emmap2 = list_compare[i]
  Name2 = list_names[i]
  c2 = list_contours[i]
  print 'Scoring...', Name2  
  if None in [Name1,Name2]: sys.exit('Scoring failed!, check input data')
  #OVR
  try:
    ccc_mask,ovr = sc.CCC_map(emmap1,emmap2,c1,c2,3)
    #print 'Percent overlap:', ovr
    if ovr < 0.0: ovr = 0.0
  except:
    print 'Exception for lccc and overlap score'
    print_exc()
    ovr = 0.0
  try: scores['overlap'].append(ovr)
  except KeyError: scores['overlap'] = [ovr]
  
  if ovr < 0.02:
    sys.exit("Maps do not overlap.")
  #SCCC
  #print 'Local correlation score: ', ccc_mask
  if ccc_mask < -1.0 or ccc_mask > 1.0:
    ccc_mask = 0.0
  try: scores['local_correlation'].append(ccc_mask)
  except KeyError: scores['local_correlation'] = [ccc_mask]
  
  #LMI
  try:
    mi_mask = sc.MI(emmap1,emmap2,c1,c2,3)
    #print 'Local Mutual information score: ', mi_mask
    if mi_mask < 0.0: mi_mask = 0.0
  except:
    print 'Exception for MI score'
    print_exc()
    mi_mask = 0.0
  try: scores['local_mutual_information'].append(mi_mask)
  except KeyError: scores['local_mutual_information'] = [mi_mask]
  
  #CCC
  try:
    ccc,ovr = sc.CCC_map(emmap1,emmap2,c1,c2)
    #print 'Correlation score: ', ccc
    if ovr < 0.0: ovr = 0.0
  except:
    print 'Exception for ccc score'
    print_exc()
    ovr = 0.0
  if ovr < 0.02:
    sys.exit("Maps do not overlap.")
  try: scores['correlation'].append(ccc)
  except KeyError: scores['correlation'] = [ccc]
    
  #NMI
  try:
    nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
    #print 'Normalized Mutual information score:', nmi
    if nmi < 0.0: nmi = 0.0
  except:
    print 'Exception for NMI score'
    print_exc()
    nmi = 0.0
  try: scores['normalised_mutual_information'].append(nmi)
  except KeyError: scores['normalised_mutual_information'] = [nmi]
  
  #CD
  try:
    chm = sc._surface_distance_score(emmap1,emmap2,c1,c2,'Minimum')
    if chm == 0.0 or chm is None:
      chm = sc._surface_distance_score(emmap1,emmap2,c1,c2,'Mean')
    #print 'Surface distance score: ', chm
    if chm < 0.0: chm = 0.0
  except:
    print 'Exception for surface distance score'
    print_exc()
    chm = 0.0
  try: scores['surface_distance'].append(chm)
  except KeyError: scores['surface_distance'] = [chm]
  
  #NV
  try:
    #nv = sc.normal_vector_score(emmap_1,emmap_2,float(c1)-(emmap_1.std()*0.05),float(c1)+(emmap_1.std()*0.05),None)
    nv = sc.normal_vector_score(emmap1,emmap2,float(c1),float(c1)+(emmap1.std()*0.05),'Minimum')
    if nv == 0.0 or nv is None: nv = sc.normal_vector_score(emmap1,emmap2,float(c1),float(c1)+(emmap1.std()*0.05))
    #print 'Normal vector score: ', nv
    if nv < 0.0: 
      nv = 0.0
  except:
    print 'Exception for NV score'
    print_exc()
    nv = 0.0
  try: scores['normal_vector'].append(nv)
  except KeyError: scores['normal_vector'] = [nv]

for i in range(len(scores['local_correlation'])):
  if i == 0:   print '{:15} {:6} {:6} {:6} {:6} {:6} {:6} {:6}'.format('Model', 'LCCC', 'LMI', 'OVR', 'CCC', 'NMI', 'SD', 'NV')
  print '{:15} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f}'.format(list_names[i], scores['local_correlation'][i], scores['local_mutual_information'][i] \
                        , scores['overlap'][i], scores['correlation'][i], scores['normalised_mutual_information'][i], \
                        scores['surface_distance'][i], scores['normal_vector'][i])
