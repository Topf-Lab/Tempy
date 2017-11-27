#===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#	  Copyright  2015 Birkbeck College University of London. 
#
#				Authors: Maya Topf, Daven Vasishtan, Arun Prasad Pandurangan,
#						Irene Farabella, Agnel-Praveen Joseph, Harpal Sahota
# 
#     This software is made available under GPL V3 license
#     http://www.gnu.org/licenses/gpl-3.0.html
#     
#     
#     Please cite your use of TEMPy in published work:
#     
#     Farabella, I., Vasishtan, D., Joseph, A.P., Pandurangan, A.P., Sahota, H. & Topf, M. (2015). J. Appl. Cryst. 48.
#
#===============================================================================

from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom
import TEMPy.Vector as Vector
import urllib
from numpy import append
import sys, os


class mmCIFParser:
    """A class to read mmCIF files either directly from the mmCIF or a structure instance from Biopython"""
    def __init__(self):
        pass

    @staticmethod
    def read_mmCIF_file(structure_id, filename,hetatm=False,water= False):
        """
        
        Read mmCIF file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of mmCIF file       
            *filename*
                name of mmCIF file
            *hetatm*
                Boolean representing whether the mmCIF file contains hetatom.
                Default and recommended is False.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.
        
        Return:
            Structure Instance
        """
        from Bio.PDB import MMCIFParser as MMCIFParserBiopy
        p=MMCIFParserBiopy()#permissive default True
        structure=p.get_structure(structure_id, filename)
        return mmCIFParser._biommCIF_strcuture_to_TEMpy(filename,structure,hetatm,water)

    @staticmethod
    def fetch_mmCIF(structure_id, filename,hetatm=False,water= False):
        
        """
        
        Fetch mmCIF file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of mmCIF file       
            *filename*
                name of mmCIF file
            *hetatm*
                Boolean representing whether the mmCIF file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.
        
        Return:
            Structure Instance
         """
        from Bio.PDB import MMCIFParser as MMCIFParserBiopy
        
        p=MMCIFParserBiopy()
        url = 'http://www.rcsb.org/pdb/files/%s.cif' % structure_id
        urllib.urlretrieve(url, filename)
        structure=p.get_structure(structure_id, filename)
        return mmCIFParser._biommCIF_strcuture_to_TEMpy(filename,structure,hetatm,water)

    @staticmethod
    def _biommCIF_strcuture_to_TEMpy(filename,structure,hetatm=False,water= False):
            #imported if and when the function is executed.
        """
        PRIVATE FUNCTION to convert to Structure Instance
        filename = name of mmCIF file
        hetatm = Boolean representing whether to add hetatm to the structure.Default and Raccomanded is False.
        water = Boolean representing whether to add water to the structure.Default and Raccomanded is False.
        """
        from Bio.PDB import MMCIFParser as MMCIFParserBiopy
        
        p=MMCIFParserBiopy()
        
        atomList = []
        hetatomList=[]
        wateratomList=[]
        footer = ''
        header = ''
        cif_code=filename.split("/")[-1]#use os.1FAT.cif
        structure_id="%s" % cif_code[:-4]
        structure=p.get_structure(structure_id, filename)
        residues = structure.get_residues()
        for res in residues:
            hetfield=res.get_id()[0]
            if hetfield[0]=="H":
                for atom in res:
                    BioPyAtom(atom)
                    hetatomList.append(BioPyAtom(atom))
            elif hetfield[0]=="W":
                for atom in res:
                    BioPyAtom(atom)
                    wateratomList.append(BioPyAtom(atom))
            else:
                for atom in res:
                    BioPyAtom(atom)
                    atomList.append(BioPyAtom(atom))
        if hetatm:
            atomList = append(atomList, hetatomList)
        if water:
            atomList = append(atomList, wateratomList)
        
        return BioPy_Structure(atomList, filename=filename, header=header, footer=footer)


class PDBParser:
    """A class to read PDB files either directly from the pdb or a structure instance from Biopython"""
    def __init__(self):
        pass

    @staticmethod
    def read_PDB_file(structure_id, filename,hetatm=False,water= False,chain=None):
        """
        
        Read PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
        
        p=PDBParserBiopy(QUIET=True)#permissive default True
#        try:
        structure=p.get_structure(structure_id, filename)
        #except AssertionError:
         #   sys.stderr.write('unknown element in BioPyton\n')
          #  sys.stderr.write(structure_id)
           # sys.exit()
        return PDBParser._bio_strcuture_to_TEMpy(filename,structure,hetatm,water)


    @staticmethod
    def read_PDB_file_BioPy(structure_id, filename,hetatm=False,water= False,chain=None):
        """
        
        Read PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
        
        p=PDBParserBiopy(QUIET=True)#permissive default True
#        try:
        structure=p.get_structure(structure_id, filename)
        #except AssertionError:
         #   sys.stderr.write('unknown element in BioPyton\n')
          #  sys.stderr.write(structure_id)
           # sys.exit()
        return structure
    
    @staticmethod
    def fetch_PDB(structure_id, filename,hetatm=False,water= False):       
        """
 
        Fetch PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
 
        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % structure_id
        p=PDBParserBiopy(QUIET=True)#permissive default True
        urllib.urlretrieve(url, filename)
        structure=p.get_structure(structure_id, filename)
        return PDBParser._bio_strcuture_to_TEMpy(filename,structure,hetatm,water)
        
    @staticmethod
    def fetch_PDB_BioPy(structure_id, filename,hetatm=False,water= False):       
        """
 
        Fetch PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
 
        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % structure_id
        p=PDBParserBiopy(QUIET=True)#permissive default True
        urllib.urlretrieve(url, filename)
        structure=p.get_structure(structure_id, filename,hetatm=hetatm,water=water)
        return structure

    @staticmethod
    def _bio_strcuture_to_TEMpy(filename,structure,hetatm=False,water= False):
            #imported if and when the function is executed.
        """
        PRIVATE FUNCTION to convert to Structure Instance
        filename = name of mmCIF file
        hetatm = Boolean representing whether to add hetatm to the structure.Default and Raccomanded is False.
        water = Boolean representing whether to add water to the structure.Default and Raccomanded is False.
        """
        #from Bio.PDB import PDBParser as PDBParserBiopy
        atomList = []
        hetatomList=[]
        wateratomList=[]
        footer = ''
        header = ''
        #pdb_code=filename.split("/")[-1]#use os.
        #p=PDBParserBiopy()#permissive default True
        #structure_id="%s" % pdb_code[:-4]
        #structure=p.get_structure(structure_id, filename)
        residues = structure.get_residues()
        for res in residues:
            hetfield=res.get_id()[0]
            if hetfield[0]=="H":
                for atom in res:
                    BioPyAtom(atom)
                    hetatomList.append(BioPyAtom(atom))
            elif hetfield[0]=="W":
                for atom in res:
                    BioPyAtom(atom)
                    wateratomList.append(BioPyAtom(atom))
            else:
                for atom in res:
                    BioPyAtom(atom)
                    atomList.append(BioPyAtom(atom))
        if hetatm:
            atomList = append(atomList, hetatomList)
        if water:
            atomList = append(atomList, wateratomList)
        
        return BioPy_Structure(atomList, filename=filename, header=header, footer=footer)
    
    @staticmethod
    def calc_SA(self,pdbfile,rsa=True,outsafile=None):
        assert os.path.isfile(pdbfile)
        if outsafile is None: outsafile = os.path.basename(pdbfile)+'_sa.out'
        #o = open(outsafile,'w')
        cmd = "~/data/packages/freesasa/freesasa-1.1/src/freesasa %s --rsa_file=%s\
         --no-log --radii=naccess"%(pdbfile,outsafile)
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
    
    @staticmethod
    def write_sasd_to_txt(sasds,pdb):

        """
        
        Outputs sasds to .txt file
        
        Arguments:
    
       *sasds*
           dictionary of sasds
       *pdb*
           .pdb file sasds were calculated on
        """

        if not os.path.exists('./Jwalk_results'):
            os.makedirs('./Jwalk_results')

        with open('./Jwalk_results/%s_crosslink_list.txt' % pdb[:-4],'w') as outf:

            outf.write(' '.join('{0:<13}'.format(col) for col in ['Index','Model','Atom1','Atom2','SASD','Euclidean Distance']))
            outf.write('\n')
            index = 1

            for xl in sasds:
                (aa1,chain1,res1)=xl[0]
                (aa2,chain2,res2)=xl[1]
                atom1 = ('%s-%d-%s-CA' % (res1,aa1,chain1) )
                atom2 = ('%s-%d-%s-CA' % (res2,aa2,chain2) )
                sasd=xl[2]
                ed=xl[3]
                outf.write(' '.join('{0:<13}'.format(col) for col in [index,pdb,atom1,atom2,sasd,ed]))
                outf.write('\n')
                index +=1

    @staticmethod
    def write_sasd_to_pdb(dens_map,sasds,pdb):
        """
            
        Outputs sasds to .pdb file
            
        Arguments:
           
           *dens_map*
               Solvent accessible surface on masked array
           *sasds*
               dictionary of sasds
           *pdb*
               .pdb file sasds were calculated on
        """
    
        if not os.path.exists('./Jwalk_results'):
            os.makedirs('./Jwalk_results')
    
        apix = dens_map.apix
        origin = dens_map.origin
        path_coord = {}
    
        for xl in sasds:
            a = []
            for (x,y,z) in sasds[xl]:
                a.append([(x*apix)+origin[0], (y*apix)+origin[1], (z*apix)+origin[2]])
            path_coord[xl] = a
    
        with open('./Jwalk_results/%s_crosslinks.pdb' % pdb[:-4],'w') as pdb:
    
            m_count = 1
            for xl in path_coord:
                (aa1,chain1,res1)=xl[0]
                (aa2,chain2,res2)=xl[1]
                sasd=xl[2]
                count = 1    
                pdb.write('MODEL %d %s%d%s-%s%d%s\n' % (m_count,res1,aa1,chain1,res2,aa2,chain2))
                m_count = m_count+1
                for (x,y,z) in path_coord[xl]:
                    p=Vector.Vector(x,y,z)
                    a=p.to_atom()
                    a.record_name = 'ATOM'
                    a.serial = count
                    a.atom_name = 'C'
                    a.alt_loc = ''
                    a.res = 'GLY'
                    a.chain = 'A'
                    a.res_no = count
                    a.icode = ''
                    a.occ = 1
                    a.temp_fac = 0
                    a.elem = 'C'
                    a.charge = ''
                    #print a.__dict__
                    #atom = BioPyAtom(a)
                    pdb.write(a.write_to_PDB())
                    count +=1
                pdb.write('END\n')
        

