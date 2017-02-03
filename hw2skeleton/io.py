import glob
import os
from .utils import Atom, Residue, ActiveSite
import numpy as np


def read_active_sites(dir):
    """
    Read in all of the active sites from the given directory.

    Input: directory
    Output: list of ActiveSite instances
    """
    files = glob.glob(dir + '/*.pdb')

    active_sites = []
    # iterate over each .pdb file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.pdb")):

        active_sites.append(read_active_site(filepath))

    print("Read in %d active sites"%len(active_sites))

    return active_sites


def read_active_site(filepath):
    """
    Read in a single active site given a PDB file

    Input: PDB file path
    Output: ActiveSite instance
    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".pdb":
        raise IOError("%s is not a PDB file"%filepath)

    active_site = ActiveSite(name[0])

    r_num = 0

    # open pdb file
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in f:
            if line[0:3] != 'TER':
                # read in an atom
                atom_type = line[13:17].strip()
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                atom = Atom(atom_type)
                atom.coords = (x_coord, y_coord, z_coord)

                residue_type = line[17:20] # amino acid
                residue_number = int(line[23:26])

                # make a new residue if needed
                if residue_number != r_num:
                    residue = Residue(residue_type, residue_number)
                    r_num = residue_number

                # add the atom to the residue
                residue.atoms.append(atom)

            else:  # I've reached a TER card
                active_site.residues.append(residue)
                
        # Categorize residues in an active site by biochemical properties, ouput a dictionary where keys are amino acids and values are the count of that aa
        Aresidues = {}
        aminoacids = ['ASP','GLU','HIS','ARG','LYS','SER','THR','ASN','GLN','CYS','PHE','TRP','TYR','GLY','ALA','VAL','LEU','ILE','MET','PRO']
        for aa in aminoacids:
            Aresidues[aa] = 0 # initialize Aresidues dictionary with aa's as keys and all values set to 0
        for res in active_site.residues:
            if res.type not in Aresidues:
                Aresidues[res.type] = 1
            else:
                Aresidues[res.type] += 1
        totalresidues = sum(Aresidues.values()) # total residues in active site


        # Categorize residues by biochemical properties, store in list where values are counts of aa's in that category. 
        #Arescategories = [acidic, basic, hydrophobic, nonpolar, polar] fractions
        Arescategories = np.array([0.0]*5)
        Arescategories[0] = Aresidues['ASP'] + Aresidues['GLU'] #'acidic'
        Arescategories[1] = Aresidues['HIS'] + Aresidues['ARG'] + Aresidues['LYS'] #'basic'
        Arescategories[2] = Aresidues['PHE'] + Aresidues['TRP'] + Aresidues['TYR'] #'hydrophobic'
        Arescategories[3] = Aresidues['GLY'] + Aresidues['ALA'] + Aresidues['VAL'] + Aresidues['LEU'] + Aresidues['ILE'] + Aresidues['MET'] + Aresidues['PRO'] #'nonpolar'
        Arescategories[4] = Aresidues['SER'] + Aresidues['THR'] + Aresidues['ASN'] + Aresidues['GLN'] + Aresidues['CYS'] #'polar'
       
        # Normalize by total number of residues in the active site so that each category will represent the fraction of amino acids in the active site in that category (e.g. 0.5 polar, 0.5 basic)
        Arescategories = Arescategories/totalresidues
            
        active_site.categories = Arescategories # add categories to the active site        
        return active_site
        
#%% 
         # Categorize residues by biochemical properties, store in dictionary where keys are categories and values are counts of aa's in that category
#        Arescategories = {}
#        Arescategories['acidic'] = Aresidues['ASP'] + Aresidues['GLU']
#        Arescategories['basic'] = Aresidues['HIS'] + Aresidues['ARG'] + Aresidues['LYS']
#        Arescategories['polar'] = Aresidues['SER'] + Aresidues['THR'] + Aresidues['ASN'] + Aresidues['GLN'] + Aresidues['CYS']
#        Arescategories['hydrophobic'] = Aresidues['PHE'] + Aresidues['TRP'] + Aresidues['TYR']
#        Arescategories['nonpolar'] = Aresidues['GLY'] + Aresidues['ALA'] + Aresidues['VAL'] + Aresidues['LEU'] + Aresidues['ILE'] + Aresidues['MET'] + Aresidues['PRO']
#        
#        # Normalize by total number of residues in the active site so that each category will represent the fraction of amino acids in the active site in that category (e.g. 0.5 polar, 0.5 basic)
#        for category in Arescategories.keys():
#            Arescategories[category] = (Arescategories[category])/totalresidues
#            
#        active_site.categories.update(Arescategories) # add categories to the active site
        #print(active_site.categories)

#    return active_site
#%%

def write_clustering(filename, clusters):
    """
    Write the clustered ActiveSite instances out to a file.

    Input: a filename and a clustering of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusters)):
        out.write("\nCluster %d\n--------------\n" % i)
        for j in range(len(clusters[i].points)):
            out.write("%s\n" % clusters[i].points[j])
           # out.write("%s\n" % clusters[i].points[j].residues)
            out.write("%s\n" % clusters[i].points[j].categories)
    out.close()


def write_mult_clusterings(filename, clusterings):
    """
    Write a series of clusterings of ActiveSite instances out to a file.

    Input: a filename and a list of clusterings of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusterings)):
        clusters = clusterings[i]

        for j in range(len(clusters)):
            out.write("\nCluster %d\n------------\n" % j)
            for k in range(len(clusters[j])):
                out.write("%s\n" % clusters[j][k])

    out.close()
