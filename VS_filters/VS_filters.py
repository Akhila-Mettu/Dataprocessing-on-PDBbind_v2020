#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
import os
from rdkit import RDConfig
from rdkit.Chem.PandasTools import SaveXlsxFromFrame




import rdkit
from rdkit import Chem
def ro5_property_estimation(mol):  #Adjusted rule of five #filters 1_5
    molecular_weight = Descriptors.ExactMolWt(mol)
    n_HBA = Descriptors.NumHAcceptors(mol)
    n_HBD = Descriptors.NumHDonors(mol)
    Logp = Descriptors.MolLogP(mol)
    Rot_bonds=Descriptors.NumRotatableBonds(mol)
    ro5_conditions = [(molecular_weight >250 and molecular_weight < 600), n_HBA < 10, n_HBD < 5, (Logp>=1 and Logp <= 5), Rot_bonds < 14]
    ro5_fulfilled = sum(ro5_conditions) == 5
    return pd.Series(
        [molecular_weight, n_HBA, n_HBD, Logp, Rot_bonds, ro5_fulfilled],
        index=["Molecular_Weight", "n_HBA", "n_HBD", "LogP",'Rot_bonds', "ro5_fulfilled"],
    )



import seaborn as sns
def img_ro5(df):
    plt.figure(figsize=(10,5))
    sns.set_style("white")
    sns.histplot(df["Molecular_Weight"],bins=50)
    plt.savefig('Molecular_Weight.png')
    sns.histplot(df["n_HBA"],bins=30)
    plt.savefig('n_HBA.png')
    sns.histplot(df["n_HBD"],bins=20)
    plt.savefig('n_HBD.png')
    sns.histplot('LogP.png',bins=30)
    plt.savefig('n_HBD.png')
    sns.histplot(df15["Rot_bonds"],bins=20)
    plt.savefig('Rot_bonds.png')
   
    return Molecular_Weight.png, n_HBA.png, n_HBD.png,Rot_bonds.png




def mol_with_atom_index(mol): #function to understand atom index in a molecule
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol


# Chem.Draw.MolsToGridImage(df1['mol'].tolist()[0:5],legends=df1['ligand_name'].tolist()[0:5], molsPerRow=5) #To look into molecules of interest in a df



def P_containing_molecules(df):
    p=[]
    for m in df['mol']:
        atoms=[x for x in m.GetAtoms()]
        atom_num=[x.GetAtomicNum() for x in atoms]
        c=atom_num.count(15)
        p.append(c)
    df['P_containing']=p        
    return df




def aliphatic_amino_count(mol): #Aliphatic amino group count
    atoms=[x for x in mol.GetAtoms()] 
    ind=[x.GetIdx() for x in atoms]
    atom_num=[x.GetAtomicNum() for x in atoms]
    atom_hyb=[x.GetHybridization() for x in atoms]
    comb=list(zip(atom_num, atom_hyb))
    Natoms=[(i[1]) for _,i in enumerate(comb) if i[0]==7]
    count=len([x for x in Natoms if x == rdkit.Chem.rdchem.HybridizationType.SP3])
    comb=list(zip(ind,atom_num))
    a=[i for i,j in comb if j==7]
    b=[i for i,j in comb if j==16]
    c=[[i+1, i-1, i+2, i-2, i+3, i-3, i+4, i-4, i+5, i-5] for i,j in comb if j==7]
    m=[]
    n=[]
    for i in b:
        for j in c:
            for k in j:
                if i==k:
                    m.append(i)
                   
                else:
                    None
    count=max(count-len(m),0)
   
    return count





# df1['ali_N']=df1['mol'].apply(aliphatic_amino_count) #how to apply the above function



def car_acids(df):
    p=Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')
    y=[]
    for i in df['mol']:
        y.append(len(i.GetSubstructMatches(p)))
    df['carboxyl_group_count']=y
    return df



def aliphatic_atoms(mol): #Gives the generator of list of aliphatic chains if present in each compounds.
    rot_atom_pairs = list(mol.GetSubstructMatches(Chem.MolFromSmarts("[R0;D2]")))
    l=[list(x) for x in rot_atom_pairs]
    f=[(i[0]) for i in l]
    import itertools
    for i, j in itertools.groupby(enumerate(f), lambda x: x[1] - x[0]):
        j = list(j)
        start = j[0][1]
        length = len(j)
        if length == 1:
            yield (start-start)
        else:
            yield ((start+length)-start)

def connect_aa(df):#connects to df
    ali_atoms=[]
    for mol in df['mol']:
        ali_atoms.append(list(aliphatic_atoms(mol)))    
    
    return ali_atoms

def Aliphatic_c(ali_atoms): #Selects the longest sidechain amongst the list of sidechains.
    ali_c=[]
    for i in ali_atoms:
        if i and max(i):
            ali_c.append(max(i))
        else:
            ali_c.append(0)
    return ali_c

def aliphatic_atom_count(df):#overall program to count the longest aliphatic chain in a given molecule.
    result=connect_aa(df)
    ali_c= Aliphatic_c(list(result))
    df['Aliphatic_chain_len']=ali_c
    return df


def chiral_center_and_ringcount(df):
    y=[]
    z=[]
    for mol in df['mol']:
        chirals_c=len(Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=False,useLegacyImplementation=True))
        num_rings= len(Chem.GetSSSR(mol))
        y.append(chirals_c)
        z.append(num_rings)
    
    df['chiral_c']=y
    df['ring_c']=z      
    
    return df


def FusedBondCount(df):#number of rings in fused ring system #cannot give accurate count for larger fused ring system(fused ring size>8 atoms).
    import math
    RingFusedAtom = Chem.MolFromSmarts('[*R2]')
    matches=[]
    
    for mol in df['mol']:    
        matches.append(math.ceil(len(mol.GetSubstructMatches(RingFusedAtom))/2))
    
    df['fused_bond_count']=matches
    return df



def count_four_fused_rings(mol): #Gives the count of fused ring system >4 as 1 
    ri = mol.GetRingInfo()    
    d=[]  
    for ring in ri.AtomRings():
        d.append(set(ring))                
    total_rings=len(d)
    fused_ring_count=0
       
    from itertools import combinations
    f=[]
    for items in list(combinations(d, 4)):
        if items[0].intersection(items[1]):
            if items[1].intersection(items[2]):
                if items[2].intersection(items[3]):
                    f.append(items[0].union(items[1], items[2], items[3]))
                    fused_ring_count=+1
                        
    return  fused_ring_count

def four_fusedring_count(df): #applying the above function on df
    r=[]
    for mol in df['mol']:
        r.append(count_four_fused_rings(mol))
    df['count_four_fused_rings']=r
    
    return df
        


def protonate_ali_amino2(mol):
    atoms=[x for x in mol.GetAtoms()] 
    ind=[x.GetIdx() for x in atoms]
    atom_num=[x.GetAtomicNum() for x in atoms]
    atom_hyb=[x.GetHybridization() for x in atoms]
    comb=list(zip(ind,atom_num))   
    a=[i for i,j in comb if j==7] #N
    b=[i for i,j in comb if j==16] #S
    c=[[i+1, i-1, i+2, i-2, i+3, i-3, i+4, i-4, i+5, i-5] for i,j in comb if j==7] 
    m=[]
    for i in b:#index of s
        for j in c:#index of sourounding
            for k in j:
                m.append(i)
                
    for at in mol.GetAtoms(): 
        if at.GetAtomicNum() == 7 and at.GetHybridization()==rdkit.Chem.rdchem.HybridizationType.SP3 and at.GetFormalCharge()==0:
            if m in c: #if s in  list of sorrounding index
                at.SetFormalCharge(0)
    
        elif at.GetAtomicNum() == 7 and at.GetHybridization()==rdkit.Chem.rdchem.HybridizationType.SP3 and at.GetFormalCharge()==0:
            at.SetFormalCharge(1)
    
        
    return mol
        
       
def deprotonation_cooh(df):
    deprotonate_cooh  =  AllChem.ReactionFromSmarts("[C:1](=[O:2])-[OH1:3]>>[C:1](=[O:2])-[O-H0:3]")
    mols_deprot  =  []
    for  m  in  df['mol']:
        m_deprot  =  deprotonate_cooh.RunReactants((m,))
        mols_deprot.append(m_deprot[0][0]  if  m_deprot  else  m)
    df['new_mol']=mols_deprot
    return df

