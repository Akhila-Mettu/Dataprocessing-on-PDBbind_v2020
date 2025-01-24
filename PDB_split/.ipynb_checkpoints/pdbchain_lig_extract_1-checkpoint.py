#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import time


# In[2]:


ligand_list=['BBY', 'UFV', '5RL', '5WD', 'ZT2', 'ZT4', 'CLQ', 'TPR', 'E7E', 'JJN', '82J', '09I', 'A7M', 'SB2','35I', '0JN', '3JG', 'HYS','GM5','ACV']


# In[3]:


len(ligand_list)


# In[4]:


ions_list=["ZN","MN", "CD", "NI", "NA", "MG", "ARS", "K", "CA", "HG"]


# In[5]:


def pdb_split_chain(input_dir, output_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith('.pdb'):
            d=filename
            with open(os.path.join(input_dir, filename)) as g:
                file_lines=g.readlines()

                for lines in file_lines:
                    if lines.split()[0]=='HETATM'  and (lines.split()[4][0]=='A' or lines.split()[4][0]=='X' or lines.split()[4][0]=='U' or lines.split()[4][0]=='H') and (lines.split()[3] in ligand_list):
                        print(filename+":A")

                        with open(os.path.join(input_dir, filename)) as g:
                            file_lines=g.readlines()
                            f= open(output_dir+filename,'a+')

                            for lines in file_lines:    

                                if lines.split()[0]=='HEADER':
                                    f.writelines(lines)
                                if lines.split()[0]=='TITLE':
                                    f.writelines(lines)                
                                if lines.split()[0]=='ATOM' and (lines.split()[4][0]=='A' or lines.split()[4][0]=='X'or lines.split()[4][0]=='U' or lines.split()[4][0]=='H'):
                                    f.writelines(lines)
                                if lines.split()[0]=='TER' and (lines.split()[3][0]=='A' or lines.split()[3][0]=='X'or lines.split()[3][0]=='U' or lines.split()[3][0]=='H'):
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and (lines.split()[3]  in ligand_list) and (lines.split()[4][0]=='A' or lines.split()[4][0]=='X' or lines.split()[4][0]=='U' or lines.split()[4][0]=='H'):
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and (lines.split()[4][0]=='A'or lines.split()[4][0]=='X' or lines.split()[4][0]=='U' or lines.split()[4][0]=='H') and lines.split()[3]  in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)


                            f.close()
                            break


                    elif lines.split()[0]=='HETATM' and lines.split()[4][0]=='B' and lines.split()[3] in ligand_list:
                        print(filename+":B")

                        with open(os.path.join(input_dir, filename)) as g:
                            file_lines=g.readlines()
                            f= open(output_dir+filename,'a+')


                            for lines in file_lines:
                                if lines.split()[0]=='HEADER':
                                    f.writelines(lines)
                                if lines.split()[0]=='TITLE':
                                    f.writelines(lines)
                                if lines.split()[0]=='ATOM' and lines.split()[4][0]=='B':
                                    f.writelines(lines)
                                if lines.split()[0]=='TER' and lines.split()[3][0]=='B':
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and lines.split()[3] in ligand_list and lines.split()[4][0]=='B':
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and lines.split()[4][0]=='B' and lines.split()[3] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break

                    elif lines.split()[0]=='HETATM' and lines.split()[4][0]=='C' and lines.split()[3] in ligand_list:
                        print(filename+":C")

                        with open(os.path.join(input_dir, filename)) as g:
                            file_lines=g.readlines()
                            f= open(output_dir+filename,'a+')


                            for lines in file_lines:
                                if lines.split()[0]=='HEADER':
                                    f.writelines(lines)
                                if lines.split()[0]=='TITLE':
                                    f.writelines(lines)
                                if lines.split()[0]=='ATOM' and lines.split()[4][0]=='C':
                                    f.writelines(lines)
                                if lines.split()[0]=='TER' and lines.split()[3][0]=='C':
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and lines.split()[3] in ligand_list and lines.split()[4][0]=='C':
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and lines.split()[4][0]=='C' and lines.split()[3] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break

                    elif lines.split()[0]=='HETATM' and lines.split()[4][0]=='D' and lines.split()[3] in ligand_list:
                        print(filename+":D")

                        with open(os.path.join(input_dir, filename)) as g:
                            file_lines=g.readlines()
                            f= open(output_dir+filename,'a+')


                            for lines in file_lines:
                                if lines.split()[0]=='HEADER':
                                    f.writelines(lines)
                                if lines.split()[0]=='TITLE':
                                    f.writelines(lines)
                                if lines.split()[0]=='ATOM' and lines.split()[4][0]=='D':
                                    f.writelines(lines)
                                if lines.split()[0]=='TER' and lines.split()[3][0]=='D':
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and lines.split()[3] in ligand_list and lines.split()[4][0]=='D':
                                    f.writelines(lines)
                                if lines.split()[0]=='HETATM' and lines.split()[4][0]=='D' and lines.split()[3] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break





# In[6]:


#20 psbs in seven seconds
t1=time.time()
pdb_split_chain(input_dir='./pdbs/', output_dir='./results_1/')
t2=time.time()
print(t2-t1)


# In[7]:


def ligand_from_pdb(input_dir,output_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith('.pdb'):
            d=filename.replace(".pdb","")

            with open(os.path.join(input_dir, filename)) as file:
                for lines in file.readlines():
                    if lines.split()[0]=='HETATM' and lines.split()[3]  in ligand_list:
                        for j in ligand_list:
                            if lines.split()[3]==j:
                                g=j

            with open(os.path.join(input_dir, filename)) as file:
                f= open(output_dir+d+"_"+g+".pdb", 'a+')
                f.writelines(f"pdb_id {d}\n")
                for lines in file.readlines():
                    if lines.split()[0]=='HETATM' and lines.split()[3]  in ligand_list:
                        f.writelines(lines)
                
                f.writelines("END")
                f.close()

            with open(os.path.join(input_dir, filename)) as file:
                g=open(output_dir+d+".pdb", "a+")
                for lines in file.readlines():
                    if not (lines.split()[0]=='HETATM' and lines.split()[3] in ligand_list):

                        g.writelines(lines)
                g.close()
            







# In[8]:


t1=time.time()
ligand_from_pdb(input_dir='./results_1/', output_dir='./results_3/')
t2=time.time()
print(t2-t1)


# In[ ]:





# In[ ]:




