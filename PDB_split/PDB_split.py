#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import time

def pdb_split_chain(input_dir, output_dir, ligand_list, ions_list):
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


def ligand_from_pdb(input_dir,output_dir, ligand_list):
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
            


def pdb_split_chain2(input_dir, output_dir, ligand_list, ions_list):
    for filename in os.listdir(input_dir):
        if filename.endswith('.pdb'):
            d=filename
            with open(os.path.join(input_dir, filename)) as g:
                file_lines=g.readlines()

                for lines in file_lines:
                    if lines.split()[0][0:6]=='HETATM' and lines.split()[3][0] =='A' and lines.split()[2][-3:]  in ligand_list:
                        print(filename+":Aa")

                        with open(os.path.join(input_dir, filename)) as g:
                            file_lines=g.readlines()
                            f= open(output_dir+filename,'a+')


                            for lines in file_lines:
                                if lines.split()[0]=='HEADER':
                                    f.writelines(lines)
                                if lines.split()[0]=='TITLE':
                                    f.writelines(lines)
                                if lines.split()[0]=='ATOM' and lines.split()[4][0]=='A':
                                    f.writelines(lines)
                                if lines.split()[0]=='TER' and lines.split()[3][0]=='A':
                                    f.writelines(lines)
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:] in ligand_list and lines.split()[3][0]=='A':
                                    f.writelines(lines)
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[3][0]=='A' and lines.split()[2] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break
                    elif lines.split()[0][0:6]=='HETATM' and lines.split()[3][0] =='B' and lines.split()[2][-3:] in ligand_list:
                        print(filename+":Bb")

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
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:] in ligand_list and lines.split()[3][0]=='B':
                                    f.writelines(lines)
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[3][0]=='B' and lines.split()[2] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break

                    elif lines.split()[0][0:6]=='HETATM' and lines.split()[3][0] =='C' and lines.split()[2][-3:] in ligand_list:
                        print(filename+":Cc")

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
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:] in ligand_list and lines.split()[3][0]=='C':
                                    f.writelines(lines)
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[3][0]=='C' and lines.split()[2] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break
                    elif lines.split()[0][0:6]=='HETATM' and lines.split()[3][0] =='D' and lines.split()[2][-3:] in ligand_list:
                        print(filename+":Dd")

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
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:] in ligand_list and lines.split()[3][0]=='C':
                                    f.writelines(lines)
                                if lines.split()[0][0:6]=='HETATM' and lines.split()[3][0]=='C' and lines.split()[2] in ions_list:
                                    f.writelines(lines)
                                if lines.split()[0]=='MASTER':
                                    f.writelines(lines)
                                if lines.split()[0]=='END':
                                    f.writelines(lines)

                            f.close() 
                            break
def ligand_from_pdb2(input_dir, output_dir, ligand_list):
    for filename in os.listdir(input_dir):
        if filename.endswith('.pdb'):
            d=filename.replace(".pdb","")

            with open(os.path.join(input_dir, filename)) as file:
                for lines in file.readlines():
                    if lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:]  in ligand_list:
                        for j in ligand_list:
                            if lines.split()[2][-3:]==j:
                                g=j

            with open(os.path.join(input_dir, filename)) as file:
                f= open(output_dir+d+"_"+g+".pdb", 'a+')
                f.writelines(f"pdb_id {d}\n")
                
                for lines in file.readlines():
                
                    
                    if lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:]  in ligand_list:
                        f.writelines(lines)
                
                f.writelines("END")
                f.close()

            with open(os.path.join(input_dir, filename)) as file:
                g=open(output_dir+d+".pdb", "a+")
                for lines in file.readlines():
                    if not (lines.split()[0][0:6]=='HETATM' and lines.split()[2][-3:] in ligand_list):

                        g.writelines(lines)
                g.close()
            
