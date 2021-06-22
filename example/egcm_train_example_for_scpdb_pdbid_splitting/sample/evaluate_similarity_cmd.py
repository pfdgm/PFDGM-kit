#!/usr/bin/env python
# coding=utf-8
import pickle 
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem,rdmolfiles,Draw
import argparse as arg
import os
import tensorflow as tf
def show_similarity(refmol,mollist,imgname=''):
    if len(mollist)>0:
        mollist=[mol for mol in mollist if mol is not None]
        cononical_smiles=[rdmolfiles.MolToSmiles(mol,canonical=True) for mol in mollist if mol is not None]
        reffp=AllChem.GetMorganFingerprint(refmol,2)
        fplist=[AllChem.GetMorganFingerprint(mol,2) for mol in mollist if mol is not None]
        similaritylist=np.array([DataStructs.TanimotoSimilarity(reffp,fp) for fp in fplist])
        simi_mat=np.array([[DataStructs.TanimotoSimilarity(fplist[i],fplist[j])  for i in range(len(fplist)) if i!=j] for j in range(len(fplist))])
        unique=float(len(list(set(cononical_smiles)))/len(cononical_smiles))
        unlike=np.min(simi_mat)
        average_simi=np.mean(simi_mat)
        order=np.argsort(-similaritylist)
        order_s=[similaritylist[i] for i in order]
        order_m=[mollist[i] for i in order]
        
        #if imgname!='' and np.max(similaritylist)>0.7 and len(set(order_s))>10:
        #    try:
        #        image=Draw.MolsToGridImage([refmol]+order_m,molsPerRow=5,subImgSize=(500,500),legends=['refmol']+['%.2f'%sim for sim in order_s])
        #        image.save(imgname)
        #    except:
        #        pass
        return similaritylist,unique,unlike,average_simi
    else:
        return []

def smileslisttomol(smileslist,ifreturnsmileslist=True):
    mollist=[]
    smileslist2=[]
    for smiles in smileslist:
        try:
            mol=rdmolfiles.MolFromSmiles(smiles)
            if mol is not None:
                mollist.append(mol)
                smileslist2.append(smiles)
        except Exception as e:
            print (smiles)
    try:
        validity=float(len(smileslist2)/len(smileslist))    
    except:
        validity=0
    return mollist,smileslist2,validity

parser=arg.ArgumentParser(description='Grep cativity from a complex pdb')
parser.add_argument('-i','--inputdataset')
parser.add_argument('-o','--outputname')
parser.add_argument('--descriptor')
args=parser.parse_args()
inputdataset=args.inputdataset
outputname=args.outputname
descriptor=args.descriptor


with open(inputdataset,'rb') as f:
    generatemolslist=pickle.load(f)
print (generatemolslist[0])
genmols=[]
similaritylist=[];validitylist=[]
#refmollist=[]
#refsmileslist=[]

if not os.path.exists('./similarity/%s'%descriptor):
    os.system('mkdir -p ./similarity/%s'%descriptor)

egcm_pocket_similaritylist=[]
dt_pocket_similaritylist=[]
egcmdislist=[]
dtdislist=[]
uniquelist=[]
unlikelist=[]
average_similist=[]
num=0
similaritydic={}
for id,moldict in enumerate(generatemolslist):
    try:
        print ('Deal with %d-----'%id)
        refmol=rdmolfiles.MolFromSmiles(moldict['refmol'])
        gmollist,tmpsmilelist,validity=smileslisttomol(moldict[descriptor+'_mol_0'][1])
        similarity,unique,unlike,average_simi=show_similarity(refmol,gmollist,imgname='./similarity/%s/%s.png'%(descriptor,moldict["name"].split('/')[-1]))
        similaritylist.append(np.max(similarity))
        uniquelist.append(unique)
        unlikelist.append(unlike)
        average_similist.append(average_simi)
        validitylist.append(validity)
    except:
        pass
np.savetxt(f'{outputname}.similarity',similaritylist,fmt='%.2f',header=outputname)
np.savetxt(f'{outputname}.validity',validitylist,fmt='%.2f',header=outputname)
np.savetxt(f'{outputname}.unique',uniquelist,fmt='%.2f',header=outputname)
np.savetxt(f'{outputname}.unlike',unlikelist,fmt='%.2f',header=outputname)
np.savetxt(f'{outputname}.average_similarity',average_similist,fmt='%.2f',header=outputname)

    
