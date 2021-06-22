#!/usr/bin/env python
# coding=utf-8
import numpy as np 
import rdkit
from rdkit import Chem
import h5py, ast, pickle
from ddc_pub import ddc_v3 as ddc
import os
from rdkit.Chem import rdmolfiles 
import argparse as arg
from tqdm import tqdm 

parser=arg.ArgumentParser(description='sample molecules')
parser.add_argument('-c','--cuda')
parser.add_argument('-m','--model_name')
parser.add_argument('-d','--dataset_name')
parser.add_argument('-o','--outputdataset')
parser.add_argument('--descriptor')
parser.add_argument('--ifnoise')
parser.add_argument('--noisenum')
parser.add_argument('-n','--samplenum')

args=parser.parse_args()
cuda=args.cuda
modelname=args.model_name
dataset_name=args.dataset_name
outputname=args.outputdataset
ifnoise=int(args.ifnoise)
descriptor=args.descriptor
noisenum=int(args.noisenum)
samplenum=int(args.samplenum)
os.environ["CUDA_VISIBLE_DEVICES"]=cuda
model=ddc.DDC(model_name=modelname)

with open(dataset_name,'rb') as f:
    pocketlist=pickle.load(f)

smileslist=[];mollist=[];
generatemollist=[]

model.batch_input_length=256

for id,complexstruc in tqdm(enumerate(pocketlist)):
    try:
        if complexstruc['rdkitsmiles'] and len(complexstruc[descriptor])>0:
            smiles=complexstruc['rdkitsmiles']
            if descriptor=='dtdescriptor':
                dp=np.array(complexstruc[descriptor])[0]
            else:
                dp=np.array(complexstruc[descriptor])
            dpsmileslist=[]
            try:
                noise=(np.random.rand(len(dp))*2-1)*noisenum/float(10)+1
                if ifnoise:
                    smiles_out,_=model.predict(latent=dp*noise,temp=0)
                else:
                    smiles_out,_=model.predict(latent=dp,temp=0)
                dpsmileslist.append(smiles_out)
            except Exception as e1:
                print (e1)
            try:
                noise=(np.random.rand(len(dp))*2-1)*noisenum/10+1
                if ifnoise:
                    smiles_out,_=model.predict_batch(latent=np.array([dp*noise]),temp=1)
                else:
                    smiles_out,_=model.predict_batch(latent=np.array([dp]),temp=1)
                dpsmileslist.append(smiles_out)
            except Exception as e3:
                print (descriptor,e3)
            generatemoldict={}
            generatemoldict["name"]=complexstruc["name"]
            generatemoldict["refmol"]=smiles
            generatemoldict[descriptor+'_mol_'+'%d'%ifnoise]=dpsmileslist
            generatemollist.append(generatemoldict)
    except Exception as e2:
        print (e2)
    #print (generatemollist)
    if id%50==0:
        with open(outputname,'wb') as f:
            pickle.dump(generatemollist,f) 
	
