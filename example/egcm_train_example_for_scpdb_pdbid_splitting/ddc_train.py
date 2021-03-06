#!/usr/bin/env python
# coding=utf-8
import numpy as np 
import rdkit
from rdkit import Chem
import h5py, ast, pickle
from ddc_pub import ddc_v3 as ddc
import os
os.environ["CUDA_VISIBLE_DEVICES"]="2"

dataset_name="../../basic_datasets/scPDB_pdbid_train_half_2021.pickle"

with open(dataset_name,'rb') as f:
    pocketlist=pickle.load(f)

mollist=[];egcmlist=[];molsmileslist=[]
for id,pocket in enumerate(pocketlist):
    try:
        if len(pocket["rdkitsmiles"])>0 and len(pocket["rdkitsmiles"])<300 and len(pocket["sortedegcm"])>0:
            mol=Chem.MolFromSmiles(pocket["rdkitsmiles"])
            if mol:
                mollist.append(mol)
                egcmlist.append(pocket["sortedegcm"])
                molsmileslist.append(pocket["rdkitsmiles"])
    except Exception as e:
        print (e,pocket["name"])

maxlen=np.max([len(molsmiles) for molsmiles in molsmileslist])
charset=''
for amol in molsmileslist:
    for achar in amol:
        if achar not in charset:
            charset+=achar
print ('Charset:',charset)
print (len(mollist),len(egcmlist),len(molsmileslist))

#print (PCBmodel.__dict__)
dataset_info={"maxlen":maxlen+20,"charset":charset,"name":"scPDB_2021_pdbid"}

x=egcmlist 
y=mollist
model=ddc.DDC(x=x,y=y,dataset_info=dataset_info,scaling=True)
model.fit(epochs=2000,lr=0.001,mini_epochs=2,model_name='scPDB_2021_pdbid',gpus=1,patience=1,checkpoint_dir="./model/",save_period=50,lr_decay=True,sch_epoch_to_start=400,sch_last_epoch=2000,sch_lr_init=1e-3,sch_lr_final=1e-6)

