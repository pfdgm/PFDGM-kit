import time
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdmolfiles

from scipy.spatial import distance 
import numpy as np
import json
import pickle
import os
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, QED
from rdkit import Chem, DataStructs
from ddc_pub import ddc_v3 as ddc
from egcm_descriptor import * 
from tqdm import tqdm 
from multiprocessing import Pool,Queue,Manager,Process

def cal_descriptor(setting,com,Dqueue):
    try:
        filename=com['name']+'/site'
        tmp=os.popen('/bin/bash -c "structconvert -imol2 %s -osd %s > tmp"'%(filename+'.mol2',filename+'.sdf'))
        site=Chem.SDMolSupplier(filename+'.sdf')[0]
        pocketinfo=CalPocketCoulombMatrix(setting,site)
        descriptor=pocketinfo['sortedegcm']
        com['sortedegcm']=descriptor
    except Exception as e:
        print (com['name'],e)
        com['sortedegcm']=[]
    Dqueue.put(com)
    return 

def receive_descriptor(filename,DQueue):
    newlist=[]
    num=0
    while True:
        num+=1
        com=DQueue.get()
        if com:
            newlist.append(com)
        else:
            break
    with open(filename,'wb') as f:
        pickle.dump(newlist,f)
    return

import argparse as arg
parser=arg.ArgumentParser(description='Grep cativity from a complex pdb')
parser.add_argument('-i','--inputdataset')
parser.add_argument('-o','--outputdataset')
parser.add_argument('-c','--ctrlinfo')
args=parser.parse_args()
inputdataset=args.inputdataset
outputdataset=args.outputdataset
ctrlinfo=args.ctrlinfo

###############Load control parameters###################
setting=PARAMS()
setting.Loadjsonparm(ctrlinfo)
###############Load control parameters###################

with open(inputdataset,'rb') as f:
    pockets=pickle.load(f)
newpockets=[]

manager=Manager()
DQueue=manager.Queue()
p=Pool(28)
resultlist=[]
filename=outputdataset
for com in pockets:
    result =p.apply_async(cal_descriptor,(setting,com,DQueue))
    resultlist.append(result)
p.close()
receive_process=Process(target=receive_descriptor,args=(filename,DQueue))
receive_process.start()
for i in tqdm(range(len(resultlist))):
    tmp=resultlist[i].get()
    if tmp:
        print (tmp)
p.terminate()
p.join()
DQueue.put(None)
receive_process.join()
print ("program finished")

