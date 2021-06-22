#!/usr/bin/env python
# coding: utf-8
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
fragmentdict={'COO':'C(O)=O','NHCO':'O=CN','ARG':'NC(N)=N','Methoxy':'C1=CN=CN1',\
              'TRP':'C1=CNC2=C1C=CC=C2','TYR':'C1=CC=C(O)C=C1','PHE':'C1=CC=CC=C1',\
              'LYS':'CN','MET':'CSC','CYS':'CS','OH':'CO'}

fragidxdict={}
for id,i in enumerate(fragmentdict.keys()):
    fragidxdict[i]=(id+1)

class Groupsetting:
    def __init__(self):
        self.GroupList=[]
        self.GroupNumCtrl=[]
        self.Dim=0
    def update(self):
        self.Dim=np.sum(self.GroupNumCtrl)
        self.GroupPT=np.zeros((len(self.GroupList),2),dtype=int)
        self.Groupidx=np.zeros(self.Dim)
        self.GroupPT[0][0]=0;self.GroupPT[0][1]=self.GroupNumCtrl[0]-1
        for j in range(self.GroupPT[0][0],self.GroupPT[0][1]+1):
            self.Groupidx[j]=fragidxdict[self.GroupList[0]]
        for i in range(1,len(self.GroupList)):
            self.GroupPT[i][0]=self.GroupPT[i-1][1]+1
            self.GroupPT[i][1]=self.GroupPT[i][0]+self.GroupNumCtrl[i]-1
            for j in range(self.GroupPT[i][0],self.GroupPT[i][1]+1):
                self.Groupidx[j]=fragidxdict[self.GroupList[i]]
        return 
        
class PARAMS:
    def __init__(self):
        self.Groupsetting=Groupsetting()

    def Loadjsonparm(self,jsonfile):
        with open(jsonfile,'r') as f:
            jsondict=json.load(f)
            if 'Groupsetting' in jsondict.keys():
                self.Groupsetting=Groupsetting()
                Loaddict2obj(jsondict['Groupsetting'],self.Groupsetting)
                self.Groupsetting.update()
def Loaddict2obj(dict,obj):
    objdict=obj.__dict__
    for i in dict.keys():
        if i not in objdict.keys():
            print ("%s is not a standard setting option!"%i)
        objdict[i]=dict[i]
    obj.__dict__=objdict

    
def CalPocketCoulombMatrix(setting,sitemol,name=''):
    coulombmatrix=np.zeros((setting.Groupsetting.Dim,setting.Groupsetting.Dim),dtype=float)
    for i in range(setting.Groupsetting.Dim):
        coulombmatrix[i][i]=0.5*setting.Groupsetting.Groupidx[i]**(2.4)
    outlist=[]
    siteenv={}
    sitecrdenv={}
    for i in setting.Groupsetting.GroupList:
        fraglist=[]
        fragsmiles=rdmolfiles.MolFromSmiles(fragmentdict[i])
        matchfrag=sitemol.GetSubstructMatches(fragsmiles)
        for atruple in matchfrag: 
            if all(j not in outlist for j in atruple):
                outlist+=list(atruple)
                fraglist.append(atruple)
        siteenv[i]=fraglist
    sitefragnumdict={}
    for i in siteenv.keys():
        sitefragnumdict[i]=len(siteenv[i])

    sitecoord=sitemol.GetConformer().GetPositions()

    for i in siteenv.keys():
        centerlist=[]
        tmpfraglist=siteenv[i]
        for atomlist in tmpfraglist:
            atomcoord=[sitecoord[atomlist[j]] for j in range(len(atomlist))]
            atommass=[sitemol.GetAtomWithIdx(atomlist[j]).GetMass() for j in range(len(atomlist))]
            masscenter=np.sum(np.array(atomcoord)*np.array(atommass).reshape(-1,1),axis=0)/np.sum(atommass)
            centerlist.append(masscenter)
        sitecrdenv[i]=centerlist

    for i in range(len(setting.Groupsetting.GroupList)):
        alist=siteenv[setting.Groupsetting.GroupList[i]]
        for j in range(len(setting.Groupsetting.GroupList)):
            blist=siteenv[setting.Groupsetting.GroupList[j]]
            for m in range(len(alist)):
                pt1=setting.Groupsetting.GroupPT[i][0]+m
                crd1=sitecrdenv[setting.Groupsetting.GroupList[i]][m]
                for n in range(len(blist)):
                    pt2=setting.Groupsetting.GroupPT[j][0]+n
                    crd2=sitecrdenv[setting.Groupsetting.GroupList[j]][n]
                    if pt1!=pt2:
                        coulombmatrix[pt1][pt2]=setting.Groupsetting.Groupidx[pt1]*setting.Groupsetting.Groupidx[pt2]\
                                                /np.sqrt(np.sum((crd1-crd2)**2))
    EGCM,EVEC=np.linalg.eig(coulombmatrix)
    EGCM=-np.sort(-EGCM)
    pocketdict={'name':name,'pocket':sitemol,'coulombmatrix':coulombmatrix,'sortedegcm':EGCM,'groupnum':sitefragnumdict}
    return pocketdict

def pcb_descriptor(mol,qsar_model=None):
    descriptors=[]
    try:
        logp  = Descriptors.MolLogP(mol)
        tpsa  = Descriptors.TPSA(mol)
        molwt = Descriptors.ExactMolWt(mol)
        hba   = rdMolDescriptors.CalcNumHBA(mol)
        hbd   = rdMolDescriptors.CalcNumHBD(mol)
        qed   = QED.qed(mol)
        # Calculate fingerprints
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=2048)
        ecfp4 = np.zeros((2048,))
        DataStructs.ConvertToNumpyArray(fp, ecfp4)
        # Predict activity and pick only the second component
        active = qsar_model.predict_proba([ecfp4])[0][1]
        descriptors=np.array([logp, tpsa, molwt, qed, hba, hbd, active])
    except Exception as e:
        print (e)
    return descriptors

