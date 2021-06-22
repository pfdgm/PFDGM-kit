from pfdgm.base  import * 
import pickle 
import rdkit 
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, QED,rdmolfiles
from rdkit import Chem, DataStructs
import h5py
import random
import argparse as arg

parser=arg.ArgumentParser(description='compare dock result')
parser.add_argument('-i','--inputdataset')
parser.add_argument('-o','--outputname')
parser.add_argument('--mollistkey')
parser.add_argument('--mode')
parser.add_argument('-n','--docknum')
parser.add_argument('--nproc')
parser.add_argument('--boxsize')
parser.add_argument('--nproc_per_job')
args=parser.parse_args()

inputdataset=args.inputdataset
outputname=args.outputname
mollistkey=args.mollistkey
mode=args.mode
boxsize=int(args.boxsize)
nproc_per_job=int(args.nproc_per_job)

docknum=int(args.docknum)
nproc=int(args.nproc)
def dock_ligand(pocket,smileslist=[],PQueue=None,mode=''):
    #path='/mnt/home/myxu/Database/'+pocket['name'].strip('/home/myxu/Workspace/Database/')
    path=pocket['name']
    print (path)
    try:
        oriligpath=path+'/_ligand.pdb'
        orilig=ligand(pdbname=oriligpath)
        orilig.cal_masscenter()
        center=orilig.masscenter
        receptorpath=path+'/_protein.pdb'
        receptor=protein(pdbname=receptorpath)
        terminalcmd('/bin/bash -c "source ~/.bashrc && pythonsh $DOCKSCRIPS/prepare_receptor4.py -r %s.pdb -o %s.pdbqt && source activate TF12"'%(receptor.name,receptor.name))
    except Exception as e:
        print ('Error :',e)
    pocket['vina_affinity%s'%mode]=[]
    num=0
    for id,smiles in enumerate(smileslist):
        try:
            ligname='_Ligand%d%s'%(id,mode)
            ligpath=path+'/'+ligname
            rlig=rdmolfiles.MolFromSmiles(smiles)
            AllChem.EmbedMolecule(rlig, randomSeed=10)
            rdmolfiles.MolToPDBFile(rlig,path+'/%s.pdb'%ligname)
            terminalcmd('/bin/bash -c "source ~/.bashrc && pythonsh $DOCKSCRIPS/prepare_ligand4.py -l %s.pdb -o %s.pdbqt -A hydrogens&& source  activate TF12"'%(ligpath,ligpath))
            with open("%s/%s_r.conf"%(path,ligname),'w') as f:
                f.write("receptor = _protein.pdbqt\n")
                f.write("ligand   = %s.pdbqt\n"%ligname)
                f.write("center_x = %.3f\n"%center[0])
                f.write("center_y = %.3f\n"%center[1])
                f.write("center_z = %.3f\n"%center[2])
                f.write("size_x   = %d\n"%boxsize)
                f.write("size_y   = %d\n"%boxsize)
                f.write("size_z   = %d\n"%boxsize)
                f.write("energy_range = 4\n")
            terminalcmd('/bin/bash -c "cd %s && vina --config %s_r.conf --out dock_%s.pdbqt --log dock_%s.log --exhaustiveness %d >/dev/null && cd -"'%(path,ligname,ligname,ligname,nproc_per_job))
            with open("%s/dock_%s.log"%(path,ligname),"r") as f:
                affinitylist=[]
                line=f.readline()
                while line:
                    if '-----+----' in line:
                        for i in range(9):
                            line=f.readline()
                            if '  '+str(i+1)+' ' in line:
                                var=line.split()
                                affinity=var[1]
                                affinitylist.append(affinity)
                    line=f.readline()
            pocket['vina_affinity%s'%mode].append((smiles,affinitylist))
            num+=1
            if num>=5:
               break
        except Exception as e:
            print (pocket['vina_affinity%s'%mode])
    if PQueue:
        PQueue.put(pocket)
        return 
    else:
        return pocket

def receive_pocket(PQueue,filename):
    pocketlist=[]
    num=0
    while True:
        pocket=PQueue.get()
        if pocket:
            pocketlist.append(pocket)
        else:
            break
        num+=1
        print ('dock finished %d'%num)
        if num%100==0:
            with open(filename,'wb') as f:
                pickle.dump(pocketlist,f)
    with open(filename,'wb') as f:
        pickle.dump(pocketlist,f)
    return 

from multiprocessing import Pool,Queue,Manager,Process
from tqdm import tqdm
manager=Manager()
pocketQueue=manager.Queue()
p=Pool(nproc)
resultlist=[]
dataset_name="/home/myxu/Workspace/Database/ChemBL/CHEMBL25_TRAIN_MOLS.h5"
chemblset=h5py.File(dataset_name,'r')
sminame="./reinvent_sample.smi"
smilist=[]
with open(sminame,'r') as f:
    for line in f:
        smilist.append(line.strip('\n'))
mols=list(chemblset["mols"][:])
with open(inputdataset,'rb') as f:
    dockpocketlist=pickle.load(f)
print (len(dockpocketlist),dockpocketlist[0].keys())
    
for pocket in dockpocketlist:
    if ('rdkitsmiles' in pocket.keys() or (mode!='random' and mode!='original')):
        if mode=='random':
            smileslist=[Chem.MolToSmiles(Chem.Mol(mol)) for mol in random.sample(mols,docknum)]
            result=p.apply_async(dock_ligand,(pocket,smileslist,pocketQueue,mode))
        elif mode=='original':
            smileslist=[pocket['rdkitsmiles']]
            result=p.apply_async(dock_ligand,(pocket,smileslist,pocketQueue,mode))
        elif mode=='reinvent':
            smileslist=[smiles for smiles in random.sample(smilist,docknum)]
            result=p.apply_async(dock_ligand,(pocket,smileslist,pocketQueue,mode))
        else :
            smileslist=[pocket[mollistkey][0]]+list(pocket[mollistkey][1])#[:docknum]
            print (smileslist)
            #print (smileslist,docknum)
            result=p.apply_async(dock_ligand,(pocket,smileslist,pocketQueue,mollistkey))
        resultlist.append(result)
p.close()

receive_process=Process(target=receive_pocket,args=(pocketQueue,outputname))
receive_process.start()
for i in tqdm(range(len(resultlist))):
    tmp=resultlist[i].get()
    print (tmp)
p.terminate()
p.join()
pocketQueue.put(None)
receive_process.join()
print ("program finished")



