# Pocket focused deep generative model （PFDGM-kit）

PFDGM was proposed by integrating the three-dimensional (3D) structural information of the protein binding pocket into the conditional RNN (cRNN) model to control the generation of 9 drug-like molecules. In this model, the composition of the protein binding pocket is effectively characterized through a coarse-grain strategy and the 3D information of the pocket can be represented by the sorted eigenvalues of the Coulomb matrix (EGCM) of the coarse-grained atoms composing the binding pocket. It has been shown that the model trained with the constraint of protein environment information has a clear tendency 16 on generating compounds with higher similarity to the original X-ray-bound ligand than the normal RNN model and also better docking scores. Our results demonstrate the potential application of the controlled generative model for the targeted molecule 18 generation and guided exploration on the drug-like chemical space.
Here is the main achievement of EGCM controlled cRNN for Code for the purposes of [De Novo Molecule Design Through Molecular Generative Model Conditioned by 3D Information of Protein Binding Sites](https://chemrxiv.org/engage/chemrxiv/article-details/60c7537a702a9bbac118c383).

### Custom Dependencies
- [DeepDrugCoder (DDC)](https://github.com/pcko1/Deep-Drug-Coder) 
- [structconvert] (contained in Schrödinger)
- [Vina] (Vina docking program)

### File structure
- [DDC]: it's the source code of a condensed version of DeepDrugCoder, the installation tutorial can be found from here(https://github.com/pcko1/Deep-Drug-Coder) 
- [example/egcm_train_example_for_scpdb_pdbid_splitting]: An example of training the EGCM controlled CRNN Model
- [example/egcm_train_example_for_scpdb_pdbid_splitting/sample]: An example of sampling molecules with EGCM controlled CRNN model
- [scripts]: the main code of PFDGM-kit
- [baisc dataset]: it is only a part of sc-PDB dataset due the file size limitation of github.

### Main code
the main code of PFDGM-kit contains following files:
- [egcm_descriptor.py]: It is the main achievement of EGCM descriptor. 
- [add_descriptor.py]: Generate EGCM descriptor for dataset. It need the egcm_descriptor.py during calculation and an input control file (input.json).
- [input.json]: it is the setting of EGCM descriptor calculations. 
- [model_train.py]: it is the training script of EGCM controlled CRNN

- [sample_mol.py]: it is the sampling script of EGCM controlled CRNN
- [evaluate_similarity_cmd.py]: it is the script for evaluation of sampled molecules
- [vinda_docking.py]: it is the docking script for sampled molecules

### details:
- [step 1]: python -u add_descriptor.py -i INPUTDATASET -o OUTPUTDATASET -c CTRLINFO
- [step_2]: python -u model_train.py -i TRAININGSET -m MODELNAME -e EPOCHS --cuda GPU_ID
- [step_3]: python -u sample_mols.py --cuda GPU_ID -m MODEL_NAME -d DATASET_NAME -o OUTPUTDATASET --descriptor DESCRIPTOR --ifnoise IFNOISE --noisenum NOISENUM -n SAMPLENUM
- [step 4]: python -u evaluate_similarity_cmd.py -i INPUTDATASET -o OUTPUTNAME --descriptor DESCRIPTOR
- [step_5]: python -u vina_docking.py  -i INPUTDATASET -o OUTPUTNAME --mollistkey MOLLISTKEY --mode MODE -n DOCKNUM, --docknum DOCKNUM --nproc NPROC --boxsize BOXSIZE --nproc_per_job NPROC_PER_JOB
