

''' Preprocessing smiles with RDKIT can be a pain, especially working with
uncurated databases so I use this script to make sure all the SMILES can 
be read
'''

def sanitize_smiles(smi):
    #Return a canonical smile representation of smi
    try:
        mol = Chem.MolFromSmiles(smi, sanitize=True)
        smi_canon = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
        return (smi_canon)
    except:
        return (False)

def preprocess_smiles(smi_file, outfile=None):
    infile = pd.read_csv(smi_file, header=None)
    infile.columns = ["Smile", "ID"]
    infile.dropna(inplace=True)
    infile.drop_duplicates(keep=False,inplace=True)
    sanitized = pd.DataFrame(data=None, columns=["Smile", "ID"])
    for mol in infile.iloc:
        #print(mol['Smile'])
        smi_canon = sanitize_smiles(mol["Smile"])
        if smi_canon:
            row = pd.DataFrame({"Smile": [smi_canon], 'ID': [mol['ID']]})
            sanitized = pd.concat([sanitized, row])
    if outfile == None:
        outfile = smi_file[:-4]
    sanitized.to_csv(outfile+"_sanitized.smi", index=False, header=False)
    
preprocess_smiles("data/Kappa_Opioid.smi")
