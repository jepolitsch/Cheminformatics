def gen3D(smi_file, outfile, conformations, max_iter=200):
    time = 0
    energy = {}
    supplier = Chem.SmilesMolSupplier(smi_file, delimiter=',',smilesColumn=0,
                                      nameColumn=1, titleLine=False)
    sdfWriter = Chem.SDWriter(str(outfile))
    for i in tqdm(range(len(supplier))):
        try:
            mol = supplier[i]
            if (mol.HasProp("_Name")):
                molName = mol.GetProp("_Name")
            molH = AllChem.AddHs(mol, addCoords=True)
            cids = AllChem.EmbedMultipleConfs(molH, numThreads=0, numConfs=int(conformations), randomSeed=51254)
            mp = AllChem.MMFFGetMoleculeProperties(molH)
            for cid in cids:
                field = AllChem.MMFFGetMoleculeForceField(molH, mp, confId = cid,)
                field.Minimize(maxIts=max_iter)
                e = field.CalcEnergy()
                molH.SetProp("_Name", (str(molName)+'.'+str(cid)))
                molH.SetProp('Energy', str(e))
                if silent == False:
                    print('{0:}\t\t{1:.3f}\t\t{2:}'.format(molName, e, cid))
                sdfWriter.write(molH, confId = cid)
                energy.update({(molName+"."+str(cid)): e})
        except:
            pass
        time += i
    print("SDF file written for {}  in {}".format(smi_file, outfile))
    sdfWriter.close()
    return(outfile, energy)
