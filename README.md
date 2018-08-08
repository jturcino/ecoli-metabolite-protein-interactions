# ecoli-metabolite-protein-interactions
[Emili lab](http://www.emililab.org/) and [Vajda lab](https://structure.bu.edu/) collaborative project examining metabolite-protein interactions in *E. coli*

## Analysis Agenda
Preliminary analysis plan (subject to change)

### 1. Metabolites
Identify as many metabolites as possible and generate 3D structures

- Pull metabolite SMILES, InChI, and InChIKey from [EcoCyc](https://ecocyc.org/) and [HMDB](http://www.hmdb.ca/)
- Find compounds with only one name and/or agreement between EcoCyc and HMDB
- **[Pending]** Identify compounds with experimental standards (need list)
- Generate `.pdb` files for compounds' 3D structure using SMILES, InChI, and/or InChIKey via [OpenBabel](http://openbabel.org/wiki/Main_Page) or [ROCS](https://www.eyesopen.com/rocs)

### 2. Proteins
Use representative PDB ID matches to model proteins

- **[[Done]](http://www.bioinf.org.uk/pdbsws/index.html)** Look up PDB IDs from UniProt IDs
- Pick representatives from PDB ID matches (e.g. *E. coli* protein only) for mapping. Criteria for selection TBD; considering docking.
- Run representative PDB ID matches in [FTMap](http://ftmap.bu.edu/login.php) utilizing script submission

### 3. Metabolite-Protein Interactions
Identify experimentally verified interactions

- Meet with [Segrè lab](https://www.bu.edu/segrelab/)
- Look up databases and other resources for metabolite-protein interactions. Dr. Vajda will do some manual validation if necessary.

### 4. Metabolite-Protein Docking
With Amanda and Marcello!

- Select validation metabolite-protein dockings to run on different docking servers
- Select best small-molecule docker, based on performance in the validation set
