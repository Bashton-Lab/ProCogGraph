# ProCogGraph
Graph based ligand-domain interaction database for exploring and mining domain-cognate ligand interactions. Powered by PDBe-graph and Neo4j.

## To Do

| Task | Status | Priority |
| ---- | ------ | -------- |
| Obtain biological ligand names from database source | not started | 2 |
| Find new way to combine duplicate biological ligands from different databases | not started | 3 |
| Add new biological ligands to database | not started | 5 |
| Write script to process PDB files into database | not started | 4 |
| Write pipline to annotated CDS from metagenomic libraries with cognate ligand binding profiles | not started | 1 |

## PDBe Graph Data

ProCogGraph utilises PDBe graph data to obtain the following information:

* Protein chains: EC number annotation
* Protein chains: SCOP domain annotation
* Protein chains: CATH domain annotation
* Protein chains: InterPro domain annotation
* Bound Molecules: Component entities

## Biological Ligands

### Sources

In ProCogGraph, we collect biological ligands from a variety of database sources including:

* [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
* [KEGG](https://www.genome.jp/kegg/)
* [ChEBI](https://www.ebi.ac.uk/chebi/)
* [Rhea](https://www.rhea-db.org/)
* [GlyTouCan](https://glytoucan.org/)

Potential databases which will be included in future versions include:

* [PDB Ligand Expo](https://ligand-expo.rcsb.org/ld-download.html)
* [UniProt](https://www.uniprot.org/)
* [Reactome](https://reactome.org/)
* [MetaCyc](https://metacyc.org/)
* [LIPID MAPS](https://www.lipidmaps.org/)
* [HMDB](http://www.hmdb.ca/)
* [BindingDB](https://www.bindingdb.org/bind/index.jsp)
* [PDB](https://www.rcsb.org/)
* [DrugBank](https://www.drugbank.ca/)
* [ChEMBL](https://www.ebi.ac.uk/chembl/)

### Obtaining Biological Ligands

The biological ligands are obtained from the databases listed above using the following steps:

1. Process the enzyme.data from ExPasY to obtain the up to data EC numbers list.
2. Search for KEGG enzyme record matching EC number. If found, extract the EC substrate codes, EC product codes, EC dbxrefs and reactions IDs.
3. For each reaction ID, search for reaction record and extract the substrate codes and product codes.
4. Combine the compound codes (from EC record and reaction record) to obtain a set of compoound IDs for each EC.
5. Search for compound ID records from KEGG and where possible obtain a SMILES string.
6. Search for compound ID records in ChEBI and obtain SMILES string where possible.
7. Search for compound ID records in PubChem (by CID) and obtain SMILES string where possible.
8. Search for glycan IDs in GlyTouCan and obtain SMILES string where possible.
9. Search for EC records in Rhea and obtain reaction SMILES string where possible. Split into component SMILES strings.
10. Combine into a single biological ligands dataframe.

## Defining Biological Ligand Similarity

### Existing methods of biological ligand similarity

In PROCOGNATE biological ligand similarity is defined as ... EXPLAIN HERE

Subsequently, the PARITY method was developed in 2018 by the Thornton group. EXPLAIN HERE

### ProCogGraph method of biological ligand similarity

In ProCogGraph, we define biological ligand similarity using the PARITY method. Various cutoffs for biological similarity are defined for this score in the literature: 2018 paper uses 0.7, subsequent 2020 paper uses as low as 0.3.

To determine the appropriate cutoff score for ProCogGraph, a receiver operating characteristic (ROC) curve was generated and the score threshold determined using Youden's J statistic. The ROC curve was generated using the following steps:

1. A manually curated set of "cognate" ligands were identified which were already bound to a PDB protein chain.
2. These were scored against their corresponding cognate ligands using the PARITY score. These constitute the positive examples in the ROC curve.
3. PDB ligands were then randomly matched to biological ligands, and scored using the PARITY score. These constitute the negative examples in the ROC curve.
4. The ROC curve was generated using the scikit-learn python package.
5. The score threshold was determined using Youden's J statistic.

The distribution of positive and negative scores is shown in the histogram below:

![Score Distribution]()

The ROC curve is shown below:

![ROC Curve]()

The score threshold was determined to be 0.55 (Youden's index = 0.92).

Potential other avenues for scoring to be explored: the Tanimoto similarity between the Morgan fingerprints of the ligands.

## Domain Ownership

### Existing methods of domain ownership

In PROCOGNATE domain ownership is defined as ...

domain which is most likely to be responsible for the binding of a given ligand. We define this as the domain which has the highest similarity to the ligand. We define similarity as the Tanimoto similarity between the ligand and the domain's cognate ligands. The cognate ligands are defined as the ligands which are bound by the domain's protein chains. The Tanimoto similarity is calculated using the Morgan fingerprint of the ligand and the Morgan fingerprints of the cognate ligands. The Morgan fingerprints are calculated using the RDKit python package.

### ProCogGraph method of domain ownership

Since the original PROCOGNATE database, several tools have been developed, and are integrated into PDBe-graph, which allow a rich description of contacts - primarily, Arpeggio.

In ProCogGraph, domain contact counts are enumerated based on annotations from Arpeggio in PDBe-graph.

. The domain which has the highest number of contacts with the ligand is defined as the domain which is most likely to be responsible for the binding of the ligand.

## Accessing the Database

The ProCogGraph database is available to access at [ProCogGraph](procoggraph.com)

## Creating the Database from source

The ProCogGraph database is created using Neo4J. The database is created using the following steps:

1. Download the latest version of [Neo4J Desktop](https://neo4j.com/download/)
2. Copy the Neo4j files from the `neo4j` folder into the `import` folder in the Neo4J Database.
3. Run the admin-import tool:
```
bin/neo4j-admin database import full --array-delimiter="|" --skip-bad-relationships --delimiter="\t" --nodes=boundMolecule=import/bound_molecules.csv --nodes=boundEntity=import/bound_entities.csv --relationships=CONTAINS=import/bm_be_rels.csv --nodes=biologicalLigand=import/biological_ligand_nodes.csv --relationships=HAS_SIMILARITY=import/bound_entity_parity_score_rels.csv.gz --nodes=ecID=import/ec_id_nodes.csv --relationships=IS_IN_EC=import/biological_ligands_ec.csv --nodes=ecClass=import/ec_nodes_class.csv --relationships=IS_IN_CLASS=import/ec_class_subclass_rel.csv --nodes=ecSubClass=import/ec_nodes_subclass.csv --relationships=IS_IN_SUBCLASS=import/ec_subclass_subsubclass_rel.csv --nodes=ecSubSubClass=import/ec_nodes_subsubclass.csv --relationships=IS_IN_SUBSUBCLASS=import/ec_subsubclass_id_rel.csv --nodes=scopDomain=import/scop_domains_nodes.csv --relationships=IS_IN_SCOP_FAMILY=import/scop_domain_family_rels.csv --nodes=scopFamily=import/scop_family_nodes.csv --relationships=IS_IN_SCOP_SUPERFAMILY=import/scop_family_superfam_rels.csv --nodes=scopSuperfamily=import/scop_superfamily_nodes.csv --relationships=IS_IN_SCOP_FOLD=import/scop_superfam_fold_rels.csv --nodes=scopFold=import/scop_fold_nodes.csv --relationships=IS_IN_SCOP_CLASS=import/scop_fold_class_rels.csv --nodes=IS_IN_SCOP_CLASS=import/scop_class_nodes.csv --nodes=cathClass=import/cath_class_nodes.csv --relationships=IS_IN_CATH_CLASS=import/cath_class_architecture_rels.csv --nodes=cathArchitecture=import/cath_architecture_nodes.csv --relationships=IS_IN_CATH_ARCHITECTURE=import/cath_architecture_topology_rels.csv --nodes=cathTopology=import/cath_topology_nodes.csv --relationships=IS_IN_CATH_TOPOLOGY=import/cath_topology_homology_rels.csv --nodes=cathHomology=import/cath_homology_nodes.csv --relationships=IS_IN_CATH_HOMOLOGY=import/cath_homology_domain_rels.csv --nodes=cathDomain=import/cath_domains_nodes.csv --nodes=interproDomain=import/interpro_domain_nodes.csv --relationships=INTERACTS_WITH_LIGAND=import/scop_domain_ligand_interactions.csv --relationships=INTERACTS_WITH_LIGAND=import/cath_domain_ligand_interactions.csv --relationships=INTERACTS_WITH_LIGAND=import/interpro_domain_ligand_interactions.csv --nodes=proteinChain=import/pdb_protein_chain_nodes.csv --nodes=entry=import/pdb_entry_nodes.csv --relationships=IS_IN_PDB=import/bound_molecules_pdb_rels.csv --relationships=IS_IN_PDB=import/pdb_protein_rels.csv --relationships=IS_IN_PROTEIN_CHAIN=import/cath_protein_rels.csv --relationships=IS_IN_PROTEIN_CHAIN=import/scop_protein_rels.csv --relationships=IS_IN_PROTEIN_CHAIN=import/interpro_protein_rels.csv --relationships=IS_IN_EC=import/protein_ec_rels.csv --overwrite-destination neo4j
```
This creates the database.

4. The database can then be started and accessed using the Neo4J Desktop application.

## Database Schema

## Citations
