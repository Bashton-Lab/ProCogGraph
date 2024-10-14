# ProCogGraph Tutorial

By following this tutorial, you should be able to familiarise yourself with the way of navigating the ProCogGraph NeoDash interface, and explore its main features. The tutorial makes use of the search tab, the second tab in the interface, here there are various boxes which can be used to search for different entities within the graph. In each case when a match is found matching hits are displayed to the right of the search box, pressing on these blue links will then take you to other tabs in the interface.

1. First, we need to log-on to ProCogGraph, point a web browser at: 
   [http://procoggraphprod.uksouth.cloudapp.azure.com:5005/](http://procoggraphprod.uksouth.cloudapp.azure.com:5005/)

2. Then use the password `procoggraph` to enter the NeoDash interface.

**Question 1**: Look up the PDB code `1LDM` *(use Search tab and entre the value into PDB Search box, press on blue 1ldm in PDB results box to the left to be taken to the second tab*).

1. What is this enzyme?

2. What is this enzymes EC number? (*See Summary box*). Look this up in the Enzyme

3. Query (in Search tab), and KEGG. 

4. What is bound to this enzyme in the PDB? (*See PDB Ligand Table*)

5. What cognate ligands are here, and which are missing? (*See cognate Ligand column in PDB Ligand Table*)

6. What is bound in their place? (*See the Het Code column, Het codes are the native PDB ligand codes*)

7. Press the Score button in the PDB Ligand table to visualise the atoms in the match, which have changed?

**Question 2**: Look up PDB `2SCU` (*use Search tab and entre the value into PDB Search box, press on blue 1ldm in PDB results box to the left to be taken to the second tab*).

- What is this enzyme?

- What is this enzymes EC number? (*See Summary box*). Look this up in the Enzyme

- Query (in Search tab), and KEGG.

- What is bound to this enzyme in the PDB? (*See PDB Ligand Table*)

- What cognate ligands are here, and which are missing? (*See cognate Ligand column in PDB Ligand Table*)

- What is bound in their place? (*See the Het Code column, Het codes are the native PDB ligand codes*)

- Press the Score button in the PDB Ligand table to visualise the atoms in the match, which have changed?

**Question 3**: Look up PDB `1DLJ`

- What is this enzyme?

- What is this enzymes EC number? (See Summary box). Look this up in the Enzyme Query (in Search tab), and KEGG. 

- What is bound to this enzyme in the PDB? (See PDB Ligand Table)

- What cognate ligands are here, and which are missing? (Lee cognate Ligand column in PDB Ligand Table)

- Which domains and residues in them share the binding of this ligand? 

HINT: *Here you want to press on the blue box, with a percentage, in the Contact % column in the CATH Domain Interactions table in the PDB view, this will display the interaction between this domain and the ligand in the Interface view to the right. Note that residues from the domain in question are blue, those from other domains are purple. And that interacting residues for this domain are listed*.

**Question 5**: Look up the 3 EC numbers from Questions 1-3 using the Search Tab and the Enzyme Search box. How many CATH Homologous Superfamilies and PDB structures are associated with each reaction?

**Question 6**: How does a domain duplication change the type of cognate ligands from PDB code `1KND` to PDB code `1BH5`?

- How are these ligands bound by domains and chains in each structure?

HINT: *You will want to view each ligand and its interactions again, press on the blue box under Contact % in the CATH Domain Interactions table in the PDB view, this will display the interaction between this domain and the ligand in the Interface view to the right. Note that residues from the domain in question are blue, those from other domains are purple. And that interacting residues for this domain are listed.  Also take a look at Figure 4 in the paper, can you view this yourself with the interface buttons. Remember to cross-ref the Het code of the ligand from the PDB Ligand Table with the PDB ligand column in the Domain Interactions table*.

**Question 7**: Using the domain search mode, compare one specialised and one generalised domain. SPECIALISED: `2.40.110.10` and GENERALISED: `3.20.20.100`

- How many cognate ligands are listed for each domain?

- Take a look at the Combinatorial Interactions, which other domains do we frequently see with the one looked up? 

HINT: *Use the search tab and enter the CATH codes given into the Domain Search box*.

**Question 8**: Find a structure that has a number of putative cognate ligands assigned. What structural features of the best scoring cognate ligand differentiate it from other potential cognates? Do all potential cognates have a similar structure? HINT: *You can work down the list in the PDB Ligand Table pressing on the Score buttons, note that if there is more than one cognate ligand then there will be separate sets of hits for each*.
