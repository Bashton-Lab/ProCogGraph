#! /usr/bin/env nextflow


// be sure to set conda.enabled = true in nextflow.config

params.sifts_file = '/raid/MattC/repos/sifts.csv'
params.data_dir = '/raid/MattC/repos/ProCogGraphData/data_files'
params.cache_dir = '/raid/MattC/repos/ProCogGraphData/cache_files'
params.publish_dir = "test_out"
params.rhea_reactions = "/raid/MattC/repos/ProCogGraphData/procoggraph_20240418/cognate_ligands/rhea_reactions.pkl"
params.neo4j_user = 'neo4j'
params.neo4j_password = 'password'
params.neo4j_bolt_uri = 'bolt://localhost:7687'

params.publish_dir = '/raid/MattC/repos/ProCogGraphData/procoggraph'
params.data_files_dir = '/raid/MattC/repos/ProCogGraphData/data_files/'
params.cache_files_dir = '/raid/MattC/repos/ProCogGraphData/cache_files/'
params.scripts_dir = '/raid/MattC/repos/ProCogGraph/snakemake/scripts/'
params.parity_threshold = 0.4
params.domain_contact_cutoff = 3
params.threads = 75

params.enzyme_dat_file = 'enzyme.dat'
params.enzyme_class_file = 'enzclass.txt'
params.pdbe_graph_yaml = 'pdbe_graph_queries.yaml'
params.glycoct_cache = 'glycoct_cache.pkl'
params.smiles_cache = 'smiles_cache.pkl'
params.csdb_linear_cache = 'csdb_linear_cache.pkl'
params.sifts_ec_mapping = 'pdb_chain_enzyme.tsv.gz'
params.scop_domains_info_file = 'dir.cla.scop.1_75.txt'
params.scop_descriptions_file = 'dir.des.scop.1_75.txt'
params.interpro_xml = 'interpro.xml.gz'
params.pfam_clans = 'clan.txt.gz'
params.pfam_clan_rels = 'clan_membership.txt.gz'
params.pfam_a_file = 'pfamA.txt.gz'
params.cddf = 'cath-domain-description-file.txt'
params.ccd_cif = 'ccd.cif'

params.pubchem_mapping = 'pubchem_substance_id_mapping.txt'
params.chebi_kegg = 'ChEBI_Results.tsv'
params.chebi_relations = 'relation.tsv'
params.chebi_names = 'chebi_names.tsv.gz'
params.brenda_ligands = 'brenda_ligands_2023_1.csv'
params.brenda_mol_dir = 'molfiles_in_2023_1'
params.brenda_json = 'brenda_2023_1.json'
params.rhea2ec = 'rhea2ec.tsv'
params.rhea_directions = 'rhea-directions.tsv'
params.rhea_reactions_smiles = 'rhea-reaction-smiles.tsv'
params.rhea_rd_dir = 'rd'
params.rhea_reactions = 'rhea_reactions.pkl'
params.parity_cache = 'cache_parity_calcs.pkl'
params.sugar_cifs_cache_dir = 'sugar_cifs'
params.cognate_ligands = '/raid/MattC/repos/ProCogGraphData/procoggraph_20240418/cognate_ligands/biological_ligands_df.pkl'

params.kegg_enzyme_cache = 'kegg_enzyme_strings.txt'
params.kegg_reaction_cache = 'kegg_reaction_strings.txt'
params.kegg_compound_cache_dir = 'kegg_compound_cache_dir'
params.glytoucan_cache = 'gtc_cache.pkl'

//maybe we can split the bio-h and sifts download to post process mmcif so that we dont unnecessarily donwload them? might help with the out of sync problem too
process DOWNLOAD_MMCIF {
    storeDir "${params.cache_dir}/mmcif"
    publishDir "${params.publish_dir}/mmcif"
    input:
    val pdb_id

    output:
    tuple val(pdb_id), path("${pdb_id}_updated.cif"), path("${pdb_id}_bio-h.cif")
    script:
    """
    wget -O "${pdb_id}_updated.cif" "https://www.ebi.ac.uk/pdbe/entry-files/${pdb_id}_updated.cif"; 
    wget -O "${pdb_id}_bio-h.cif" "https://www.ebi.ac.uk/pdbe/model-server/v1/${pdb_id}/full?encoding=cif&data_source=pdb-h"
    """
}

process PROCESS_MMCIF {
    publishDir "${params.publish_dir}/process_struct"
    errorStrategy 'ignore'
    cache 'lenient'
    input:
        tuple( val(pdb_id), path(updated_cif) )
        
    output:
        tuple val(pdb_id), path("${pdb_id}_bound_entity_info.pkl"), path("${pdb_id}_arpeggio.csv")
    script:
    """
    python3 ${workflow.projectDir}/bin/process_pdb_structure.py --cif ${updated_cif} --pdb_id ${pdb_id}
    """

}

process RUN_ARPEGGIO {
    cache 'lenient'
    conda '/raid/MattC/repos/envs/arpeggio-env.yaml'
    publishDir "${params.publish_dir}/arpeggio"
    input:
        tuple val(pdb_id), path(arpeggio_selections), path(bio_h_cif)

    output:
        tuple val(pdb_id), path("${pdb_id}_bio-h.json")

    script:
    """
    pdbe-arpeggio -sf ${arpeggio_selections} ${bio_h_cif}
    """

}

process PROCESS_CONTACTS {
    cache 'lenient'
    publishDir "${params.publish_dir}/process_contacts"
    input:
        tuple val(pdb_id), path(updated_cif), path(bound_entity_pickle), path(arpeggio_json)
    output:
        path("${pdb_id}_bound_entity_contacts.tsv")
    script:
    """
    python3 ${workflow.projectDir}/bin/process_pdb_contacts.py --bound_entity_pickle ${bound_entity_pickle} --cif ${updated_cif} --contacts ${arpeggio_json} --pdb_id ${pdb_id}
    """
}

process PROCESS_ALL_CONTACTS {
    publishDir "${params.publish_dir}/contacts"
    cache 'lenient'
    input:
        path combined_contacts
    output:
        path("bound_entities_to_score.pkl"), emit: bound_entities
        path("cath_pdb_residue_interactions.csv.gz"), emit: cath
        path("scop_pdb_residue_interactions.csv.gz"), emit: scop
        path ("pfam_pdb_residue_interactions.csv.gz"), emit: pfam
        path ("interpro_pdb_residue_interactions.csv.gz"), emit: interpro
        path ("scop2b_pdb_residue_interactions.csv.gz"), emit: scop2b
        
    script:
    """
    python3 ${workflow.projectDir}/bin/process_all_pdb_contacts.py --contacts ${combined_contacts} --ccd_cif "${params.data_dir}/${params.ccd_cif}" --pfam_a_file "${params.data_dir}/${params.pfam_a_file}" --pfam_clan_rels "${params.data_dir}/${params.pfam_clan_rels}" --pfam_clans "${params.data_dir}/${params.pfam_clans}" --scop_domains_info_file "${params.data_dir}/${params.scop_domains_info_file}" --scop_descriptions_file "${params.data_dir}/${params.scop_descriptions_file}" --interpro_xml "${params.data_dir}/${params.interpro_xml}" --cddf "${params.data_dir}/${params.cddf}" --glycoct_cache "${params.cache_dir}/${params.glycoct_cache}" --smiles_cache "${params.cache_dir}/${params.smiles_cache}" --csdb_linear_cache "${params.cache_dir}/${params.csdb_linear_cache}" --enzyme_dat_file "${params.data_dir}/${params.enzyme_dat_file}" --enzyme_class_file "${params.data_dir}/${params.enzyme_class_file}" --sifts_ec_mapping "${params.data_dir}/${params.sifts_ec_mapping}"
    """

}

process GET_COGNATE_LIGANDS {
   storeDir "${params.cache_dir}/cognate_ligands"
   publishDir "${params.publish_dir}/cognate"
   input:
       path enzyme_dat_file
       path enzyme_class_file
   output:
   path "biological_ligands_df.pkl"
   script:
   """
   python3 ${workflow.projectDir}/bin/get_ec_information.py --ec_dat ${enzyme_dat_file} --enzyme_class_file ${enzyme_class_file} --pubchem ${params.data_dir}/${params.pubchem_mapping} --chebi ${params.data_dir}/${params.chebi_kegg} --chebi_relations ${params.data_dir}/${params.chebi_relations} --rhea_reactions ${params.rhea_reactions} --kegg_enzyme_string ${params.cache_dir}/${params.kegg_enzyme_cache} --kegg_reaction_string ${params.cache_dir}/${params.kegg_reaction_cache} --smiles_cache ${params.cache_dir}/${params.smiles_cache} --csdb_cache ${params.cache_dir}/${params.csdb_linear_cache} --compound_cache_dir ${params.cache_dir}/${params.kegg_compound_cache_dir} --gtc_cache ${params.cache_dir}/${params.glytoucan_cache} --brenda_mol_dir ${params.data_dir}/${params.brenda_mol_dir} --brenda_json ${params.data_dir}/${params.brenda_json} --brenda_ligands ${params.data_dir}/${params.brenda_ligands}
   """
}

// rule preprocess_rhea:
//     input:
//         rhea2ec = config["data_files_dir"] + config["rhea2ec"],
//         rhea_directions = config["data_files_dir"] + config["rhea_directions"],
//         rd_dir = config["data_files_dir"] + config["rhea_rd_dir"],
//         chebi_names = config["data_files_dir"] + config["chebi_names"],
//         script = config["scripts_dir"] + "preprocess_rhea.py"
//     output:
//         rhea_reactions = output_dir + "/cognate_ligands/rhea_reactions.pkl"
//     params:
//         outdir = output_dir + "/cognate_ligands"
//     threads: 1
//     shell:
//         "python3 {input.script} --rhea_ec_mapping {input.rhea2ec} --rhea_reaction_directions {input.rhea_directions} --rd_dir {input.rd_dir} --chebi_names {input.chebi_names} --outdir {params.outdir}"



process SCORE_LIGANDS {
    cache 'lenient'
    cpus "${params.threads}"
    publishDir "${params.publish_dir}/scores"
    input:
        path bound_entities_to_score
        //path cognate_ligands
    
    output:
        path "all_parity_calcs.pkl", emit: all_parity_calcs
        path "cache_parity_calcs.pkl", emit: parity_cache
    
    script:
    """
    python3 ${workflow.projectDir}/bin/get_pdb_parity.py --cache ${params.cache_dir}/${params.parity_cache}  --processed_ligands_file ${bound_entities_to_score} --cognate_ligands_file ${params.cognate_ligands} --threshold ${params.parity_threshold} --threads ${task.cpus}
    """
}

process PRODUCE_NEO4J_FILES {
    cache 'lenient'
    publishDir "${params.publish_dir}/neo4j"
    input:
        path all_parity_calcs
        path cognate_ligands
        path bound_entities
        path cath_domain_ownership
        path scop_domain_ownership
        path interpro_domain_ownership
        path pfam_domain_ownership
    output:
        path "ec_id_nodes.csv.gz"
        path "ec_nodes_class.csv.gz"
        path "ec_nodes_subclass.csv.gz"
        path "ec_nodes_subsubclass.csv.gz"
        path "ec_class_subclass_rel.csv.gz"
        path "ec_subclass_subsubclass_rel.csv.gz"
        path "ec_subsubclass_id_rel.csv.gz"
        path "cognate_ligand_nodes.csv.gz"
        path "cognate_ligands_ec.csv.gz"
        path "pdb_entry_nodes.csv.gz"
        path "pdb_ec_rels.csv.gz"
        path "scop_domains_nodes.csv.gz"
        path "cath_domains_nodes.csv.gz"
        path "interpro_domain_nodes.csv.gz"
        path "pfam_domains_nodes.csv.gz"
        path "scop_family_nodes.csv.gz"
        path "scop_superfamily_nodes.csv.gz"
        path "scop_fold_nodes.csv.gz"
        path "scop_class_nodes.csv.gz"
        path "scop_domain_family_rels.csv.gz"
        path "scop_family_superfam_rels.csv.gz"
        path "scop_superfam_fold_rels.csv.gz"
        path "scop_fold_class_rels.csv.gz"
        path "cath_class_nodes.csv.gz"
        path "cath_architecture_nodes.csv.gz"
        path "cath_topology_nodes.csv.gz"
        path "cath_homologous_superfamily_nodes.csv.gz"
        path "cath_class_architecture_rels.csv.gz"
        path "cath_architecture_topology_rels.csv.gz"
        path "cath_topology_homology_rels.csv.gz"
        path "cath_homologous_superfamily_domain_rels.csv.gz"
        path "pfam_clans.csv.gz"
        path "pfam_clan_rels.csv.gz"
        path "bound_descriptors.csv.gz"
        path "bound_entities.csv.gz"
        path "be_bd_rels.csv.gz"
        path "be_pdb_rels.csv.gz"
        path "bound_entity_parity_score_rels.csv.gz"
        path "cath_domain_ligand_interactions.csv.gz"
        path "scop_domain_ligand_interactions.csv.gz"
        path "pfam_domain_ligand_interactions.csv.gz"
        path "interpro_domain_ligand_interactions.csv.gz"
        path "cath_pdb_rels.csv.gz"
        path "scop_pdb_rels.csv.gz"
        path "pfam_pdb_rels.csv.gz"
        path "interpro_pdb_rels.csv.gz"
        path "procoggraph_node.csv.gz"
    script:
    """
    python3 ${workflow.projectDir}/bin/produce_neo4j_files.py --enzyme_dat_file ${params.data_dir}/${params.enzyme_dat_file}  --enzyme_class_file ${params.data_dir}/${params.enzyme_class_file} --cognate_ligands ${cognate_ligands} --cath_domain_ownership ${cath_domain_ownership} --scop_domain_ownership ${scop_domain_ownership} --interpro_domain_ownership ${interpro_domain_ownership} --pfam_domain_ownership ${pfam_domain_ownership} --bound_entities ${bound_entities} --parity_calcs ${all_parity_calcs} --parity_threshold ${params.parity_threshold} --rhea2ec ${params.data_dir}/${params.rhea2ec} --rhea_dir ${params.data_dir}/${params.rhea_directions} --rhea_reaction_smiles ${params.data_dir}/${params.rhea_reactions_smiles}
    """
}

workflow {
    cognate_ligands = GET_COGNATE_LIGANDS(Channel.fromPath("${params.data_dir}/${params.enzyme_dat_file}"), Channel.fromPath("${params.data_dir}/${params.enzyme_class_file}"))
    pdb_ids = Channel.fromPath(params.sifts_file) | splitCsv(header:true) | map { row -> row.pdb_id }
    download_mmcif = DOWNLOAD_MMCIF(pdb_ids)
    process_mmcif = PROCESS_MMCIF( download_mmcif
        .map { all_out -> [all_out[0], all_out[1]] } )
    arpeggio = RUN_ARPEGGIO(
        process_mmcif
            .map { all_out -> [all_out[0], all_out[2]] }
        .join( 
        download_mmcif
            .map { all_out -> [all_out[0], all_out[2]] } ))
    contacts = PROCESS_CONTACTS(
        download_mmcif.map { all_out -> [all_out[0], all_out[1]] } 
        .join(
        process_mmcif.map { all_out -> [all_out[0], all_out[1]] })
        .join(
        arpeggio.map { all_out -> [all_out[0], all_out[1]] } ) )
    collected = contacts.collectFile(name: 'combined_files.tsv', keepHeader: true, skip: 1)
    all_contacts = PROCESS_ALL_CONTACTS( collected )
    score_ligands = SCORE_LIGANDS( all_contacts.bound_entities )
    produce_neo4j_files = PRODUCE_NEO4J_FILES( score_ligands.all_parity_calcs, cognate_ligands , all_contacts.bound_entities, all_contacts.cath, all_contacts.scop, all_contacts.interpro, all_contacts.pfam)
}
