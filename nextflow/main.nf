#! /usr/bin/env nextflow

//maybe we can split the bio-h and sifts download to post process mmcif so that we dont unnecessarily donwload them? might help with the out of sync problem too
process DOWNLOAD_MMCIF {
    errorStrategy 'ignore' //ignore errors - modelserver may return 404 if assembly not available.

    storeDir "${params.cache_dir}/mmcif"
    publishDir "${params.publish_dir}/mmcif"
    input:
        tuple ( val(pdb_id), val(assembly_id) )

    output:
    tuple val(pdb_id), val(assembly_id), path("${pdb_id}_updated.cif.gz"), path("${pdb_id}_bio-h.cif")
    script:
    """
    wget -O "${pdb_id}_updated.cif" "https://www.ebi.ac.uk/pdbe/entry-files/${pdb_id}_updated.cif"; gzip "${pdb_id}_updated.cif" ;
    wget -O "${pdb_id}_bio-h.cif" "https://www.ebi.ac.uk/pdbe/model-server/v1/${pdb_id}/full?encoding=cif&data_source=pdb-h"
    """
}

process PROCESS_MMCIF {
    publishDir "${params.publish_dir}/process_struct"
    errorStrategy 'ignore'
    cache 'lenient'
    input:
        tuple( val(pdb_id), val(assembly_id), path(updated_cif) )
        
    output:
        tuple val(pdb_id), val(assembly_id), path("${pdb_id}_bound_entity_info.pkl"), path("${pdb_id}_arpeggio.csv")
    script:
    """
    python3 ${workflow.projectDir}/bin/process_pdb_structure.py --cif ${updated_cif} --pdb_id ${pdb_id} --assembly_id ${assembly_id}
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
    errorStrategy 'ignore'
    input:
        tuple val(pdb_id), val(assembly_id), path(updated_cif), path(bound_entity_pickle), path(arpeggio_json), val(domain_contact_cutoff)
    output:
        path("${pdb_id}_bound_entity_contacts.tsv")
    script:
    """
    python3 ${workflow.projectDir}/bin/process_pdb_contacts.py --bound_entity_pickle ${bound_entity_pickle} --cif ${updated_cif} --contacts ${arpeggio_json} --pdb_id ${pdb_id} --assembly_id ${assembly_id} --domain_contact_cutoff ${domain_contact_cutoff}
    """
}

process PROCESS_ALL_CONTACTS {
    publishDir "${params.publish_dir}/contacts"
    cache 'lenient'
    input:
        path combined_contacts
        path ccd_cif
        path pfam_a
        path pfam_clan_rels
        path pfam_clans
        path scop_domains_info
        path scop_domains_description
        path interpro_xml
        path cath_domains_descriptions
        path glycoct_cache
        path smiles_cache
        path csdb_linear_cache
        path enzyme_dat_file
        path enzyme_class_file
        path sifts_ec_mapping


    output:
        path("bound_entities_to_score.pkl"), emit: bound_entities
        path("cath_pdb_residue_interactions.csv.gz"), emit: cath
        path("scop_pdb_residue_interactions.csv.gz"), emit: scop
        path ("pfam_pdb_residue_interactions.csv.gz"), emit: pfam
        path ("interpro_pdb_residue_interactions.csv.gz"), emit: interpro
        path ("scop2b_pdb_residue_interactions.csv.gz"), emit: scop2b
        
    script:
    """
    python3 ${workflow.projectDir}/bin/process_all_pdb_contacts.py --contacts ${combined_contacts} --ccd_cif ${ccd_cif} --pfam_a_file ${pfam_a} --pfam_clan_rels ${pfam_clan_rels} --pfam_clans ${pfam_clans} --scop_domains_info_file ${scop_domains_info} --scop_descriptions_file ${scop_domains_description} --interpro_xml ${interpro_xml} --cddf ${cath_domains_descriptions} --glycoct_cache ${glycoct_cache} --smiles_cache ${smiles_cache} --csdb_linear_cache ${csdb_linear_cache} --enzyme_dat_file ${enzyme_dat_file} --enzyme_class_file ${enzyme_class_file} --sifts_ec_mapping ${sifts_ec_mapping}
    """

}

process GET_COGNATE_LIGANDS {
   storeDir "${params.cache_dir}/cognate_ligands"
   publishDir "${params.publish_dir}/cognate"
   input:
       path enzyme_dat_file
       path enzyme_class_file
       path pubchem_mapping
       path chebi_kegg_mapping
       path chebi_relations
       path rhea_reactions
       path kegg_enzyme_string
       path kegg_reaction_string
       path smiles_cache
       path csdb_linear_cache
       path compound_cache_dir
       path gtc_cache
       path brenda_mol_dir
       path brenda_json
       path brenda_ligands

   output:
   path "biological_ligands_df.pkl"
   script:
   """
   python3 ${workflow.projectDir}/bin/get_ec_information.py --ec_dat ${enzyme_dat_file} --enzyme_class_file ${enzyme_class_file} --pubchem ${pubchem_mapping} --chebi ${chebi_kegg_mapping} --chebi_relations ${chebi_relations} --rhea_reactions ${rhea_reactions} --kegg_enzyme_string ${kegg_enzyme_string} --kegg_reaction_string ${kegg_reaction_string} --smiles_cache ${smiles_cache} --csdb_cache ${csdb_linear_cache} --compound_cache_dir ${compound_cache_dir} --gtc_cache ${gtc_cache} --brenda_mol_dir ${brenda_mol_dir} --brenda_json ${brenda_json} --brenda_ligands ${brenda_ligands}
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
        path cognate_ligands
        path parity_cache
        val parity_threshold
    
    output:
        path "all_parity_calcs.pkl", emit: all_parity_calcs
        path "cache_parity_calcs.pkl", emit: parity_cache
    
    script:
    """
    python3 ${workflow.projectDir}/bin/get_pdb_parity.py --cache ${parity_cache} --processed_ligands_file ${bound_entities_to_score} --cognate_ligands_file ${cognate_ligands} --threshold ${parity_threshold} --threads ${task.cpus}
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
        path enzyme_dat
        path enzyme_class
        val parity_threshold
        path rhea2ec
        path rhea_directions
        path rhea_reaction_smiles

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
    python3 ${workflow.projectDir}/bin/produce_neo4j_files.py --enzyme_dat_file ${enzyme_dat} --enzyme_class_file ${enzyme_class} --cognate_ligands ${cognate_ligands} --cath_domain_ownership ${cath_domain_ownership} --scop_domain_ownership ${scop_domain_ownership} --interpro_domain_ownership ${interpro_domain_ownership} --pfam_domain_ownership ${pfam_domain_ownership} --bound_entities ${bound_entities} --parity_calcs ${all_parity_calcs} --parity_threshold ${parity_threshold} --rhea2ec ${rhea2ec} --rhea_dir ${rhea_directions} --rhea_reaction_smiles ${rhea_reaction_smiles}
    """
}

process PROCESS_SIFTS {
    input:
        path sifts_file
        path assemblies_file
    output:
        path "sifts_assemblies.csv"
    script:
    """
    python3 ${workflow.projectDir}/bin/process_sifts.py --sifts_file ${sifts_file} --assemblies_file ${assemblies_file}
    """
}

workflow {
    cognate_ligands = GET_COGNATE_LIGANDS(Channel.fromPath("${params.enzyme_dat_file}"), Channel.fromPath("${params.enzyme_class_file}"), Channel.fromPath("${params.pubchem_mapping}"), Channel.fromPath("${params.chebi_kegg}"), Channel.fromPath("${params.chebi_relations}"), Channel.fromPath("${params.rhea_reactions}"), Channel.fromPath("${params.kegg_enzyme_cache}"), Channel.fromPath("${params.kegg_reaction_cache}"), Channel.fromPath("${params.smiles_cache}"), Channel.fromPath("${params.csdb_linear_cache}"), Channel.fromPath("${params.kegg_compound_cache_dir}", type: "dir"), Channel.fromPath("${params.glytoucan_cache}"), Channel.fromPath("${params.brenda_mol_dir}", type: "dir"), Channel.fromPath("${params.brenda_json}"), Channel.fromPath("${params.brenda_ligands}"))
    enzyme_structures = PROCESS_SIFTS(Channel.fromPath("${params.sifts_file}"), Channel.fromPath("${params.assemblies_file}"))
    pdb_ids = enzyme_structures | splitCsv(header:true) | map { row -> [row.PDB, row.ASSEMBLY_ID] } | randomSample( 10, 42 )
    download_mmcif = DOWNLOAD_MMCIF( pdb_ids
            .map { all_out -> [all_out[0], all_out[1]] })
    process_mmcif = PROCESS_MMCIF( download_mmcif
        .map { all_out -> [all_out[0], all_out[1], all_out[2]] } )
    arpeggio = RUN_ARPEGGIO(
        process_mmcif
            .map { all_out -> [all_out[0], all_out[3]] }
        .join( 
        download_mmcif
            .map { all_out -> [all_out[0], all_out[3]] } ))
    contacts = PROCESS_CONTACTS(
        download_mmcif.map { all_out -> [all_out[0], all_out[1], all_out[2]] } 
        .join(
        process_mmcif.map { all_out -> [all_out[0], all_out[2]] })
        .join(
        arpeggio.map { all_out -> [all_out[0], all_out[1]] } )
        .combine( Channel.from(params.domain_contact_cutoff) ) )
    collected = contacts.collectFile(name: 'combined_files.tsv', keepHeader: true, skip: 1)
    all_contacts = PROCESS_ALL_CONTACTS( collected, Channel.fromPath("${params.ccd_cif}"), Channel.fromPath("${params.pfam_a_file}"), Channel.fromPath("${params.pfam_clan_rels}"), Channel.fromPath("${params.pfam_clans}"), Channel.fromPath("${params.scop_domains_info_file}"), Channel.fromPath("${params.scop_descriptions_file}"), Channel.fromPath("${params.interpro_xml}"), Channel.fromPath("${params.cddf}"), Channel.fromPath("${params.glycoct_cache}"), Channel.fromPath("${params.smiles_cache}"), Channel.fromPath("${params.csdb_linear_cache}"), Channel.fromPath("${params.enzyme_dat_file}"), Channel.fromPath("${params.enzyme_class_file}"), Channel.fromPath("${params.sifts_file}") )
    score_ligands = SCORE_LIGANDS( all_contacts.bound_entities, cognate_ligands, Channel.fromPath(params.parity_cache), Channel.from(params.parity_threshold) )
    produce_neo4j_files = PRODUCE_NEO4J_FILES( score_ligands.all_parity_calcs, cognate_ligands , all_contacts.bound_entities, all_contacts.cath, all_contacts.scop, all_contacts.interpro, all_contacts.pfam, Channel.fromPath("${params.enzyme_dat_file}"), Channel.fromPath("${params.enzyme_class_file}"), Channel.from(params.parity_threshold), Channel.fromPath("${params.rhea2ec}"), Channel.fromPath("${params.rhea_directions}"), Channel.fromPath("${params.rhea_reactions_smiles}") )
}

