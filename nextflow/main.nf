#! /usr/bin/env nextflow

process PROCESS_MMCIF {
    label 'medmem' //gives 8gb and single core (max observed usage is 8gb)
    publishDir "${params.publish_dir}/process_struct", mode: 'copy'
    errorStrategy { task.exitStatus in 120..122 ? 'ignore' : 'terminate' }
    cache 'lenient'
    input:
        path manifest
        
    output:
        path("process_mmcif_log.txt"), emit: log
        path("bio_h_cif_*.csv"), emit: arpeggio_batch
    script:
    """
    python3 ${workflow.projectDir}/bin/process_pdb_structure.py --manifest ${manifest} --threads ${task.cpus}
    """

}

process RUN_ARPEGGIO {
    label 'arpeggio' //gives 18gb mem for arpeggio in slurm submission (offering overhead for complex structures)
    cache 'lenient'
    conda '/raid/MattC/repos/envs/arpeggio-env.yaml'
    publishDir "${params.publish_dir}/arpeggio", mode: 'copy'
    input:
        path arpeggio_batch

    output:
        path("*_bio-h.json")

    script:
    """
    ${workflow.projectDir}/bin/run_arpeggio.sh ${arpeggio_batch}
    """

}

process PROCESS_CONTACTS {
    label 'medmem' //gives 12gb and single core (max observed usage is below 12gb)
    cache 'lenient'
    publishDir "${params.publish_dir}/process_contacts", mode: 'copy'
    errorStrategy { task.exitStatus in 124..127 ? 'ignore' : 'terminate' }
    input:
        path combined_json
        path manifest
        val domain_contact_cutoff
    output:
        path "*_bound_entity_contacts.tsv", emit: contacts
        path "process_contacts_log.txt", emit: log
    script:
    """
    sed '\$ s/,\$//' ${combined_json} > combined_json_formatted.json
    echo "{ \$(cat combined_json_formatted.json) }" > combined_json_formatted.json
    python3 ${workflow.projectDir}/bin/process_pdb_contacts.py --contacts_json combined_json_formatted.json --manifest ${manifest} --domain_contact_cutoff ${domain_contact_cutoff}
    """
}

process PROCESS_ALL_CONTACTS {
    label 'largecpu_largemem'
    publishDir "${params.publish_dir}/contacts", mode: 'copy'
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
        path cath_names
        path cddf
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
    python3 ${workflow.projectDir}/bin/process_all_pdb_contacts.py --contacts ${combined_contacts} --ccd_cif ${ccd_cif} --pfam_a_file ${pfam_a} --pfam_clan_rels ${pfam_clan_rels} --pfam_clans ${pfam_clans} --scop_domains_info_file ${scop_domains_info} --scop_descriptions_file ${scop_domains_description} --interpro_xml ${interpro_xml} --cath_names ${cath_names} --cddf ${cddf} --glycoct_cache ${glycoct_cache} --smiles_cache ${smiles_cache} --csdb_linear_cache ${csdb_linear_cache} --enzyme_dat_file ${enzyme_dat_file} --enzyme_class_file ${enzyme_class_file} --sifts_ec_mapping ${sifts_ec_mapping}
    """

}

process SCORE_LIGANDS {
    label 'largecpu_largemem'
    cache 'lenient'
    publishDir "${params.publish_dir}/scores", mode: 'copy'
    input:
        path bound_entities_to_score
        path cognate_ligands
        path parity_cache
        val parity_threshold
    
    output:
        path "all_parity_calcs.pkl", emit: all_parity_calcs
        path "cache_parity_calcs_new.pkl", emit: parity_cache
    
    script:
    """
    python3 ${workflow.projectDir}/bin/get_pdb_parity.py --cache ${parity_cache} --processed_ligands_file ${bound_entities_to_score} --cognate_ligands_file ${cognate_ligands} --threshold ${parity_threshold} --threads ${task.cpus}
    """
}

process PRODUCE_NEO4J_FILES {
    label 'largecpu_largemem'
    cache 'lenient'
    publishDir "${params.publish_dir}/neo4j", mode: 'copy'
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

workflow {
    processed_struct_manifest = PROCESS_MMCIF( Channel.fromPath(params.manifest) )
    arpeggio = RUN_ARPEGGIO( processed_struct_manifest.arpeggio_manifest ).collectFile(name: 'combined_contacts.json', storeDir: "${params.publish_dir}/arpeggio", cache: true )
    contacts = PROCESS_CONTACTS( arpeggio, processed_struct_manifest.updated_manifest, params.domain_contact_cutoff )
    collected_contacts = contacts.contacts.collectFile(name: 'combined_contacts.tsv', storeDir: "${params.publish_dir}/contacts", cache: true)
    all_contacts = PROCESS_ALL_CONTACTS( collected_contacts, Channel.fromPath("${params.ccd_cif}"), Channel.fromPath("${params.pfam_a_file}"), Channel.fromPath("${params.pfam_clan_rels}"), Channel.fromPath("${params.pfam_clans}"), Channel.fromPath("${params.scop_domains_info_file}"), Channel.fromPath("${params.scop_descriptions_file}"), Channel.fromPath("${params.interpro_xml}"), Channel.fromPath("${params.cath_names}"), Channel.fromPath("${params.cddf}"), Channel.fromPath("${params.glycoct_cache}"), Channel.fromPath("${params.smiles_cache}"), Channel.fromPath("${params.csdb_linear_cache}"), Channel.fromPath("${params.enzyme_dat_file}"), Channel.fromPath("${params.enzyme_class_file}"), Channel.fromPath("${params.sifts_file}") )
    score_ligands = SCORE_LIGANDS( all_contacts.bound_entities, Channel.fromPath(params.cognate_ligands), Channel.fromPath(params.parity_cache), Channel.from(params.parity_threshold) )
    produce_neo4j_files = PRODUCE_NEO4J_FILES( score_ligands.all_parity_calcs, Channel.fromPath(params.cognate_ligands) , all_contacts.bound_entities, all_contacts.cath, all_contacts.scop, all_contacts.interpro, all_contacts.pfam, Channel.fromPath("${params.enzyme_dat_file}"), Channel.fromPath("${params.enzyme_class_file}"), Channel.from(params.parity_threshold), Channel.fromPath("${params.rhea2ec}"), Channel.fromPath("${params.rhea_directions}"), Channel.fromPath("${params.rhea_reactions_smiles}") )
}