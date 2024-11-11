/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { MULTIQC                } from '../modules/local/multiqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_multiplesequencealign_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Local subworkflows
//
include { STATS                  } from '../subworkflows/local/stats'
include { CREATE_TCOFFEETEMPLATE } from '../modules/local/create_tcoffee_template'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR                          } from '../modules/nf-core/untar/main'
include { CSVTK_JOIN as MERGE_STATS_EVAL } from '../modules/nf-core/csvtk/join/main.nf'
include { PIGZ_COMPRESS                  } from '../modules/nf-core/pigz/compress/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT CLASS-MODULES MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MSA_ALIGNMENT } from '../subworkflows/mirpedrol/msa_alignment/main'
include { MSA_GUIDETREE } from '../subworkflows/mirpedrol/msa_guidetree/main'
include { MSA_TREEALIGN } from '../subworkflows/mirpedrol/msa_treealign/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MULTIPLESEQUENCEALIGN {

    take:
    ch_input    // channel: [ meta, path(sequence.fasta), path(reference.fasta), path(pdb_structures.tar.gz), path(templates.txt) ]
    ch_versions // channel: [ path(versions.yml) ]

    main:
    def ch_multiqc_files             = Channel.empty()
    def stats_summary                = Channel.empty()

    ch_input
        .map {
            meta, fasta, ref, str, template ->
                [ meta, file(fasta) ]
        }
        .set { ch_seqs }

    ch_input
        .filter { input -> input[2].size() > 0}
        .map {
            meta, fasta, ref, str, template ->
                [ meta, file(ref) ]
        }
        .set { ch_refs }

    ch_input
        .filter { input -> input[4].size() > 0}
        .map {
            meta, fasta, ref, str, template ->
                [ meta, file(template) ]
        }
        .set { ch_templates }

    ch_input
        .map {
            meta, fasta, ref, str, template ->
                [ meta, str ]
        }
        .filter { input -> input[1].size() > 0 }
        .set { ch_structures }

    // ----------------
    // STRUCTURES
    // ----------------
    // Structures are taken from a directory of PDB files.
    // If the directory is compressed, it is uncompressed first.
    ch_structures
        .branch { structures ->
            compressed:   structures[1].endsWith('.tar.gz')
            uncompressed: true
        }
        .set { ch_structures }

    UNTAR (ch_structures.compressed)
        .untar
        .mix(ch_structures.uncompressed)
        .map {
            meta,dir ->
                [ meta,file(dir).listFiles().collect() ]
        }
        .set { ch_structures }

    // ----------------
    // TEMPLATES
    // ----------------
    // If a family does not present a template but structures are provided, create one.
    ch_structures
        .join(ch_templates, by:0, remainder:true)
        .set { ch_structures_template }
    ch_structures_template
        .branch { it ->
            template: it[2] != null
            no_template: true
        }
        .set { ch_structures_branched }

    // Create the new templates and merge them with the existing templates
    CREATE_TCOFFEETEMPLATE (
        ch_structures_branched.no_template
            .map {
                meta,structures,template ->
                    [ meta, structures ]
            }
    )
    def new_templates = CREATE_TCOFFEETEMPLATE.out.template
    ch_structures_branched.template
        .map {
            meta,structures,template ->
                [ meta, template ]
        }
        .set { forced_templates }

    def ch_templates_merged = forced_templates.mix(new_templates)

    // Merge the structures and templates channels, ready for the alignment
    ch_structures_template = ch_templates_merged.combine(ch_structures, by:0)

    //
    // Compute summary statistics about the input sequences
    //
    if (!params.skip_stats) {
        STATS (
            ch_seqs,
            ch_structures
        )
        ch_versions   = ch_versions.mix(STATS.out.versions)
        stats_summary = stats_summary.mix(STATS.out.stats_summary)
    }

    def msa_alignment = Channel.empty()

    if (params.guidetree && params.treealign) {
        //
        // Compute tree
        //
        MSA_GUIDETREE (
            ch_seqs
        )
        ch_versions = ch_versions.mix(MSA_GUIDETREE.out.versions)

        // Prepare channels for treealign to make sure the correct tree is used for the respective alignment
        ch_seqs
            .combine(MSA_GUIDETREE.out.tree, by:0)
            .set { ch_seqs_trees }
        ch_seqs_trees
            .multiMap { meta, seq, tree ->
                sequences: [meta, seq]
                trees: [meta, tree]
            }
            .set { ch_seqs_trees_multi }

        //
        // Align with a given tree
        //
        MSA_TREEALIGN (
            ch_seqs_trees_multi.sequences,
            ch_seqs_trees_multi.trees
        )
        ch_versions = ch_versions.mix(MSA_TREEALIGN.out.versions)
        msa_alignment.mix(MSA_TREEALIGN.out.alignment)
    }

    if (params.alignment) {
        //
        // Align
        //
        MSA_ALIGNMENT (
            ch_seqs
        )
        ch_versions = ch_versions.mix(MSA_ALIGNMENT.out.versions)
        msa_alignment.mix(MSA_ALIGNMENT.out.alignment)
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    def multiqc_out = Channel.empty()
    if (!params.skip_multiqc && (!params.skip_stats || !params.skip_eval)) {

        def ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        def ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        def ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            []
        )
        multiqc_out = MULTIQC.out.report.toList()
    }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
    multiqc  = multiqc_out // channel: [ path(multiqc_report.html) ]
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


