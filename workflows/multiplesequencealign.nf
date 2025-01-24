/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
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
include { EVALUATE               } from '../subworkflows/local/evaluate'
include { TEMPLATES              } from '../subworkflows/local/templates'
include { PREPROCESS             } from '../subworkflows/local/preprocess'
include { VISUALIZATION          } from '../subworkflows/local/visualization'


//
// MODULE: local modules
//
include { PREPARE_MULTIQC } from '../modules/local/prepare_multiqc'
include { PREPARE_SHINY   } from '../modules/local/prepare_shiny'

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
include { FASTAVALIDATOR                 } from '../modules/nf-core/fastavalidator/main'

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
    ch_input    // channel: [ meta, path(sequence.fasta), path(reference.fasta), path(dependency_files.tar.gz), path(templates.txt) ]
    ch_tools    // channel: [ meta_guidetree_treealign, meta_alignment ]

    main:
    ch_multiqc_files             = Channel.empty()
    ch_multiqc_table             = Channel.empty()
    evaluation_summary           = Channel.empty()
    stats_summary                = Channel.empty()
    stats_and_evaluation_summary = Channel.empty()
    ch_refs                      = Channel.empty()
    ch_templates                 = Channel.empty()
    ch_optional_data             = Channel.empty()

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_input
        .filter { it[1].size() > 0}
        .map {
            meta, fasta, ref, str, template ->
                [ meta, file(fasta) ]
        }
        .set { ch_seqs }

    ch_input
        .filter { it[2].size() > 0}
        .map {
            meta, fasta, ref, str, template ->
                [ meta, file(ref) ]
        }
        .set { ch_refs }

    ch_input
        .filter { it[4].size() > 0}
        .map {
            meta, fasta, ref, str, template ->
                [ meta, file(template) ]
        }
        .set { ch_templates }

    // ----------------
    // DEPENDENCY FILES
    // ----------------

    /*
    * We currently support 2 ways of reading in the optional_data:
    * 1. Provide a folder containing the optional_data via the `optional_data_dir` parameter
    * 2. Provide the dependency files directly in the input samplesheet
    */

    // If the optional_data folder is provided, use it to identify the optional_data based on sequence IDs
    if(params.optional_data_dir){

        // Identify the sequence IDs from the input fasta file(s)
        ch_seqs.splitFasta(record: [ id: true ] )
            .map { id, seq_id -> [ seq_id, id ] }
            .set { ch_seqs_split }

        // if compressed, uncompress the optional_data folder
        if(params.optional_data_dir.endsWith('.tar.gz')){

            optional_data_dir = Channel.fromPath(params.optional_data_dir)
                                        .map { it -> [[id: it.baseName],it] }

            UNTAR (optional_data_dir)
                .untar
                .map { meta, dir -> [ file(dir).listFiles() ] }
                .flatten()
                .set{ optional_data_to_be_mapped }
            ch_versions = ch_versions.mix(UNTAR.out.versions)

        }
        // otherwise, directly use the optional_data within the folder
        else {
            optional_data_to_be_mapped = Channel.fromPath(params.optional_data_dir+"/**")
        }

        // Map the optional_data to the sequence IDs
        optional_data_to_be_mapped
            .map { it -> [ [ id: it.baseName ], it ] }
            .combine(ch_seqs_split, by: 0)
            .map { dep_id, dep, fasta_id -> [ fasta_id, dep ] }
            .groupTuple(by: 0)
            .set { ch_optional_data }
    } else {

        // otherwise, use the dependency files provided in the input samplesheet
        ch_input
            .map {
                meta, fasta, ref, str, template ->
                    [ meta, str ]
            }
            .filter { it[1].size() > 0 }
            .set { ch_optional_data }

        // Dependency files are taken from a directory.
        // If the directory is compressed, it is uncompressed first.
        ch_optional_data
            .branch {
                compressed:   it[1].endsWith('.tar.gz')
                uncompressed: true
            }
            .set { ch_optional_data }

        UNTAR (ch_optional_data.compressed)
            .untar
            .mix(ch_optional_data.uncompressed)
            .map {
                meta,dir ->
                    [ meta,file(dir).listFiles().collect() ]
            }
            .set { ch_optional_data }
        ch_versions   = ch_versions.mix(UNTAR.out.versions)
    }

    //
    // VALIDATE AND PREPROCESS INPUT FILES
    //

    FASTAVALIDATOR(ch_seqs)
    ch_versions = ch_versions.mix(FASTAVALIDATOR.out.versions)

    if (!params.skip_preprocessing) {
        PREPROCESS(ch_optional_data)
        ch_optional_data = PREPROCESS.out.preprocessed_optionaldata
        ch_versions      = ch_versions.mix(PREPROCESS.out.versions)
    }


    //
    // TEMPLATES
    //
    TEMPLATES (
        ch_optional_data,
        ch_templates,
        "${params.templates_suffix}"
    )
    ch_optional_data_template = TEMPLATES.out.optional_data_template

    //
    // Compute summary statistics about the input sequences
    //
    if (!params.skip_stats) {
        STATS (
            ch_seqs,
            ch_optional_data
        )
        ch_versions   = ch_versions.mix(STATS.out.versions)
        stats_summary = stats_summary.mix(STATS.out.stats_summary)
    }

    ch_seqs
        .combine(ch_tools)
        // Add tools and arguments to the meta
        .multiMap {
            meta, fasta, meta_guidetree_treealign, meta_alignment ->
                guidetree: [ meta + ["guidetree":meta_guidetree_treealign.guidetree, "args_guidetree":meta_guidetree_treealign.args_guidetree, "args_guidetree_clean":meta_guidetree_treealign.args_guidetree_clean], fasta, meta_guidetree_treealign.guidetree]
                alignment: [ meta + ["alignment":meta_alignment.alignment, "args_alignment":meta_alignment.args_alignment, "args_alignment_clean":meta_alignment.args_alignment_clean], fasta, meta_alignment.alignment]
        }
        .set { ch_fasta_tools }

    ch_fasta_tools.guidetree
        .filter{ it -> it[0].guidetree }
        .unique()
        .set { ch_fasta_guidetree }
    ch_fasta_tools.alignment
        .filter{ it -> it[0].alignment }
        .unique()
        .set { ch_fasta_alignment }

    ch_fasta_guidetree.dump( tag: 'ch_fasta_guidetree' )
    ch_fasta_alignment.dump( tag: 'ch_fasta_alignment' )

    //
    // Compute tree
    //
    MSA_GUIDETREE (ch_fasta_guidetree)
    ch_versions = ch_versions.mix(MSA_GUIDETREE.out.versions)

    ch_seqs
        .map { meta, fasta ->
            [ meta.id, meta, fasta ]
        }
        .combine(
                MSA_GUIDETREE.out.tree
                    .map { meta, tree ->
                        [ meta.id, meta, tree ]
                    }
                , by:0
        ) // combine by meta ID
        .map { meta_id, meta_fasta, fasta, meta_tree, tree ->
            [ meta_tree.guidetree, meta_tree, fasta, tree ]
        }
        .combine(
            ch_tools
                .map { meta_guidetree_treealign, meta_alignment ->
                    [ meta_guidetree_treealign.guidetree, meta_guidetree_treealign ]
                }
            , by: 0
        ) // combine by guidetree
        .map {
            guidetree, meta, fasta, tree, meta_guidetree_treealign ->
                [meta + ["treealign":meta_guidetree_treealign.treealign, "args_treealign":meta_guidetree_treealign.args_treealign, "args_treealign_clean":meta_guidetree_treealign.args_treealign_lean], fasta, tree, meta_guidetree_treealign.treealign]
        }
        .multiMap {
            meta, fasta, tree, treealign ->
                fastas: [ meta, fasta, treealign ]
                trees: [ meta, tree, treealign ]
        }
        .set { ch_tree_treealign }

    ch_alignment_output = Channel.empty()

    //
    // Align with a given tree
    //
    MSA_TREEALIGN (
        ch_tree_treealign.fastas,
        ch_tree_treealign.trees
    )
    ch_versions = ch_versions.mix(MSA_TREEALIGN.out.versions)
    ch_alignment_output = ch_alignment_output.mix(MSA_TREEALIGN.out.alignment)


    //
    // Align
    //
    compress_during_align = !(params.skip_compression || (!params.skip_eval || params.build_consensus))

    MSA_ALIGNMENT ( ch_fasta_alignment )
    ch_versions = ch_versions.mix(MSA_ALIGNMENT.out.versions)
    ch_alignment_output = ch_alignment_output.mix(MSA_ALIGNMENT.out.alignment)

    if (!params.skip_compression && !compress_during_align) {
        PIGZ_COMPRESS (ch_alignment_output)
        ch_versions = ch_versions.mix(PIGZ_COMPRESS.out.versions)
    }

    //
    // Evaluate the quality of the alignment
    //
    if (!params.skip_eval) {
        EVALUATE (ch_alignment_output, ch_refs, ch_optional_data_template)
        ch_versions        = ch_versions.mix(EVALUATE.out.versions)
        evaluation_summary = evaluation_summary.mix(EVALUATE.out.eval_summary)
    }

    //
    // Combine stats and evaluation reports into a single CSV
    //
    if (!params.skip_stats || !params.skip_eval) {
        stats_summary_csv = stats_summary.map{ meta, csv -> csv }
        eval_summary_csv  = evaluation_summary.map{ meta, csv -> csv }
        stats_summary_csv.mix(eval_summary_csv)
                        .collect()
                        .map {
                            csvs ->
                                [ [ id:"summary_stats_eval" ], csvs ]
                        }
                        .set { stats_and_evaluation }
        MERGE_STATS_EVAL (stats_and_evaluation)
        stats_and_evaluation_summary = MERGE_STATS_EVAL.out.csv
        ch_versions                  = ch_versions.mix(MERGE_STATS_EVAL.out.versions)
    }

    //
    // MODULE: Shiny
    //
    if (!params.skip_shiny) {
        shiny_app = Channel.fromPath(params.shiny_app)
        PREPARE_SHINY (stats_and_evaluation_summary, shiny_app)
        ch_versions = ch_versions.mix(PREPARE_SHINY.out.versions)
    }


    if (!params.skip_visualisation) {
        VISUALIZATION (
            ch_alignment_output,
            MSA_GUIDETREE.out.tree,
            ch_optional_data
        )
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'multiplesequencealign_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params                        = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    PREPARE_MULTIQC (stats_and_evaluation_summary)
    ch_multiqc_files                      = ch_multiqc_files.mix(PREPARE_MULTIQC.out.multiqc_table.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


