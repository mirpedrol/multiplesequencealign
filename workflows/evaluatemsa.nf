/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES
include { MULTIQC                } from '../modules/local/multiqc'
include { PREPARE_MULTIQC } from '../modules/local/prepare_multiqc'
include { PREPARE_SHINY   } from '../modules/local/prepare_shiny'

//SUBWORKFLOWS
include { EVALUATE               } from '../subworkflows/local/evaluate'

// FUNCTIONS
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

workflow EVALUATEMSA {

    take:
    msa_alignment
    ch_refs
    ch_structures_template
    stats_summary
    ch_versions

    main:
    evaluation_summary           = Channel.empty()
    stats_and_evaluation_summary = Channel.empty()
    ch_multiqc_table             = Channel.empty()

    //
    // Evaluate the quality of the alignment
    //
    if (!params.skip_eval) {
        EVALUATE (msa_alignment, ch_refs, ch_structures_template)
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
        ch_shiny_stats = PREPARE_SHINY.out.data.toList()
        ch_versions = ch_versions.mix(PREPARE_SHINY.out.versions)
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
    multiqc_out = Channel.empty()
    if (!params.skip_multiqc && (!params.skip_stats || !params.skip_eval)) {
        ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

        PREPARE_MULTIQC (stats_and_evaluation_summary)
        ch_multiqc_table = ch_multiqc_table.mix(PREPARE_MULTIQC.out.multiqc_table.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_multiqc_table
        )
        multiqc_out = MULTIQC.out.report.toList()
    }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
    multiqc  = multiqc_out // channel: [ path(multiqc_report.html) ]
}
