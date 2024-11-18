#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/multiplesequencealign
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/multiplesequencealign
    Website: https://nf-co.re/multiplesequencealign
    Slack  : https://nfcore.slack.com/channels/multiplesequencealign
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIPLESEQUENCEALIGN   } from './workflows/multiplesequencealign'
include { EVALUATEMSA             } from './workflows/evaluatemsa'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_multiplesequencealign_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_multiplesequencealign_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_MULTIPLESEQUENCEALIGN {

    take:
    samplesheet // channel: samplesheet read in from --input
    outdir

    main:
    def ch_versions = Channel.empty()

    //
    // WORKFLOW: Run pipeline
    //
    MULTIPLESEQUENCEALIGN (
        samplesheet,
        ch_versions,
        outdir
    )

    emit:
    multiqc_report =  MULTIPLESEQUENCEALIGN.out.multiqc

}

workflow NFCORE_EVALUATEMSA {

    take:
    evaluate_samplesheet // channel: [ /path/to/file.csv ]
    stats_summary        // channel: [ meta, /path/to/file.csv ]
    outdir

    main:
    def ch_versions = Channel.empty()

    //
    // WORKFLOW: Run evaluation pipelines
    //
    EVALUATEMSA (
        evaluate_samplesheet,
        stats_summary,
        ch_versions,
        outdir
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    if (params.evaluate) {
        // WORKFLOW: Run evaluation workflow
        NFCORE_EVALUATEMSA (
            "${params.outdir}/downstream_samplesheets/evaluation.csv",
            "${params.outdir}/downstream_samplesheets/stats.csv",
            params.outdir
        )
    } else {
        //
        // WORKFLOW: Run main workflow
        //
        NFCORE_MULTIPLESEQUENCEALIGN (
            PIPELINE_INITIALISATION.out.samplesheet,
            params.outdir
        )
    }

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_MULTIPLESEQUENCEALIGN.out.multiqc_report,
        "${params.outdir}/shiny_app",
        "${params.outdir}/pipeline_info",
        params.shiny_trace_mode,
        params.evaluate
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
