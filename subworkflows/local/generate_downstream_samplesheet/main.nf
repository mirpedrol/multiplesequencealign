/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: GENERATE A DOWNSTREAM SAMPLESHEET FOR EVALUATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList      } from 'plugin/nf-schema'

workflow SAMPLESHEET_EVALUATION {
    take:
    ch_msa        // channel: [ meta, /path/to/file.aln ]
    ch_references // channel: [ meta, /path/to/file.aln ]
    ch_structures // channel: [ meta, /path/to/file.pdb ]
    outdir        // params.outdir

    main:
    // Try reading an existing samplesheet
    def samplesheet = file("${outdir}/downstream_samplesheets/evaluation.csv")
    def ch_existing_samplesheet = Channel.empty()
    if (samplesheet.exists()) {
        ch_existing_samplesheet = Channel.fromList(samplesheetToList(samplesheet, "${projectDir}/assets/schema_evaluate.json"))
    }
    // Create a channel with the new values for the samplesheet
    def ch_info_for_samplesheet = ch_msa
        .join(ch_references, by: 0, remainder: true)
        .join(ch_structures, by: 0, remainder: true)
        .map { meta, msa, reference, structure ->
            [id: meta.id, msa: msa, reference: reference, structure: structure]
        }
    // Join both channels
    ch_existing_samplesheet
        .mix(ch_info_for_samplesheet)
        .unique()
        .set { ch_list_for_samplesheet }

    channelToSamplesheet(ch_list_for_samplesheet, "${outdir}/downstream_samplesheets/evaluation")
}

workflow SAMPLESHEET_STATS {
    take:
    stats_summary // channel: [ meta, /path/to/file.csv ]
    outdir

    main:
    // Try reading an existing samplesheet
    def samplesheet = file("${outdir}/downstream_samplesheets/stats.csv")
    def ch_existing_samplesheet = Channel.empty()
    if (samplesheet.exists()) {
        ch_existing_samplesheet = Channel.fromList(samplesheetToList(samplesheet, "${projectDir}/assets/schema_stats.json"))
    }
    // Create a channel with the new values for the samplesheet
    def ch_info_for_samplesheet = stats_summary
        .map { meta, csv ->
            [id: meta.id, stats: csv]
        }
    // Join both channels
    ch_existing_samplesheet
        .mix(ch_info_for_samplesheet)
        .unique()
        .set { ch_list_for_samplesheet }

    channelToSamplesheet(ch_list_for_samplesheet, "${outdir}/downstream_samplesheets/stats")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW CALLING SPECIFIC SAMPLESHEET GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_DOWNSTREAM_SAMPLESHEETS {
    take:
    evaluation_msa        // channel: [ meta, /path/to/file.aln ]
    evaluation_references // channel: [ meta, /path/to/file.aln ]
    evaluation_structures // channel: [ meta, /path/to/file.pdb ]
    stats_summary         // channel: [ meta, /path/to/file.csv ]
    outdir                // params.outdir

    main:
    SAMPLESHEET_EVALUATION(
        evaluation_msa,
        evaluation_references,
        evaluation_structures,
        outdir
    )
    SAMPLESHEET_STATS(
        stats_summary,
        outdir
    )
}

// Input can be any channel with a dictionary
def channelToSamplesheet(ch_list_for_samplesheet, path) {
    ch_list_for_samplesheet
        .first()
        .map { it -> it.keySet().join(",") }
        .concat(ch_list_for_samplesheet.map { it -> it.values().join(",").replace("null", "") })
        .collectFile(
            name: "${path}.csv",
            newLine: true,
            sort: false
        )
}
