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
include { STATS                            } from '../subworkflows/local/stats'
include { CREATE_TCOFFEETEMPLATE           } from '../modules/local/create_tcoffee_template'
include { GENERATE_DOWNSTREAM_SAMPLESHEETS } from '../subworkflows/local/generate_downstream_samplesheet/main'
include { EXTRACT_STRUCTURES               } from '../subworkflows/local/extract_structures/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR                          } from '../modules/nf-core/untar/main'
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
    ch_input       // channel: [ meta, path(sequence.fasta), path(reference.fasta), path(pdb_structures.tar.gz), path(templates.txt) ]
    ch_versions    // channel: [ path(versions.yml) ]
    outdir         // params.outdir
    alignment      // params.alignment
    alignment_args // params.alignment_args
    guidetree      // params.guidetree
    guidetree_args // params.guidetree_args
    treealign      // params.treealign
    treealign_args // params.treealign_args

    main:
    def ch_multiqc_files             = Channel.empty()
    def stats_summary                = Channel.empty()

    ch_input
        .map {meta, fasta, ref, structure, template ->
            def alignment_clean = alignment ? alignment.replace("/", "-") : ""
            def alignment_args_clean = alignment_args ? alignment_args.toString().trim().replace("  ", " ").replace(" ", "-").replaceAll("==", "-").replaceAll("\\s+", "") : ""
            def guidetree_clean = guidetree ? guidetree.replace("/", "-") : ""
            def guidetree_args_clean = guidetree_args ? guidetree_args.toString().trim().replace("  ", " ").replace(" ", "-").replaceAll("==", "-").replaceAll("\\s+", "") : ""
            def treealign_clean = treealign ? treealign.replace("/", "-") : ""
            def treealign_args_clean = treealign_args ? treealign_args.toString().trim().replace("  ", " ").replace(" ", "-").replaceAll("==", "-").replaceAll("\\s+", "") : ""
            [
                [
                    "id": meta.id,
                    "alignment": alignment_clean, "alignment_args": alignment_args_clean,
                    "guidetree": guidetree_clean, "guidetree_args": guidetree_args_clean,
                    "treealign": treealign_clean, "treealign_args": treealign_args_clean
                ],
                fasta, ref, structure, template
            ]
        }
        .multiMap { meta, fasta, ref, structure, template ->
            seqs: [ meta, fasta ]
            refs: [ meta, ref ]
            structures: [ meta, structure ]
            templates: [ meta, template ]
        }
        .set { ch_input_multi }

    ch_input_multi.refs
        .filter { meta, ref ->
            ref.size() > 0
        }
        .set { ch_refs }
    ch_input_multi.structures
        .filter { meta, structure ->
            structure.size() > 0
        }
        .set { ch_structures_tar }
    ch_input_multi.templates
        .filter { meta, template ->
            template.size() > 0
        }
        .set { ch_templates }

    // ----------------
    // STRUCTURES
    // ----------------
    EXTRACT_STRUCTURES(ch_structures_tar)
    ch_structures = EXTRACT_STRUCTURES.out.structures

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
            ch_input_multi.seqs,
            ch_structures
        )
        ch_versions   = ch_versions.mix(STATS.out.versions)
        stats_summary = stats_summary.mix(STATS.out.stats_summary)
    }

    def msa_alignment = Channel.empty()

    if (guidetree && treealign) {
        //
        // Compute tree
        //
        MSA_GUIDETREE (
            ch_input_multi.seqs
        )
        ch_versions = ch_versions.mix(MSA_GUIDETREE.out.versions)

        // Prepare channels for treealign to make sure the correct tree is used for the respective alignment
        ch_input_multi.seqs
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
        msa_alignment = msa_alignment.mix(MSA_TREEALIGN.out.alignment)
    }

    if (alignment) {
        //
        // Align
        //
        MSA_ALIGNMENT (
            ch_input_multi.seqs
        )
        ch_versions = ch_versions.mix(MSA_ALIGNMENT.out.versions)
        msa_alignment = msa_alignment.mix(MSA_ALIGNMENT.out.alignment)
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // SUBWORKFLOW: Generate samplesheets for downstream workflows
    //
    GENERATE_DOWNSTREAM_SAMPLESHEETS (
        msa_alignment,
        ch_refs,
        ch_structures_tar,
        stats_summary,
        outdir
    )

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
        ch_multiqc_files                          = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                          = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                          = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

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


