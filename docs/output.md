# nf-core/multiplesequencealign: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. **Input files summary**: (Optional) computation of summary statistics on the input fasta file, such as the average sequence similarity across the input sequences, their length, etc. Skip by `--skip_stats` as a parameter.

2. **Guide Tree**: (Optional) Renders a guide tree. This is meant to try different combinations of guide trees and assembly tools, if you are interested in a standard alignment procedure you can ignore this step.
3. **Align**: aligns the sequences.
4. **Evaluate**: (Optional) The obtained alignments are evaluated with different metrics: Sum Of Pairs (SoP), Total Column score (TC), iRMSD, Total Consistency Score (TCS), etc. Skip by passing `--skip_eval` as a parameter.
5. **Report**: Reports about the collected information of the runs are reported in a shiny app and a summary table in multiqc. Skip by passing `--skip_shiny` and `--skip_multiqc`.

## Input files summary

Statistics about the input files are collected and summarized into a final csv file.

<details markdown="1">
<summary>Output files</summary>

- `summary/stats/`
  - `complete_summary_stats.csv`: csv file containing the summary for all the statistics computed on the input file.
  - `complete_summary_stats_with_trace.csv`: csv file containing the content of complete_summary_stats merged with the information of the trace file. This will not be produced if `-resume` is used.
  - `sequences/`
    - `seqstats/*_seqstats.csv`: file containing the sequence input length for each sequence in the family defined by the file name. If `--calc_seq_stats` is specified.
    - `perc_sim/*_txt`: file containing the pairwise sequence similarity for all input sequences. If `--calc_sim` is specified.
  - `structures/` - `plddt/*_full_plddt.csv`: file containing the plddt of the structures for each sequence in the input file. If `--extract_plddt` is specified.
  </details>

## Trees

If you explicitly specifified (via the toolsheet) to compute guidetrees to be used by the MSA tool, those are stored in the **trees** directory.

<details markdown="1">
<summary>Output files</summary>

- `trees/`
  - `*/*.dnd`: guide tree files.

</details>

## Alignment

All MSA computed are stored in the **alignment** directory.

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `*/*.fa`: each subdirectory is named after the sample id. It contains all the alignments computed on it. The filename contains all the informations of the input file used and the tool.
    The file naming convention is:
    {Input*file}*{Tree}\_args-{Tree_args}\_{MSA}\_args-{MSA_args}.aln

</details>

## Evaluation

Files with the summary of the computed evaluation statistics are stored in the **evaluation** directory.

<details markdown="1">
<summary>Output files</summary>

- `evaluation/`
  - `tcoffee_irmsd/`: directory containing the files with the complete iRMSD files. If `--calc_irmsd` is specified.
  - `tcoffee_tcs/`: directory containing the files with the complete TCS files. If `--calc_tcs` is specified.
  - `complete_summary_eval.csv`: csv file containing the summary of all evaluation metrics for each input file.
  </details>

## Shiny App

A Shiny app is created to explore interactively your results. It can be found in the **shiny_app** folder.

A shiny app is prepared to visualize the summary statistics and evaluation of the produced alignments (skip with `--skip_shiny`).
To run the Shiny app use the following commands from the results directory:
`cd shiny_app`
`./run.sh`

<details markdown="1">
<summary>Output files</summary>

- `shiny_app/`
  - `run.sh`: executable to start the shiny app.
  - `*.py*`: shiny app files.
  - `*.csv`: csv file used by shiny app.
  - `trace.txt`: trace file used by shiny app.
  </details>

Be aware that you have to have [shiny](https://shiny.posit.co/py/) installed to access this feature.

### MultiQC

You can find the MultiQC report in the `multiqc` folder.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. A table with all the collected statistics and evaluted metrics is reported as well as all the versions of the tools used for the computation.

### Pipeline information

Extra information about the pipeline execution are stored in the **pipeline_info** folder.

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
