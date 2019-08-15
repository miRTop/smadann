#!/usr/bin/env nextflow


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --csv                         CSV file with first column read file and second column protocol
      --protocol                    Library preparation protocol. Default: "illumina". Can be set as "illumina", "nextflex", "qiaseq" or "cats"

    Trimming options
      --three_prime_adapter         3’ Adapter to trim. Default: None
      --min_length [int]            Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18
      --clip_R1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1
      --three_prime_clip_R1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed

    Other options
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


if (params.csv){
    Channel
        .fromPath(params.csv, checkIfExists: true)
        .splitCsv(header:false)
        .map { row -> [file(row[0]), row[1]] }
        .ifEmpty { exit 1, "params.csv was empty - no input files supplied" }
        .set{raw_reads_trimgalore}
} else if (params.reads) {
    Channel
            .fromPath( params.reads )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
            .set {raw_reads_trimgalore}
}  

ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    label 'process_low'
    tag "$reads"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
     set file(reads), val(protocol) from raw_reads_trimgalore

    output:
    file '*.gz' 
    file '*trimming_report.txt' into trimgalore_results
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    // Define regular variables so that they can be overwritten
    clip_R1 = params.clip_R1
    three_prime_clip_R1 = params.three_prime_clip_R1
    three_prime_adapter = params.three_prime_adapter

    // Presets
    if (protocol == "illumina"){
        clip_R1 = 0
        three_prime_clip_R1 = 0
        three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
    } else if (protocol == "nebnext"){
        clip_R1 = 0
        three_prime_clip_R1 = 0
        three_prime_adapter = "AGATCGGAAGAGCACACGTCT"
    } else if (protocol == "nextflex"){
        clip_R1 = 4
        three_prime_clip_R1 = 4
        three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
    } else if (protocol == "qiaseq"){
        clip_R1 = 0
        three_prime_clip_R1 = 0
        three_prime_adapter = "AACTGTAGGCACCATCAAT"
    } else if (protocol == "cats"){
        clip_R1 = 3
        three_prime_clip_R1 = 0
        // three_prime_adapter = "GATCGGAAGAGCACACGTCTG"
        three_prime_adapter = "AAAAAAAA"
    } else {
        //custom protocol 
        clip_R1 = params.clip_R1
        three_prime_clip_R1 = params.three_prime_clip_R1
        three_prime_adapter = params.three_prime_adapter
        protocol = params.protocol
    }
    tg_length = "--length 15"
    c_r1 = clip_R1 > 0 ? "--clip_R1 ${clip_R1}" : ''
    tpc_r1 = three_prime_clip_R1 > 0 ? "--three_prime_clip_R1 ${three_prime_clip_R1}" : ''
    """
    trim_galore --adapter ${three_prime_adapter} $tg_length $c_r1 $tpc_r1 --max_length 40 --gzip $reads --fastqc
    """
}


/*
 * STEP 8 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.collect()
    file ('trim_galore/*') from trimgalore_results.collect()
    file multiqc_config from ch_multiqc_config

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc . --config $multiqc_config -f -m cutadapt -m fastqc 
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/smadann] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/smadann] FAILED: $workflow.runName"
    }
}
