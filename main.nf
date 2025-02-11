/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This path is part of 'RNA-Toy'.
 *
 *   RNA-Toy is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNA-Toy is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNA-Toy.  If not, see <http://www.gnu.org/licenses/>.
 */
 
/* 
 * Proof of concept Nextflow based RNA-Seq pipeline
 * 
 * Authors:
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Emilio Palumbo <emiliopalumbo@gmail.com> 
 */ 

 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = 'results'


/*
* This is written to the .nextflow.log.X file in the directory that the nextflow command is called
* from (*not* $baseDir)
*/
log.info """\
         R N A T O Y   P I P E L I N E    
         =============================
         genome: ${params.genome}
         annot : ${params.annot}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()

/*
 * the reference genome file
 * Seems like file() is not a function, only file(). Can't find that info anywhere in the docs!
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)
 
/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 * see https://www.nextflow.io/docs/latest/channel.html?highlight=fromfilepairs#fromfilepairs
 * Q: What does this do if one/some of the files are missing their pair?
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs } 
 
/*
 * Step 1. Builds the genome index required by the mapping process
 * Q: How can you make this check whether the index file exists and skip the 
 * indexing process if so?
 */
process buildIndex {
    tag "$genome_file.baseName"

    input:
    path genome from genome_file
     
    output:
    path 'genome.index*' into genome_index
       
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}
 
/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 * The "set pair_id" call, derived from fromFilePairs,seems to be how the sample id is kept track of when 
 * there are multiple input files per sample. The syntax is kind of weird but nicer than Snakemake!
 *
 * Q: Why do you have to mv the tophat output into the current dir?
 *    I assume because that is where the process looks for the output file, and possibly it can
 *    *only* look in the current dir and not in any other. It's also possible this was fixed with 
 *    the new path qualifier.
 * 
 * Q: Why isn't the genome or genome index file specified in the tophat2 options?
 * A: It is the base file name of both files that is specified.
 * Q: What is the "genome.index" specification in the tophat2 options?
 * A: That is the actual base file name of the genome + idx files. How does this get into the current
 *    dir? It must be placed there by the "input: file index from genome_index" call.
 *    Mechanically, the genome_index channel contains the "genome.index*" files. Calling
 *    the file from the channel in the input must write the file to the working directory of
 *    the process that calls it.
 * Tophat docs: 
 * Usage: tophat [options]* <genome_index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]
 * Q: What is the "set pair_id" language do?
 * A: set is an operator (see https://www.nextflow.io/docs/latest/operator.html#set) that
 *    assigns a channel to a variable name. In this case, it creates a channel called pair_id
 *    that comes from the first field of the channel created by fromFilePairs. The syntax is weird though...
 * Q: Is pair_id a channel? Or some other kind of variable? 
 * Q: Where is task.cpus specified?
 * A: It can be manually specified by the "cpus X" statement in a process, otherwise it appears to
 *    be automatically determined by available resources. If it is set to more than available cores,
 *    and error will be thrown.
 */
process mapping {
    tag "$pair_id" // tag = User provided identifier associated this task.
     
    input:
    path genome from genome_file 
    path annot from annotation_file
    path index from genome_index
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, "accepted_hits.bam" into bam
 
    """
    tophat2 -p ${task.cpus} --GTF $annot genome.index $reads
    mv tophat_out/accepted_hits.bam .
    """
}
  
/*
 * Step 3. Assembles the transcript by using the "cufflinks" tool
 * Q: Does publishDir require the directory to exist? Or will it make it?
 * A: publishDir *will* create the directory for you (great!)
 * See https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir
 */
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'  // copy the output files to this directory,
                                            // rather than only storing in process' tmp dir
       

    input:
    path annot from annotation_file
    set pair_id, file(bam_file) from bam
     
    output:
    set pair_id, file('transcript_*.gtf') into transcripts
 
    """
    cufflinks --no-update-check -q -p $task.cpus -G $annot $bam_file
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
