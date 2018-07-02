"""
Project: cider
File: fusions

Generate reads in FASTQ format for 1,2 or 3 genes.


 -g,--gene-model <arg>                   *Required* Path to gene model
                                         file in refFlat format
 -h,--help                               print usage info
 -v,--version                            Display version info

==================================================================
Output
==================================================================
 -f,--fasta-output <arg>                 File name of FASTA output
 -t,--text-output <arg>                  File name of text output
 -r,--reference <arg>                    Path to indexed reference genome
                                         fasta file (.fai). Required for
                                         FASTA output.

==================================================================
Fusion Types
==================================================================
 -j,--tri-fusion <arg>                   Number of fusions with three
                                         genes
 -n,--fusions <arg>                      Number of fusions to generate
                                         using two randomly selected genes
 -s,--self-fusion <arg>                  Number of self-fusions (fusions
                                         with single gene)
 -x,--read-through <arg>                 Number of read through fusions
 -y,--intra-chrom <arg>                  Number of intra-chromosome
                                         fusions (fusions within single
                                         chrom)

==================================================================
Fusion Options
==================================================================
 -a,--auto-correct-orientation           Auto correct orientation of
                                         fusion sequence if genes are
                                         located on different strands
 -c,--cds-only                           Only include CDS exons
 -d,--out-of-frame                       Allow fusions outside of reading
                                         frames. By default the reading
                                         frame is preserved
 -e,--keep-exon-boundary                 Generate fusion breaks on exon
                                         boundaries only
 -u,--foreign-insertion-length <arg>     Maxium length of randomly
                                         generated sequence to insert
                                         between fusion breakpoints
 -w,--foreign-insertion-perecent <arg>   Percent of fusions to insert
                                         foreign sequence between fusion
                                         breakpoints

==================================================================
Gene Selection
==================================================================
 -1,--gene1 <arg>                        Filter for gene1
 -2,--gene2 <arg>                        Filter for gene2
 -3,--gene3 <arg>                        Filter for gene3
 -b,--background-reads <arg>             Path to BAM file containing
                                         background reads. Genes will be
                                         selected for fusions according to
                                         the read profile of the
                                         background reads.
 -k,--rpkm-cutoff <arg>                  RPKM cutoff when using background
                                         BAM file. Genes below the cutoff
                                         will be ignored
 -l,--limit <arg>                        Limit all fusions to specific
                                         geneId, transcriptId, or chrom
 -m,--gene-selection-method <arg>        Method to use when selecting
                                         genes for fusions:
                                         uniform|empirical|binned
 -p,--threads <arg>                      Number of threads to spawn when
                                         processing background BAM file

==================================================================
Convert GTF/GFF gene model to refFlat format for use with Fusim
==================================================================
 -i,--gtf <arg>                          Input GTF file for conversion
 -o,--output <arg>                       Output refFlat file for
                                         conversion
 -z,--convert                            Convert GTF/GFF to refFlat
                                         (genePred) format


Contact:
---------
salturki@gmail.com

16:13 
2/21/15
"""

import os
import click
from ampsim import download
import pysam
import gzip
from itertools import combinations
from ampsim.utils import create_manifest

# DETANGO = CONFIG.get('spring', 'detango')

# CS_indexer = CONFIG.get('chimerascan', 'chimerascan_indexer')
# CS_scanner = CONFIG.get('chimerascan', 'chimerascan_run')
# CS_html_table = CONFIG.get('chimerascan', 'chimerascan_html_table')
# CS_bowtie_index = CONFIG.get('chimerascan', 'bowtie_index')
# CS_bowtie = CONFIG.get('chimerascan', 'bowtie')
# CS_bowtie_builder = CONFIG.get('chimerascan', 'bowtie_builder')
# CS_bowtie_folder = CONFIG.get('chimerascan', 'bowtie_folder')


# Get the resources required to run various steps
hg19 = download.Hg19()
samtools = download.Samtools()
art = download.Art()
bedtools = download.Bedtools()
vt = download.Vt()
bcftools = download.Bcftools()
fusim = download.Fusim()
gene_model = download.RefSeq()


FLANKING_SIZE = 200  # add flanking seq to the gene boundaries when generating normal ref fasta file


def generate_normal_reference(genes, output, filename):
    """
    For each gene in genes, get the largest transcript boundaries and add flanking 200bp on each side to generate
    a fasta reference file
    :param genes: List of 1,2 or 3 genes
    :return: A path to a temporary normal reference file
    """
    genome_fasta = pysam.Fastafile(hg19.full_path)
    holder = dict(zip(genes, [{} for x in genes if x is not None and x != '']))
    o_path = os.path.join(output, '%s_normal_ref.fasta' % filename)
    o = open(o_path, 'w')
    f = open(gene_model.full_path, 'r')
    for line in f:
        line = line.strip().split("\t")
        gene = line[0]
        transcript = line[1]
        chrom = line[2]
        start = min([int(line[4]), int(line[5])])
        end = max([int(line[4]), int(line[5])])
        diff = end - start
        coding_starts = line[9]
        coding_ends = line[10]
        if gene in genes:
            if gene not in holder:
                holder[gene] = {}

            holder[gene][diff] = (chrom, start, end, transcript, coding_starts, coding_ends)

    for gene in holder:
        chrom, start, end, transcript, coding_starts, coding_ends = holder[gene][max(holder[gene], key=holder[gene].get)]
        coding_exons = zip([int(x) for x in coding_starts.split(",") if x != ''],
                           [int(x) for x in coding_ends.split(",") if x != ''])
        for start, end in sorted(coding_exons):
            start -= FLANKING_SIZE
            end += FLANKING_SIZE
            o.write(">%s_%d_%d_%s\n" % (chrom, start, end, transcript))
            fasta_seq = genome_fasta.fetch(chrom, start, end)
            o.write(fasta_seq + "\n")
    o.close()

    return o_path


def generate_fusion_reference(genes, output, frameness, keep_exon, add_insertion, filename):
    """
    """
    ke = fr = ins = ''
    if keep_exon:
        ke = "--keep-exon-boundary"

    if frameness:
        fr = "--out-of-frame"

    if add_insertion:
        ins = '--foreign-insertion-perecent 1 --foreign-insertion-length 30'

    o_fasta_path = os.path.join(output, '%s_ref_normal.fasta' % filename)
    o_text_path = os.path.join(output, '%s_fusion_truth_set.txt' % filename)
    genes_arguments = []
    for idx, gene in enumerate(genes):
        genes_arguments.append('--gene%d=%s' % (idx + 1, gene))

    cmd = 'java -jar %s --gene-model=%s --fusions=1 --reference=%s %s --fasta-output=%s --text-output=%s %s %s %s ' % \
          (fusim.full_path, gene_model.full_path, hg19.full_path,
           " ".join(genes_arguments), o_fasta_path, o_text_path, ke, fr, ins)
    print cmd
    os.system(cmd)
    return o_fasta_path


def generate_fastq(ref_path, output_dir, f_prefix, coverage):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    o_path = os.path.join(output_dir, f_prefix)
    cmd = '%s -i %s -o %s  -l 150 -f %2.f -p -m 400 -s 10 --rndSeed 9999' % (art.full_path, ref_path, o_path, coverage)
    print cmd
    os.system(cmd)
    return o_path


def merge_fastq(normal_fastq, fusion_fastq, output, filename):
    r1_path = os.path.join(output, "%s_fusion1.fastq.gz" % filename)
    r2_path = os.path.join(output, "%s_fusion2.fastq.gz" % filename)
    r1 = gzip.open(r1_path, 'w')
    r2 = gzip.open(r2_path, 'w')
    for path in [normal_fastq, fusion_fastq]:
        f1 = open(path + "1.fq", 'r').read()
        f2 = open(path + "2.fq", 'r').read()
        r1.write(f1)
        r2.write(f2)
    r1.close()
    r2.close()
    return r1_path, r2_path


def generate_multiple_fusions():
    holder = {}
    genes = ["ALK", "EML4"]#, "EML4"] #,"TMPRSS2", "ERG"]
    keep_exons = [0]  #, 1]
    frameness = [0]   #, 1]
    add_insertion = [0]  #, 1]
    fusion_fraction = [0.4]  #[0.2, 0.1, 0.05, 0.01]
    gene_pairs = combinations(genes, 2)
    for gp in gene_pairs:
        for ex in keep_exons:
            for fr in frameness:
                for ff in fusion_fraction:
                    for ins in add_insertion:
                        key = (gp[0], gp[1], ex, fr, ff, ins)
                        holder[key] = {"gene1": gp[0],
                                       "gene2": gp[1],
                                       "keep_exon": ex,
                                       "frameness": fr,
                                       "fraction": ff,
                                       "add_insertion": ins,
                                       }
                    print key
    print len(holder)
    return holder


def generate_manifest(fastq1, fastq2):
    pass


# @click.group(invoke_without_command=True)
# @click.option('--gene1', '-1', help='Select gene 1')
# @click.option('--gene2', '-2', help='Select gene 2')
# @click.option('--gene3', '-3', help='Select gene 3')
# @click.option('--out-of-frame', '-f', is_flag=True, default=False, help='Allow fusions outside of reading frames')
# @click.option('--keep-exon-boundary', '-e', is_flag=True, default=False,
#               help='Generate fusion breaks on exon boundaries only')
# @click.option('--output', '-o', help='Output directory')
def cli(gene1, gene2, gene3, frameness, keep_exon, fusion_fraction,
        add_insertion, total_coverage, output, common_filename):
    """[Simulator] Fusion generator."""
    normal_coverage = total_coverage * (1. - fusion_fraction)
    fusion_coverage = total_coverage * fusion_fraction
    normal_ref = generate_normal_reference([gene1, gene2, gene3], output, common_filename)
    fusion_ref = generate_fusion_reference([gene1, gene2, gene3], output,
                                           keep_exon, frameness, add_insertion, common_filename)
    normal_fastq = generate_fastq(normal_ref, output, 'normal', normal_coverage)
    fusion_fastq = generate_fastq(fusion_ref, output, 'fusion', fusion_coverage)
    merged1, merged2 = merge_fastq(normal_fastq, fusion_fastq, output, common_filename)
    # chimerascan_bedpe = run_chimerascan(merged1, merged2, output)
    # print chimerascan_bedpe
    # generate_manifest(merged1, merged2)
    # run_detango(merged1, merged2, output)
    return merged1, merged2


if __name__ == '__main__':
    # cli ()
    total_coverage = 200
    candidates = generate_multiple_fusions()

    o = open('/Users/saeed/Desktop/temp/fastq_list.txt', 'w')
    output = '/Users/saeed/Desktop/temp/syn_fusions/fusim/'
    fastq_pairs = {}
    # o = open('/PHShome/sha13/tmp_files/syn_fusions/fastq_list.txt', 'w')
    # output = '/PHShome/sha13/tmp_files/syn_fusions'

    for key in candidates:
        gene3 = None
        gene1, gene2, keep_exon, frameness, fusion_fraction, add_insertion = key
        common_filename = "%s_%s_ke%s_fr%d_ff%.2f_ins%d" % key

        # create a folder for the results
        output = os.path.join(output, common_filename)
        os.system('rm -rf %s' % output)
        os.system('mkdir -p %s' % output)

        print "Processing ", key
        f1, f2 = cli(gene1, gene2, gene3,
                     frameness, keep_exon, fusion_fraction, add_insertion,
                     total_coverage, output, common_filename)
        fastq_pairs[common_filename] = {'somatic': [f1, f2]}
        o.write("\t".join([f1, f2]) + "\n")
    o.close()
    fusion_manifest_path = create_manifest(os.path.join(output, 'manifests'),
                                           pipeline_name='solid_fusion',
                                           paths=fastq_pairs)
    output = os.path.dirname(output)
    os.system('for i in `find %s -name *.aln`; do echo $i; rm $i; done' % output)
    os.system('for i in `find %s -name *.fq`; do echo $i; rm $i; done' % output)
    print "compressing all files to archive.tar.gz"
    os.system('cd %s && tar -zcvf archive.tar.gz *.*' % output)
    print(fusion_manifest_path)
    print("cider solid_fusion -m {} -n solid_fusion_run1 -b".format(fusion_manifest_path))
