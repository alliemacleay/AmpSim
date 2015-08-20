"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import subprocess as sp

import click
import pysam
from pybedtools import BedTool

from ampsim.utils import create_dir, bgzip, create_manifest
from ampsim.tools.variants import get_synthetic_variants


picard = '/Users/saeed/bioApps/picrad_htsjdk/picard.jar'
bam_surgeon_addsnv = "/Users/saeed/bioApps/bamsurgeon/addsnv.py"
bam_surgeon_addindel = "/Users/saeed/bioApps/bamsurgeon/addindel.py"

flanking_size = 0  # to add on both ends of each target in order to generate NGS reads


def generate_target_reference(genome_path, bed_path, output_dir):
    """
    Generate a fasta file to hold the reference sequence for all regions in a BED file
    :param genome_path: path to the genome reference file
    :param bed_path:
    :return: path to reference
    """
    ref_dir = create_dir(os.path.join(output_dir, 'reference'))
    target_genome_path = os.path.join(ref_dir, 'targets_ref.fasta')
    ref_genome = pysam.FastaFile(genome_path)
    target_genome = open(target_genome_path, 'w')
    bed = BedTool(bed_path)
    for record in bed:
        chrom, start, end = record[:3]
        start = int(start) + flanking_size
        end = int(end) + flanking_size

        contig_name = ">" + "_".join([chrom, str(start), str(end)])
        target_genome.write(contig_name + "\n")
        target_genome.write(ref_genome.fetch(chrom, start, end) + "\n")
    ref_genome.close()
    target_genome.close()
    return target_genome_path


def run_art(art_illumina, target_genome_path, coverage, read_length, output_dir):
    """
    Generate NGS reads (pair-end, amplicon), compress with gzip and return paths to R1 and R2
    :param art_illumina: path to art_illumina bin
    :param target_genome_path: Path to genome reference
    :param coverage: The coverage folds
    :param read_length: The length of reads.
    :param output_dir: The output directory
    :return: paths to fq1 and fq2
    """
    file_basename = os.path.basename(target_genome_path).split('.')[0]  # normal or somatic
    # The following line needs a '/' at the end of the path for ART otherwise it writes the files to the parent dir
    fastq_dir = create_dir(os.path.join(output_dir, 'fastq')) + '/normal'
    cmd = '%s --noALN --paired --in %s --len %d --fcov %d --mflen 150 --sdev 1 --rndSeed 0 --out %s' \
          % (art_illumina, target_genome_path, read_length, coverage, fastq_dir)
    sp.check_call(cmd, shell=True)
    fastq1_path = bgzip(fastq_dir + "1.fq")
    fastq2_path = bgzip(fastq_dir + "2.fq")
    return fastq1_path, fastq2_path


def generate_bam(ref_path, fastq1_path, fastq2_path, wd):
    print("Generating NGS BAM (normal) ..")
    if not os.path.isdir(os.path.join(wd, "bam")):
        os.mkdir(os.path.join(wd, "bam"))
    output_path = os.path.join(wd, 'bam', 'normal')
    sp.check_call('bwa mem -a -M -O 6 %s %s %s > %s.sam' %
                  (ref_path, fastq1_path, fastq2_path, output_path), shell=True)
    print("SortSam (normal) ..")
    sp.check_call('java -jar %s SortSam INPUT=%s.sam OUTPUT=%s.bam SORT_ORDER=coordinate' %
                  (picard, output_path, output_path), shell=True)
    # print("AddOrReplaceReadGroups (normal) ..")
    # sp.check_call('java -jar %s AddOrReplaceReadGroups INPUT=%s.bam OUTPUT=%s.ordered.bam SORT_ORDER=coordinate '
    #               'RGID=normal RGLB=normal RGPL=Illumina RGSM=normal RGPU=normal ' %
    #               (picard, output_path, output_path), shell=True)

    print("Indexing the BAM file")
    sp.check_call('samtools index -n %s %s.ordered' % (output_path, output_path), shell=True)
    bam_path = output_path + ".ordered.bam"
    os.remove(output_path + ".bam")
    os.remove(output_path + ".sam")
    return bam_path


def bam_2_fastq(bam_path, wd):
    filename = os.path.basename(bam_path)[:-4]
    fastq1 = os.path.join(wd, "fastq", filename + "1.fq")
    fastq2 = os.path.join(wd, "fastq", filename + "2.fq")
    cmd = 'java -jar %s SamToFastq INPUT=%s FASTQ=%s SECOND_END_FASTQ=%s' % (picard, bam_path, fastq1, fastq2)
    sp.check_call(cmd, shell=True)
    return bgzip(fastq1), bgzip(fastq2)


def save_bam_surgeon_variants(snvs, indels, all, wd):
    create_dir(os.path.join(wd, 'bamsurgeon'))
    snvs_path = os.path.join(wd, 'bamsurgeon', 'snvs.txt')
    indels_path = os.path.join(wd, 'bamsurgeon', 'indels.txt')
    all_path = os.path.join(wd, 'bamsurgeon', 'all.txt')
    open(snvs_path, 'w').write(snvs)
    open(indels_path, 'w').write(indels)
    open(all_path, 'w').write(all)
    return snvs_path, indels_path


def run_bam_surgeon(normal_bam_path, genome_ref, surgeon_snv_path, surgeon_indel_path, wd):
    print("Running BAM surgeon for SNVs variants ..")
    tmp_dir = os.path.join(wd, "bamsurgeon_tmp", "logs")
    if not os.path.isdir(tmp_dir):
        os.system('mkdir -p %s' % tmp_dir)
    somatic_bam_snv_path = normal_bam_path.replace('normal.ordered.', 'somatic_SNV.')
    cmd = '%s --aligner mem --picardjar %s --procs 4 --varfile %s --bamfile %s --reference %s --outbam %s --tmpdir %s' %\
          (bam_surgeon_addsnv, picard, surgeon_snv_path,
           normal_bam_path, genome_ref, somatic_bam_snv_path, tmp_dir)
    print(cmd)
    sp.check_call(cmd, shell=True)

    somatic_bam_path = somatic_bam_snv_path.replace('somatic_SNV.', 'somatic_SNV_ord.')
    cmd = 'java -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=coordinate' %\
          (picard, somatic_bam_snv_path, somatic_bam_path)
    print(cmd)
    sp.check_call(cmd, shell=True)
    sp.check_call("samtools index %s" % somatic_bam_path, shell=True)

    print("Running BAM surgeon for INDELs variants ..")
    somatic_bam_snv_indel_path = somatic_bam_snv_path.replace('somatic_SNV.', 'somatic_SNV_INDEL.')
    somatic_bam_snv_indel_sorted_path = somatic_bam_snv_indel_path.replace('somatic_SNV_INDEL.', 'somatic.')
    cmd = '%s --aligner mem --picardjar %s --procs 4 --varfile %s --bamfile %s --reference %s --outbam %s --tmpdir %s' %\
          (bam_surgeon_addindel, picard, surgeon_indel_path,
           somatic_bam_path, genome_ref, somatic_bam_snv_indel_path, tmp_dir)
    print(cmd)
    sp.check_call(cmd, shell=True)
    cmd = 'java -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=coordinate' %\
          (picard, somatic_bam_snv_indel_path, somatic_bam_snv_indel_sorted_path)
    print(cmd)
    sp.check_call(cmd, shell=True)
    sp.check_call("samtools index %s" % somatic_bam_snv_indel_sorted_path, shell=True)

    os.remove(somatic_bam_snv_path.replace('somatic_SNV.', 'somatic_SNV_ord.'))
    os.remove(somatic_bam_snv_path)
    os.remove(somatic_bam_snv_indel_path)
    os.system('rm -rf {}'.format(tmp_dir))
    return somatic_bam_snv_indel_sorted_path


@click.group(invoke_without_command=True)
@click.option('--bed-path', '-b', help='Path to target regions in BED format')
@click.option('--output-dir', '-o', help='Path to output directory')
@click.option('--genome', '-g', help='Path to fasta genome reference')
@click.option('--art-illumina', '-a', help='Path to art_illumina software')
@click.option('--read-length', '-r', default=100, help='Read length')
@click.option('--coverage', '-c', default=200, help='Coverage')
def cli(bed_path=None, output_dir=None, genome=None, art_illumina=None, read_length=None, coverage=None):

    target_genome_path = generate_target_reference(genome, bed_path, output_dir)

    normal_fastq1_path, normal_fastq2_path = run_art(art_illumina, target_genome_path, coverage,
                                                     read_length, output_dir)

    normal_bam_path = generate_bam(genome, normal_fastq1_path, normal_fastq2_path, output_dir)

    snvs, indels, all_variants = get_synthetic_variants(bed_path, genome)

    snvs_path, indels_path = save_bam_surgeon_variants(snvs, indels, all_variants, output_dir)

    somatic_bam_path = run_bam_surgeon(normal_bam_path, genome, snvs_path, indels_path, output_dir)

    somatic_fastq1_path, somatic_fastq2_path = bam_2_fastq(somatic_bam_path, output_dir)

    snapshot_manifest_path = create_manifest(os.path.join(output_dir, 'manifests'), pipeline_name='snapshot',
                                             paths=[somatic_fastq1_path, somatic_fastq2_path])

    # onek_manifest_path = create_manifest(os.path.join(output_dir, 'manifests'), pipeline_name='onek',
    #                                      paths=[normal_fastq1_path, normal_fastq2_path,
    #                                             somatic_fastq1_path, somatic_fastq2_path])

    print(snapshot_manifest_path)
    # print(onek_manifest_path)

    os.system('source ~/.virtualenvs/cider2/bin/activate && cider snapshot -n tmp -m {} -f -r '
              .format(snapshot_manifest_path))


if __name__ == '__main__':
    cli()

