"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import subprocess as sp

import click
import pysam
import gzip
import csv
from collections import namedtuple

from ampsim.utils import create_dir, bgzip, create_manifest, vcf_parser, strip_chrom, get_var_type_size
from ampsim.tools.variants import get_synthetic_variants, get_ref
from ampsim import download
from ampsim.tools import venner



def get_resources():
    hg19 = download.Hg19()
    samtools = download.Samtools()
    art = download.Art()
    bedtools = download.Bedtools()
    vt = download.Vt()
    bcftools = download.Bcftools()
    ArtificialFastqGenerator_path = os.getenv("ARTFQ", '/Users/saeed/Downloads/ArtificialFastqGenerator/ArtificialFastqGenerator.jar')
    return hg19, samtools, art, bedtools, vt, bcftools, ArtificialFastqGenerator_path

# Get the resources required to run various steps
hg19, samtools, art, bedtools, vt, bcftools, ArtificalFastqGenerator_path = get_resources()

def edit_reference(target_start, target_seq, var_type, var_pos, var_ref, var_alt):
    """
    Given a sequence from normal human genome at a target region, edit it to add the variant sequence
    :param target_start:
    :param target_seq:
    :param var_type:
    :param var_pos:
    :param var_seq:
    :return:
    """
    target_seq = list(target_seq)
    relative_position = abs(target_start - var_pos)  # the position as 0-index in python data structures
    if var_type == 'SNV':

        new_sequence = target_seq
        new_sequence[relative_position] = var_alt
    elif var_type == 'DEL':
        relative_position = abs(target_start - var_pos)
        del_start = relative_position + 1
        del_end = del_start + len(var_ref) - 1  # take 1 away since var_ref include the a single normal ref base
        new_sequence = target_seq[:del_start] + target_seq[del_end:]
    elif var_type == 'INS':
        relative_position = abs(target_start - var_pos)
        # var_alt[1:] start from 1 index to ignore the REF nucleotide
        new_sequence = target_seq[:relative_position + 1] + list(var_alt[1:]) + target_seq[relative_position + 1:]
    return ''.join(new_sequence)


def generate_normal_somatic_references(variants, output_dir):
    """
    Generate a fasta file to hold the reference sequence for all regions in a BED file
    :param variants: list of target regions with their internal synthetic variants to be added to the normal reference
    :param output_dir: Path to save output for this method
    :return: path to reference
    """
    ref_dir = create_dir(os.path.join(output_dir, 'reference'))
    normal_ref_path = os.path.join(ref_dir, 'normal_ref.fasta')
    somatic_ref_path = os.path.join(ref_dir, 'somatic_ref.fasta')
    original_genome = pysam.FastaFile(hg19.full_path)
    normal_ref_file = open(normal_ref_path, 'w')
    somatic_ref_file = open(somatic_ref_path, 'w')
    for target in variants:
        chrom, target_start, target_end = target
        variant = variants[target]

        contig_name = ">" + "_".join([chrom, str(target_start), str(target_end)])

        normal_ref_file.write(contig_name + "\n")
        normal_seq = get_ref(original_genome, chrom, target_start, target_end)
        normal_ref_file.write(normal_seq + "\n")

        somatic_ref_file.write(contig_name + "\n")
        somatic_seq = edit_reference(target_start, normal_seq, variant.var_type, variant.start, variant.ref, variant.alt)
        somatic_ref_file.write(somatic_seq + "\n")
    original_genome.close()
    normal_ref_file.close()
    return normal_ref_path, somatic_ref_path


def save_truth_set(variants, wd):
    """
    Write the variants to a truth set file
    :param variants:
    :return: Path to the bed file
    """
    wd = create_dir(os.path.join(wd, 'truth_set'))
    targets_path = os.path.join(wd, 'targets.bed')
    bed_path = os.path.join(wd, 'truth_set_variants.txt')
    vcf_path = os.path.join(wd, 'truth_set_variants.before.vcf')
    with open(bed_path, 'w') as bed_file:
        with open(vcf_path, 'w') as vcf_file:
            with open(targets_path, 'w') as targets_file:
                vcf_file.write('##fileformat=VCFv4.1' + "\n")
                vcf_file.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	truth_set' + "\n")
                for target in sorted(variants):
                    variant = variants[target]
                    variant_end = variant.end if variant.end != variant.start else variant.start + 1  # 1 bp for snp

                    bed_file.write(strip_chrom(variant.chrom) + "\t" + "\t".join(
                        [str(variant.start),
                         str(variant_end),
                         variant.ref,
                         str(variant.alt),
                         variant.var_type]) + "\n")

                    targets_file.write("\t".join([strip_chrom(variant.chrom),
                                                  str(variant.start - 100),
                                                  str(variant.start + 100),
                                                  variant.ref, variant.alt]) + "\n")

                    vcf_file.write("\t".join([strip_chrom(variant.chrom), str(variant.start), '.', variant.ref,
                                              variant.alt, '.', '.', 'empty', 'DP', '200']) + "\n")
    return bed_path, vcf_path


def convert_sam_fastq(sam_path):
    """
    Convert sam file to fastq R1 and R2 (used when we need the *_errFree.sam files)
    :param sam_path:
    :return:
    """
    fq1_path = sam_path[:-4] + '1.fq'
    fq2_path = sam_path[:-4] + '2.fq'
    bam_path = sam_path[:-3] + 'bam'
    cmd = '{0} view -Sb {1} > {2}'.format(samtools.full_path, sam_path, bam_path)  # sam 2 bam
    sp.check_call(cmd, shell=True)
    cmd = '{0} bamtofastq -i {1} -fq {2} -fq2 {3}'.format(bedtools.full_path, bam_path, fq1_path, fq2_path)  # bam 2 fq
    sp.check_call(cmd, shell=True)
    return fq1_path, fq2_path


def run_CuReSim(name, ref_path, coverage, read_length, output_dir, errors):
    #/Users/saeed/Downloads/CuReSim1.3/CuReSim.jar
    pass


def run_ArtificialFastqGenerator(name, ref_path, coverage, read_length, output_dir, errors):
    """
    Generate NGS reads (pair-end, amplicon), compress with gzip and return paths to R1 and R2
    :param ref_path: The path to the reference file (could be partial reference)
    :param name: The base line name for output files. (e.g. normal , somatic)
    :param coverage: The coverage folds
    :param read_length: The length of reads.
    :param output_dir: The output directory
    :param errors: Add errors to NGS reads
    :return: paths to fq1 and fq2
    """
    fastq_dir = create_dir(os.path.join(output_dir, 'fastq', name)) + '/{}_'.format(name)
    test_fastq1 = "/Users/saeed/Downloads/ArtificialFastqGenerator/test1.fastq"
    test_fastq2 = "/Users/saeed/Downloads/ArtificialFastqGenerator/test2.fastq"
    if errors:
        cmd = 'java -jar {0} -R {1} -RL {2} -O {3} -CMP {4} -S ">" -E ">" -URQS false -SE true'\
            .format(ArtificialFastqGenerator_path, ref_path, read_length, fastq_dir, coverage)
        print(cmd)
        sp.check_call(cmd, shell=True)
    else:
        cmd = 'java -jar {0} -R {1} -RL {2} -O {3} -CMP {4} -F1 {5} -F2 {6} -S ">1" -E ">" -URQS false -SE true'\
            .format(ArtificialFastqGenerator_path, ref_path, read_length, fastq_dir, coverage, test_fastq1, test_fastq2)
        print(cmd)
        sp.check_call(cmd, shell=True)
    fastq1_path = bgzip(fastq_dir + ".1.fastq")
    fastq2_path = bgzip(fastq_dir + ".2.fastq")

    return fastq1_path, fastq2_path


def run_art(name, ref_path, coverage, read_length, output_dir, errors):
    """
    Generate NGS reads (pair-end, amplicon), compress with gzip and return paths to R1 and R2
    :param ref_path: The path to the reference file (could be partial reference)
    :param name: The base line name for output files. (e.g. normal , somatic)
    :param coverage: The coverage folds
    :param read_length: The length of reads.
    :param output_dir: The output directory
    :param errors: Add errors to NGS reads
    :return: paths to fq1 and fq2
    """
    # The following line needs a '/' at the end of the path for ART otherwise it writes the files to the parent dir
    fastq_dir = create_dir(os.path.join(output_dir, 'fastq', name)) + '/{}_'.format(name)  # --rndSeed 0
    cmd = '%s --errfree --noALN --paired --in %s --len %d --fcov %d --mflen 200 --sdev 0 --out %s -d %s' \
              % (art.full_path, ref_path, read_length, coverage, fastq_dir, name)
    print(cmd)
    sp.check_call(cmd, shell=True)
    # if _platform in ['linux']:
    #     bsub_cmd = 'bsub -q medium -oo o.log "{}"'.format(cmd)
    #     sp.check_call(bsub_cmd, shell=True)
    if errors:
        fastq1_path = bgzip(fastq_dir + "1.fq")
        fastq2_path = bgzip(fastq_dir + "2.fq")
    else:
        sam_path = fastq_dir + "_errFree.sam"
        cleaned_sam_path = check_cigar(sam_path)
        fastq1_path, fastq2_path = convert_sam_fastq(cleaned_sam_path)
        fastq1_path = bgzip(fastq1_path)
        fastq2_path = bgzip(fastq2_path)
    return fastq1_path, fastq2_path


def merge_fastq_files(normal_fq1, normal_fq2, somatic_fq1, somatic_fq2):
    for n_fq_path, s_fq_path in [(normal_fq1, somatic_fq1), (normal_fq2, somatic_fq2)]:
        with gzip.open(n_fq_path, 'r') as n_fq:
            with gzip.open(s_fq_path, 'a') as s_fq:
                for line in n_fq:
                    s_fq.write(line)
    return somatic_fq1, somatic_fq2


def match_targets_with_normalized_variants(targets, variants):
    """
    After loading variants from normalized variants, we need to create a dictionary of all targets regions with
    their new variants as values
    :param targets: dictionary of BED target regions
    :param variants: dictionary of all variants from normalized VCF file
    :return: dictionary where keys are (chrom, start,end) of targets and values are matched variants
    """
    holder = {}
    for variant in variants:
        for target in targets:
            t_chrom, t_start, t_end = target
            v_chrom, v_start, v_end = variant
            if t_chrom == v_chrom and t_start <= v_start <= t_end:
                holder[(t_chrom, t_start, t_end)] = variants[variant]
                break
    return holder


def normalize_variants(targets, variant_path):
    """
    Run 'vt normalize' tool to normalize the synthetic variants
    :param variant_path: BED file of synthetic variants
    :return: path to the output vcf file
    """
    output_vcf_path = variant_path[:-10] + 'normalized.vcf'
    cmd = '{0} normalize {1} -r {2} > {3}'.format(vt.full_path, variant_path, hg19.full_path, output_vcf_path)
    sp.check_call(cmd, shell=True)
    variants = {}
    Variant = namedtuple('Variant', 'chrom start end ref alt var_type var_size')
    for record in vcf_parser(output_vcf_path):
        key, chrom, start, end, ref, alt, var_type, var_size = record
        chrom = strip_chrom(chrom)
        key = (chrom, start, end)
        variants[key] = Variant(chrom, start, end, ref, alt.sequence, var_type, var_size)
    variants = match_targets_with_normalized_variants(targets, variants)
    return output_vcf_path, variants


def check_cigar(sam_path):
    """
    ART simulation has a bug when generating *_errFree.sam files where CIGAR is not equal the actual sequence which is
    also different from the fq files. This function ignore those records.
    :param sam_path:
    :return: Path to new cleaned file
    """
    output_path = sam_path[:-4] + '_cleaned.sam'
    output_file = open(output_path, 'w')
    reads_with_issues = []
    with open(sam_path, 'r') as samfile:
        for line in samfile:
            if line.startswith('@'):
                continue
            line = line.split('\t')
            seq = line[9]
            cigar = int(line[5].replace("=", ''))
            if len(seq) != cigar:
                reads_with_issues.append(line[0])
                print("Error in read length and cigar. Ignoring this read:", line)
        samfile.seek(0)  # restart
        for line in samfile:
            if line.startswith('@'):
                output_file.write(line)
                continue
            line = line.split('\t')
            if line[0] not in reads_with_issues:
                output_file.write('\t'.join(line))
    output_file.close()
    # os.remove(sam_path)
    return output_path


def parse_user_variant_list(variant_list_path, flanking):
    """
    Parse the user variant file (chrom, start, end, ref, alt)
    :param variant_list_path:
    :param flanking: Integer
    :return: Dictionary of variants where keys are tuples: (chrom, start_target, target_end)
    """
    variants = {}
    Variant = namedtuple('Variant', 'chrom start end ref alt var_type var_size')
    with open(variant_list_path, 'rU') as variants_file:
        reader = csv.reader(variants_file, delimiter='\t')
        for row in reader:
            var_type, var_size = get_var_type_size(row[3], row[4])
            target_start = int(row[1]) - flanking
            target_end = int(row[2]) + flanking
            variants[(row[0], target_start, target_end)] = Variant(chrom=row[0], start=int(row[1]), end=int(row[2]),
                                                                   ref=row[3], alt=row[4],
                                                                   var_type=var_type, var_size=var_size)
    return variants


@click.group(invoke_without_command=True)
@click.option('--bed-path', '-b', help='Path to target regions in BED format')
@click.option('--output-dir', '-o', help='Path to output directory')
@click.option('--read-length', '-l', default=150, help='Read length (default=150)')
@click.option('--coverage', '-c', default=200, help='Coverage (default=200)')
@click.option('--errors/--no-errors', default=True, help='Adding errors to NGS reads')
@click.option('--flanking', '-f', default=200, help='Flanking size (default=100)')
@click.option('--indel-size', '-i', default="1,20", help='Range of indel size as "min,max" (default=1,15)')
@click.option('--fractions', '-r', default="0.2", #0.005,0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
              help='Fraction of alternative read')
@click.option('--variant-types', '-v', default="SNV,DEL,INS", help='SNV,INS or DEL (default=SNV)')
@click.option('--variant-list', '-n', help='A BED file of variant to be simulated (chrom, start, end, ref, alt)')
def cli(bed_path=None, variant_path=None, output_dir=None, read_length=None, coverage=None, errors=None,
        flanking=None, fractions=None, indel_size=None, variant_types=None, variant_list=None):
    fractions = fractions.split(",")
    if variant_list is None:
        # if the user supply the target regions but not the variants
        min_indel, max_indel = indel_size.split(",")
        indel_sizes = range(int(min_indel), int(max_indel))
        var_types = variant_types.split(",")
        variants = get_synthetic_variants(read_length, bed_path, hg19.full_path, flanking, var_types, indel_sizes)
    else:
        # if the user provided a list of variants, parse them then load them into a dictionary
        variants = parse_user_variant_list(variant_list, flanking)

    fastq_pairs = {}

    variant_bed_path, variant_vcf_path = save_truth_set(variants, output_dir)
    variant_vcf_path, variants = normalize_variants(variants, variant_vcf_path)

    # Get reference sequence based on the input target regions defined in BED file.
    normal_ref_path, somatic_ref_path = generate_normal_somatic_references(variants, output_dir)

    # write new normal ref fasta files

    normal_fq1_path, normal_fq2_path = run_art('normal', normal_ref_path, coverage, read_length,
                                               output_dir, errors)

    # Generate fq1/fq2 for various fraction combination
    for somatic_fraction in fractions:
        somatic_fraction = float(somatic_fraction)
        somatic_coverage = float(coverage) * float(somatic_fraction)
        normal_coverage = float(coverage) * (1. - somatic_fraction) if somatic_fraction != 1. else 1  # at least 1 read

        normal_fq1_path, normal_fq2_path = run_art('normal_%s' % str(somatic_fraction), normal_ref_path,
                                                   normal_coverage, read_length, output_dir, errors)

        somatic_fq1_path, somatic_fq2_path = run_art('somatic_%s' % str(somatic_fraction),
                                                     somatic_ref_path, somatic_coverage, read_length,
                                                     output_dir, errors)

        somatic_fq1_path, somatic_fq2_path = merge_fastq_files(normal_fq1_path, normal_fq2_path,
                                                               somatic_fq1_path, somatic_fq2_path)

        fastq_pairs[str(somatic_fraction)] = {'somatic': {'fastq': [somatic_fq1_path, somatic_fq2_path]},
                                              'normal': {'fastq': [normal_fq1_path, normal_fq2_path]}
                                              }

    snapshot_manifest_path = create_manifest(os.path.join(output_dir, 'manifests'), pipeline_name='snapshot',
                                             paths=fastq_pairs)

    onek_manifest_path = create_manifest(os.path.join(output_dir, 'manifests'), pipeline_name='onek',
                                         paths=fastq_pairs)

    print(snapshot_manifest_path)
    print(onek_manifest_path)

    # os.system('source ~/.virtualenvs/cider2/bin/activate && cider snapshot -n snvs_test -m {} -f '
    #           '-b /Users/saeed/Dropbox/Git/LMM_CID/CID/AmpSim/ampsim/tests/data/200.bed --aligner bwa'
    #           .format(snapshot_manifest_path))

    # os.system('source ~/.virtualenvs/cider2/bin/activate && cider snapshot -n snvs_freeBayes -m {} -f '
    #           '-b /PHShome/sha13/projects/cid/300.bed'
    #           .format(snapshot_manifest_path))

if __name__ == '__main__':
    cli()

