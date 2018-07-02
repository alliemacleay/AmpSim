"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import subprocess as sp
from collections import namedtuple

import logging
import click
import pysam

from ampsim.utils import vcf_parser, bgzip, strip_chrom
from ampsim.tools.generate_ngs import save_truth_set, normalize_variants, create_manifest
from ampsim import download

logger = logging.getLogger(__name__)

padding_size = 1000  # the minimum length to separate two modified variants

samtools = download.Samtools()
bedtools = download.Bedtools()
picard = download.Picard()
bamutil = download.Bamutil()


def load_truth(vcf_path, bam_path):
    """
    Return list of variants to be halved in the BAM file. Only SNV will be considered and they must be separated by
    at least padding_size
    :param vcf_path: Path
    :return: Dictionary  of chrom:pos:ref,alt
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    Variant = namedtuple('Variant', 'chrom start end ref alt')
    variants = {}
    last_position = 0  # the position of last loaded variant in the dictionary. Keep 1000 bp between them
    for record in vcf_parser(vcf_path):
        key, chrom, start, end, ref, alt, var_type, var_size = record
        # if start != 30257560:
        #     continue
        start = int(start)
        if not count_read_coverage(bam, chrom, start, start)[0] >= 40:
            continue  # skip variants that are not covered by at least 20 reads

        if var_type != 'SNV' or abs(start - last_position) < padding_size:
            continue
        variants[(chrom, start, end)] = Variant(chrom, start, end, str(ref), str(alt))
        last_position = start
    return variants


def sort_index(bam_path):
    print('Sorting and index ..')
    new_sorted_path = bam_path[:-4] + '_sorted.bam'
    sp.call(['samtools', 'sort', bam_path, new_sorted_path[:-4]])
    sp.call(['samtools', 'index', new_sorted_path])
    return new_sorted_path


def bam_2_fastq(bam_path, alt_count):
    """
    Convert bam file to fastq R1 and R2 (used when we need the *_errFree.sam files)
    :param bam_path: Path to BAM file
    :return:
    """
    fq1_path = bam_path[:-4] + '_{}_1.fastq'.format(alt_count)
    fq2_path = bam_path[:-4] + '_{}_2.fastq'.format(alt_count)
    unmapped_path = bam_path[:-4] + '_{}_unmapped.fastq'.format(alt_count)
    sorted_bam = bam_path[:-4] + '_{}_sorted.bam'.format(alt_count)

    cmd = '{0} sort {1} {2}'.format(samtools.full_path, bam_path, sorted_bam[:-4])

    sp.check_call(cmd, shell=True)

    # cmd = '{0} bam2FastQ --in {1} --firstOut {2} --secondOut {3} --unpairedOut {4}'\
    #     .format(bamutil.full_path, sorted_bam, fq1_path, fq2_path, unmapped_path)
    cmd = 'java -jar {0} SamToFastq INPUT={1} FASTQ={2} SECOND_END_FASTQ={3} UNPAIRED_FASTQ={4}'\
        .format(picard.full_path, bam_path, fq1_path, fq2_path, unmapped_path)
    # cmd = '{0} bamtofastq -i {1} -fq {2} -fq2 {3}'
    # .format(bedtools.full_path, bam_path, fq1_path, fq2_path)  # bam 2 fq
    print(cmd)
    sp.check_call(cmd, shell=True)
    return bgzip(fq1_path), bgzip(fq2_path)


def read_is_high_quality(var_start, read, mate_read):
    """
    Check if the read passes few QC checks like MQ > mpq, no duplicate, mapped
    :param var_start: The variant position (integer)
    :param read: A pysam read object
    :return: Boolean
    """
    base_idx = abs((read.reference_start + 1) - var_start)
    if base_idx > read.rlen - 1:
        return False  # the variant is not in this read (the read is at its boundary)
    if read.is_duplicate:
        return False
    if read.is_qcfail:
        return False
    if read.is_unmapped:
        return False
    if read.mapq < 20:
        return False
    if read.mate_is_unmapped:  # when the mate is unmapped it will break the mate(read) method below. Ignore these.
        return False
    if not read.is_paired:   # when the mate is missing it will break the mate(read) method below. Ignore these.
        return False
    if var_start in read.positions and var_start in mate_read.positions:  # pair-read overlaps and inflating counts
        return False
    return True


def edit_reads(reads, start, ref, alt, alt_count):
    Ns = {}
    new_reads = []
    total = len(reads)

    if total % 2 == 0:  # if the number of reads  is even, get the half
        alt_counts = total / 2
    else:
        reads.pop()  # remove last element
        alt_counts = len(reads) / 2

    for counter in range(0, len(reads)):
        read, mate_read = reads[counter]
        base_idx = abs((read.reference_start + 1) - start)

        if counter < alt_counts + alt_count:  # alt_count can be 0 to get 50/50 or (+5 or -5) to tip the balance off
            allele = str(alt)
        else:
            allele = str(ref)

        if allele not in Ns:
            Ns[allele] = 0
        Ns[allele] += 1
        read.reference_id -= 1  # illumina reference 1 is MT so we need to covert it to our reference
        mate_read.reference_id -= 1
        new_seq = list(read.seq)
        new_seq[base_idx] = allele
        new_seq = ''.join(new_seq)
        new_qual = read.query_qualities  # Very important since assigning to seq will invalidate any quality scores.
        read.seq = new_seq
        read.query_qualities = new_qual
        new_reads.append((read, mate_read))

    print(start, ref, alt, Ns, total, alt_counts)
    return new_reads


def edit_bam(bam_path, variants, alt_count, output_dir):
    """
    Go through variants and edit the BAM file to make the alternative allele ratio 50/50 (or add / minus 5 reads away
    from 50/50 ratio).
    :param bam_path:
    :param variants:
    :param alt_count: Number of reads to add to the 50/50 ratio to tip it off.
    :return:
    """
    print('Editing variant reads to make it 50/50 allele ratio ..')
    output_path = os.path.join(output_dir, bam_path.split("/")[-1][:-4] + '_{}_edited.bam'.format(alt_count))
    template_bam = pysam.AlignmentFile('/Users/saeed/BaseSpace/NA12878/chr21/0.5.cid_0.A0_P0.tumor.all_chroms.marked_dup.bam', 'rb')
    in_bam = pysam.AlignmentFile(bam_path, 'rb')
    cid_header = template_bam.header
    illumina_header = in_bam.header
    cid_header['RG'] = illumina_header['RG']
    del cid_header['PG']
    out_bam = pysam.AlignmentFile(output_path, 'wb', header=cid_header)

    for key in variants:
        chrom, start, end, ref, alt = variants[key]
        reads = []
        processed_reads_names = {}

        for read in in_bam.fetch(region='%s:%d-%d' % (chrom, start, start)):
            # if the read has no mate , skip it since we can't add it later in the FASTQ (bedtools will complain) and
            # this will so this will disturb the 50/50 without us knowing
            try:
                mate_read = in_bam.mate(read)
                if mate_read.qname in processed_reads_names:
                    continue
            except:
                # No mate is found. Move on to the next read.
                continue
            # Remove reads that are likely to be ignored by the caller so we need to keep the 50/50 ratio holds
            if read_is_high_quality(start, read, mate_read):
                if read.qname not in processed_reads_names:
                    reads.append((read, mate_read))
                    processed_reads_names[read.qname] = ''
                else:
                    continue  # either the read or its mate has been processed for this variant
            else:
                continue
        edited_reads = edit_reads(reads, start, ref, alt, alt_count)
        for read, mate_read in edited_reads:
            out_bam.write(read)
            out_bam.write(mate_read)

    in_bam.close()
    out_bam.close()
    output_path = sort_index(output_path)
    return output_path


def count_read_coverage(bam, chrom, start, end):
    """ calculate coverage of aligned reads over region
    """

    coverage = []
    start = int(start)
    end = int(end)
    for i in range(end-start+1):
        coverage.append(0.0)

    i = 0
    if chrom in bam.references:
        for pcol in bam.pileup(region="{0}:{1}-{2}".format(chrom, start, end)):
            n = 0
            if pcol.pos >= start and pcol.pos <= end:
                for read in pcol.pileups:
                    if read.alignment.mapq >= 0 and not read.alignment.is_duplicate:
                        n += 1
                coverage[i] = n
                i += 1

    return coverage


@click.group(invoke_without_command=True)
@click.option('--bam-path', '-b', type=click.Path(exists=True), help='Path to sorted and indexed BAM file')
@click.option('--vcf-path', '-v', type=click.Path(exists=True), help='Path to variants in a VCF format')
@click.option('--output-dir', '-o', help='Path save the output files')
def cli(bam_path=None, vcf_path=None, output_dir=None):
    """
    For SNV in truth set (VCF), alter the BAM file to make the REF:ALT ration 50:50. Write new BAM file and List of
    edited variants.
    :param bam_path: Path to a bam file
    :param vcf_path: Path to a truth set in a VCF file
    :return: None
    """
    variants = load_truth(vcf_path, bam_path)
    variant_bed_path, variant_vcf_path = save_truth_set(variants, output_dir)

    fastq_holder = {}

    alt_count_dict = {'0.50': 0, "0.51": 2, "0.49": -2}

    for alt_count_name, alt_count_value in alt_count_dict.items():
        edited_bam_path = edit_bam(bam_path, variants, alt_count_value, output_dir)
        fq1_path, fq2_path = bam_2_fastq(edited_bam_path, alt_count_name)
        fastq_holder[alt_count_name] = {'somatic':
                                            {'fastq': [fq1_path, fq2_path],
                                             'bam': edited_bam_path}
                                        }

    snapshot_manifest_path = create_manifest(os.path.join(output_dir, 'manifests'),
                                             pipeline_name='snapshot',
                                             paths=fastq_holder)

    # variant won't match because of the chr prefix so leave this step at the end.

    variant_vcf_path, variants = normalize_variants(variants, variant_vcf_path)

    print(snapshot_manifest_path)

if __name__ == '__main__':
    cli()
