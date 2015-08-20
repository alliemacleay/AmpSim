"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

from collections import namedtuple
from random import choice, randint, seed
from pybedtools import BedTool
import pysam


indel_sizes = [1, 2, 3]
indel_sizes += range(4, 31, 3)  # 51
dna = ['A', 'C', 'T', 'G']
fractions = [0.2]  # [0.11, 0.15]  # 0.005, 0.01, 0.03, 0.05, 0.07, 0.09,
var_types = ['SNV', 'DEL', 'INS']

seed(102)


def load_bed(path):
    """
    parse bed file
    :param path:
    :return: dictioanary
    """
    targets = []
    bed = BedTool(path)
    for record in bed:
        targets.append(record[:3])
    return targets


def generate_insertion_seq(size):
    """
    Generate DNA insertion with size=size
    :param size:
    :return:
    """
    holder = ''
    while len(holder) <= size:
        holder += choice(dna)
    return holder


def get_ref(genome, chrom, start, end):
    """
    Fetch genome reference given chrom, start, end
    :param genome:
    :param chrom:
    :param start:
    :param end:
    :return: DNA sequence
    """
    return genome.fetch(chrom, int(start), int(end))


def get_alt(ref, var_size, var_type):
    """
    Generate alternative allele from SNV and INS but not DEL.
    :param var_type:  SNV or INS
    :return: dna string
    """
    if var_type == 'SNV':
        nucleotides = [x for x in dna if x != ref]
        return choice(nucleotides)
    elif var_type == 'INS':
        return generate_insertion_seq(var_size)


def get_bamsurgeon_output(variants):
    snvs_output = ''
    indels_output = ''
    merged_output = '\t'.join(['#chrom', 'start', 'end', 'ref', 'alt', 'fraction', 'var_type']) + "\n"
    for variant in variants:
        merged_output += '\t'.join([str(x) for x in variant]) + "\n"
        if variant.var_type == 'SNV':
            snvs_output += '\t'.join([str(x) for x in[variant.chrom, variant.start, variant.end,
                                      variant.fraction, variant.alt]]) + "\n"
        elif variant.var_type == 'DEL':
            indels_output += '\t'.join([str(x) for x in[variant.chrom, variant.start, variant.end,
                                        variant.fraction, variant.var_type, '']]) + "\n"
        elif variant.var_type == 'INS':
            indels_output += '\t'.join([str(x) for x in [variant.chrom, variant.start, variant.end,
                                        variant.fraction, variant.var_type, variant.alt]]) + "\n"

    return snvs_output, indels_output, merged_output


def generate_variants(genome, targets):
    """
    Add random variant SNV, DEL or INS with random fraction (one per target region)
    :param targets: list of regions (chrom, start,end
    :return: dictionary of variants
    """
    variants = []
    Variant = namedtuple('Variant', 'chrom start end ref alt fraction var_type')
    for target in targets:
        var_type = choice(var_types)
        var_size = choice(indel_sizes)  # used from indels
        fraction = choice(fractions)
        chrom, start, end = target
        # padding with var_size to ensure we don't create a variant that extend outside the target
        position = abs(int(start) - int(end)) / 2 + int(start)  # randint(int(start) + var_size, int(end) - var_size)
        if var_type == 'SNV':
            ref = get_ref(genome, chrom, position, position + 1)
            alt = get_alt(ref, var_size, var_type)
            start = end = position + 1
        elif var_type == 'INS':
            ref = get_ref(genome, chrom, position, position + 1)
            alt = ref + get_alt(ref, var_size, var_type)
            start = position + 1
            end = position + var_size
        elif var_type == 'DEL':
            ref = get_ref(genome, chrom, position, position + var_size)
            alt = ref[0]  # the first nucleotides
            start = position + 1
            end = position + var_size
        variants.append(Variant(chrom, start, end, ref, alt, fraction, var_type))
    return variants


def get_synthetic_variants(targets_bed, genome_ref):
    """
    The is the main method. Generate synthetic snv and indel based one per target region, where the fraction
    and variant type SNV,INS or DEL is random values from fractions, indel_sizes and fractions defined at the top
    of this module.
    :param targets_bed:
    :param genome_ref:
    :return: String for Bam_surgeon format (SNVS and INDELs and merged as our truth set)
    """
    genome = pysam.FastaFile(genome_ref)
    targets = load_bed(targets_bed)
    variants = generate_variants(genome, targets)
    snvs, indels, merged = get_bamsurgeon_output(variants)
    return snvs, indels, merged

