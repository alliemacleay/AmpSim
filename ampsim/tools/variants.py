"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

from collections import namedtuple
from random import choice, seed
from pybedtools import BedTool
import pysam
from bx.intervals.intersection import Intersecter, Interval  # pylint: disable=E0611


dna = ['A', 'C', 'T', 'G']


seed(102)


def load_bed(read_length, path, flanking):
    """
    parse bed file and ignore overlapping regions
    :param path: BED file path
    :param flanking: Integer to add to each target's start and end
    :return: dictionary of targets where keys are chrom,start,end
    """
    targets = {}
    unique_targets = []
    bed = BedTool(path)
    for record in bed:
        chrom, start, end = record[0], int(record[1]), int(record[2])
        start -= flanking
        end += flanking
        if int(end) < int(start):   # if this target on minus strand flip the start/end
            start, end = end, start
        if chrom not in targets:
            targets[chrom] = Intersecter()
        if targets[chrom].find(start, end) == [] and abs(start - end) >= (read_length * 2):  # no overlaps yet
            targets[chrom].add_interval(Interval(start, end))
            unique_targets.append([chrom, start, end])

    return unique_targets


def get_ref(genome, chrom, start, end):
    """
    Fetch genome reference given chrom, start, end
    :param genome: pysam fasta file
    :param chrom:
    :param start:
    :param end:
    :return: DNA sequence
    """
    return genome.fetch(chrom, int(start) - 1, int(end) - 1)


def get_alt(ref, var_size, var_type):
    """
    Generate alternative allele from SNV and INS but not DEL.
    :param var_type:  SNV or INS
    :return: A string of nucleotides
    """
    if var_type == 'SNV':
        nucleotides = [x for x in dna if x != ref]
        return choice(nucleotides)
    elif var_type == 'INS':
        nucleotides = ''
        while len(nucleotides) < var_size:
            nucleotides += choice(dna)
        return nucleotides


def generate_variant(genome, target, var_types, indel_sizes):
    """
    Add random variant SNV, DEL or INS with random fraction (one per target region, in the middle)
    :param genome: pysam fasta genome object
    :param target: A tuple of chrom, start,end
    :param var_types: SNV, INS, DEL (later will add CNV_DUP, CNV_DEL)
    :param indel_sizes: Integer of specific size other wise for indels otherwise it chose randomly from allowed indel sizes
    :return: Variant class
    """
    Variant = namedtuple('Variant', 'chrom start end ref alt var_type var_size')

    var_type = choice(var_types)
    var_size = choice(indel_sizes)

    chrom, start, end = target
    # padding with var_size to ensure we don't create a variant that extend outside the target
    position = (abs(int(start) - int(end)) / 2) + int(start)  # // is floor division
    if var_type == 'SNV':
        var_size = 1
        ref = get_ref(genome, chrom, position, position + 1)
        alt = get_alt(ref, var_size, var_type)
        start = end = position
    elif var_type == 'INS':
        ref = get_ref(genome, chrom, position, position + 1)
        alt = ref + get_alt(ref, var_size, var_type)
        start = position
        end = position + var_size
    elif var_type == 'DEL':
        ref = get_ref(genome, chrom, position, position + var_size + 1)
        alt = ref[0]  # the first nucleotides
        start = position
        end = position + var_size
    return Variant(chrom, start, end, ref, alt, var_type, var_size)


def get_synthetic_variants(read_length, targets_bed, genome_ref, flanking, var_types, indel_sizes):
    """
    The is the main method. Generate synthetic snv and indel based one per target region, where the fraction
    and variant type SNV,INS or DEL is random values from fractions, indel_sizes and fractions defined at the top
    of this module.
    :param targets_bed:
    :param genome_ref:
    :param flanking: Integer of number of bp to add on each side of the traget region to ensure well coverage and avoid
    small target that covered by the same reads (block of reads with the same start and end).
    :param var_types: SNV,INS,DEL
    :param indel_sizes: List of allowed sizes for indels
    :return: String for Bam_surgeon format (SNVS and INDELs and merged as our truth set)
    """
    genome = pysam.FastaFile(genome_ref)
    targets = load_bed(read_length, targets_bed, flanking)
    variants = {}
    for target in targets:
        variants[tuple(target)] = generate_variant(genome, target, var_types, indel_sizes)
    return variants

