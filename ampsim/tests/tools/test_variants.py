"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

from nose.tools import assert_equal

import random

import pysam

from ampsim.tools import generate_ngs
from ampsim.tools import variants
from ampsim.download import Hg19

random.seed(102)  # to make sure testing cases are the same every time we run the tests

hg19 = Hg19()
genome = pysam.Fastafile(hg19.full_path)


def test_target():
    target = ('1', 101185921, 101186231)
    target_start = target[1]
    target_seq = variants.get_ref(genome, *target)
    variant = variants.generate_variant(genome, target)
    somatic_seq = generate_ngs.edit_reference(target_start, target_seq, 'SNV', variant.start, variant.ref, variant.alt)


def test_snv_synthesis():
    """Unit test synthesizing SNV in a reference genome"""
    chrom = '1'
    position = 11184622
    expected_ref = 'T'
    expected_alt = 'G'
    normal_seq = 'CATCA'
    somatic_seq = 'CAGCA'  # T > G in the middle

    # test fetching a genome
    observed_ref = variants.get_ref(genome, chrom, position, position + 1)
    assert_equal(expected_ref, observed_ref)

    # testing adding to a read
    observed_seq = generate_ngs.edit_reference(position - 2, normal_seq, 'SNV', position, expected_ref, expected_alt)
    assert_equal(somatic_seq, observed_seq)

    # test ref / alt meeting VCF specifications


def test_deletion_synthesis():
    """Unit test synthesizing DEL in a reference genome"""
    chrom = '1'
    position = 11184622
    var_size = 1
    expected_ref = 'TC'
    expected_alt = 'T'
    normal_seq = 'CATCA'
    somatic_seq = 'CATA'  # TCACA > T in the middle

    # test fetching a genome
    observed_ref = variants.get_ref(genome, chrom, position, position + var_size + 1)
    assert_equal(expected_ref, observed_ref)

    # testing adding to a read
    observed_seq = generate_ngs.edit_reference(position - 2, normal_seq, 'DEL', position, expected_ref, expected_alt)
    assert_equal(somatic_seq, observed_seq)

    # test variant synthesis
    target = ('1', 11184622, 11184632)
    variant = variants.generate_variant(genome, target, var_size=2, var_type='DEL')
    assert_equal(variant.ref, 'CGC')
    assert_equal(variant.alt, 'C')

    target = ('1', 101197913, 101198300)
    variant = variants.generate_variant(genome, target, var_size=16, var_type='DEL')
    normal_seq = variants.get_ref(genome, *target)
    observed_seq = generate_ngs.edit_reference(target[1], normal_seq, 'DEL', variant.start, variant.ref, variant.alt)
    # assert_equal(variant.ref, 'CGC')
    # assert_equal(variant.alt, 'C')


def test_insertion_synthesis():
    """Unit test synthesizing INS in a reference genome"""
    chrom = '1'
    position = 11184622
    var_seq = 'A'
    expected_ref = 'T'
    expected_alt = 'T' + var_seq
    normal_seq = 'CATCA'
    somatic_seq = 'CATACA'  # T > TACTT in the middle

    # test fetching a genome
    observed_ref = variants.get_ref(genome, chrom, position, position + 1)
    assert_equal(expected_ref, observed_ref)

    # testing adding to a read
    observed_seq = generate_ngs.edit_reference(position - 2, normal_seq, 'INS', position, expected_ref, expected_alt)
    assert_equal(somatic_seq, observed_seq)

    # test variant synthesis
    target = ('1', 11184622, 11184632)
    variant = variants.generate_variant(genome, target, var_size=2, var_type='INS')
    assert_equal(variant.ref, 'C')
    assert_equal(variant.alt, 'CAT')




