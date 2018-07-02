"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import click
import vcf
import pysam
import matplotlib.pyplot as plt

from ampsim.utils import strip_chrom

plt.style.use('ggplot')


def plot(data):
    plt.figure(figsize=(7, 5), dpi=80)
    plt.xlabel("allele fractions", fontsize=12)
    plt.ylabel("Number of variants", fontsize=12)  # with error profile #  error-free
    data = list(filter(lambda x: 0.45 < x < 0.55, data))
    plt.hist(sorted(data), histtype='step', bins=50)
    plt.tight_layout()
    plt.show()


@click.group(invoke_without_command=True)
@click.option('--vcf-path', '-v', type=click.Path(exists=True), help='Path to variants in a VCF format')
@click.option('--bam-path', '-b', type=click.Path(exists=True), help='Path to variants in a BAM format')
def cli(vcf_path=None, bam_path=None):
    """
    Plot allele fraction
    :param vcf_path: Path VCF file
    :param bam_path: Path BAM file
    :return: None
    """
    data = []
    bam = pysam.AlignmentFile(bam_path, 'rb')
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    for record in vcf_reader:
        if not record.is_snp:
            continue
        data.append(record.samples[0].data.VF)
        if record.samples[0].data.VF == 0.5:
            AD = record.samples[0].data.AD
            counts = bam.count_coverage(record.CHROM, record.POS - 1, record.POS, quality_threshold=0)
            print(record, AD, dict(zip('ACGT', [int(x.tolist()[0]) for x in counts])))
        # break
        # print(record, record.samples[0].data.VF)
        # print(dir(record))
        # break
    # plot(data)

if __name__ == '__main__':
    cli()
