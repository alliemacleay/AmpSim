"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import re
import os
from collections import namedtuple
import click
import pandas as pd
from pybedtools import BedTool
import matplotlib.pyplot as plt

from ampsim.utils import safe_divide, vcf_parser, get_vcfs

pd.set_option('display.expand_frame_repr', False)  # to print the dataframe without breaking the row into multiple lines
plt.style.use('ggplot')


def calculate_target_size(path):
    """
    Given a bed for target regions, calculate the total size of the regions (bp)
    :param path: To  a bed file
    :return: Integer size
    """
    size = 0
    bed = BedTool(path)
    for target in bed:
        size += abs(target.start - target.end)
    return size


def load_variants(truth_path, vcf_paths, allowed_callers):
    """
    Load variant from the Truth set and all VCF files.
    :param truth_path: BED file of the truth variant sets
    :param vcf_paths: List of paths to the VCF files
    :return: dataframe of variants and a list of caller column names
    """
    columns = ['key', 'chrom', 'start', 'end', 'ref', 'alt', 'var_type', 'var_size', 'is_truth', 'fraction']
    callers = {caller: False for caller in [os.path.basename(path).split(".")[-3].split('_')[0]
                                            for path in vcf_paths] if caller in allowed_callers}
    columns.extend(callers)
    Variant = namedtuple('Variant', columns)
    truth = {}

    for vcf_path in vcf_paths:
        caller = os.path.basename(vcf_path).split(".")[-3].split('_')[0]
        if caller not in allowed_callers:
            continue

        fraction = re.findall("\d+\.\d+", os.path.basename(vcf_path))[0]
        # fraction = os.path.basename(vcf_path).split("_")[0][:-4]
        # load variant from the truth set
        for record in vcf_parser(truth_path):
            key, chrom, start, end, ref, alt, var_type, var_size = record
            key = "_".join([chrom, str(start), ref, alt.sequence])
            if fraction not in truth:
                truth[fraction] = {}
            if key not in truth[fraction]:
                truth[fraction][key] = Variant(key, chrom, start, end, ref, alt, var_type, int(var_size),
                                               is_truth=True, fraction=fraction, **callers)._asdict()

        for record in vcf_parser(vcf_path):
            key, chrom, start, end, ref, alt, var_type, var_size = record
            if key in truth[fraction]:
                truth[fraction][key][caller] = True
            else:
                truth[fraction][key] = Variant(key, chrom, start, end, ref, alt, var_type, abs(var_size),
                                               is_truth=False, fraction=fraction, **callers)._asdict()
                truth[fraction][key][caller] = True

    holder = [truth[fraction][key] for fraction in truth for key in truth[fraction]]
    df = pd.DataFrame.from_records(data=holder, columns=columns)
    return df, callers


def calculate_true_negative(target_size, var_type, caller, df):
    """
    The true negative are all nucleotide bases that don't overlap with truth set variants or false positive variants
    :param target_size: Integer
    :param caller: String caller name
    :param df: dataframe
    :return: Integer
    """
    expected_variant_size = sum(df.loc[(df['var_type'].isin(var_type)) &
                                       (df['is_truth'] == True)]['var_size'].tolist())
    observed_fp_size = sum(df.loc[(df['var_type'].isin(var_type)) &
                                  (df['is_truth'] == False) &
                                  (df[caller] == True)]['var_size'].tolist())
    return target_size - expected_variant_size - observed_fp_size


def caller_performance(df, target_size, callers, fractions=None, var_type=['SNV', 'DEL', 'INS', 'CNV_DEL', 'CNV_DUP']):
    """
    Return the FP, FN, TP, TN, PPV, NPV, specificity, sensitivity, recall, precision
    :param df: The variant including truth set
    :param var_type: list of variant types snv, ins, del, cnv_del, cnv_dup
    :param callers: list of caller column names in the dataframe
    :param fractions: list of fractions to filter in
    :return: tp, tn, fp, fn, total, sen, spe, per, acc
    """
    columns = ['caller', 'fraction', 'tp', 'tn', 'fp', 'fn', 'total', 'sensitivity', 'specificity',
               'ppv', 'npv', 'accuracy', 'precision', 'recall', 'f1']
    performance_df = pd.DataFrame(columns=columns)

    if fractions:
        df = df.loc[df['fraction'].isin(fractions)]
        fraction_name = ",".join([str(x) for x in fractions])
    else:
        fraction_name = 'all'

    for caller in callers:
        if caller not in df.columns:
            continue
        record = {k: 0. for k in columns}
        record['caller'] = caller
        record['fraction'] = fraction_name
        record['tp'] = float(df.loc[(df['var_type'].isin(var_type)) & (df['is_truth'] == True) &
                                    (df[caller] == True)].shape[0])
        record['tn'] = calculate_true_negative(target_size, var_type, caller, df)
        record['fp'] = float(df.loc[(df['var_type'].isin(var_type)) & (df['is_truth'] == False) &
                                    (df[caller] == True)].shape[0])
        record['fn'] = float(df.loc[(df['var_type'].isin(var_type)) & (df['is_truth'] == True) &
                                    (df[caller] == False)].shape[0])
        record['total'] = record['tp'] + record['tn'] + record['fp'] + record['fn']
        record['sensitivity'] = safe_divide(record['tp'], (record['tp'] + record['fn']))
        record['specificity'] = safe_divide(record['tn'], (record['tn'] + record['fp']))
        record['ppv'] = safe_divide(record['tp'], (record['tp'] + record['fp']))
        record['npv'] = safe_divide(record['tn'], (record['tn'] + record['fn']))
        record['accuracy'] = safe_divide(record['tp'] + record['tn'], record['total'])
        record['precision'] = safe_divide(record['tp'], record['tp'] + record['fp'])
        record['recall'] = safe_divide(record['tp'], record['tp'] + record['fn'])
        record['f1'] = safe_divide(2 * (record['precision'] * record['recall']),
                                   (record['precision'] + record['recall']))
        performance_df.loc[caller] = pd.Series(record)

    return performance_df


def plot_performance_per_bin(df, callers, var_types, targets_size, output_png):
    """
    Plot the sensitivity and specificity for each caller at every fraction bin and save it to the output_file[:-3] + png
    :param df: Dataframe of all variants
    :param callers: list of callers (each will have its own line)
    :param var_types:  list of variants to include
    :param targets_size: Total size of the targets
    :param output_png: Path to save the plot
    :return: None
    """
    results = {}
    fractions = sorted(list(set(([float(x) for x in df['fraction'].tolist()]))))
    for caller in callers:
        results[caller] = {'sensitivity': {}, 'specificity': {}}
        for fraction in fractions:
            caller_performance_df = caller_performance(df, targets_size, [caller], [str(fraction)], var_types)
            results[caller]['sensitivity'][fraction] = caller_performance_df.loc[caller]['sensitivity']
            results[caller]['specificity'][fraction] = caller_performance_df.loc[caller]['specificity']
    results_df = pd.DataFrame.from_dict(data=results, orient='index')

    # plotting
    plt.figure(figsize=(7, 5), dpi=80)
    plt.title("Callers sensitivity {0} variants ({1})".format(output_png.split("/")[-1][:-4], ",".join(var_types)),
              fontsize=14)
    plt.xlabel("ALT fractions", fontsize=12)
    plt.ylabel("Sensitivity", fontsize=12)  # with error profile #  error-free

    for caller in callers:
        fractions = [fraction for fraction in fractions if fraction <= 1.0]
        y_values = [results_df.loc[caller]['sensitivity'][fraction] for fraction in fractions]
        if caller == 'ensemble':
            linestyle = '--'
        else:
            linestyle = '-'
        plt.plot(fractions, y_values, linestyle=linestyle, linewidth=2, label=caller)
    plt.legend(fontsize=11, loc='best')
    plt.tight_layout()
    plt.savefig(output_png)
    plt.show()


@click.group(invoke_without_command=True)
@click.option('--truth-path', '-t', help='Path to the variant truth set')
@click.option('--target-bed', '-b', help='Path to teget bed')
@click.option('--cider-path', '-c', help='Path to cider output to get all vcf files at once')
@click.option('--vcf-paths', '-v', multiple=True, help='Path to the VCF paths (multiple files)')
@click.option('--var-types', '-r', default='SNV,INS,DEL', help='List SNV,DEL or INS (comma separated). Default is SNV.')
@click.option('--output-png', '-o', help='Path to save the results')
def cli(truth_path=None, target_bed=None, cider_path=None, vcf_paths=None, var_types=None, output_png=None):
    var_types = var_types.split(",")
    allowed_callers = ['lofreq', 'gatk', 'mutect', 'mutect2', 'ensembler']#, 'hotspotter']
    targets_size = calculate_target_size(target_bed)
    vcf_paths = get_vcfs(cider_path)  # TODO remove this in production

    df, callers = load_variants(truth_path, vcf_paths, allowed_callers)
    print(df.to_csv(sep='\t'))

    fractions = sorted(list(set(([float(x) for x in df['fraction'].tolist()]))))
    for fraction in fractions:
        results_df = caller_performance(df, targets_size, allowed_callers, fractions=[str(fraction)], var_type=var_types)
        print(results_df.to_csv(sep='\t'))

    # results_df = caller_performance(df, targets_size, allowed_callers, var_type=var_types)
    # print(results_df.to_csv(sep='\t'))

    plot_performance_per_bin(df, callers, var_types, targets_size, output_png)

if __name__ == '__main__':
    cli()
