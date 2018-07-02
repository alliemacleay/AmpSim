"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import fnmatch
import subprocess as sp
import errno
import vcf


def create_dir(path):
    """
    Delete a path to directory and create it again.
    :param path:
    :return: path
    """
    if os.path.isdir(path):
        os.system('rm -fr {}'.format(path))
    os.system('mkdir -p {}'.format(path))
    return path


def mkdir_p(path):
    """
    Behave like 'mkdir -p' in unix shell.
    :param path:
    :return: path
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    return path


def strip_chrom(chrom):
    """
    Remove Chr or chr from the chromosome id
    :param chrom: String
    :return: String
    """
    if 'chr' in chrom.lower():
        return chrom[3:]
    else:
        return chrom


def which(program):
    """
    Simulate which command
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    :param program: String
    :return: Path to program or None if it doesn't exist
    """
    def is_exe(fpath):
        """check if file exists"""
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath = os.path.split(program)[0]
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def bgzip(path):
    """
    Run bgzip and return the file path plus the new .gz extension
    :param path: Input file path
    :return: Output file path
    """
    bgzip_path = path + '.gz'
    bgzip_tool = os.getenv('BGZIP', 'bgzip')
    sp.check_call([bgzip_tool, '-f', path])
    return bgzip_path


def create_manifest(output_dir, pipeline_name, paths):
    """
    Generate a manifest file for snapshot or onek pipeline
    :param output_dir:
    :param pipeline_name:
    :param paths:
    :return: New manifest files directory
    """
    if not os.path.isdir(output_dir):
        create_dir(output_dir)
    new_manifest_path = os.path.join(output_dir, pipeline_name + '.manifest')
    header = """##fileformat=SOMATIC,1.0
IN_PIPELINE	OUT_PIPELINE	BATCH_ID	CID_ID	NAE_ID	SPECIMEN_ID	SPECIMEN_TYPE	P5_BARCODE	P7_BARCODE	RUN_FOLDER	R1	R2	R3	R4	BAM	VCF"""
    with open(new_manifest_path, 'w') as fileobj:
        fileobj.write(header + "\n")
        for idx, fraction in enumerate(paths.keys()):
            somatic_fq1 = paths[fraction]['somatic'].get('fastq', '')[0]
            somatic_fq2 = paths[fraction]['somatic'].get('fastq', '')[1]
            somatic_bam = paths[fraction]['somatic'].get('bam', '')
            somatic_vcf = paths[fraction]['somatic'].get('vcf', '')
            if pipeline_name == 'onek':
                normal_fq1 = paths[fraction]['normal'].get('fastq', '')[0]
                normal_fq2 = paths[fraction]['normal'].get('fastq', '')[1]
                normal_bam = paths[fraction]['normal'].get('bam', '')
                normal_vcf = paths[fraction]['normal'].get('vcf', '')
                # make a copy of the normal and save it with its somatic pair
                line = [pipeline_name, '', str(fraction), 'cid_{}'.format(idx), 'nae_{}'.format(idx),
                        'specimen_{}'.format(fraction), 'normal', 'A{}'.format(idx), 'P{}'.format(idx), '',
                        normal_fq1, normal_fq2, '', '',
                       normal_bam, normal_vcf]
                fileobj.write("\t".join(line) + "\n")
            line = [pipeline_name, '', str(fraction), 'cid_{}'.format(idx), 'nae_{}'.format(idx),
                    'specimen_{}'.format(fraction), 'tumor', 'A{}'.format(idx), 'P{}'.format(idx), '',
                    somatic_fq1, somatic_fq2, '', '',
                   somatic_bam, somatic_vcf]
            fileobj.write("\t".join(line) + "\n")
    return new_manifest_path


def safe_divide(nm, dm):
    """
    To avoid dividing by zero
    :param nm: Nominator
    :param dm: Dominator
    :return: Float
    """
    if dm == 0:
        return 0.
    else:
        return nm / float(dm)


def get_var_type_size(ref, alt):
    """
    Calculate the variant type SNV, INS or DEL and length
    :param ref: String
    :param alt: String
    :return: var_type and var_size
    """
    var_size = len(alt) - len(ref)
    if var_size == 0:
        var_type = 'SNV'
    elif var_size > 0:
        var_type = 'INS'
    else:
        var_type = 'DEL'
    return var_type, var_size


def vcf_parser(vcf_path):
    """
    Parse VCF records and return a generator
    :param vcf_path: Path to the vcf file
    :return: list
    """
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    for record in vcf_reader:
        if record.FILTER:
            continue
        if not (record.is_snp or record.is_deletion):
            continue  # ignore cnv and sv variants
        chrom, start, ref, alt = record.CHROM, record.POS, record.REF, record.ALT[0]
        if alt is None:
            continue
        key = "_".join([str(x) for x in [chrom, start, ref, alt]])
        var_type, var_size = get_var_type_size(ref, alt)
        end = start + abs(var_size)
        yield key, chrom, start, end, ref, alt, var_type, var_size


def get_vcfs(path):
    vcf_paths = []
    dirs = ['CALLING_MUTECT2',
            'CALLING_MUTECT',
            'CALLING_LOFREQ',
            # '12.VARSCAN_CALL',
            'CALLING_HOTSPOTTER',
            'CALLING_GATK',
            'VCF_ENSEMBLER']
    pattern = '*.vcf.gz'
    for d in dirs:
        dir_path = os.path.join(path, d)
        files = [os.path.join(dir_path, filename) for filename in os.listdir(dir_path)]
        vcf_paths.extend(fnmatch.filter(files, pattern))
    return vcf_paths
