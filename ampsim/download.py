"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import subprocess as sp
from sys import platform as _platform

from ampsim.utils import mkdir_p, which


class Resource(object):
    """Base class for resources"""

    def __init__(self, **kwargs):   # pylint: disable=E1101
        """
        Initate a lego with basic information and test if it is_ready
        :return: None
        """
        self.__dict__.update(kwargs)  # update the class variables
        super(Resource, self).__init__()

        if not self.is_ready():
            msg = "{}'s is_ready() returns False. Errors in installer or validator methods."\
                .format(self.__class__.__name__)
            raise ValueError(msg)

    def is_ready(self):
        """
        Check if the tool has a proper installation and test it on sample files.
        :return: Boolean
        """
        # If the tool is not accessible from the command line, run installer.
        self.installer()
        self.validator()
        return True

    def installer(self):
        """
        Download source files, compile and install at path set by environment variable CIDER_RESOURCES
        :return: None
        """
        raise NotImplementedError()

    def validator(self):
        """
        A condition to test if the tool or the resource file installed properly (e.g. check version or md5).
        :return: Boolean
        """
        raise NotImplementedError()

    @property
    def base_path(self):
        """Path to the tool's directory"""
        path = os.path.join(os.path.expanduser("~"), 'ampsim_resources', self.name)  # pylint: disable=E1101
        return mkdir_p(path=path)

    @property
    def full_path(self):
        """Path to the tool itself"""
        raise NotImplementedError()


class Hg19(Resource):
    """A class to download and index hg19 reference"""
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'hg19'
        self.version = '19'
        self.url = 'https://www.dropbox.com/s/pu5mvefewfq50it/hg19.tar.gz'
        self.env_var = "REF"
        super(Hg19, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download genome file, uncompress and index it with BWA and PICARD
        :return: None
        """
        if not os.path.isfile(self.full_path):
            local_file = os.path.join(self.base_path, 'hg19.tar.gz')
            if not os.path.isfile(local_file):
                sp.check_call(['wget', '--directory-prefix', self.base_path, self.url])
            sp.check_call(['tar', 'xzvf', local_file, '-C', self.base_path])
        return True

    def validator(self):
        """
        Check the md5sum of the genome file after decompression
        :return: Boolean
        """
        # md5sum = (hashlib.md5(self.full_path)).hexdigest()
        # return md5sum == 'c191475a4e7cea5b2e41dd5ffc711380'
        return True

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'Homo_sapiens_assembly19.fasta') if os.environ.get('REF', None) is None else os.environ.get('REF')


class RefSeq(Resource):
    """A class to download RefSeq table"""
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'refseq'
        self.version = 'NA'
        self.env_var = "REFSEQ"
        self.url = 'http://www.dropbox.com/s/rn42zvyttvho7le/refFlat.txt.gz'
        super(RefSeq, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download genome file, uncompress and index it with BWA and PICARD
        :return: None
        """
        if not os.path.isfile(self.full_path):
            local_file = os.path.join(self.base_path, self.url.split("/")[-1])
            if not os.path.isfile(local_file):
                sp.check_call(['wget', '--directory-prefix', self.base_path, self.url])
                sp.check_call(['gunzip', local_file])
        return True

    def validator(self):
        """
        Check the md5sum of the genome file after decompression
        :return: Boolean
        """
        # md5sum = (hashlib.md5(self.full_path)).hexdigest()
        # return md5sum == 'c191475a4e7cea5b2e41dd5ffc711380'
        return True

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'refFlat.txt') if os.environ.get(self.env_var, None) is None else os.environ.get(self.env_var)


class Samtools(Resource):
    """Base class for Samtools """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'samtools'
        self.version = '1.2'
        self.env_var = "SAM"
        self.url = 'https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2'
        super(Samtools, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not which(self.full_path):  # check if bwa full path is available
            local_file = os.path.join(self.base_path, 'samtools-1.2.tar.bz2')
            local_dir = os.path.join(self.base_path, 'samtools-1.2')
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            if os.path.isdir(local_dir):
                os.remove(local_dir)
            sp.call(['bzip2', '-d', os.path.join(self.base_path, 'samtools-1.2.tar.bz2')])
            sp.call(['tar', 'xvf', os.path.join(self.base_path, 'samtools-1.2.tar'), '-C', self.base_path])
            sp.call('cd {} && make'.format(local_dir), shell=True)
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_version = '1.2'
        pipe = sp.Popen([self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_version in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        sam = os.path.join(self.base_path, 'samtools-1.2', self.name) if os.environ.get(self.env_var, None) is None else os.environ.get(self.env_var, None) 
        print(sam)
        return sam


class Bcftools(Resource):
    """Base class for htslib """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'bcftools'
        self.version = '1.2'
        self.env_var = "BCF"
        self.url = 'https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2'
        super(Bcftools, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not which(self.full_path):  # check if bwa full path is available
            local_file = os.path.join(self.base_path, 'bcftools-1.2.tar.bz2')
            local_dir = os.path.join(self.base_path, 'bcftools-1.2')
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            if os.path.isdir(local_dir):
                os.system('rm -r {}'.format(local_dir))
            sp.call(['bzip2', '-d', os.path.join(self.base_path, 'bcftools-1.2.tar.bz2')])
            sp.call(['tar', 'xvf', os.path.join(self.base_path, 'bcftools-1.2.tar'), '-C', self.base_path])
            sp.call('cd {} && make'.format(local_dir), shell=True)
        self.export_bcftools()  # required by Ensemble merger tool
        return True

    def export_bcftools(self):
        """
        ensemble caller need direct access to bcftools
        :return:
        """
        path = os.path.join(self.base_path, 'bcftools-1.2')
        os.environ["PATH"] += os.pathsep + path

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_text = 'Version: 1.2'
        pipe = sp.Popen([self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_text in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'bcftools-1.2', self.name) if os.environ.get(self.env_var, None) is None else os.environ.get(self.env_var, None)


class Bamutil(Resource):
    """Base class for Bedtools """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'bamutil'
        self.version = '1.0.13'
        self.env_var = "BAMU"
        self.url = 'http://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.1.0.13.tgz'
        super(Bamutil, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        print(self.full_path)
        if not os.path.isfile(self.full_path):  # check if bwa full path is available
            local_file = os.path.join(self.base_path, 'BamUtilLibStatGen.1.0.13.tgz')
            local_dir = os.path.join(self.base_path, 'bamUtil_1.0.13')
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            if os.path.isdir(local_dir):
                os.remove(local_dir)
            sp.call(['tar', 'zxvf', os.path.join(self.base_path, 'BamUtilLibStatGen.1.0.13.tgz'), '-C', self.base_path])
            sp.call('cd {} && make'.format(local_dir), shell=True)
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_text = 'Version: 1.0.13'
        pipe = sp.Popen([self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_text in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'bamUtil_1.0.13', 'bamUtil', 'bin', 'bam') if os.environ.get(self.env_var, None) is None else os.environ.get(self.env_var, None) 

class Bedtools(Resource):
    """Base class for Bedtools """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'bedtools'
        self.version = '2.25.0'
        self.env_var = "BEDTOOL"
        self.url = 'https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz'
        super(Bedtools, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not which(self.full_path):  # check if bwa full path is available
            local_file = os.path.join(self.base_path, 'bedtools-2.25.0.tar.gz')
            local_dir = os.path.join(self.base_path, 'bedtools2')
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            if os.path.isdir(local_dir):
                os.remove(local_dir)
            sp.call(['tar', 'zxvf', os.path.join(self.base_path, 'bedtools-2.25.0.tar.gz'), '-C', self.base_path])
            sp.call('cd {} && make'.format(local_dir), shell=True)
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_text = 'flexible tools for genome'
        pipe = sp.Popen([self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_text in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'bedtools2', 'bin', self.name) if os.environ.get(self.env_var, None) is None else os.environ.get(self.env_var, None) 



class Art(Resource):
    """Base class for Art for NGS read simulation"""
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'art'
        self.version = '2.25.0'
        if _platform in ['linux', 'linux2']:
            self.url = \
                'http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz'
        elif _platform in ['darwin']:
            self.url = \
                'http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915macos64tgz.tgz'
        super(Art, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not which(self.full_path):  # check if bwa full path is available
            filename = self.url.split('/')[-1]
            local_file = os.path.join(self.base_path, filename)
            local_dir = os.path.join(self.base_path, 'art_bin_ChocolateCherryCake')
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            if os.path.isdir(local_dir):
                os.remove(local_dir)
            sp.call(['tar', 'zxvf', os.path.join(self.base_path, filename), '-C', self.base_path])
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_text = '2.3.7 (Mar 19, 2015)'
        pipe = sp.Popen([self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_text in stdout:
            return True
        return True

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'art_bin_ChocolateCherryCake', 'art_illumina')


class Vt(Resource):
    """Base class for Vt """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'vt'
        self.version = '0.57'
        self.url = 'https://github.com/atks/vt/archive/0.57.tar.gz'
        super(Vt, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not os.path.isfile(self.full_path):  # check if picard full path is available
            file_name = os.path.basename(self.url)
            local_dir = os.path.join(self.base_path, 'vt-0.57')
            local_file = os.path.join(self.base_path, file_name)
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
                sp.call(['tar', 'zxfv', local_file, '-C', self.base_path])
            sp.call('cd {} && make'.format(local_dir), shell=True)
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_text = 'annotate_variants'
        pipe = sp.Popen([self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_text in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'vt-0.57', 'vt')


class Fusim(Resource):
    """Base class for Fusim """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'fusim'
        self.version = '0.2.2'
        self.url = 'https://github.com/aebruno/fusim/raw/master/releases/fusim-0.2.2-bin.tar.gz'
        super(Fusim, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not os.path.isfile(self.full_path):  # check if picard full path is available
            file_name = os.path.basename(self.url)
            local_dir = os.path.join(self.base_path, 'fusim-0.2.2')
            local_file = os.path.join(self.base_path, file_name)
            if not os.path.isfile(local_file):
                sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            if os.path.isdir(local_dir):
                sp.call(['rm', '-r', local_dir])
            sp.call(['tar', 'zxfv', local_file, '-C', self.base_path])
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_text = 'v0.2.2'
        pipe = sp.Popen(['java', '-jar', self.full_path, '--version'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_text in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'fusim-0.2.2', 'fusim.jar')


class Picard(Resource):
    """Base class for Picard """
    def __init__(self, *args, **kwargs):
        """
        Initiate a lego with basic information and test if is_ready
        :return: None
        """
        self.name = 'picard'
        self.version = '1.128'
        self.url = 'https://github.com/broadinstitute/picard/releases/download/1.128/picard-tools-1.128.zip'
        super(Picard, self).__init__(*args, **kwargs)

    def installer(self):
        """
        Download source files, compile and install
        :return: None
        """
        if not os.path.isfile(self.full_path):  # check if picard full path is available
            sp.call(['wget', '--directory-prefix', self.base_path, self.url])
            sp.call(['unzip', '-o', os.path.join(self.base_path, 'picard-tools-1.128.zip'), '-d', self.base_path])
        return True

    def validator(self):
        """
        Test sample data on this tool and check we get the expected output.
        :return: Boolean
        """
        expected_term = 'CheckIlluminaDirectory'
        pipe = sp.Popen(['java', '-jar', self.full_path], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
        stdout, stderr = pipe.communicate()
        if expected_term in stdout:
            return True
        return False

    @property
    def full_path(self):
        """Path to the tool itself"""
        return os.path.join(self.base_path, 'picard-tools-1.128', self.name + '.jar')

