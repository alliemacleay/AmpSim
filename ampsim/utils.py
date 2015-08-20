"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import subprocess as sp
from pkg_resources import resource_filename


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


def bgzip(path):
    bgzip_path = path + '.gz'
    sp.check_call('bgzip -f {}'.format(path), shell=True)
    return bgzip_path


def create_manifest(output_dir, pipeline_name, paths):
    """
    Generate a manifest file for snapshot or onek pipeline
    :param output_dir:
    :param pipeline_name:
    :param paths:
    :return: New manifest files directory
    """
    create_dir(output_dir)
    raw_manifest_path = resource_filename('ampsim.resources', '{}_template.manifest'.format(pipeline_name))
    new_manifest_path = os.path.join(output_dir, pipeline_name + '.manifest')
    manifest_text = open(raw_manifest_path, 'r').read().format(*paths)
    with open(new_manifest_path, 'w') as fileobj:
        fileobj.write(manifest_text)
    return new_manifest_path


