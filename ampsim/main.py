"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import


import os
import sys
from datetime import date
import pkg_resources

import click

from ampsim.version import __version__

###################################
# Set global logging configurations
###################################
import logging
from logging.config import dictConfig

logging_config = dict(
    version=1,
    formatters={
        'general': {'format': '%(levelname)s: %(asctime)s: %(message)s',
                    'datefmt': '%Y-%m-%d %H:%M:%S'}
    },
    handlers={
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'general',
            'level': os.environ.get('AMPSIM_LOGGER', 'DEBUG'),
        },
        'file': {
            'class': 'logging.handlers.RotatingFileHandler',
            'level': os.environ.get('AMPSIM_LOGGER', 'DEBUG'),
            'formatter': 'general',
            'filename': pkg_resources.resource_filename('ampsim', 'ampsim.log'),
            'maxBytes': "10000000",  # rotate the file when size reaches 10 Mb
            'backupCount': "1",  # backup one file only e.g. cider.log.1
        },
    },
    loggers={
        'cider': {
            'handlers': ['file', 'console'],
            'level': os.environ.get('AMPSIM_LOGGER', 'DEBUG')
        }
    }
)
dictConfig(logging_config)
logger = logging.getLogger(__name__)
###################################

msg = 'AmpSim is a command-line to simulate SNVs and INDELs for Amplicon NGS data in target regions' \
      'Current version is v{0} ({1})'.format(__version__, date.today().year)


folders_of_interest = [
    os.path.join(os.path.join(os.path.dirname(__file__)), 'tools'),
]

cli_files = {}


def has_cli_method(script_path):
    """
    Check if a script has a cli() method in order to add it to the main
    :param script_path: to a python script inside Cider packages
    :return: Boolean
    """
    file_obj = open(script_path, 'r').read()
    if "cli()" in file_obj:
        return True
    else:
        return False


class MyCLI(click.MultiCommand):
    """CLI main class"""
    def list_commands(self, ctx):
        rv_obj = []
        for folder in folders_of_interest:
            rv_part = []
            for filename in os.listdir(folder):
                if filename.endswith('.py') and not filename.startswith("__init__"):
                    if not has_cli_method(os.path.join(folder, filename)):
                        continue
                    rv_part.append(filename[:-3])
                    cli_files[filename[:-3]] = folder
            rv_part.sort()
            rv_obj.extend(rv_part)  # to sort pipelines then helpers instead of mixing them when help message is printed
        return rv_obj

    def get_command(self, ctx, name):
        ns_obj = {}
        if not cli_files:
            self.list_commands(ctx)
        try:
            file_path = os.path.join(cli_files[name], name + '.py')
        except ValueError, err:
            sys.stderr.write('ERROR: Unknown command %s\n' % str(err))
            sys.exit(1)
        with open(file_path) as file_obj:
            code = compile(file_obj.read(), file_path, 'exec')
            eval(code, ns_obj, ns_obj)  # pylint: disable=W0123
        return ns_obj['cli']


cli = MyCLI(help=msg)

if __name__ == "__main__":
    cli()
