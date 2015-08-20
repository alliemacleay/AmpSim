"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import


from setuptools import setup, find_packages
exec(open('ampsim/version.py').read())


def load_requirements():
    f = open('requirements.txt')
    return [l.strip(' ') for l in f]

setup(
    name='Ampsim',
    version=__version__,
    py_modules=['cider'],
    packages=find_packages(),
    include_package_data=True,
    url='https://bitbucket.org/mghcid/ampsim',
    author='Saeed Al Turki',
    author_email='salturki@gmail.com',
    license='MIT',
    install_requires=load_requirements(),
    entry_points='''
        [console_scripts]
        ampsim=ampsim.main:cli
    ''',
)
