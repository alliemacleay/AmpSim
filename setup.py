"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import

from setuptools import setup, find_packages
import subprocess as sp
exec(open('ampsim/version.py').read())

def is_downloadable(version):
    """ Can the package be installed from bitbucket? """
    version = version.replace('v', '').replace('V', '')
    try:
        int(version[0])
    except ValueError:
        return False
    return True


def get_pip_version(module):
    """ Return installed version from pip freeze """
    resp = sp.Popen(['pip', 'freeze'], stdout=sp.PIPE)
    try:
        ver = sp.check_output(['grep', module], stdin=resp.stdout)
        resp.wait()
        ver = ver.strip().split('=')[-1].replace('v', '').replace('V', '')
    except sp.CalledProcessError:
        ver = None
    return ver


def install_from_git(url):
    """ Update packages from url """
    simple_name = url.split('=')[-1]
    if '@' in url:
        exp_version = url.split('@')[-1].split('#')[0].replace('v', '').replace('V', '')
        version = get_pip_version(simple_name)
        if is_downloadable(exp_version) and version != exp_version:
            sp.call(['pip', 'uninstall', '-y', simple_name])
            sp.call(['pip', 'install', url])
    else:
        print('Check that {} is properly installed'.format(simple_name))


def load_requirements():
    f = open('requirements.txt')
    pkgs = []
    for l in f:
        if l.startswith('-e') or '://' in l:
            install_from_git(l.strip())
        else:
            pkgs.append(l.strip())
    return pkgs

print(load_requirements())
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
