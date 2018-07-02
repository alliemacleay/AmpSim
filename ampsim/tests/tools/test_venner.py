"""
(c) MGH Center for Integrated Diagnostics
"""
from __future__ import print_function
from __future__ import absolute_import

from pkg_resources import resource_filename
from nose.tools import assert_equal
from ampsim.tools import venner


def test_calculate_target_size():
    """Unit test calculate_target_size method"""
    path = resource_filename('ampsim.tests.data', 'targets.bed')
    assert_equal(venner.calculate_target_size(path), 2418)
