"""Automatically find variant plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
from inspect import getmembers, isfunction

from nose.plugins.skip import SkipTest
from nose.tools import nottest

import mitty.lib


@nottest
def test_wrapper(func):
  func[1]()


@nottest
def plugin_has_no_tests(_):
  raise SkipTest('No tests')


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def self_test_all_found_plugins():
  """Plugin self test"""
  for name, module in mitty.lib.discover_all_reads_plugins():
    model = mitty.lib.load_reads_plugin(name)
    tests = [v for v in getmembers(model, isfunction) if v[0].startswith('test')]
    if len(tests) == 0:
      plugin_has_no_tests.description = name + ' (read plugin) self test(s)'
      yield plugin_has_no_tests, None
    else:
      for test in tests:
        test_wrapper.description = name + ' (read plugin) self test(s): ' + (test[1].func_doc or test[1].__name__)
        # We can't ensure that a dev will provide us with a function doc, so we use the name if can't find a doc string
        yield test_wrapper, test