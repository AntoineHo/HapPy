from unittest import TestCase
from happy.Hap import test

class TestConsole(TestCase):
    def test_basic(self):
        s = test()
        self.assertTrue(s == True)
