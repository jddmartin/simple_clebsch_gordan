"""Tests of cg.py
"""

from __future__ import division, print_function

import sys, math, unittest, time, numpy.testing
from math import sqrt

import cg

class cg_test(unittest.TestCase):
    def runTest(self):
        # test cases are taken from Particle Data Group table:
        # http://pdg.lbl.gov/2015/reviews/rpp2014-rev-clebsch-gordan-coefs.pdf
        test_cases = (
            (1/2, 1/2, 1/2, -1/2, 0, 0, sqrt(1/2)),
            (2, 1/2, 1, -1/2, 3/2, 1/2, sqrt(3/5)),
            (2, 1/2, -2, 1/2, 1.5, -3/2, -sqrt(4/5)),
            (1, 1/2, 0, -1/2, 3/2, -1/2, sqrt(2/3)),
            (3/2, 1/2, 3/2, -1/2, 2, 1, sqrt(1/4)),
            (3/2, 1, 3/2, 0, 5/2, 3/2, sqrt(2/5)),
            (2, 1, 2, -1, 3, 1, sqrt(1/15)),
            (1, 1, 0, 0, 0, 0, -sqrt(1/3)),
            (3/2, 1, -3/2, -1, 5/2, -5/2, 1),
            (2, 3/2, 2, -3/2, 5/2, 1/2, sqrt(6/35)),
            (3/2, 3/2, -1/2, 1/2, 0, 0, sqrt(1/4)),
            (3/2, 3/2, -3/2, 1/2, 1, -1, sqrt(3/10)),
            (2, 3/2, 0, -1/2, 1/2, -1/2, -sqrt(1/5)),
            (2, 2, 2, -2, 4, 0, sqrt(1/70)),
            (2, 2, 1, -2, 4, -1, sqrt(1/14)),
            (2, 2, -2, 1, 3, -1, -sqrt(3/10)),
            (2, 2, -1, -1, 3, -2, 0),
            )
            
        for case in test_cases:
            j1, j2, m1, m2, j, m, test_cg = case
            my_cg = cg.cg(j1,j2,m1,m2,j,m)
            print("< %g %g ; %g %g | %g %g > = " % case[:-1])
            print("    (calc), (expected) %.16g %.16g" % 
                  (my_cg, test_cg))
            numpy.testing.assert_approx_equal(test_cg, my_cg, 14)
            # test permutation symmetry:
            perm_cg = cg.cg(j2,j1,m2,m1,j,m)
            numpy.testing.assert_approx_equal(
                (-1)**(j-j1-j2)*test_cg,
                perm_cg, 14)
            
if __name__ == "__main__":
    unittest.main()
