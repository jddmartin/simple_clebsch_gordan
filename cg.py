"""A simple numerical Clebsch-Gordan calculator using recursion relations

See: "Notes on the numerical calculation of Clebsch-Gordan coefficients"
for the relevant equations referenced in the code.

See "test_cg.py" for testing.

Note that this is not an efficient implementation.  Memoization could be
used to improve efficiency by avoiding repeated calls to determine the
same coefficients.  Also, the recursive approach used here will be more
efficient in a language that performs "tail call elimination" (e.g. Haskell).

J. D. D. Martin, Fall 2015.
"""

import math

def cg(j1, j2, m1, m2, j, m):
    """Compute the Clebsch-Gordan coefficient: <j1 j2; m1, m2 | j, m>"""
    if ((m1+m2) != m): return 0.0
    if ((j < abs(j1-j2)) or (j > j1+j2)): return 0.0
    if not (_valid_pair(j1, m1) and _valid_pair(j2, m2) and _valid_pair(j, m)):
        return 0.0
    if ((m1+m2) == j):
        if m1 == j1:
            return _cg_special_case(j1, j2, j)
        else: # Use Eq. 4
            return (-_c(1.0,j2,m2-1)/_c(1.0,j1,m1) *
                    cg(j1, j2, m1+1, m2-1, j, j))
    else: # Use Eq. 3
        # for those about to recurse ... we salute you
        return 1.0/_c(-1.0,j,m+1)*(
            _c(-1.0,j1,m1+1)*cg(j1,j2,m1+1,m2,j,m+1) +
            _c(-1.0,j2,m2+1)*cg(j1,j2,m1,m2+1,j,m+1)
            )
    
def _cg_special_case(j1, j2, j):
    """Compute <j1 j2; j1, j-j1 | j, j>"""
    def f(m1): return (_c(1,j2,j-m1)/_c(1,j1,m1-1))**2 # Eq. 10
    def s(m1):
        if ((abs(m1-1) > j1) or (abs(j-(m1-1)) > j2)):
            return 1.0
        else: # Eq. 11
            return 1.0+f(m1)*s(m1-1)
    return 1.0/math.sqrt(s(j1)) # Eq. 12

def _c(sign,j,m):
    """Helper function, Eq. 2"""
    return math.sqrt(j*(j+1)-m*(m+sign))

def _valid_pair(j, m):
    """Check that j and m are valid angular momentum quantum numbers"""
    if ((j < 0.0) or (not (2.0*j).is_integer())): return False
    if (not (2.0*m).is_integer()): return False
    if ((abs(m) > j)): return False
    if (((2*j) % 2) != ((2*abs(m)) % 2)): return False
    return True

