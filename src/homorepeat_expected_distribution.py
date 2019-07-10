"""
Calculate the expected distribution of homorepeats in sequence, modelling character distribution with bernoulli trials.
"""

import logging
LOG = logging.getLogger(__name__)


def check_parameters(n, k, p):
    if p < 0 or p > 1:
        raise Exception("The input probability p: {} is not sensible.".format(p))
    if k > n or n < 1 or k < 1:
        raise Exception("The input parameters n: {}, k: {} are not sensible.".format(n,k))


def expected_frequency_runs(k, p):
    """
    Calculate the expected frequency of single character runs of length `k` in an unbounded sequence,
     given the probability of the character `p`.

    Similar to Proposition 2.1 (part I, expected value) of
    Makri, F. S., & Psillakis, Z. M. (2011). On success runs of a fixed length in Bernoulli sequences: Exact and asymptotic results. Computers & Mathematics with Applications, 61(4), 761–772. http://doi.org/10.1016/j.camwa.2010.12.023

    """
    return ((1-p) ** 2) * (p ** k)


def expected_number_runs(n, k, p):
    """
    Calculate the expected number of single character runs of length `k` in a sequence of length
    `n`, given the probability of the character `p`.

    Implements Proposition 2.1 (part I, expected value) of
    Makri, F. S., & Psillakis, Z. M. (2011). On success runs of a fixed length in Bernoulli sequences: Exact and asymptotic results. Computers & Mathematics with Applications, 61(4), 761–772. http://doi.org/10.1016/j.camwa.2010.12.023

    """

    check_parameters(n, k, p)

    if n == k:
        # Fix me
        return p ** k
    elif n >= k + 1:
        return (1-p) * (p ** k) * (2 + (n-k-1)*(1-p))
    else:
        raise Exception("The input parameters n: {}, k: {}, p: {} are not sensible".format(n,k,p))


def test_expectation():

    assert expected_number_runs(1,1,0.5) == 0.5
    assert expected_number_runs(1,1,0.3) == 0.3
    assert expected_number_runs(2,2,0.5) == 0.25
    assert expected_number_runs(2,1,0.5) == 0.5
    assert expected_number_runs(3,1,0.5) == 5/8


def variance_number_runs(n, k, p):
    """
    Calculate the variance of the number of single character runs of length `k` in a sequence of length
    `n`, given the probability of the character `p`.

    Implements Proposition 2.1 (part II, variance) of
    Makri, F. S., & Psillakis, Z. M. (2011). On success runs of a fixed length in Bernoulli sequences: Exact and asymptotic results. Computers & Mathematics with Applications, 61(4), 761–772. http://doi.org/10.1016/j.camwa.2010.12.023

    """

    check_parameters(n, k, p)

    q = 1 - p

    if n == k:
        # Fix me
        return (p**k) * (1-(p**k))
    elif n <= 2*k or n == 2*k+1:
        v1 = 2*q(p**k) -4*(q**2)(p**(2*k)) +(n-k-1)*(q**2)(p**k) -(n-k-1)*2*(q**4)(p**(2*k)) -4*(n-k-1)*(q**3)*(p**(2*k))
        if n <= 2*k:
            return v1
        else:
            return v1 + 2*q(p**(2*k))
    elif n >= 2*k+2:
        return 2*q(p**k) +2(q**2)*(p**(2*k)) +(n-k-1)*(q**2)*(p**k) + 2*(n-4*k-4)*(q**3)*(p**(2*k)) -((n-k-1)**2 -(n-2*k-2)*(n-2*k-3)) * (q**4)*(p**(2*k))
    else:
        raise Exception("The input parameters n: {}, k: {}, p: {} are not sensible".format(n,k,p))


def test_variance():

    assert variance_number_runs(1,1,0.5) == 0.25


def expected_number_runs_cumulated_over_ns(ns, ps, k_max=10):
    """
    Calculate the expected number of single character runs for multiple lengths `k` in multiple sequences of lengths
    `ns`, given multiple character probabilities `ps`.

    `k_max` is the maximum k for which the expected number of single character runs is calculated.

    E.g.:
    ps = {"A": 0.5, "B": 0.3, "C": 0.2}
    ns = [1, 2, 3, 2, 17]

    """

    unique_n = set(ns)
    largest_possible_k = max(unique_n)
    result = {aa: {k: 0 for k in range(1, min(k_max,largest_possible_k)+1)} for aa in ps.keys()}

    for aa, p in ps.items():
        print(aa)
        for n in unique_n:
            n_count = ns.count(n)
            for k in range(1, min(k_max,n)+1):
                e = expected_number_runs(n, k, p)
                result[aa][k] += n_count * e

    return result


if __name__ == '__main__':

    p_q = 0.05
    n = 10
    k = 4
    print(expected_number_runs(n, k, p_q))

    ps = {"A": 0.5, "B": 0.3, "C": 0.2}
    ns = [1,2,3,4,5,6,7,8,9,10,2]
    print(expected_number_runs_cumulated_over_ns(ns, ps))

    test_expectation()
    test_variance()