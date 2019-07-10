### How many homorepeats of type $a$ with $n$ repeat units do we expect in a sequence of length $s$?

This problem is described by sequential runs of successes in a Bernoulli trial. Exact solutions to
the expected value and variance of the number of runs of a given length $n$ in a bounded sequence
of length $s$ are derived in, e.g., [Makri, F. S., & Psillakis, Z. M. (2011). On success runs of a fixed length in Bernoulli sequences: Exact and asymptotic results. Computers & Mathematics with Applications, 61(4), 761–772. http://doi.org/10.1016/j.camwa.2010.12.023
](http://www.sciencedirect.com/science/article/pii/S0898122110009284). We implemented the derived expressions in Python3.  

The probability of amino acid $a$ is equated to the probability of a success, and the expected values and variances can
be derived for all sequences or subsequences of different lengths in the sequence set.
This calculation is repeated for every amino acid in the sequence set.

#### Calculating expected homorepeat run lengths for a set of regions.

In short: Just add up the expected values for all regions.


#### Our derivations

On a sidenote, we derived the above described expected values, disregarding the effect of the bound on the region.
Thus, these results are a) of little use compared to the above results and b) only valid in good approximation when $n$ << $s$.

(Todo: Better definitions. Expected values are calculated for random variables, not parameters.)

$a$:          An an amino acid. $a \in [A,C,D,...]$  
$p(a)$:       The probability distribution of amino acids in the sequence.  
$s$:          The length of the entire sequence.  
$l$:          The repeat unit length. Here, we are only looking at homorepeats. Therefore, $l = 1$.  
$n$:          The number of repeat units. For homorepeats, $n$ equals the length of the entire repeat region.  
$N$:          Random variable, describes the number of repeat units in a homorepeat obtained by drawing amino acids independently from the given distribution.  
$E[f(a|s)]$   The expected number of amino acids $a$ in a sequence of length s.  
$p_a(n)$      The probability distribution for different $n$ of tandem repeats of type $a$.  
$E_a[N]$      The expected number of repeat units for a tandem repeat of type $a$.  
$E_a[C]$      The expected count of tandem repeats of type $a$ in the sequence.  

We derive:

$E[f(a|s)] = p(a) \cdot s$

$p_a(n) = p(a)^n \cdot (1 - p(a))$

$E_a[N] = \sum_{n=1}^{\inf} n \cdot p(a)^n \cdot (1-p(a))$

$E_a[C] = \frac{E[f(a|s)]}{E_a[N]}$

And finally, what we are really interested in: How many homorepeats of type $a$ with $n$ repeat units to expect in a sequence of length $s$:

$E_a[c(n=k)] = p_a(n) \cdot E_a[C] = \frac{p(a)^2 \cdot s}{\sum_{n=1}^{\inf} n \cdot p(a)^n \cdot (1-p(a))}$
