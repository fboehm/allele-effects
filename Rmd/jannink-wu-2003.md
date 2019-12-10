Methods from JANNINK and WU (2003)
================
Frederick J. Boehm
12/6/2019

Last modified: 2019-12-06 11:17:07.

## Overview

We write code below to implement Jannink and Wu’s methods (JANNINK and
WU 2003).

They use a scalar Metropolis-Hastings procedure.

Since we have no \(X\), per their notation, in our model, we omit it
from the MCMC.

In Jannink and Wu’s notation, \(Q\) is our HMM-derived allele
probabilities matrix.

\(C\) is the matrix that does the “collapsing”, from 8 (for CC and DO)
alleles to \(l\) alleles.

\(a\) is the allele effects vector.

Thus, for MCMC, we take the procedure of page 135 of JANNINK and WU
(2003) and modify it to:

1.  update Q (allele probabilities matrix)
2.  update allele effects a
3.  update residual variance
4.  update C given l
5.  update l, number of alleles (between 2 and 8 for DO & CC)

## update QTL inheritance matrix Q

JANNINK and WU (2003) cite a 1994 paper by RC Jansen. I believe this is
a mis-specified citation. The article that they cite doesn’t contain any
discussion of full conditionals. Unless they’re citing it for another
reason…

## update QTL allelic effects a

## update family means b and residual variance s2

## update QTL allelic configuration conditional on allelic number, C| l;

## update the number of alleles l

<div id="refs" class="references">

<div id="ref-jannink2003estimating">

JANNINK, JEAN-LUC, and XIAO-LIN WU. 2003. “Estimating Allelic Number and
Identity in State of Qtls in Interconnected Families.” *Genetics
Research* 81 (2). Cambridge University Press: 133–44.

</div>

</div>
