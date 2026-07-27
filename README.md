# two-distance-sets-over-finite-fields
This code is supplementary material for the paper "Two-distance sets over finite fields" by Jozsef Solymosi and Chi Hoi Yip.

This is an exact SageMath verification of Proposition 6.5.

Run with

    sage d6_finite_check.sage

The program follows Section 6.4 directly:
  1. generate all graphs on eight vertices up to isomorphism with nauty;
  2. use exact resultants over Z[t] to reduce to finitely many
     (Gamma, p, lambda);
  3. construct the auxiliary compatibility graphs over GF(p);
  4. use Sage's exact Cliquer interface to compute their clique numbers.

All arithmetic is exact.  The program stops if the resultant reduction
leaves a characteristic-independent set of 20 candidates or if an
auxiliary graph contains a clique of order 20.
