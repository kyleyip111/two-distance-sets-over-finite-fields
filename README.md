# Two-distance sets over finite fields

This repository contains the SageMath code accompanying the paper  
**Two-distance sets over finite fields**.

The computation verifies the finite graph problem arising in the exceptional
case of dimension 6. In particular, it rules out a 28-point two-distance set
in dimension 6 over every finite field of odd characteristic.

## Requirements

The computation uses SageMath, including its interfaces to:

- `nauty`, to generate graphs on eight vertices up to isomorphism;
- `Cliquer`, to compute clique numbers exactly.

All polynomial, resultant, finite-field, and graph computations are exact. No
floating-point arithmetic is used.

## Running the computation

From the repository directory, run

```bash
sage d6_finite_check.sage
```

## What the program does

A hypothetical 28-point two-distance set in dimension 6 leads to a finite
problem involving:

- a graph on eight vertices;
- an odd prime `p`;
- a scalar `lambda` over the prime field;
- a collection of binary vectors of length eight satisfying certain quadratic,
  affine, and pairwise compatibility conditions.

The program verifies that no such system exists.

The computation proceeds in four stages.

1. **Generate eight-vertex graphs.**  
   All simple graphs on eight vertices are generated up to isomorphism using
   `nauty`.

2. **Reduce to finitely many characteristics.**  
   Exact integer polynomials and resultants are used to produce a finite list
   of possible odd primes. This avoids imposing any arbitrary upper bound on
   the characteristic.

3. **Construct admissible binary vectors.**  
   For each surviving triple `(Gamma, p, lambda)`, the program checks the
   required quadratic, affine, and isotropic conditions over the corresponding
   prime field.

4. **Solve the compatibility problem.**  
   The admissible binary vectors are used as the vertices of an auxiliary
   graph. Two vectors are adjacent when they satisfy the required pairwise
   compatibility condition. A hypothetical 28-point configuration would
   produce a clique of size 20. The clique number of every auxiliary graph is
   computed exactly using `Cliquer`.

## Output

A complete run gives:

- 12,346 isomorphism classes of graphs on eight vertices;
- 13,809 relevant polynomial factors;
- 60 exceptional odd primes, the largest being 1367;
- 17,347 candidate triples `(Gamma, p, lambda)`;
- 10,867 nonsingular triples;
- 4,075 auxiliary compatibility graphs.

The largest clique number among all 4,075 auxiliary graphs is **17**.

Since a hypothetical 28-point configuration would require a clique of size 20,
no such configuration exists. This completes the computational part of the
dimension-6 argument.

## Reproducibility

The paper cites a fixed commit of this repository so that the exact version of
the code used in the proof remains permanently identifiable.

For the mathematical reduction leading to this computation, see the
dimension-6 section of the paper, especially the subsection **Exact
verification of the finite problem**.
