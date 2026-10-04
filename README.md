# two-distance-sets-over-finite-fields
This code is supplementary material for the paper Two-distance sets over finite fields by József Solymosi and Chi Hoi Yip.  
It gives the exact SageMath verification used in Proposition 4.5, in the subsection “Exact verification of the finite problem.”  

Run with
sage d6_finite_check.sage

The computation:
1. generates all graphs on eight vertices up to isomorphism using nauty;
2. uses exact resultants over \(\mathbb Z[t]\) to reduce to finitely many triples \((\Gamma,p,\lambda)\);
3. constructs the corresponding compatibility graphs over \(\mathbb F_p\);
4. computes their clique numbers exactly using Cliquer.

The expected output is recorded in output.txt. In particular, the program constructs \(4{,}075\) auxiliary graphs and finds maximum clique number \(17\), ruling out the required clique of order \(20\).
