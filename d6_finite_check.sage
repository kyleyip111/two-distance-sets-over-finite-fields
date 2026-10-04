"""
Exact verification for Proposition 4.5 in Section 4.3.

The notation follows the manuscript:

    Delta_Gamma(t) = det(C_Gamma - t I_8),
    S_Gamma(t)     = 1^T adj(C_Gamma - t I_8) 1,
    F_{Gamma,b}(t) = b^T adj(C_Gamma - t I_8) b + t Delta_Gamma(t),
    G_{Gamma,b}(t) = b^T adj(C_Gamma - t I_8) 1 - Delta_Gamma(t).

For a retained triple (Gamma,p,lambda), the code checks equations
(4.5), (4.7), and (4.8), and then uses equation (4.6) as the
compatibility relation in the final auxiliary graph.

All computations are exact.  The assertions at the end reproduce the
numerical counts quoted in the proof of Proposition 4.5.
"""

from collections import Counter

from sage.all import (
    GF,
    Graph,
    PolynomialRing,
    ZZ,
    gcd,
    graphs,
    identity_matrix,
    matrix,
    vector,
)


# Size of the principal block in Proposition 4.4.
N = 8
TARGET = 20
VERTICES = list(range(N))
ONE_ZZ = vector(ZZ, [1] * N)
BINARY_ZZ = [
    vector(ZZ, [(mask >> i) & 1 for i in range(N)])
    for mask in range(1 << N)
]

# Numerical data reported in Section 4.3.
EXPECTED = {
    "graph_classes": 12346,
    "relevant_pairs": 13809,
    "maximum_Bh": 16,
    "exceptional_primes": 60,
    "largest_exceptional_prime": 1367,
    "raw_triples": 17347,
    "singular_triples": 6480,
    "nonsingular_triples": 10867,
    "auxiliary_graphs": 4075,
    "cases_by_prime": {3: 2455, 5: 1562, 7: 54, 19: 2, 47: 2},
    "maximum_clique": 17,
}

Zt = PolynomialRing(ZZ, "t")
t = Zt.gen()
I_Zt = identity_matrix(Zt, N)
ONE_Zt = vector(Zt, ONE_ZZ)


def check_reported_value(name, value):
    """Check a numerical value against the value quoted in the paper."""
    expected = EXPECTED[name]
    if value != expected:
        raise RuntimeError(
            f"Reported count mismatch for {name}: got {value}, expected {expected}"
        )


def polynomial_data(G):
    """Return Delta_Gamma, S_Gamma, and all 256 pairs (F_{Gamma,b},G_{Gamma,b})."""
    C = G.adjacency_matrix(vertices=VERTICES, base_ring=ZZ)
    M = matrix(Zt, C) - t * I_Zt
    Delta = M.det()
    Adj = M.adjugate()
    Adj_one = Adj * ONE_Zt
    S = ONE_Zt.dot_product(Adj_one)

    FG = []
    for b0 in BINARY_ZZ:
        b = vector(Zt, b0)
        Adj_b = Adj * b
        F = b.dot_product(Adj_b) + t * Delta
        G_b = b.dot_product(Adj_one) - Delta
        FG.append((F, G_b))
    return Delta, S, FG


def resultant_reduction(graph_list):
    """
    Carry out the resultant reduction from the proof of Proposition 4.5.

    Returns the nonsingular triples (graph_index,p,lambda) and a summary of
    the exact counts appearing before the final clique computation.
    """
    raw_triples = set()
    nonsingular_triples = set()
    exceptional_primes = set()
    relevant_pairs = 0
    maximum_Bh = 0

    for graph_index, G in enumerate(graph_list):
        Delta, S, FG = polynomial_data(G)

        # The content of S_Gamma divides its leading coefficient -8, hence it
        # has no odd prime divisor.  Sage may include nonunit integer content
        # among the factors over ZZ; degree-zero factors are not polynomial
        # factors h in the manuscript.
        for h, _multiplicity in S.factor():
            if h.degree() == 0:
                continue

            # A relevant root lambda is not 0 or -1, and h must not divide
            # Delta_Gamma over Q[t].  Since h is primitive and irreducible,
            # resultant(h,Delta_Gamma)=0 is equivalent to h | Delta_Gamma.
            if h(0) == 0 or h(-1) == 0 or h.resultant(Delta) == 0:
                continue

            relevant_pairs += 1
            Bh = 0
            prime_count = Counter()

            for F, G_b in FG:
                rF = ZZ(h.resultant(F))
                rG = ZZ(h.resultant(G_b))

                # This is precisely b in B_h: both divisibilities hold over
                # Z[t], equivalently both resultants vanish over Q.
                if rF == 0 and rG == 0:
                    Bh += 1
                    continue

                # For b not in B_h, r_{h,b} is the positive gcd of the two
                # resultants.  A prime p can support b only if p | r_{h,b}.
                r = abs(gcd(rF, rG))
                for p0 in r.prime_divisors():
                    p = int(p0)
                    if p != 2:
                        prime_count[p] += 1

            maximum_Bh = max(maximum_Bh, Bh)
            if Bh >= TARGET:
                raise RuntimeError(
                    f"|B_h| >= {TARGET} for graph {graph_index} and h = {h}"
                )

            # This is exactly the necessary filter (4.13):
            # |B_h| + #{b not in B_h : p | r_{h,b}} >= 20.
            primes = sorted(
                p for p, count in prime_count.items() if Bh + count >= TARGET
            )
            exceptional_primes.update(primes)

            for p in primes:
                # Direct evaluation is simpler than factoring h modulo p.
                # range(1,p-1) excludes lambda=0 and lambda=-1.
                for lam in range(1, p - 1):
                    if h(lam) % p != 0:
                        continue
                    triple = (graph_index, p, lam)
                    raw_triples.add(triple)
                    if Delta(lam) % p != 0:
                        nonsingular_triples.add(triple)

    singular_count = len(raw_triples - nonsingular_triples)
    summary = {
        "relevant_pairs": relevant_pairs,
        "maximum_Bh": maximum_Bh,
        "exceptional_primes": len(exceptional_primes),
        "largest_exceptional_prime": max(exceptional_primes),
        "raw_triples": len(raw_triples),
        "singular_triples": singular_count,
        "nonsingular_triples": len(nonsingular_triples),
    }

    print("relevant pairs (Gamma,h):", summary["relevant_pairs"])
    print("maximum |B_h|:", summary["maximum_Bh"])
    print("exceptional odd primes:", summary["exceptional_primes"])
    print("largest exceptional prime:", summary["largest_exceptional_prime"])
    print("triples before nonsingularity:", summary["raw_triples"])
    print("singular triples discarded:", summary["singular_triples"])
    print("nonsingular triples:", summary["nonsingular_triples"])

    return sorted(nonsingular_triples), summary


def admissible_vectors(G, p, lam_integer):
    """
    Check (4.8), then return R and all binary vectors satisfying (4.5),(4.7).

    Here B=C-lambda I_8 and R=B^{-1}.  Returning None means that B is
    singular or that equation (4.8) fails for the triple.
    """
    Fp = GF(p)
    lam = Fp(lam_integer)
    C = G.adjacency_matrix(vertices=VERTICES, base_ring=Fp)
    B = C - lam * identity_matrix(Fp, N)
    if B.det() == 0:
        return None

    R = B.inverse()
    one = vector(Fp, [1] * N)

    # Equation (4.8): 1^T R 1 = 0.
    if one.dot_product(R * one) != 0:
        return None

    candidates = []
    for b0 in BINARY_ZZ:
        b = vector(Fp, b0)

        # Equation (4.5): b^T R b = -lambda.
        if b.dot_product(R * b) != -lam:
            continue

        # Equation (4.7): b^T R 1 = 1.
        if b.dot_product(R * one) != 1:
            continue

        candidates.append(b)

    return R, candidates


def compatibility_graph(R, candidates):
    """Return the auxiliary graph in which adjacency is equation (4.6)."""
    H = Graph()
    H.add_vertices(range(len(candidates)))
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            value = candidates[i].dot_product(R * candidates[j])
            if value == 0 or value == 1:  # equation (4.6)
                H.add_edge(i, j)
    return H


def finite_check(graph_list, triples):
    """Run the final exact clique computation for all nonsingular triples."""
    cases_by_prime = Counter()
    maximum_candidates = 0
    maximum_clique = 0

    for graph_index, p, lam in triples:
        data = admissible_vectors(graph_list[graph_index], p, lam)
        if data is None:
            continue

        R, candidates = data
        if len(candidates) < TARGET:
            continue

        H = compatibility_graph(R, candidates)
        omega = int(H.clique_number(algorithm="Cliquer"))
        if omega >= TARGET:
            raise RuntimeError(
                f"{TARGET}-clique found for graph {graph_index}, "
                f"p={p}, lambda={lam}"
            )

        cases_by_prime[p] += 1
        maximum_candidates = max(maximum_candidates, len(candidates))
        maximum_clique = max(maximum_clique, omega)

    cases_by_prime = dict(sorted(cases_by_prime.items()))
    summary = {
        "auxiliary_graphs": sum(cases_by_prime.values()),
        "cases_by_prime": cases_by_prime,
        "maximum_candidates": maximum_candidates,
        "maximum_clique": maximum_clique,
    }

    print("auxiliary graphs:", summary["auxiliary_graphs"])
    print("cases by characteristic:", summary["cases_by_prime"])
    print("maximum number of admissible vectors:", summary["maximum_candidates"])
    print("maximum clique number:", summary["maximum_clique"])

    return summary


def main():
    graph_list = list(graphs.nauty_geng("8"))
    print("8-vertex graph classes:", len(graph_list))
    check_reported_value("graph_classes", len(graph_list))

    triples, reduction_summary = resultant_reduction(graph_list)
    for name in (
        "relevant_pairs",
        "maximum_Bh",
        "exceptional_primes",
        "largest_exceptional_prime",
        "raw_triples",
        "singular_triples",
        "nonsingular_triples",
    ):
        check_reported_value(name, reduction_summary[name])

    final_summary = finite_check(graph_list, triples)
    check_reported_value("auxiliary_graphs", final_summary["auxiliary_graphs"])
    check_reported_value("cases_by_prime", final_summary["cases_by_prime"])
    check_reported_value("maximum_clique", final_summary["maximum_clique"])

    print("Proposition 4.5 verified.")


if __name__ == "__main__":
    main()
