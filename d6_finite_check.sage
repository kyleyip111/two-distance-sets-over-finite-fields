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

N = 8
VERTICES = list(range(N))
ONE_ZZ = vector(ZZ, [1] * N)
BINARY_ZZ = [
    vector(ZZ, [(mask >> i) & 1 for i in range(N)])
    for mask in range(1 << N)
]

Zt = PolynomialRing(ZZ, "t")
t = Zt.gen()
I_Zt = identity_matrix(Zt, N)
ONE_Zt = vector(Zt, ONE_ZZ)


def polynomial_data(G):
    """Return Delta, S, and the 256 pairs (F_b,G_b) in ZZ[t]."""
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
    """Return the nonsingular triples (graph index,p,lambda)."""
    raw_triples = set()
    nonsingular_triples = set()
    exceptional_primes = set()
    relevant_factors = 0
    maximum_Bh = 0

    for graph_index, G in enumerate(graph_list):
        Delta, S, FG = polynomial_data(G)

        # The content divides the leading coefficient -8, so it has no odd
        # prime divisor.  Sage includes nonunit integer content among the
        # factors of a polynomial over ZZ; those degree-zero factors must not
        # be treated as irreducible polynomial factors h.
        for h, _multiplicity in S.factor():
            if h.degree() == 0:
                continue

            # A relevant root is not 0 or -1, and h must not divide Delta over
            # QQ[t].  The resultant test is valid even when h is nonmonic.
            if h(0) == 0 or h(-1) == 0 or h.resultant(Delta) == 0:
                continue

            relevant_factors += 1
            Bh = 0
            prime_count = Counter()

            for F, G_b in FG:
                rF = ZZ(h.resultant(F))
                rG = ZZ(h.resultant(G_b))

                # Since h is irreducible over QQ, both resultants vanish iff
                # h divides both F and G_b over QQ[t].
                if rF == 0 and rG == 0:
                    Bh += 1
                    continue

                r = abs(gcd(rF, rG))
                for p0 in r.prime_divisors():
                    p = int(p0)
                    if p != 2:
                        prime_count[p] += 1

            maximum_Bh = max(maximum_Bh, Bh)
            if Bh >= 20:
                raise RuntimeError(
                    f"|B_h| >= 20 for graph {graph_index} and h = {h}"
                )

            # This is exactly the necessary filter (4.13)
            # [eq:d6-exceptional-filter].
            primes = sorted(
                p for p, count in prime_count.items() if Bh + count >= 20
            )
            exceptional_primes.update(primes)

            for p in primes:
                # Direct evaluation is simpler than factoring h modulo p.
                for lam in range(1, p - 1):  # excludes 0 and -1
                    if h(lam) % p != 0:
                        continue
                    triple = (graph_index, p, lam)
                    raw_triples.add(triple)
                    if Delta(lam) % p != 0:
                        nonsingular_triples.add(triple)

    singular_count = len(raw_triples - nonsingular_triples)
    print("relevant factors:", relevant_factors)
    print("maximum |B_h|:", maximum_Bh)
    print("exceptional odd primes:", len(exceptional_primes))
    print("largest exceptional prime:", max(exceptional_primes))
    print("witness-root triples before nonsingularity:", len(raw_triples))
    print("singular triples discarded:", singular_count)
    print("nonsingular triples:", len(nonsingular_triples))

    return sorted(nonsingular_triples)


def admissible_vectors(G, p, lam_integer):
    """Return R and the vectors satisfying (4.5), (4.7), and (4.8)."""
    Fp = GF(p)
    lam = Fp(lam_integer)
    C = G.adjacency_matrix(vertices=VERTICES, base_ring=Fp)
    B = C - lam * identity_matrix(Fp, N)
    if B.det() == 0:
        return None

    R = B.inverse()
    one = vector(Fp, [1] * N)
    if one.dot_product(R * one) != 0:                 # equation (4.8) [eq:d6-isotropic]
        return None

    candidates = []
    for b0 in BINARY_ZZ:
        b = vector(Fp, b0)
        if b.dot_product(R * b) != -lam:              # equation (4.5) [eq:d6-norm]
            continue
        if b.dot_product(R * one) != 1:               # equation (4.7) [eq:d6-main]
            continue
        candidates.append(b)
    return R, candidates


def compatibility_graph(R, candidates):
    """Return the auxiliary graph whose edges encode (4.6) [eq:d6-compat]."""
    H = Graph()
    H.add_vertices(range(len(candidates)))
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            value = candidates[i].dot_product(R * candidates[j])
            if value == 0 or value == 1:
                H.add_edge(i, j)
    return H


def finite_check(graph_list, triples):
    """Compute every auxiliary clique number using exact Cliquer."""
    cases_by_prime = Counter()
    maximum_candidates = 0
    maximum_clique = 0

    for graph_index, p, lam in triples:
        data = admissible_vectors(graph_list[graph_index], p, lam)
        if data is None:
            continue
        R, candidates = data
        if len(candidates) < 20:
            continue

        H = compatibility_graph(R, candidates)
        omega = int(H.clique_number(algorithm="Cliquer"))
        if omega >= 20:
            raise RuntimeError(
                f"20-clique found for graph {graph_index}, p={p}, lambda={lam}"
            )

        cases_by_prime[p] += 1
        maximum_candidates = max(maximum_candidates, len(candidates))
        maximum_clique = max(maximum_clique, omega)

    cases_by_prime = dict(sorted(cases_by_prime.items()))
    print("auxiliary graphs:", sum(cases_by_prime.values()))
    print("cases by characteristic:", cases_by_prime)
    print("maximum number of admissible vectors:", maximum_candidates)
    print("maximum clique number:", maximum_clique)


def main():
    graph_list = list(graphs.nauty_geng("8"))
    print("8-vertex graph classes:", len(graph_list))

    triples = resultant_reduction(graph_list)
    finite_check(graph_list, triples)
    print("Proposition 4.5 verified.")


if __name__ == "__main__":
    main()
