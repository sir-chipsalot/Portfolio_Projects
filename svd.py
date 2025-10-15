import numpy as np

# ==================================
# Part A: SVD I made
# ==================================

def my_svd(A, tol=1e-12, raise_if_singular=True):
    """
    Do SVD of A (N x M) using eigenvalues/eigenvectors only.

    Plan:
      - If N >= M: use B = A^T A (M x M).  B is symmetric.
        eigen(B) gives eigenvalues (s^2) and eigenvectors (columns of V).
        Then U = (A V_i) / s_i.
      - Else (N < M): use C = A A^T to get U first, then V = (A^T U_i) / s_i.

    Returns:
      U (N x N), S (N x M with diag = singular values), V (M x M),
      cond2 (2-norm condition number), and A_inv if square & full rank else None.
    """
    A = np.array(A, dtype=float, copy=True)
    if A.ndim != 2:
        raise ValueError("A must be 2D")
    N, M = A.shape

    if N >= M:
        ATA = A.T @ A
        evals, V = np.linalg.eigh(ATA)              # real, sorted ascending
        order = np.argsort(evals)[::-1]             # want biggest first
        evals = evals[order]
        V = V[:, order]
        s = np.sqrt(np.clip(evals, 0.0, None))

        # make U by pushing V through A and dividing by s
        U = np.zeros((N, N))
        r = min(N, M)
        for i in range(r):
            if s[i] > tol:
                U[:, i] = (A @ V[:, i]) / s[i]
            else:
                U[:, i] = 0.0
        # QR to clean up columns / complete basis
        U, _ = np.linalg.qr(U)
        # polish V too
        V, _ = np.linalg.qr(V)
    else:
        AAT = A @ A.T
        evals, U = np.linalg.eigh(AAT)
        order = np.argsort(evals)[::-1]
        evals = evals[order]
        U = U[:, order]
        s = np.sqrt(np.clip(evals, 0.0, None))
        r = min(N, M)

        V = np.zeros((M, M))
        for i in range(r):
            if s[i] > tol:
                V[:, i] = (A.T @ U[:, i]) / s[i]
            else:
                V[:, i] = 0.0
        V, _ = np.linalg.qr(V)

    # build rectangular S with the singular values on its diagonal
    S = np.zeros((N, M))
    m = min(N, M)
    S[np.arange(m), np.arange(m)] = s[:m]

    # condition number: biggest sigma / smallest non-tiny sigma
    nz = s[s > tol]
    cond2 = float('inf') if nz.size == 0 else float(nz.max() / nz.min())

    # try to build inverse if square and full rank
    A_inv = None
    if N == M:
        full_rank = (nz.size == m) and (nz.min() > tol)
        if full_rank:
            S_inv = np.zeros((M, N))
            S_inv[np.arange(m), np.arange(m)] = 1.0 / s[:m]
            A_inv = V @ S_inv @ U.T
        else:
            if raise_if_singular:
                raise np.linalg.LinAlgError("Matrix looks singular (tiny sigma); no inverse.")

    return U, S, V, cond2, A_inv


def check_my_svd_once():
    """Run SVD on one random matrix and print a few sanity checks."""
    np.random.seed(42)
    A = np.random.rand(4, 4)

    print(f"The A matrix is: {A}\n")

    U, S, V, cond2, A_inv = my_svd(A)

    A_back = U @ S @ V.T
    print("rebuild error:", np.linalg.norm(A - A_back))

    s_np = np.linalg.svd(A, compute_uv=False)
    s_me = np.diag(S)[:min(A.shape)]
    print("my sigmas:", np.round(s_me, 6))
    print("np sigmas:", np.round(np.sort(s_np)[::-1], 6))
    print("cond2 mine:", cond2, " vs np:", float(s_np.max()/s_np.min()))


def compare_svd_quick():
    """Another small comparison against numpy — prints numbers only."""
    np.random.seed(0)
    A = np.random.randn(5, 5) + np.eye(5)
    U, S, V, cond2, Ainv = my_svd(A)
    s_np = np.linalg.svd(A, compute_uv=False)
    print("my singular values:", np.round(np.diag(S)[:min(A.shape)], 8))
    print("np singular values:", np.round(np.sort(s_np)[::-1], 8))
    print("my cond2:", cond2, " np cond2:", float(s_np.max()/s_np.min()))
    print("||A - USV^T||:", np.linalg.norm(A - U @ S @ V.T))
    print("||Ainv - inv(A)||:", np.linalg.norm(Ainv - np.linalg.inv(A)))


def make_pinv_with_my_svd(A, tol=1e-12):
    """Build A^+ (pseudoinverse) using my_svd only."""
    U, S, V, cond2, _ = my_svd(A, tol=tol, raise_if_singular=False)
    N, M = A.shape
    m = min(N, M)
    s = np.diag(S)[:m]

    S_plus = np.zeros((M, N))  # size M x N so that A^+ = V S^+ U^T
    for i in range(m):
        if s[i] > tol:
            S_plus[i, i] = 1.0 / s[i]
    A_pinv = V @ S_plus @ U.T
    return A_pinv, cond2, s


# ==============================
# Part B: little spring-chain toy
# ==============================

def make_k_matrix(N, k_list, ends=("fixed", "free")):
    """
    Build stiffness matrix K for a straight chain of N masses.
    ends can be ("fixed"|"free", "fixed"|"free").

    If an end is "fixed" we add a spring to the wall there.
    Between each neighbor mass we add an internal spring.
    """
    left, right = ends
    K = np.zeros((N, N), float)
    springs = []  # list of (a, b, k); node -1 means ground (u=0)
    idx = 0

    if left == "fixed":
        kL = float(k_list[idx]); idx += 1
        springs.append((-1, 0, kL))
        K[0, 0] += kL

    for i in range(N-1):
        k = float(k_list[idx]); idx += 1
        springs.append((i, i+1, k))
        K[i, i]     += k
        K[i+1, i+1] += k
        K[i, i+1]   -= k
        K[i+1, i]   -= k

    if right == "fixed":
        kR = float(k_list[idx]); idx += 1
        springs.append((N-1, -1, kR))
        K[N-1, N-1] += kR

    expected = (1 if left == "fixed" else 0) + (N-1) + (1 if right == "fixed" else 0)
    if len(k_list) != expected:
        raise ValueError(f"k_list needs {expected} values for ends={ends}, got {len(k_list)}")

    return K, springs


def solve_springs_with_my_svd(N, k_list, masses, ends=("fixed","free"), f=None, gravity=None, area=None, tol=1e-12):
    """
    Solve K u = f using my SVD (via pseudoinverse). Works even if K is singular.

    - masses: used only if gravity is given; we add masses*gravity to f
    - area:   if provided, we compute stress = spring_force / area
    """
    masses = np.array(masses, float).reshape(-1)
    if masses.size != N:
        raise ValueError("masses must have length N")

    if f is None:
        f = np.zeros(N)
    f = np.array(f, float).reshape(-1)
    if f.size != N:
        raise ValueError("f must have length N")

    if gravity is not None:
        f = f + masses * float(gravity)

    K, springs = make_k_matrix(N, k_list, ends=ends)
    print("Stiffness matrix K:\n", K)
    print("f (after gravity if used):\n", f)

    K_pinv, cond2, svals = make_pinv_with_my_svd(K, tol=tol)
    u = K_pinv @ f

    # elongations per spring (u_b - u_a). ground (-1) has u=0
    elong = []
    for a, b, k in springs:
        ua = 0.0 if a == -1 else u[a]
        ub = 0.0 if b == -1 else u[b]
        elong.append(ub - ua)
    elong = np.array(elong)

    spring_forces = np.array([k * e for (_, _, k), e in zip(springs, elong)])
    spring_stress = None if area is None else spring_forces / float(area)

    # print a bit of spectral info so I can see conditioning
    print("\nSVD info for K:\n")

    print("singular values (desc):", np.round(np.sort(svals)[::-1], 8))
    # also show eigenvalues (nice for a symmetric K)
    evals = np.linalg.eigvals(K)
    evals_sorted = np.sort(np.real_if_close(evals))[::-1]
    print("eigenvalues (desc):", np.round(evals_sorted, 8))
    print("cond2(K):", cond2)

    note = ""
    if ends == ("free", "free"):
        note = (
        "When both ends are free, the whole chain can slide left or right without stretching any springs."
        "This means the stiffness matrix (K) can't be inverted — it's singular. "
        "If the total external force isn't zero, the system can't balance exactly, "
        "so the pseudoinverse finds the closest possible fit. "
        "If the total external force is zero, there are infinite positions that work, "
        "and the pseudoinverse picks the one with the smallest overall movement. "
        "\n"
        )

    return {
        "K": K,
        "u": u,
        "springs": springs,
        "elongations": elong,
        "spring_forces": spring_forces,
        "spring_stress": spring_stress,
        "cond2": cond2,
        "singular_values": svals,
        "note": note,
    }


# ===================
# Demos / try-it runs
# ===================

def show_all_spring_setups():
    print("SPRING-MASS CHAIN USING MY SVD")
    print("-" * 40)

    print("one fixed, one free")
    N = 3
    out = solve_springs_with_my_svd(
        N, k_list=[100.0, 200.0, 300.0],  # k_left, k01, k12
        masses=[1.0, 1.0, 1.0], ends=("fixed","free"), f=[0.0, 10.0, 0.0])
    print("displacements:", np.round(out["u"], 6))
    print("spring forces:", np.round(out["spring_forces"], 6), '\n')

    print("both ends fixed")
    out = solve_springs_with_my_svd(
        N, k_list=[100.0, 200.0, 300.0, 400.0],  # k_left, k01, k12, k_right
        masses=[1.0, 1.0, 1.0], ends=("fixed","fixed"), f=[0.0, 10.0, 0.0])
    print("displacements:", np.round(out["u"], 6))
    print("spring forces:", np.round(out["spring_forces"], 6), '\n')

    print("both ends free (singular K)")
    out = solve_springs_with_my_svd(
        N, k_list=[200.0, 300.0],  # k01, k12
        masses=[1.0, 1.0, 1.0], ends=("free","free"), f=[5.0, -5.0, 0.0])
    print("displacements (min-norm):", np.round(out["u"], 6))
    print("spring forces:", np.round(out["spring_forces"], 6), '\n')
    if out["note"]:
        print("Discussion on finindings:\n" + out["note"], '\n')


def main():
    print("PART A: my_svd quick checks")
    print("-" * 40, "\n")

    print("My checks against np:\n\n")
    check_my_svd_once()

    print("\nmore checks vs numpy")
    print("-" * 40, "\n")
    compare_svd_quick()

    print("\n\nPART B: spring chain runs")
    print("-" * 40, "\n")
    show_all_spring_setups()


if __name__ == "__main__":
    main()
