# SVD and Spring-Mass Chain Solver

This project contains two main parts written in Python:

1. A **Singular Value Decomposition (SVD)** routine that computes the decomposition of any matrix using only eigenvalue and eigenvector calculations (no direct SVD library calls).
2. A **spring–mass chain solver** that uses the same SVD function to solve static equilibrium problems of spring systems under different boundary conditions.

The goal of this project is to understand how SVD works internally, see how it can be used to solve systems like \( Ku = f \), and explore what happens when the stiffness matrix becomes singular.

---

## Part A: Custom SVD Function

### Overview
The function `my_svd(A, tol=1e-12, raise_if_singular=True)` calculates the SVD of a general \( N \times M \) matrix using only eigenvalue/vector routines (`np.linalg.eigh`).

It returns:
- `U`: left singular vectors  
- `S`: rectangular diagonal matrix of singular values  
- `V`: right singular vectors  
- `cond2`: 2-norm condition number (`max(s)/min(s)`)  
- `A_inv`: matrix inverse (if square and non-singular)

If the matrix is singular and `raise_if_singular=True`, it raises an error instead of computing an inverse.

### Features
- Builds \( A^T A \) or \( A A^T \) depending on matrix size  
- Sorts eigenvalues in descending order to align with singular values  
- Uses QR factorization to ensure orthogonality of `U` and `V`  
- Computes matrix inverse using \( A^{-1} = V S^{-1} U^T \)  
- Calculates the condition number manually from singular values  

### Comparison
Helper functions are included to compare results to NumPy’s built-in version:
- `check_my_svd_once()` – tests reconstruction and condition number  
- `compare_svd_quick()` – quick side-by-side comparison for random matrices  

---

## Part B: Spring–Mass Chain Solver

### Overview
The second part applies the custom SVD solver to analyze spring–mass systems governed by \( Ku = f \), where:
- \( K \) is the stiffness matrix,  
- \( u \) is the displacement vector,  
- \( f \) is the external force vector.

### Inputs
- Number of masses (`N`)  
- Spring constants (`k_list`)  
- Masses (`masses`)  
- Boundary conditions (`ends`) — can be `'fixed'` or `'free'`  
- Optional external forces, gravity, and cross-sectional area (for stresses)

### Outputs
For each configuration, the solver computes:
- Equilibrium displacements (`u`)  
- Spring elongations  
- Internal spring forces  
- Optional stresses  
- Singular values, eigenvalues, and condition number of \( K \)

### Boundary Conditions Tested
1. **One fixed end, one free end** — system is stable and invertible.  
2. **Two fixed ends** — system is also invertible, with a lower condition number.  
3. **Two free ends** — stiffness matrix \( K \) becomes singular because the entire chain can move without stretching.  
   - If total external force ≠ 0 → no exact equilibrium; the pseudoinverse gives the closest fit.  
   - If total external force = 0 → infinitely many equilibria; the pseudoinverse picks the one with the smallest overall movement.

---

## How to Run

1. Save the file (for example, `svd_spring_chain.py`).
2. Run in terminal:
   ```bash
   python3 svd_spring_chain.py

