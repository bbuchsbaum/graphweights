# ZoomOut Algorithm: Technical Implementation Details

## Overview

The ZoomOut algorithm performs spectral upsampling for efficient correspondence refinement between two graphs, based on the paper "ZoomOut: Spectral Upsampling for Efficient Shape Correspondence" by Melzi et al. (2019).

## Core Mathematical Framework

### Input
- `W_M`: Symmetric `nM × nM` adjacency matrix of source graph
- `W_N`: Symmetric `nN × nN` adjacency matrix of target graph  
- `C0`: Initial `k0 × k0` functional map
- `k_max`: Maximum spectral dimension to upsample to

### Output
- `C_final`: Refined `k_max × k_max` functional map
- `T_final`: Pointwise correspondence vector of length `nM`
- `Phi_M`, `Phi_N`: Spectral bases (eigenvectors)

## Algorithm Steps

### Phase 1: Spectral Decomposition

#### Laplacian Construction
- **Unnormalized**: `L = D - W` where `D = diag(rowSums(W))`
- **Normalized**: `L_sym = I - D^(-1/2) W D^(-1/2)`

#### Eigendecomposition
- Compute `k_max + 1` smallest eigenpairs using RSpectra::eigs_sym
- Extract first `k_max` non-constant eigenvectors:
  - `Phi_M_full = eigM$vec[, 2:(k_max+1)]`
  - `Phi_N_full = eigN$vec[, 2:(k_max+1)]`

### Phase 2: Iterative Upsampling

The algorithm iteratively refines the functional map through two alternating steps:

#### Step 1: Functional Map to Pointwise Map (Equation 2)

For current spectral dimension `k_cur`:

1. **Spectral Embedding Transform**:
   ```
   Z = Phi_N_k × C_cur^T
   ```
   where `Phi_N_k` are first `k_cur` eigenvectors of target graph

2. **Nearest Neighbor Search**:
   For each source point `p`, find target point `q` that minimizes:
   ```
   T(p) = argmin_q ||Phi_M_k(p,:) - Z(q,:)||_2
   ```
   
   **Implementation**: Use RANN::nn2 with Euclidean distance
   - `searchtype = "standard"` (kd-tree)
   - `eps = 0.0` for exact search, `eps = 0.1` for approximate

#### Step 2: Pointwise Map to Functional Map (Equation 1)

Update functional map to larger dimension `k_next = min(k_cur + step, k_max)`:

```
C_new = Phi_M_next^T × Phi_N_next[T_map, :]
```

where:
- `Phi_M_next`: First `k_next` eigenvectors of source graph
- `Phi_N_next[T_map, :]`: Target eigenvectors indexed by pointwise map
- This assumes orthonormal eigenvectors (avoids pseudoinverse computation)

## Key Implementation Details

### Eigendecomposition Parameters
- Use `ncv = min(max(nM, nN), 2*k_max + 20)` Lanczos vectors
- For normalized Laplacian, implement matrix-vector product function:
  ```
  Lsym_mult(x) = x - inv_sqrt_d .* (W * (inv_sqrt_d .* x))
  ```
- Handle isolated nodes: `inv_sqrt_d[d==0] = 0`

### Nearest Neighbor Search
- **Query**: Source spectral coordinates `Phi_M_k` (nM × k_cur)
- **Data**: Transformed target coordinates `Z` (nN × k_cur) 
- **Distance**: Euclidean L2 norm
- **Result**: Index vector `T_map` of length nM

### Functional Map Update
- **Input validation**: Check `T_map` indices are in range [1, nN]
- **Invalid mappings**: Set to NA, exclude from computation
- **Edge case**: If no valid mappings, return zero matrix
- **Computation**: `C_new = PhiM_valid^T * PhiN_T_valid`

### Error Handling
- **Eigendecomposition failure**: Return empty matrices, issue warnings
- **Insufficient eigenvectors**: Check `ncol(eigvec) >= k_max + 1`
- **NN search failure**: Return NA mappings with warning
- **Invalid correspondences**: Filter out invalid indices

### Reproducibility
- Set `seed` before eigen-decomposition and NN search
- Both RSpectra (ARPACK) and RANN have stochastic components

## Algorithmic Complexity

### Time Complexity
- **Eigendecomposition**: O(k_max × N × sparsity) per graph
- **NN Search per iteration**: O(nM × log(nN) × k_cur) for kd-tree
- **Functional map update**: O(k_next^2 × nM)
- **Total iterations**: ⌈(k_max - k0) / step⌉

### Space Complexity
- **Eigenvector storage**: O(N × k_max) per graph
- **Spectral coordinates**: O(max(nM, nN) × k_max)
- **Adjacency matrices**: O(N × sparsity)

## Convergence Properties

The algorithm refines correspondences by:
1. **Spectral upsampling**: Increasing spectral resolution from k0 to k_max
2. **Correspondence refinement**: Each iteration improves pointwise matching
3. **Functional map consistency**: Maintains spectral domain relationships

The iterative process typically converges to more accurate correspondences as the spectral dimension increases, leveraging higher-frequency spectral information for finer correspondence details.