# ZoomOut: Core Algorithm from Original Paper

## Overview
ZoomOut is a spectral upsampling method for efficient shape correspondence that refines functional maps through iterative two-step procedure.

## Core Mathematical Framework

### Functional Map Definition (Equation 1)
```
C = Φ_M^+ Π Φ_N
```
where:
- `Φ_M^+` is the Moore-Penrose pseudo-inverse of source eigenvectors
- `Π` is the permutation matrix representing pointwise correspondence
- `Φ_N` are the target eigenvectors

**Orthonormal Case**: When eigenfunctions are orthonormal w.r.t. area-weighted inner product (Φ_M^T A_M Φ_M = I):
```
C = Φ_M^T A_M Π Φ_N
```

### Pointwise Map Recovery (Equation 2)
```
T(p) = argmin_q ||C(Φ_N(q))^T - (Φ_M(p))^T||_2  ∀p ∈ M
```
where:
- `T(p)` maps source point `p` to target point `q`
- `Φ_M(p)` is the p-th row of source eigenvector matrix
- `Φ_N(q)` is the q-th row of target eigenvector matrix
- Implemented via nearest-neighbor search in k-dimensional spectral space

## Core Algorithm

### Basic Two-Step Procedure

**Input**: 
- Initial functional map `C_0` of size `k_0 × k_0`
- Source eigenvectors `Φ_M`
- Target eigenvectors `Φ_N`

**Steps**:
1. **Functional to Pointwise**: Convert k×k functional map to pointwise map T via Equation (2)
2. **Pointwise to Functional**: Convert pointwise map T to (k+1)×(k+1) functional map:
   ```
   C_1 = (Φ_M^{k+1})^+ A_M Π Φ_N^{k+1}
   ```

**Iteration**: Repeat steps (1) and (2) to obtain progressively larger maps:
`C_0, C_1, C_2, ..., C_n` until reaching desired size `k_max`

### Detailed Algorithm (From Theoretical Derivation)

**Initialization**:
- Given `k = k_0` and initial `C_0` of size `k_0 × k_0`

**Iterative Steps**:
1. **Step 1**: Compute `argmin_Π ||Π Φ_k^N C_k^T - Φ_k^M||_F^2`
   - This reduces to nearest-neighbor search (Equation 2)
   - Results in permutation matrix Π

2. **Step 2**: Set `k = k + 1` and compute `C_k = (Φ_k^M)^+ Π Φ_k^N`
   - Update functional map to larger dimension

3. **Step 3**: Repeat until `k = k_max`

## Theoretical Foundation

### Optimization Problem (Equation 3)
The algorithm minimizes the energy:
```
E(C) = Σ_k ||C_k^T C_k - I_k||_F^2
```
where `C_k` is the principal k×k submatrix of C, promoting orthonormal principal submatrices.

### Half-Quadratic Splitting
The optimization decouples into two subproblems:
```
min_Π ||(Φ_k^M)^+ Π Φ_k^N C_k^T - I_k||_F^2     (6)
min_C_k ||C_k - (Φ_k^M)^+ Π Φ_k^N||_F^2         (7)
```

### Regularized Formulation (Equation 8)
With regularizer R(Π), problem (6) becomes:
```
min_Π ||Π Φ_k^N C_k^T - Φ_k^M||_F^2
```

## Implementation Details

### Nearest Neighbor Search
- Query: Source spectral coordinates `Φ_M(p)` 
- Data: Transformed target coordinates `C(Φ_N(q))^T`
- Distance: Euclidean L2 norm
- Result: Pointwise correspondence T

### Acceleration Strategies

1. **Larger Step Size**: Instead of k→k+1, use k→k+step for faster upsampling
2. **Sub-sampling**: Perform refinement on subset of points, convert to dense map once
3. **Adaptive Step Size**: Use different step sizes for source vs target dimensions

### Area Weighting
- Eigenfunctions assumed orthonormal w.r.t. area-weighted inner product
- Area matrix A_M used in theoretical formulation: `Φ_M^T A_M Φ_M = I`
- In practice, may use standard orthonormal eigenvectors

## Key Properties

### Convergence
- Method promotes complete isometries through energy in Equation (3)
- Progressive upsampling extracts higher-frequency correspondence information
- Avoids local minima common in fixed-dimension methods like ICP

### Robustness
- Works with small, noisy initial maps (as few as 2×2)
- Robust to deviations from perfect isometry
- Stable w.r.t. mesh resolution and triangulation

### Computational Complexity
- Each iteration: O(n log n) for NN search + O(k^2 n) for functional map update
- Much faster than competing methods (100-500× speedup reported)
- Scalable with shape complexity

## Comparison to Related Methods

### Differences from ICP
- ICP refines in fixed spectral dimension → gets trapped in local minima
- ZoomOut progressively increases dimension → avoids local minima
- ZoomOut promotes full isometries through energy formulation

### Differences from Other Functional Map Methods
- Standard pipeline: compute large basis → optimize → convert to pointwise
- ZoomOut: start small → iteratively upsample → maintains accuracy throughout

### Point-to-Point Recovery
- Uses original functional maps recovery (Equation 2)
- Differs from deblurring methods [Ezuz and Ben-Chen 2017]
- More robust nearest-neighbor based approach