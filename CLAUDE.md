# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

The `neighborweights` (formerly `graphweights` - directory name differs from package name) is an R package for constructing adjacency matrices based on spatial and feature similarity between data points. It provides various methods to create spatial and feature-weighted adjacency matrices for analyzing complex data relationships.

## Package Structure

This is a standard R package with the following key components:

- **R/**: Core R implementation files
- **src/**: C++ implementations using Rcpp/RcppArmadillo for performance-critical functions
- **man/**: Auto-generated documentation via roxygen2
- **tests/testthat/**: Unit tests using the testthat framework
- **DESCRIPTION**: Package metadata, dependencies, and configuration

## Key Architecture Components

### Core Graph Types
- `neighbor_graph`: S3 class wrapping igraph objects or sparse matrices
- `class_graph`: Specialized graph for class-based relationships

### Main Functional Areas

1. **Spatial Weights** (`R/spatial_weights.R`): Functions for spatial adjacency matrices
2. **Neighbor Graph Construction** (`R/neighbor_graph.R`): Core graph object creation and manipulation
3. **KNN Weights** (`R/knn_weights.R`): K-nearest neighbor based adjacency matrices
4. **Diffusion and Kernels** (`R/diffusion.R`): Heat kernel and diffusion-based similarity measures

### C++ Components
- Uses RcppArmadillo for linear algebra operations
- Custom implementations in `src/` for performance-critical algorithms
- Union-Find data structures and Fenwick trees for efficient graph operations

## Development Commands

### Building and Testing
```bash
# Install dependencies and build package
R -e "devtools::install_deps()"
R -e "devtools::build()"

# Run tests
R -e "devtools::test()"

# Generate documentation
R -e "devtools::document()"

# Check package
R -e "devtools::check()"
```

### C++ Compilation
- Uses `src/Makevars` with `PKG_CXXFLAGS += -DARMA_64BIT_WORD`
- Depends on RcppArmadillo for matrix operations
- Compiled objects are in `src/` directory

## Key Dependencies

- **Matrix**: Sparse matrix operations
- **Rcpp/RcppArmadillo**: C++ integration and linear algebra
- **igraph**: Graph data structures and algorithms  
- **nabor, RcppHNSW**: Nearest neighbor search algorithms
- **rflann**: Fast approximate nearest neighbor library
- **assertthat**: Input validation

## Testing Strategy

Tests are located in `tests/testthat/` and focus on:
- Spectral properties preservation in graph operations
- Matrix structure validation (symmetry, non-negativity)
- Equivalence testing between different implementations

## Important Implementation Notes

- All adjacency matrices should be symmetric and non-negative
- The package supports both dense and sparse matrix formats
- C++ implementations often provide significant performance improvements
- Deterministic testing modes available for reproducible results