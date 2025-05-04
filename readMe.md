# Spectral Clustering - Final Project

## Overview
This project implements **Unnormalized Spectral Clustering** using both Python and C.  
It includes:
- Construction of the Weighted Adjacency Matrix (WAM)
- Diagonal Degree Matrix (DDG)
- Graph Laplacian (GL)
- Jacobi eigenvalue algorithm
- Full spectral clustering with K-means++ initialization.

The project combines C for computation-intensive operations and Python for user interaction.

---

## Files
- `spkmeans.py`: Python interface to run clustering tasks and communicate with C extensions.
- `spkmeans.c`: C implementation of the core algorithms (WAM, DDG, GL, Jacobi).
- `spkmeansmodule.c`: Python C-API wrapper to expose C functions to Python.
- `spkmeans.h`: Header file for function declarations.
- `setup.py`: Python setup file for building the C extension.
- `Makefile`: Script to compile the C executable (`spkmeans`).

---

## Installation

### Build C executable
make

### Build Python C Extension
python3 setup.py build_ext --inplace

---

## Usage

### Python Script (`spkmeans.py`)
```bash
python3 spkmeans.py [k] goal input_file.txt
- `k`: (optional) number of clusters. If omitted, it will be chosen using the eigengap heuristic.
- `goal`: One of the following:
  - `spk`: Full spectral clustering
  - `wam`: Compute Weighted Adjacency Matrix
  - `ddg`: Compute Diagonal Degree Matrix
  - `gl`: Compute Graph Laplacian
  - `jacobi`: Compute eigenvalues and eigenvectors
- `input_file.txt`: Input data file.

#### Example
python3 spkmeans.py 3 spk input_1.txt

---

### C Executable (`spkmeans`)
./spkmeans goal input_file.txt
- `goal`: One of `wam`, `ddg`, `gl`, or `jacobi`.
- `input_file.txt`: Input data file.

#### Example
./spkmeans gl input_1.txt

---

## Output Format
- **spk**:
  - First line: Indices of initial centroids.
  - Following lines: Final centroids, each on a new line.
- **jacobi**:
  - First line: Eigenvalues.
  - Following lines: Corresponding eigenvectors (columns).
- **wam/ddg/gl**:
  - Matrix rows separated by newlines, values comma-separated.

*All outputs are formatted to 4 decimal places.*

---

## Assumptions
- All inputs are valid and correctly formatted.
- No error handling for invalid arguments (except general error print).
- K-means uses `ε = 0` and `max_iter = 300`.
- Random seed is fixed at `np.random.seed(0)` for reproducibility.
- Floating-point anomalies (e.g., -0.0000) are handled gracefully.

---

## Error Handling
If any error occurs (e.g., memory allocation issues), the program will print:
An Error Has Occurred

and terminate.

---

## References
- Ulrike Von Luxburg, *A tutorial on spectral clustering*, *Statistics and Computing*, 17(4):395–416, 2007.


