# SteinerQLinkDancer: an algorithm to compute q-Analog Code for S2[2,3,19]

This repo contains code for generating q-analogs of Steiner systems S2[2,3,19] over the finite field GF(2,19). Full report on the results and the code can be found in the q_analog_report.pdf file.

## Contents

- `main.cpp`: Main driver code 
- `SS2_generator.cpp`: Generates 2D subspaces
- `SS3_generator.cpp`: Generates 3D subspaces
- `kramer_generator.cpp`: Generates Kramer-Mesner matrices
- `dancing_links.cpp`: Implements Dancing Links algorithm
- `matrix.cpp`: Matrix representation for Kramer-Mesner

## Overview

The code generates q-analogs for the Steiner system S2[2,3,19] over GF(2,19) following these steps:

1. Generate all 2D subspaces under the normalizer of the Singer subgroup. Group into orbits.
2. Generate 3D subspaces. Group into orbits containing 7 different 2D orbits.  
3. Build partial Kramer-Mesner matrices with random rows from each 3D orbit.
4. Use Dancing Links to find an exact cover. If no solution, find maximal partial solutions.
5. Combine partial solutions to get largest overall solution.

Key techniques:

- Represent subspaces using generator elements 
- Use orbits to reduce search space
- Probabilistic construction of Kramer-Mesner matrix
- Dancing Links algorithm for exact cover
- Greedy search for maximal partial solutions

## Usage

The main driver code is in `main.cpp`. It will:

- Generate SS2 and SS3 orbits
- Build Kramer-Mesner matrices
- Run Dancing Links solver 
- Find best partial solutions

To use:

```
g++ main.cpp SS2_generator.cpp SS3_generator.cpp kramer_generator.cpp dancing_links.cpp matrix.cpp -o qanalog
./qanalog
```

The code is parallelized using OpenMP. Set OMP_NUM_THREADS environment variable to control number of threads.

## References

Main reference papers:

- Braun, M., Etzion, T., Östergård, P. R., Vardy, A., & Wassermann, A. (2013). Existence of q-analogs of Steiner systems. Forum of Mathematics, Pi, 1.
- Knuth, D. E. (2000). Dancing links. arXiv preprint cs/0011047.

## Authors

- Andrew Elashkin
- Andy Berger
