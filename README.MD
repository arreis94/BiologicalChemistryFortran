# Fortran subroutines and programs for biological chemical calculations 
 ![Fortran](https://fortran-lang.org/en/_static/fortran-logo-256x256.png)
## Features

### Examples
- In the data folder there are example datas for Water, Hydrogen Peroxide and Methanol molecules.
###  Subroutines
- cmass.f90
    Calculates carbon atom mass.
- crossp.f90
    Cross conjugation for p-orbitals.
- cshift.f90
    C-13 chemical shifts.
- distance.f90
    Calculatets the chemical bond distance's normal vector.
- integral.f90
    Calculates integrals with Simpson's rule and Trapezoidal Rule.
    Note: 
    - Trapezoidal Rule is a rule that evaluates the area under the curves by dividing the total area into smaller trapezoids rather than using rectangles. This integration works by approximating the region under the graph of a function as a trapezoid, and it calculates the area.
    - Simpson's Rule is a numerical method that approximates the value of a definite integral by using quadratic functions. 
- pi.f90
    Calculates pi.
- rotconst.f90
    Calculates the rotational constants.
    `Bv = Be + αe(v + ½)` 

### Calculations on hetero atoms
- Balaban-het.f90
    Balaban index is a convex function in double stars, the Balaban index (also called J index).
- fi_het.f90
    Calculates the Fidelity Index.
- getmdist_het.f90
    Floyd-Warshall algorithm for generating the shortest distance matrix.
- harary_het.f90
    The Harary index is defined as the sum of reciprocals of distances between all pairs of vertices of a connected graph.
- randic_het.f90
    Calculates the Randić (connectivity index) topological index.
- vdegree_het.f90
    Vertex-degree vector.
- wiener_het.f90
    Calculates the sum of the lengths of the shortest paths between all pairs of vertices in the chemical graph representing the non-hydrogen atoms in the molecule (Wiener index).
- xu_het.f90
    Xu index (Ren,1998)