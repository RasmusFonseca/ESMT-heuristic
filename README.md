# Heuristic for the Euclidean Steiner Minimal Tree in any dimension

The Euclidean Steiner Minimal Tree (ESTP) problem seeks a network of minimal total edge length spanning a set of n terminal points while allowing for the insertion of additional points (Steiner points) to decrease the overall length of the network. 

This software uses a heuristic to find solutions to ESTP for probleminstances of any dimension and almost any size (easily solves for n>10000). For a detailed description of this method, or if used in published research, please see 

* [Euclidean Steiner Tree Heuristic in d-Space](http://dimacs11.cs.princeton.edu/workshop/OlsenLorenzenFonsecaWinter.pdf). A.E. Olsen, S.S. Lorenzen, R. Fonseca and P. Winter.


# Compiling

```
$ cd src
$ make
```

The executable depends on having the qdelaunay executable from [qhull](http://www.qhull.org) in the systems PATH. Theres an easy-to-follow explanation [in the wiki](http://github.com/RasmusFonseca/ESMT-heuristic/wiki/qdelaunay).

# Usage


```
 esmt-heuristic esmt [options] <points>
 esmt-heuristic test esmt [options] <points>
 esmt-heuristic test sgh [options] <points>
 esmt-heuristic counttest <dimension> <no of points> <seed>

 options:
   -v           Verbose
   -npo         Disable post optimisation
   -nsc         Disable sub-graph concatenation
   -sct         Redo concatenation, adding other, non-covered FSTs
   -alg <name>  Sub-graph algorithm:
                  NO  = Numerical optimisation
                  RNO = Restricted numerical optimisation
                  SP  = Simplex partitioning
                  Default NO
   -s           Seed for random generation of point sets.
  # esmt only
   -pt          Print resulting tree
   -val         Validate the resulting tree
  # test only
   -nd          Do not include Delaunay tesselation in time measurement.
   -cs          Collect extra statistics (no of simplices, etc.).
   -out <path>  Output file for test results (CSV format).
   -lt <sec>    Loop time

 points:
   -in <file name> <set name>          Read set <set name>
                                       from file <file name>
   -g <set name> <dim> <no of points>  Generate a point set with
                                       the given number of points in
                                       the given dimension.
                                       Possible set names:
                                           random, sausage, simplex, grid
 points (test only):
   -ina <file name>              Read all sets from file <file name>
   -gn <set name> <dim> <no of points> <no of tests>
                                 Generate the given number of
                                 point set with the given
                                 number of points in the given
                                 dimension.
```

# Examples
```
$ cd src
$ ./esmt-heuristic esmt -in ../data/hard_instance/sausage_n16_d3.stp sausage_n16_d3
Done!
  |MST| = 24.4949
  |SMT| = 20.0392
  Ratio = 0.818097

$ # Time computation of SMT for 10000 random points in 4 dimensions.
$ time ./esmt-heuristic esmt -g random 4 10000  # 10000 random points in 4 dimensions
Done!
  |MST| = 140464
  |SMT| = 132175
  Ratio = 0.94099

real	0m17.459s
user	0m17.244s
sys	0m0.158s

$ # 10 random points in 4 dimensions. Print Steiner tree (-pt)
$ ./esmt-heuristic esmt -g random 4 10 -pt
## RESULT ##
# Terminals
  0 (66.1285, -69.5768, 96.2026, 98.1041)
  1 (-7.09594, 82.8511, 40.5581, 71.8561)
  2 (-86.1443, -75.8665, 9.50698, 31.8577)
  3 (-55.0567, 66.9595, -6.85308, -54.9845)
  4 (71.0011, -4.20561, -15.7684, -45.3275)
  5 (-79.4777, -27.8382, 82.7372, 46.9397)
  6 (29.3678, -99.8046, -43.8857, 42.4799)
  7 (28.3656, 71.4593, 33.4137, -0.593379)
  8 (79.9148, 76.8051, 51.8403, 2.00252)
  9 (-58.4282, 85.9589, -67.618, 59.0912)

# Steiner points
  10 (-58.1107, -54.9113, 38.3886, 40.2187)
  11 (16.5889, -53.6101, 5.33687, 31.3853)
  12 (32.7492, 0.53987, -0.666045, -11.2906)
  13 (4.65199, -56.4612, 28.3315, 44.1951)
  14 (5.21723, 48.6864, 10.5045, -9.33232)
  15 (-9.57395, 76.191, 13.364, 44.8511)
  16 (11.5441, 63.2484, 20.2351, 4.90537)

# Edges
  (7 8)
  (2 10)
  (5 10)
  (6 11)
  (4 12)
  (11 12)
  (0 13)
  (11 13)
  (10 13)
  (3 14)
  (12 14)
  (1 15)
  (9 15)
  (14 16)
  (15 16)
  (7 16)

# |MST|: 990.473
# |SMT|: 917.778
# Steiner ratio: 0.926606
```
