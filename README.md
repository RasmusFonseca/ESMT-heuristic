# Heuristic for the Euclidean Steiner Minimal Tree in any dimension

The Euclidean Steiner Minimal Tree (ESTP) problem seeks a network of minimal total edge length spanning a set of n terminal points while allowing for the insertion of additional points (Steiner points) to decrease the overall length of the network. 

This software uses a heuristic to find solutions to ESTP for probleminstances of any dimension and almost any size (easily solves for n>10000). For a detailed description of this method, or if used in published research, please see 

* [Euclidean Steiner Tree Heuristic in d-Space](http://dimacs11.cs.princeton.edu/workshop/OlsenLorenzenFonsecaWinter.pdf). A.E. Olsen, S.S. Lorenzen, R. Fonseca and P. Winter.


# Compiling

```
$ cd src
$ make
```


# Usage


```
 esmt-heuristic esmt [options] <points>
 esmt-heuristic test esmt [options] <points>
 esmt-heuristic test sgh [options] <points>

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

# Example
```
$ cd src
$ ./esmt-heuristic esmt -in ../data/hard_instance/sausage_n16_d3.stp sausage_n16_d3
```
