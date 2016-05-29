nfiles-ndisks is a program that I developed during my BSc, while following a course related with Intelligent Systems. It uses a genetic algorithm while trying to find the optimum solution in a constraint-based problem. Two different solutions are provided. The first one concerns fitting a random number of files with a random size each, in a hard disk (which also has a random size). The second problem is similar to the first, but uses a random number of disks instead of just a single one.

*Warning:* Do not focus on the C++ code of nfiles-ndisks because it needs serious refactoring. Instead, I'd suggest you to focus on the genetic algorithm part and how GAlib (and other similar genetic algorithm libraries) can be used for solving constraint-based problems.

The following Figure shows a successful execution (meaning that a solution was found) of nfiles-ndisks for the N disk problem (the results are written by default in the file `output.dat`).

![example](http://i67.tinypic.com/f2klmt.png)
