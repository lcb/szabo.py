# szabo.py
This is a Python re-implementation of the famous [szabo.f](http://www.ccl.net/cca/software/SOURCES/FORTRAN/szabo/index.html) program.
It demonstrates on a simple H-He+ molecule the SCF procedure, which is an important core part of ab-initio quantum chemical computations. For more information I highly recommend to read the book where the original code is from: Modern Quantum Chemistry - Introduction to Advanced Electronic Structure Theory, A. Szabo and N. S. Ostlund, Dover Publications, 1989 ([amazon](http://www.amazon.de/dp/0486691861) [scribd](https://www.scribd.com/book/271592819/Modern-Quantum-Chemistry-Introduction-to-Advanced-Electronic-Structure-Theory)).

## Changes
The basic idea was to stick to the original Fortran code as close as possible. Due to differencies between Fortran and Python, and developments in coding style recommendations, some changes were introduced on purpose, most notably:
* Arrays in Python use zero-based numbering, whereas Fortran uses index 1. To emulate one-based numbering in Python dictionaries were used to represent matrices.
* Eliminated use of goto statements (by using if/then/return constructs).
* Eliminated use of global variables (by modifying method signatures).
* Declaration of local variables was moved right before definition.
* Built-in Fortran method calls (DEXP, DSQRT) were replaced with corresponding Python implementations (math.exp, math.sqrt)
* Added references to equations in the book.
* Added comments explaining complex code fragments.

Newly introduced comments are in lowercase.

## Other re-implementations
of szabo.f on the internet:
* [Julia](https://github.com/SamChill/hartree-fock)
* [Mathematica](http://inside.mines.edu/~mlusk/H2_Roothaan_N_orb.nb)
* [Extended Python version](https://joshuagoings.wordpress.com/2013/04/24/hartree-fock-self-consistent-field-procedure/)