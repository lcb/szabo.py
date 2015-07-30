#!/usr/bin/env python
import math
import sys

"""
MINIMAL BASIS STO-3G CALCULATION ON HEH+

THIS IS A LITTLE DUMMY MAIN PROGRAM WHICH CALLS HFCALC

APPENDIX B: TWO-ELECTRON SELF-CONSISTENT-FIELD PROGRAM
OF MODERN QUANTUM CHEMISTRY by
Attila Szabo and Neil S. Ostlund
Ed. 2nd (1989) Dover Publications INC.

Labourly Typed by Michael Zitolo (Feb., 2005)
Edited and Compiled by Michael Zitolo and Xihua Chen

Cleaned up and debugged again by Andrew Long (2012) 
                  and Daniele (kalium) Dondi (2013)
Converted to Python by Lorenz Blum (2015)
"""

# Global constants/functions
PI = 3.1415926535898
S12 = 0.0
T11 = 0.0
T12 = 0.0
T22 = 0.0
V11A = 0.0
V12A = 0.0
V22A = 0.0
V11B = 0.0
V12B = 0.0
V22B = 0.0
V1111 = 0.0
V2111 = 0.0
V2121 = 0.0
V2211 = 0.0
V2221 = 0.0
V2222 = 0.0
H = {1: {}, 2: {}}
TT = {1: {1: {1: {}, 2: {}},
          2: {1: {}, 2: {}}, },
      2: {1: {1: {}, 2: {}},
          2: {1: {}, 2: {}}, }, }
G = {1: {}, 2: {}}
P = {1: {}, 2: {}}
F = {1: {}, 2: {}}
X = {1: {}, 2: {}}
s = {1: {}, 2: {}}
XT = {1: {}, 2: {}}
# ##########################

def MAIN():
    IOP = 2
    N = 3
    R = 1.4632
    ZETA1 = 2.0925
    ZETA2 = 1.24
    ZA = 2.0
    ZB = 1.0
    HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB)


def HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB):
    """
    DOES A HARTREE-FOCK CALCULATION FOR A TWO-ELECTRON DIATOMIC
    USING THE 1S MINIMAL STO-NG BASIS SET
    MINIMAL BASIS SET HAS BASIS FUNCTIONS 1 AND 2 ON NUCLEI A AND B

    IOP=0 NO PRINTING WHATSOEVER (TO OPTIMIZE EXPONENTS, SAY)
    IOP=1 PRINT ONLY CONVERGED RESULTS
    IOP=2 PRINT EVERY ITERATION
    N STO-NG CALCULATION (N=1,2 OR 3)
    R BONDLENGTH (AU)
    ZETA1 SLATER ORBITAL EXPONENT (FUNCTION 1)
    ZETA2 SLATER ORBITAL EXPONENT (FUNCTION 2)
    ZA ATOMIC NUMBER (ATOM A)
    ZB ATOMIC NUMBER (ATOM B)
    """
    if IOP > 0:
        print "STO-%sG FOR ATOMIC NUMBERS %d and %d" % (N, ZA, ZB)
    # CALCULATE ALL THE ONE AND TWO ELECTRON INTEGRALS
    INTGRL(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
    # BE INEFFICIENT AND PUT ALL INTEGRALS IN PRETTY ARRAYS
    COLECT(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
    # PERFORM THE SCF CALCULATION
    SCF(IOP, N, R, ZETA1, ZETA2, ZA, ZB)


def INTGRL(IOP, N, R, ZETA1, ZETA2, ZA, ZB):
    """
    CALCULATES ALL THE BASIC INTEGRALS NEEDED FOR SCF CALCULATION
    """
    A1 = [0] * N
    A2 = [0] * N
    D1 = [0] * N
    D2 = [0] * N
    # THESE ARE THE CONTRACTION COEFFICIENTS AND EXPONENTS FOR
    # A NORMALIZED SLATER ORBITAL WITH EXPONENT 1.0 IN TERMS OF
    # NORMALIZED 1S PRIMITIVE GAUSSIANS
    # COEF 1st dim: i-th gauss coefficient , 2nd dim: basis set (STO-NG)
    COEF = {0: {1: 1.0, 2: 0.678914, 3: 0.444635},
            1: {1: None, 2: 0.430129, 3: 0.535328},
            2: {1: None, 2: None, 3: 0.154329}}
    EXPON = {0: {1: 0.270950, 2: 0.151623, 3: 0.109818},
             1: {1: None, 2: 0.851819, 3: 0.405771},
             2: {1: None, 2: None, 3: 2.22766}}
    R2 = R * R
    # SCALE THE EXPONENTS (A) OF PRIMITIVE GAUSSIANS
    # INCLUDE NORMALIZATION IN CONTRACTION COEFFICIENTS (D)
    for I in range(N):
        A1[I] = EXPON[I][N] * (ZETA1 ** 2)
        D1[I] = COEF[I][N] * ((2.0 * A1[I] / PI) ** 0.75)
        A2[I] = EXPON[I][N] * (ZETA2 ** 2)
        D2[I] = COEF[I][N] * ((2.0 * A2[I] / PI) ** 0.75)
    # D AND A ARE NOW THE CONTRACTION COEFFICIENTS AND EXPONENTS
    # IN TERMS OF UNNORMALIZED PRIMITIVE GAUSSIANS
    global S12, T11, T12, T22, V11A, V12A, V22A, V11B, V12B, V22B, V1111, V2111, V2121, V2211, V2221, V2222
    # CALCULATE ONE-ELECTRON INTEGRALS
    # CENTER A IS FIRST ATOM, CETER B IS SECOND ATOM
    # ORIGIN IS ON CENTER A
    # V12A = OFF-DIAGONAL NUCLEAR ATTRACTION TO CENTER A, ETC.
    for I in range(N):
        for J in range(N):
            # RAP2 = SQUARED DISTANCE BETWEEN CENTER A AND CENTER P, ETC.
            RAP = A2[J] * R / (A1[I] + A2[J])
            RAP2 = RAP ** 2
            RBP2 = (R - RAP) ** 2
            S12 += S(A1[I], A2[J], R2) * D1[I] * D2[J]
            T11 += T(A1[I], A1[J], 0.0) * D1[I] * D1[J]
            T12 += T(A1[I], A2[J], R2) * D1[I] * D2[J]
            T22 += T(A2[I], A2[J], 0.0) * D2[I] * D2[J]
            V11A += V(A1[I], A1[J], 0.0, 0.0, ZA) * D1[I] * D1[J]
            V12A += V(A1[I], A2[J], R2, RAP2, ZA) * D1[I] * D2[J]
            V22A += V(A2[I], A2[J], 0.0, R2, ZA) * D2[I] * D2[J]
            V11B += V(A1[I], A1[J], 0.0, R2, ZB) * D1[I] * D1[J]
            V12B += V(A1[I], A2[J], R2, RBP2, ZB) * D1[I] * D2[J]
            V22B += V(A2[I], A2[J], 0.0, 0.0, ZB) * D2[I] * D2[J]

    for I in range(N):
        for J in range(N):
            for K in range(N):
                for L in range(N):
                    # CALCULATE TWO-ELECTRON INTEGRALS
                    RAP = A2[I] * R / (A2[I] + A1[J])
                    RBP = R - RAP
                    RAQ = A2[K] * R / (A2[K] + A1[L])
                    RBQ = R - RAQ
                    RPQ = RAP - RAQ
                    RAP2 = RAP * RAP
                    RBP2 = RBP * RBP
                    RAQ2 = RAQ * RAQ
                    RBQ2 = RBQ * RBQ
                    RPQ2 = RPQ * RPQ
                    V1111 += TWOE(A1[I], A1[J], A1[K], A1[L], 0.0, 0.0, 0.0) * D1[I] * D1[J] * D1[K] * D1[L]
                    V2111 += TWOE(A2[I], A1[J], A1[K], A1[L], R2, 0.0, RAP2) * D2[I] * D1[J] * D1[K] * D1[L]
                    V2121 += TWOE(A2[I], A1[J], A2[K], A1[L], R2, R2, RPQ2) * D2[I] * D1[J] * D2[K] * D1[L]
                    V2211 += TWOE(A2[I], A2[J], A1[K], A1[L], 0.0, 0.0, R2) * D2[I] * D2[J] * D1[K] * D1[L]
                    V2221 += TWOE(A2[I], A2[J], A2[K], A1[L], 0.0, R2, RBQ2) * D2[I] * D2[J] * D2[K] * D1[L]
                    V2222 += TWOE(A2[I], A2[J], A2[K], A2[L], 0.0, 0.0, 0.0) * D2[I] * D2[J] * D2[K] * D2[L]

    if IOP == 0:
        return

    print "R\tZETA1\tZETA2\tS12\tT11"
    print R, "\t", ZETA1, "\t", ZETA2, "\t", S12, "\t", T11
    print
    print "T12\tT22\tV11A\tV12A\tV22A"
    print T12, "\t", T22, "\t", V11A, "\t", V12A, "\t", V22A
    print
    print "V11B\tV12B\tV22B\tV1111\tV2111"
    print V11B, "\t", V12B, "\t", V22B, "\t", V1111, "\t", V2111
    print
    print "V2121\tV2211\tV2221\tV2222"
    print V2121, "\t", V2211, "\t", V2221, "\t", V2222
    print


def F0(ARG):
    """
    CALCULATES THE F FUNCTION
    FO ONLY (S-TYPE ORBITALS)
    Equation (A.32)
    """
    if ARG < 1.0E-6:
        # ASYMPTOTIC VALUE FOR SMALL ARGUMENTS
        return 1.0 - ARG / 3.0
    # F0 IN TERMS OF THE ERROR FUNCTION
    F0 = math.sqrt(PI / ARG) * DERFOTHER(math.sqrt(ARG)) / 2.0
    return F0


def DERFOTHER(ARG):
    """
    CALCULATES THE ERROR FUNCTION ACCORDING TO A RATIONAL
    APPROXIMATION FROM M. ARBRAMOWITZ AND I.A. STEGUN,
    HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER.
    ABSOLUTE ERROR IS LESS THAN 1.5*10**(-7)
    CAN BE REPLACED BY A BUILT-IN FUNCTION ON SOME MACHINES
    """
    A = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429]
    P = 0.327591
    T = 1.0 / (1.0 + P * ARG)
    TN = T
    POLY = A[0] * TN
    for I in [1, 2, 3, 4]:
        TN *= T
        POLY += A[I] * TN
    DERFOTHER = 1.0 - POLY * math.exp(-ARG * ARG)
    return DERFOTHER


def S(A, B, RAB2):
    """
    CALCULATES OVERLAPS FOR UN-NORMALIZED PRIMITIVES
    """
    S = (PI / (A + B)) ** 1.5 * math.exp(-A * B * RAB2 / (A + B))
    return S


def T(A, B, RAB2):
    """
    Kinetic energy integral A.10
    CALCULATES KINETIC ENERGY INTEGRALS FOR UN-NORMALIZED PRIMITIVES
    Equation (A.10)
    """
    T = A * B / (A + B) * (3.0 - 2.0 * A * B * RAB2 / (A + B)) * (PI / (A + B)) ** 1.5 * math.exp(
        -A * B * RAB2 / (A + B))
    return T


def V(A, B, RAB2, RCP2, ZC):
    """
    CALCULATES UN-NORMALIZED NUCLEAR ATTRACTION INTEGRALS
    Equation (A.33)
    """
    V = 2.0 * PI / (A + B) * F0((A + B) * RCP2) * math.exp(-A * B * RAB2 / (A + B))
    V = -V * ZC
    return V


def TWOE(A, B, C, D, RAB2, RCD2, RPQ2):
    """
    CALCULATES TWO-ELECTRON INTEGRALS FOR UN-NORMALIZED PRIMITIVES
    A,B,C,D ARE THE EXPONENTS ALPHA, BETA, ETC.
    RAB2 EQUALS SQUARED DISTANCE BETWEEN CENTER A AND CENTER B, ETC.
    Equation (A.41)
    """
    return 2.0 * (PI ** 2.5 ) / ((A + B) * (C + D) * math.sqrt(A + B + C + D)) * F0(
        (A + B) * (C + D) * RPQ2 / (A + B + C + D)) * math.exp(-A * B * RAB2 / (A + B) - C * D * RCD2 / (C + D))


def COLECT(IOP, N, R, ZETA1, ZETA2, ZA, ZB):
    """
    THIS TAKES THE BASIC INTEGRALS FROM COMMON AND ASSEMBLES THE
    RELEVENT MATRICES, THAT IS S,H,X,XT, AND TWO-ELECTRON INTEGRALS
    """
    # FORM CORE HAMILTONIAN
    global H
    H[1][1] = T11 + V11A + V11B
    H[1][2] = T12 + V12A + V12B
    H[2][1] = H[1][2]
    H[2][2] = T22 + V22A + V22B
    # FORM OVERLAP MATRIX
    global s
    s[1][1] = 1.0
    s[1][2] = S12
    s[2][1] = s[1][2]
    s[2][2] = 1.0
    # USE CANONICAL ORTHOGONALIZATION
    global X
    X[1][1] = 1.0 / math.sqrt(2.0 * (1.0 + S12))
    X[2][1] = X[1][1]
    X[1][2] = 1.0 / math.sqrt(2.0 * (1.0 - S12))
    X[2][2] = -X[1][2]
    # TRANSPOSE OF TRANSFORMATION MATRIX
    global XT
    XT[1][1] = X[1][1]
    XT[1][2] = X[2][1]
    XT[2][1] = X[1][2]
    XT[2][2] = X[2][2]
    # MATRIX OF TWO-ELECTRON INTEGRALS
    global TT
    TT[1][1][1][1] = V1111
    TT[2][1][1][1] = V2111
    TT[1][2][1][1] = V2111
    TT[1][1][2][1] = V2111
    TT[1][1][1][2] = V2111
    TT[2][1][2][1] = V2121
    TT[1][2][2][1] = V2121
    TT[2][1][1][2] = V2121
    TT[1][2][1][2] = V2121
    TT[2][2][1][1] = V2211
    TT[1][1][2][2] = V2211
    TT[2][2][2][1] = V2221
    TT[2][2][1][2] = V2221
    TT[2][1][2][2] = V2221
    TT[1][2][2][2] = V2221
    TT[2][2][2][2] = V2222
    if IOP == 0:
        return
    MATOUT(s, 2, 2, 2, 2, "S")
    MATOUT(X, 2, 2, 2, 2, "X")
    MATOUT(H, 2, 2, 2, 2, "H")
    for I in [1, 2]:
        for J in [1, 2]:
            for K in [1, 2]:
                for L in [1, 2]:
                    print "(", I, J, K, L, ")", TT[I][J][K][L]


def SCF(IOP, N, R, ZETA1, ZETA2, ZA, ZB):
    """
    PERFORMS THE SCF ITERATIONS
    """
    # CONVERGENCE CRITERION FOR DENSITY MATRIX
    CRIT = 1.0E-4
    # MAXIMUM NUMBER OF ITERATIONS
    MAXIT = 25
    # ITERATION NUMBER
    ITER = 0
    # USE CORE-HAMILTONIAN FOR INITIAL GUESS AT F, I.E. (P=0)
    global G, P
    for I in [1, 2]:
        for J in [1, 2]:
            P[I][J] = 0.0

    if IOP == 2:
        MATOUT(P, 2, 2, 2, 2, "P")

    while True:
        ITER += 1
        if IOP == 2:
            print "START OF ITERATION NUMBER =", ITER

        # FORM TWO-ELECTRON PART OF FOCK MATRIX FROM P
        FORMG()
        if IOP == 2:
            MATOUT(G, 2, 2, 2, 2, "G")

        # ADD CORE HAMILTONIAN TO GET FOCK MATRIX
        global F
        for I in [1, 2]:
            for J in [1, 2]:
                F[I][J] = H[I][J] + G[I][J]

        # CALCULATE ELECTRONIC ENERGY
        EN = 0.0
        for I in [1, 2]:
            for J in [1, 2]:
                EN += 0.5 * P[I][J] * (H[I][J] + F[I][J])

        if IOP == 2:
            MATOUT(F, 2, 2, 2, 2, "F")
            print "ELECTRONIC ENERGY = ", EN

        # TRANSFORM FOCK MATRIX USING G FOR TEMPORARY STORAGE
        global X, XT
        G = MULT(F, X)
        FPRIME = MULT(XT, G)
        # DIAGONALIZE TRANSFORMED FOCK MATRIX
        CPRIME, E = DIAG(FPRIME)
        # TRANSFORM EIGENVECTORS TO GET MATRIX C
        C = MULT(X,CPRIME)
        # FORM NEW DENSITY MATRIX
        OLDP = {1:{},2:{}}
        for I in [1,2]:
            for J in [1,2]:
                # SAVE PRESENT DENSITY MATRIX
                # BEFORE CREATING NEW ONE
                OLDP[I][J]=P[I][J]
                P[I][J]=0.0
                for K in [1]:
                    P[I][J]=P[I][J]+2.0 *C[I][K]*C[J][K]


        if IOP == 2:
            MATOUT(FPRIME,2,2,2,2,"F'  ")
            MATOUT(CPRIME,2,2,2,2,"C'  ")
            MATOUT(E,2,2,2,2,'E   ')
            MATOUT(C,2,2,2,2,'C   ')
            MATOUT(P,2,2,2,2,'P   ')
        # CALCULATE DELTA
        DELTA=0.0
        for I in [1,2]:
            for J in [1,2]:
                DELTA+=(P[I][J]-OLDP[I][J])**2

        DELTA=math.sqrt(DELTA/4.0 )
        if IOP > 0:
            print "DELTA(CONVERGENCE OF DENSITY MATRIX) =",DELTA
        # CHECK FOR CONVERGENCE
        if DELTA <= CRIT:
            break

        # NOT YET CONVERGED
        # TEST FOR MAXIMUM NUMBER OF ITERATIONS
        # IF MAXIMUM NUMBER NOT YET REACHED
        # GO BACK FOR ANOTHER ITERATION
        if ITER >= MAXIT:
            break

    # CALCULATION CONVERGED IF IT GOT HERE
    # ADD NUCLEAR REPULSION TO GET TOTAL ENERGY
    ENT=EN+ZA*ZB/R
    if IOP > 0:
        print "CALCULATION CONVERGED"
        print "ELECTRONIC ENERGY =",EN
        print "TOTAL ENERGY =",ENT

    if IOP == 1:
        # PRINT OUT THE FINAL RESULTS IF
        # HAVE NOT DONE SO ALREADY
        MATOUT(G,2,2,2,2,"G"   )
        MATOUT(F,2,2,2,2,"F"   )
        MATOUT(E,2,2,2,2,"E"   )
        MATOUT(C,2,2,2,2,"C"   )
        MATOUT(P,2,2,2,2,"P"   )
    # PS MATRIX HAS MULLIKEN POPULATIONS
    global S
    OLDP = MULT(P,s)
    if IOP > 0:
      MATOUT(OLDP,2,2,2,2,"PS" )

def FORMG():
    """
    CALCULATES THE G MATRIX FROM THE DENSITY MATRIX
    AND TWO-ELECTRON INTEGRALS
    """
    global G, P, TT
    for I in [1, 2]:
        for J in [1, 2]:
            G[I][J] = 0.0
            for K in [1, 2]:
                for L in [1, 2]:
                    G[I][J] += P[K][L] * (TT[I][J][K][L] - 0.5 * TT[I][L][K][J])


def DIAG(F):
    """
    DIAGONALIZES F TO GIVE EIGENVECTORS IN C AND EIGENVALUES IN E
    THETA IS THE ANGLE DESCRIBING SOLUTION
    """
    if (abs(F[1][1] - F[2][2]) > 1.0E-20):
        # SOLUTION FOR HETERONUCLEAR DIATOMIC
        THETA = 0.5 * math.atan(2.0 * F[1][2] / (F[1][1] - F[2][2]))
    else:
        # HERE IS SYMMETRY DETERMINED SOLUTION (HOMONUCLEAR DIATOMIC)
        THETA = PI / 4.0


    C = {1: {}, 2: {}}
    E = {1: {}, 2: {}}
    C[1][1] = math.cos(THETA)
    C[2][1] = math.sin(THETA)
    C[1][2] = math.sin(THETA)
    C[2][2] = -math.cos(THETA)
    E[1][1] = F[1][1] * math.cos(THETA) ** 2 + F[2][2] * math.sin(THETA) ** 2 + F[1][2] * math.sin(2.0 * THETA)
    E[1][1] = F[1][1] * math.cos(THETA) ** 2 + F[2][2] * math.sin(THETA) ** 2 + F[1][2] * math.sin(2.0 * THETA)
    E[2][2] = F[2][2] * math.cos(THETA) ** 2 + F[1][1] * math.sin(THETA) ** 2 - F[1][2] * math.sin(2.0 * THETA)
    E[2][1] = 0.0
    E[1][2] = 0.0
    # ORDER EIGENVALUES AND EIGENVECTORS
    if ( E[2][2] > E[1][1] ):
        return C,E
    TEMP = E[2][2]
    E[2][2] = E[1][1]
    E[1][1] = TEMP

    TEMP = C[1][2]
    C[1][2] = C[1][1]
    C[1][1] = TEMP

    TEMP = C[2][2]
    C[2][2] = C[2][1]
    C[2][1] = TEMP
    return C, E


def MULT(A, B):
    """
    MULTIPLIES TWO SQUARE MATRICES A AND B TO GET C
    """
    C = {1: {}, 2: {}}
    for I in range(1, len(A) + 1):
        for J in range(1, len(A) + 1):
            C[I][J] = 0.0
            for K in range(1, len(A) + 1):
                C[I][J] += A[I][K] * B[K][J]
    return C


def MATOUT(A, IM, IN, M, N, LABEL):
    """
    PRINT MATRICES OF SIZE M BY N
    """
    print "\nTHE %s ARRAY" % LABEL
    line = ""
    for i in A:
        line += "\t\t%d" % i
    print line
    for i in A:
        line = str(i)
        for j in A[i]:
            line += "\t%.6f" % A[i][j]
        print line
    print

# In python functions can only be called when defined *before*, thus this MAIN()
if __name__ == "__main__":
    MAIN()