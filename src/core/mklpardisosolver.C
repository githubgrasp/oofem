/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "mklpardisosolver.h"

#include "compcol.h"
#include "symcompcol.h"
#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "error.h"
#include "classfactory.h"
#include "convergedreason.h"

#include <mkl.h>

namespace oofem {
REGISTER_SparseLinSolver(MKLPardisoSolver, ST_MKLPardiso);

MKLPardisoSolver :: MKLPardisoSolver(Domain *d, EngngModel *m) : SparseLinearSystemNM(d, m) { }

MKLPardisoSolver :: ~MKLPardisoSolver() { }

ConvergedReason MKLPardisoSolver :: solve(SparseMtrx &A, FloatArray &b, FloatArray &x)
{
    MKL_INT neqs = b.giveSize();
    x.resize(neqs);

    MKL_INT mtype = -2;        // Real and structurally symmetric matrix
    CompCol *mat = dynamic_cast< SymCompCol * >(&A);
    if ( !mat ) {
        mtype = 11;        // Real unsymmetric matrix
        mat = dynamic_cast< CompCol * >(&A);
        if ( !mat ) {
            OOFEM_ERROR("CompCol matrix needed for Pardiso solver");
        }
    }

    // Pardiso's CGS-implementation can't handle b = 0.
    if ( b.computeSquaredNorm() == 0 ) {
        return CR_CONVERGED;
    }

    MKL_INT *ia = mat->giveColPtr().givePointer();
    MKL_INT *ja = mat->giveRowIndex().givePointer();
    double *a = mat->giveValues().givePointer();

    Timer timer;
    timer.startTimer();

    // RHS and solution vectors.
    MKL_INT nrhs = 1;          // Number of right hand sides.

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *pt[64]; 

    // Pardiso control parameters.
    MKL_INT iparm[64]; ///@todo pardisoinit seems to write outside this array
    MKL_INT maxfct, mnum, phase, error, msglvl;

    double ddum;           // Double dummy
    MKL_INT idum;         /* Integer dummy. */
    
    // Setup Pardiso control parameters
    /* -------------------------------------------------------------------- */
    error = 0;
    for ( int i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    // Settings are here:
    // https://software.intel.com/en-us/articles/pardiso-parameter-table#table2

    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */ ///@todo I have no idea if this is suitable value. Examples use 2. / Mikael
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */ 
    iparm[34] = 1;        /* C-based zero-based indexing (only in MKL) */
    ///@todo This is not included in the table of options for some reason!

    maxfct = 1;         // Maximum number of numerical factorizations.
    mnum   = 1;         // Which factorization to use.
    msglvl = 0;         // Print statistical information
    error  = 0;         // Initialize error flag
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for ( int i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    phase = 11; 

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, 
            &neqs, a, ia, ja,
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);   // FACTORIZATION!

    if ( error != 0 ) {
        OOFEM_WARNING("Error during symbolic factorization: %d", error);
        return CR_FAILED;
    }
    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */    
    phase = 22;

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &neqs,
        (void*)a, (int*)ia, (int*)ja,
        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    if ( error != 0 ) {
        OOFEM_WARNING("ERROR during numerical factorization: %d", error);
        return CR_FAILED;
    }
    OOFEM_LOG_DEBUG("Factorization completed ...\n");
    OOFEM_LOG_DEBUG ("Number of factorization MFLOPS = %d\n", iparm[18]);

    /* -------------------------------------------------------------------- */    
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */    
    phase = 33;
    
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &neqs,
        (void*)a, (int*)ia, (int*)ja,
        &idum, &nrhs, iparm, &msglvl, (void*)b.givePointer(), (void*)x.givePointer(), &error);

    if ( error != 0 ) {
        OOFEM_WARNING("ERROR during solution: %d, iparm(20) = %d", error, iparm[20-1]);
        return CR_FAILED;
    }

    OOFEM_LOG_DEBUG("Solve completed ... \n");

    /* -------------------------------------------------------------------- */    
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */    
    phase = -1;                 /* Release internal memory. */

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &neqs, &ddum, &idum, &idum, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    timer.stopTimer();
    OOFEM_LOG_INFO( "MKLPardisoSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    return CR_CONVERGED;
}

#if 0
///@todo Parallel mode of this.
ConvergedReason MKLPardisoSolver :: solve(SparseMtrx &A, FloatMatrix &B, FloatMatrix &X)
{
    ///@todo Write subfunction for this as to not repeat everything. Should be easy to add. It is very important to add this support(!!!!!)
}
#endif
} // end namespace oofem
