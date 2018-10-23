/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __MAST_gcmma_optimization_interface_h__
#define __MAST_gcmma_optimization_interface_h__

// MAST includes
#include "optimization/optimization_interface.h"

extern "C" {
    extern void raasta_(int *M, int *N,
                        double *RAA0, double *RAA,
                        double *XMIN, double *XMAX,
                        double *DF0DX, double *DFDX);
    extern void asympg_(int *ITER, int *M, int *N,
                        double *ALBEFA, double *GHINIT,
                        double *GHDECR, double *GHINCR,
                        double *XVAL, double *XMIN, double *XMAX,
                        double *XOLD1, double *XOLD2,
                        double *XLOW, double *XUPP,
                        double *ALFA, double *BETA);
    extern void mmasug_(int *ITER, int *M, int *N, double *GEPS, int *IYFREE,
                        double *XVAL, double *XMMA,
                        double *XMIN, double *XMAX,
                        double *XLOW, double *XUPP,
                        double *ALFA, double *BETA,
                        double *A, double *B, double *C, double *Y, double *Z,
                        double *RAA0, double *RAA, double *ULAM,
                        double *F0VAL, double *FVAL,
                        double *F0APP, double *FAPP,
                        double *FMAX, double *DF0DX, double *DFDX,
                        double *P, double *Q, double *P0, double *Q0,
                        double *UU, double *GRADF, double *DSRCH, double *HESSF);
    extern void conser_(int *M, int *ICONSE,
                        double *GEPS, double *F0NEW, double *F0APP,
                        double *FNEW, double *FAPP);
    extern void raaupd_(int *M, int *N, double *GEPS,
                        double *XMMA, double *XVAL,
                        double *XMIN, double *XMAX,
                        double *XLOW, double *XUPP,
                        double *F0NEW, double *FNEW,
                        double *F0APP, double *FAPP,
                        double *RAA0, double *RAA);
    extern void xupdat_(int *N, int *ITER,
                        double *XMMA, double *XVAL,
                        double *XOLD1, double *XOLD2);
    extern void fupdat_(int *M, double *F0NEW, double *FNEW,
                        double *F0VAL, double *FVAL);
}


namespace MAST {
    
    class GCMMAOptimizationInterface: public MAST::OptimizationInterface {
        
    public:
        
        GCMMAOptimizationInterface();
        
        virtual ~GCMMAOptimizationInterface()
        { }
        
        virtual void optimize();
        
        
    protected:
        
        void _output_iteration_data(unsigned int i,
                                    const std::vector<Real>& XVAL,
                                    const std::vector<Real>& XMIN,
                                    const std::vector<Real>& XMAX,
                                    const std::vector<Real>& XLOW,
                                    const std::vector<Real>& XUPP,
                                    const std::vector<Real>& ALFA,
                                    const std::vector<Real>& BETA);
    };
}





#endif // __MAST_gcmma_optimization_interface_h__
