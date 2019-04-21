/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

// MAST includes
#include "aeroelasticity/time_domain_flutter_solver.h"
#include "aeroelasticity/time_domain_flutter_solution.h"
#include "aeroelasticity/time_domain_flutter_root.h"
#include "aeroelasticity/time_domain_flutter_root_crossover.h"
#include "elasticity/structural_fluid_interaction_assembly.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "numerics/lapack_dggev_interface.h"
#include "base/parameter.h"
#include "base/nonlinear_system.h"


MAST::TimeDomainFlutterSolver::TimeDomainFlutterSolver():
MAST::FlutterSolverBase(),
_velocity_param(nullptr),
_V_range(),
_n_V_divs(0.) {
    
}



MAST::TimeDomainFlutterSolver::~TimeDomainFlutterSolver() {
    
    this->clear();
}



void
MAST::TimeDomainFlutterSolver::clear() {
    
    this->clear_solutions();
    
    _velocity_param   = nullptr;
    _V_range          = std::pair<Real, Real>(0.,0.);
    _n_V_divs         = 0;
    
    MAST::FlutterSolverBase::clear();
}




void
MAST::TimeDomainFlutterSolver::
initialize(MAST::Parameter&                             velocity_param,
           Real                                         V_lower,
           Real                                         V_upper,
           unsigned int                                 n_V_divs,
           std::vector<libMesh::NumericVector<Real> *>& basis) {

    _velocity_param = &velocity_param;
    _V_range.first  = V_lower;
    _V_range.second = V_upper;
    _n_V_divs       = n_V_divs;
    
    MAST::FlutterSolverBase::initialize(basis);
}




void
MAST::TimeDomainFlutterSolver::clear_solutions() {
    
    std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
    _flutter_solutions.begin();
    
    for ( ; it != _flutter_solutions.end(); it++)
        delete it->second;
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::iterator cross_it =
    _flutter_crossovers.begin();
    
    for ( ; cross_it != _flutter_crossovers.end(); cross_it++)
        delete cross_it->second;
    
    _flutter_solutions.clear();
    _flutter_crossovers.clear();
}




unsigned int
MAST::TimeDomainFlutterSolver::n_roots_found() const {
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int n = 0;
    for ( ; it!=end; it++)
        if (it->second->root) // a valid root pointer has been assigned
            n++;
    
    return n;
}





const MAST::FlutterRootBase&
MAST::TimeDomainFlutterSolver::get_root(const unsigned int n) const {
    
    libmesh_assert(n < n_roots_found());
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int root_num = 0;
    for ( ; it!=end; it++) {
        if (it->second->root) { // a valid root pointer has been assigned
            if (root_num == n)  // root num matches the one being requested
                return *(it->second->root);
            else
                root_num++;
        }
    }
    
    libmesh_assert(false); // should not get here
}





std::pair<bool, MAST::FlutterRootBase*>
MAST::TimeDomainFlutterSolver::find_next_root(const Real g_tol,
                                              const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    while ( it != end)
    {
        MAST::FlutterRootCrossoverBase* cross = it->second;
        
        if (!cross->root) {
            const unsigned int root_num = cross->root_num;
            std::pair<bool, MAST::FlutterSolutionBase*> sol;
            // first try the Newton search. If that fails, then
            // try the bisection search
            /*sol = newton_search(*cross->crossover_solutions.second,
                                root_num,
                                g_tol,
                                n_bisection_iters);
            if (!sol.first)*/
            sol =   _bisection_search(cross->crossover_solutions,
                                      root_num, g_tol, n_bisection_iters);
            
            cross->root = &(sol.second->get_root(root_num));
            
            // now, remove this entry from the _flutter_crossover points and
            // reinsert it with the actual critical velocity
            _flutter_crossovers.erase(it);
            std::pair<Real, MAST::FlutterRootCrossoverBase*>
            val(cross->root->V, cross);
            _flutter_crossovers.insert(val);
            return std::pair<bool, MAST::FlutterRootBase*> (true, cross->root);
        }
        
        it++;
    }
    
    // if it gets here, no new root was found
    return std::pair<bool, MAST::FlutterRootBase*> (false, nullptr);
}



std::pair<bool,  MAST::FlutterRootBase*>
MAST::TimeDomainFlutterSolver::find_critical_root(const Real g_tol,
                                                  const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    if (it == end) // no potential cross-over points were identified
        return std::pair<bool, MAST::FlutterRootBase*> (false, nullptr);
    
    // it is possible that once the root has been found, its velocity end up
    // putting it at a higher velocity in the map, so we need to check if
    // the critical root has changed
    while (!it->second->root)
    {
        MAST::FlutterRootCrossoverBase* cross = it->second;
        
        if (!cross->root) {
            
            const unsigned int root_num = cross->root_num;
            std::pair<bool, MAST::FlutterSolutionBase*> sol;
            // first try the Newton search. If that fails, then
            // try the bisection search
            /*sol = newton_search(*cross->crossover_solutions.second,
                                root_num,
                                g_tol,
                                n_bisection_iters);
            if (!sol.first)*/
            sol =   _bisection_search(cross->crossover_solutions,
                                     root_num, g_tol, n_bisection_iters);
            
            cross->root = &(sol.second->get_root(root_num));
            
            // now, remove this entry from the _flutter_crossover points and
            // reinsert it with the actual critical velocity
            _flutter_crossovers.erase(it);
            std::pair<Real, MAST::FlutterRootCrossoverBase*>
            val(cross->root->V, cross);
            _flutter_crossovers.insert(val);
        }
        
        // update the iterator to make sure that this is updated
        it = _flutter_crossovers.begin(); end = _flutter_crossovers.end();
    }
    
    // if it gets here, then the root was successfully found
    return std::pair<bool, MAST::FlutterRootBase*> (true, it->second->root);
}





std::pair<bool,  MAST::FlutterRootBase*>
MAST::TimeDomainFlutterSolver::
analyze_and_find_critical_root_without_tracking(const Real g_tol,
                                                const unsigned int n_bisection_iters) {

    // make sure that there are no solutions already stored
    libmesh_assert(! _flutter_solutions.size());
    libmesh_assert(!_flutter_crossovers.size());
    
    
    //
    // start with the previous velocity and increment till a single
    // root cross-over is available. This is done without need for
    // mode tracking. Presently, the algorithm expects that the number
    // of unstable roots at the first velocity be zero.
    //
    Real
    lower_V  = _V_range.first,
    upper_V  = _V_range.second,
    lower_g  = 0.,
    upper_g  = 0.,
    dV       =  (_V_range.second - _V_range.first)/_n_V_divs,
    new_V    = 0.;
    
    const MAST::FlutterRootBase
    *lower_root = nullptr,
    *upper_root = nullptr;
    
    std::pair<unsigned int, unsigned int>
    bracket_n_unstable_roots(0, 0);
    
    // first analyze at both the ends of the bracket
    MAST::TimeDomainFlutterSolution
    *sol                            = _analyze(lower_V).release();
    lower_root                      = sol->get_critical_root(g_tol);
    bracket_n_unstable_roots.first  = sol->n_unstable_roots_in_upper_complex_half(g_tol);
    if (_output)
        sol->print(*_output);

    // presently the algorithm requires that the first velocity has no unstable
    // roots
    if (lower_V > 0)
        libmesh_assert(!bracket_n_unstable_roots.first);
    
    // add the solution to this solver
    bool if_success =
    _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                              (lower_V, sol)).second;
    
    libmesh_assert(if_success);
    

    // now increment velocity till atleast one root is found to be critical
    bool
    cont  = true;

    while (cont) {
        
        // add a new analysis point
        upper_V                         = lower_V + dV;

        sol                             = _analyze(upper_V, sol).release();
        bracket_n_unstable_roots.second = sol->n_unstable_roots_in_upper_complex_half(g_tol);
        upper_root                      = sol->get_critical_root(g_tol);
        if (_output)
            sol->print(*_output);

        // add the solution to this solver
        if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (upper_V, sol)).second;
        
        libmesh_assert(if_success);

        // check if any new roots were found
        if (bracket_n_unstable_roots.second ||
            upper_V >= _V_range.second) {
            
            // we have found a pair of roots that surround a cross-over point
            cont = false;
        }
        else {
            
            // increment the velocity
            lower_V                         = upper_V;
            upper_V                        += dV;
            bracket_n_unstable_roots.first  = bracket_n_unstable_roots.second;
            bracket_n_unstable_roots.second = 0;
            lower_root                      = upper_root;
        }
    }
    
    std::pair<bool, MAST::FlutterRootBase*> rval(false, nullptr);
    
    // if no critical roots were found, the return with false
    if (!bracket_n_unstable_roots.second)
        return rval;
    
    // next, refine the bracket till a root is found till the desired accuracy
    unsigned int n_iters = 0;
    
    while (n_iters < n_bisection_iters) {
        
        // get the damping values from the two critical roots
        lower_g    =   lower_root->root.real();
        upper_g    =   upper_root->root.real();
        
        // if the lower velocity is 0, then we calculate V using upper V
        // as the reference.
        //if (lower_V > 0.)
        //    new_V      =
        //    lower_V + (upper_V-lower_V)/(upper_g-lower_g)*(0.-lower_g);    // using lower V as reference
        //else
        // a dampig factor of 0.1 is added to slow down the reduction. This is
        // done because the g value tends to increase very quickly after the
        // collision of roots.
        new_V      = upper_V +
        0.5*(upper_V-lower_V)/(upper_g-lower_g)*(1.e-4-upper_g);  // using upper V as reference
        
        sol        = _analyze(new_V, sol).release();
        if (_output)
            sol->print(*_output);
        
        // get the critical root from here
        const MAST::FlutterRootBase* root = sol->get_critical_root(g_tol);

        // add the solution to this solver
        if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (new_V, sol)).second;
        
        libmesh_assert(if_success);
        
        
        // check the new damping value
        if ((root->root.real() > 0.) &&  // only positively unstable roots will be used.
            (root->root.real() <= g_tol)) {
            
            rval.first  = true;
            rval.second = sol->get_critical_root(g_tol);
            return  rval;
        }
        
        // update the V value
        if (root->root.real() < 0.) {
            
            lower_V    = new_V;
            lower_g    = root->root.real();
            lower_root = sol->get_critical_root(g_tol);
        }
        else {
            
            upper_V    = new_V;
            upper_g    = root->root.real();
            upper_root = sol->get_critical_root(g_tol);
        }
        
        n_iters++;
    }
    
    // return false, along with the latest sol
    rval.first    = false;
    rval.second   = sol->get_critical_root(g_tol);
    
    return rval;
}




void
MAST::TimeDomainFlutterSolver::scan_for_roots() {
    
    // if the initial scanning has not been done, then do it now
    if (!_flutter_solutions.size()) {
        // march from the upper limit to the lower to find the roots
        Real current_V = _V_range.first,
        delta_V = (_V_range.second - _V_range.first)/_n_V_divs;
        
        std::vector<Real> V_vals(_n_V_divs+1);
        for (unsigned int i=0; i<_n_V_divs+1; i++) {
            V_vals[i] = current_V;
            current_V += delta_V;
        }
        V_vals[_n_V_divs] = _V_range.second; // to get around finite-precision arithmetic
        
        MAST::FlutterSolutionBase* prev_sol = nullptr;
        for (unsigned int i=0; i<_n_V_divs+1; i++) {
            current_V = V_vals[i];
            std::unique_ptr<MAST::TimeDomainFlutterSolution>
            sol = _analyze(current_V, prev_sol);
            
            prev_sol = sol.get();
            
            if (_output)
                sol->print(*_output);
            
            // add the solution to this solver
            bool if_success =
            _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                      (current_V, sol.release())).second;
            
            libmesh_assert(if_success);
        }
        
        _identify_crossover_points();
    }
}





void
MAST::TimeDomainFlutterSolver::print_sorted_roots() {
    
    // write only if the output is set
    if (!_output)
        return;
    
    std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
    sol_it = _flutter_solutions.begin(),
    sol_end = _flutter_solutions.end();
    libmesh_assert(sol_it != sol_end); // solutions should have been evaluated
    
    unsigned int nvals = sol_it->second->n_roots();
    libmesh_assert(nvals); // should not be zero
    
    // each root is written separately one after another
    for (unsigned int i=0; i<nvals; i++)
    {
        // print the headers
        *_output
        << "** Root # "
        << std::setw(5) << i << " **" << std::endl
        << std::setw(15) << "V_ref"
        << std::setw(15) << "Re"
        << std::setw(15) << "Im" << std::endl;
        
        // update the iterator for this analysis
        sol_it = _flutter_solutions.begin();
        
        // write the data from all solutions
        for ( ; sol_it != sol_end; sol_it++)
        {
            const MAST::FlutterRootBase& root =
            sol_it->second->get_root(i);
            
            *_output
            << std::setw(15) << root.V
            << std::setw(15) << std::real(root.root)
            << std::setw(15) << std::imag(root.root) << std::endl;
        }
        *_output << std::endl << std::endl;
    }
    
    
    // write the roots identified using iterative search technique
    std::streamsize prec = _output->precision();
    
    unsigned int nroots = this->n_roots_found();
    *_output << std::endl
    << "n critical roots identified: " << nroots << std::endl;
    for (unsigned int i=0; i<nroots; i++)
    {
        const MAST::FlutterRootBase& root = this->get_root(i);
        *_output
        << "** Root : " << std::setw(5) << i << " **" << std::endl
        << "V      = " << std::setw(15) << std::setprecision(15) << root.V << std::endl
        << "g      = " << std::setw(15) << std::real(root.root) << std::endl
        << "omega  = " << std::setw(15) << std::imag(root.root) << std::endl
        << std::setprecision((int) prec) // set the precision to the default value
        << "Modal Participation : " << std::endl ;
        for (unsigned int j=0; j<nvals; j++)
            *_output
            << "(" << std::setw(5) << j << "): "
            << std::setw(10) << root.modal_participation(j)
            << std::setw(3)  << " ";
        *_output << std::endl << std::endl;
    }
    
    
}


void
MAST::TimeDomainFlutterSolver::print_crossover_points() {
    
    // print only if the output was specified
    if (!_output)
        return;
    
    *_output << "n crossover points found: "
    << std::setw(5) << _flutter_crossovers.size() << std::endl;
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    unsigned int i=0;
    
    for ( ; it != end; it++) {
        *_output << "** Point : " << std::setw(5) << i << " **" << std::endl;
        it->second->print(*_output);
        *_output << std::endl;
        i++;
    }
}





std::pair<bool, MAST::FlutterSolutionBase*>
MAST::TimeDomainFlutterSolver::
_bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                  MAST::FlutterSolutionBase*>& ref_sol_range,
                  const unsigned int root_num,
                  const Real g_tol,
                  const unsigned int max_iters) {
    
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real
    lower_V  = ref_sol_range.first->ref_val(),
    lower_g  = ref_sol_range.first->get_root(root_num).root.real(),
    upper_V  = ref_sol_range.second->ref_val(),
    upper_g  = ref_sol_range.second->get_root(root_num).root.real(),
    new_V    = 0.;
    unsigned int n_iters = 0;
    
    MAST::FlutterSolutionBase* new_sol = nullptr;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, nullptr);
    
    while (n_iters < max_iters) {
        
        new_V    = lower_V +
        (upper_V-lower_V)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        
        new_sol  = _analyze(new_V, ref_sol_range.first).release();
        
        if (_output)
            new_sol->print(*_output);
        
        // add the solution to this solver
        bool if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (new_V, new_sol)).second;
        
        libmesh_assert(if_success);
        
        const MAST::FlutterRootBase& root = new_sol->get_root(root_num);
        
        // check if the new damping value
        if (fabs(root.root.real()) <= g_tol) {
            
            rval.first = true;
            rval.second = new_sol;
            return  rval;
        }
        
        // update the V value
        if (root.root.real() < 0.) {
            
            lower_V = new_V;
            lower_g = root.root.real();
        }
        else {
            
            upper_V = new_V;
            upper_g = root.root.real();
        }
        
        n_iters++;
    }
    
    // return false, along with the latest sol
    rval.first = false;
    rval.second = new_sol;
    
    return rval;
}




std::unique_ptr<MAST::TimeDomainFlutterSolution>
MAST::TimeDomainFlutterSolver::_analyze(const Real v_ref,
                                       const MAST::FlutterSolutionBase* prev_sol) {
    
    libMesh::out
    << " ====================================================" << std::endl
    << "Eigensolution" << std::endl
    << "   V_ref = " << std::setw(10) << v_ref << std::endl;
    
    RealMatrixX
    A,
    B;
    
    // initialize the matrices for the structure.
    _initialize_matrices(v_ref, A, B);
    
    MAST::LAPACK_DGGEV ges;
    ges.compute(A, B);
    ges.scale_eigenvectors_to_identity_innerproduct();
    
    MAST::TimeDomainFlutterSolution* root = new MAST::TimeDomainFlutterSolution;
    root->init(*this, v_ref, ges);
    if (prev_sol)
        root->sort(*prev_sol);
    
    libMesh::out
    << "Finished Eigensolution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return std::unique_ptr<MAST::TimeDomainFlutterSolution> (root);
}




void
MAST::TimeDomainFlutterSolver::_initialize_matrices(Real U_inf,
                                                    RealMatrixX &A,
                                                    RealMatrixX &B) {
    
    // now create the matrices for first-order model
    // original equations are
    //    M x_ddot + C x_dot + K x = q_dyn (A0 x + A1 x_dot)
    //  defining x1 = x_dot, the first order governing equations become
    //    [  I  0 ]   x_dot  =  [  0  I ]   x
    //    [  0  M ]  x1_dot  =  [ K1  C1]  x1
    //    with K1  = -K + q_dyn A0, and
    //         C1  = -C + q_dyn A1.
    //
    //   Then, the eigenproblem is formulated as
    //   A y = p B y, where
    //   y  = {x^T x1^T}^T;
    //   A  = [  0  I  ]
    //        [-K1 -C1 ], and
    //   B  = [  I  0 ]
    //        [  0  M ].
    //

    
    // set the velocity value in the parameter that was provided
    (*_velocity_param) = U_inf;
    

    // if the steady solver object is provided, then solve for the
    // steady state using this velocity
    if (_steady_solver) {
        libMesh::out
        << "***  Performing Steady State Solve ***" << std::endl;
        
        _steady_solver->solve();
    }
    
    
    const unsigned int n = (unsigned int)_basis_vectors->size();

    RealMatrixX
    m      =  RealMatrixX::Zero(n, n),
    c      =  RealMatrixX::Zero(n, n),
    k      =  RealMatrixX::Zero(n, n);

    
    // now prepare a map of the quantities and ask the assembly object to
    // calculate the quantities of interest.
    std::map<MAST::StructuralQuantityType, RealMatrixX*> qty_map;
    qty_map[MAST::MASS]       = &m;
    qty_map[MAST::DAMPING]    = &c;
    qty_map[MAST::STIFFNESS]  = &k;
    
    
    _assembly->assemble_reduced_order_quantity(*_basis_vectors,
                                               qty_map);
    
    
    // put the matrices back in the system matrices
    A.setZero(2*n, 2*n);
    B.setZero(2*n, 2*n);
    
    
    B.topLeftCorner(n, n)      = RealMatrixX::Identity(n,n );
    B.bottomRightCorner(n, n)  = m;
    
    
    A.topRightCorner(n, n)     = RealMatrixX::Identity(n, n);
    A.bottomLeftCorner(n, n)   = -k;
    A.bottomRightCorner(n, n)  = -c;
}






void
MAST::TimeDomainFlutterSolver::
_initialize_matrix_sensitivity_for_param(const MAST::FunctionBase& f,
                                         const libMesh::NumericVector<Real>& dXdp,
                                         Real U_inf,
                                         RealMatrixX& A,
                                         RealMatrixX& B) {
    
    // now create the matrices for first-order model
    // original equations are
    //    M x_ddot + C x_dot + K x = q_dyn (A0 x + A1 x_dot)
    //  defining x1 = x_dot, the first order governing equations become
    //    [  I  0 ]   x_dot  =  [  0  I ]   x
    //    [  0  M ]  x1_dot  =  [ K1  C1]  x1
    //    with K1  = -K + q_dyn A0, and
    //         C1  = -C + q_dyn A1.
    //
    //   Then, the eigenproblem is formulated as
    //   A y = p B y, where
    //   y  = {x^T x1^T}^T;
    //   A  = [  0  I  ]
    //        [-K1 -C1 ], and
    //   B  = [  I  0 ]
    //        [  0  M ].
    //
    
    const unsigned int n = (unsigned int)_basis_vectors->size();
    
    RealMatrixX
    m      =  RealMatrixX::Zero(n, n),
    c      =  RealMatrixX::Zero(n, n),
    k      =  RealMatrixX::Zero(n, n);
    
    
    // now prepare a map of the quantities and ask the assembly object to
    // calculate the quantities of interest.
    std::map<MAST::StructuralQuantityType, RealMatrixX*> qty_map;
    qty_map[MAST::MASS]       = &m;
    qty_map[MAST::DAMPING]    = &c;
    qty_map[MAST::STIFFNESS]  = &k;
    
    
    // set the velocity value in the parameter that was provided
    (*_velocity_param) = U_inf;
    
    _assembly->set_base_solution(dXdp, true);
    _assembly->assemble_reduced_order_quantity_sensitivity
    (f, *_basis_vectors, qty_map);
    _assembly->clear_base_solution(true);
    
    
    // put the matrices back in the system matrices
    A.setZero(2*n, 2*n);
    B.setZero(2*n, 2*n);
    
    
    B.topLeftCorner(n, n)      = RealMatrixX::Identity(n,n );
    B.bottomRightCorner(n, n)  = m;
    
    
    A.topRightCorner(n, n)     = RealMatrixX::Identity(n, n);
    A.bottomLeftCorner(n, n)   = -k;
    A.bottomRightCorner(n, n)  = -c;
}





void
MAST::TimeDomainFlutterSolver::_identify_crossover_points() {
    
    /////////////////////////////////////////////////////////////////////
    // only non-negative frequencies (>=0) are used, since a time-domain
    // solver provides complex-conjugate roots.
    /////////////////////////////////////////////////////////////////////
    
    // if the initial scanning has not been done, then do it now
    const Real tol = 1.0e-5;
    
    const unsigned int nvals = _flutter_solutions.begin()->second->n_roots();
    // make sure that the solution has been generated
    libmesh_assert(nvals);
    
    //
    // for some cases some roots trail along the g=0 axis
    // and should not be considered as flutter. These are simply
    // modes where aerodynamics do not provide any damping.
    // These modes will not have a damping more than tolerance
    //
    
    std::vector<bool> modes_to_neglect(nvals);
    std::fill(modes_to_neglect.begin(),
              modes_to_neglect.end(), false);
    
    // look for the max g val for a mode, which will indicate if the
    // mode is undamped or not
    for (unsigned int i=0; i<nvals; i++) {
        
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
        sol_it    = _flutter_solutions.begin(),
        sol_end   = _flutter_solutions.end();
        
        Real
        max_g_val = 0.,
        max_w_val = -1.e100,
        g_val = 0.,
        w_val = 0.;
        
        
        for ( ; sol_it!=sol_end; sol_it++) {
            
            g_val = fabs(sol_it->second->get_root(i).root.real());
            w_val = sol_it->second->get_root(i).root.imag();
            
            // maximum damping
            if (g_val > max_g_val)
                max_g_val = g_val;
            
            if (w_val > max_w_val)
                max_w_val = w_val;
        }

        // check the maximum damping seen for this mode
        if (max_g_val < tol || max_w_val < 0.)
            modes_to_neglect[i] = true;
    }
    
    // now look for oscillatory roots crossover points in increasing
    // order of V
    for (unsigned int i=0; i<nvals; i++) {
        
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
        sol_it      = _flutter_solutions.begin(), // first of the pair
        sol_itp1    = _flutter_solutions.begin(), // first of the pair
        sol_end     = _flutter_solutions.end();
        
        if (sol_it == sol_end)
            return;
        
        sol_itp1++; // increment for the next pair of results
        bool if_process = false;
        
        while (sol_itp1 != sol_end) {
            
            if_process =
            (// atleast one |g| > tol
             (fabs(sol_it->second->ref_val()) > tol ||
              fabs(sol_itp1->second->ref_val()) > tol) &&
             // both |g| < max_g
             //(fabs(sol_it->second->get_root(i).root.real()) < max_allowable_g &&
             // fabs(sol_itp1->second->get_root(i).root.real()) < max_allowable_g) &&
             // if the mode has been identified to be trailing along g =0,
             // neglect it
             !modes_to_neglect[i]);
            

            if (if_process) {
                
                // look for the flutter roots
                MAST::FlutterSolutionBase
                *lower = sol_it->second,
                *upper = sol_itp1->second;
                
                if ((lower->get_root(i).root.real() <= 0.) &&
                    (upper->get_root(i).root.real() > 0.)) {
                    
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::TimeDomainFlutterRootCrossover;
                    cross->crossover_solutions.first  = lower; // -ve g
                    cross->crossover_solutions.second = upper; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( lower->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
                else if ((lower->get_root(i).root.real() > 0.) &&
                         (upper->get_root(i).root.real() <= 0.)) {
                    
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::TimeDomainFlutterRootCrossover;
                    cross->crossover_solutions.first  = upper; // -ve g
                    cross->crossover_solutions.second = lower; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( upper->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
            }
            
            // increment the pointers for next pair of roots
            sol_it++;
            sol_itp1++;
        }
        
        // now check to see if the root started out as critical at
        // the lowest V value.
        sol_it     = _flutter_solutions.begin();
        sol_itp1   = _flutter_solutions.begin();
        sol_itp1++;
        
        if_process =
        (// atleast one |g| > tol
         (fabs(sol_it->second->ref_val()) > tol ||
          fabs(sol_itp1->second->ref_val()) > tol) &&
         // both |g| < max_g
         //(fabs(sol_it->second->get_root(i).root.real()) < max_allowable_g &&
         // fabs(sol_itp1->second->get_root(i).root.real()) < max_allowable_g) &&
         // if the mode has been identified to be trailing along g =0,
         // neglect it
         !modes_to_neglect[i]);
        
        Real g_val = sol_it->second->get_root(i).root.real();
        
        if (if_process &&
            g_val > 0 /*&& g_val < max_allowable_g*/) {
            
            MAST::FlutterRootCrossoverBase* cross =
            new MAST::TimeDomainFlutterRootCrossover;
            // here, both roots for crossover are set to be the same
            // note that this will not work with bisection search, since
            // it needs a bracket to begin with.
            // However, the Newton search should be able to find the
            // critical root location.
            cross->crossover_solutions.first  = sol_it->second;
            cross->crossover_solutions.second = sol_it->second;
            cross->root_num = i;
            std::pair<Real, MAST::FlutterRootCrossoverBase*>
            val( sol_it->second->get_root(i).V, cross);
            _flutter_crossovers.insert(val);
        }
        
    }
}




void
MAST::TimeDomainFlutterSolver::
calculate_sensitivity(MAST::FlutterRootBase& root,
                      const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>* dXdp,
                      libMesh::NumericVector<Real>* dXdV) {
    
    
    
    libMesh::out
    << " ====================================================" << std::endl
    << "Flutter Sensitivity Solution" << std::endl
    << "   V_ref = " << std::setw(10) << root.V << std::endl;
    
    Complex
    eig              = root.root,
    deig_dp          = 0.,
    deig_dV          = 0.,
    den              = 0.;
    
    // get the sensitivity of the matrices
    RealMatrixX
    mat_A,
    mat_B,
    mat_A_sens,
    mat_B_sens;
    
    ComplexVectorX v;
    
    // initialize the baseline matrices
    _initialize_matrices(root.V, mat_A, mat_B);
    
    // if the sensitivity of the solution was provided, then use that.
    // otherwise pass a zero vector
    libMesh::NumericVector<Real>* sol_sens = dXdp;
    std::unique_ptr<libMesh::NumericVector<Real> > zero_sol_sens;
    if (!dXdp) {
        zero_sol_sens.reset(_assembly->system().solution->zero_clone().release());
        sol_sens = zero_sol_sens.get();
    }
    else
        sol_sens = dXdp;
    
    // calculate the eigenproblem sensitivity
    _initialize_matrix_sensitivity_for_param(f,
                                             *sol_sens,
                                             root.V,
                                             mat_A_sens,
                                             mat_B_sens);
    
    
    // the eigenproblem is     y^T A x - lambda y^T B x = 0
    // therefore, the denominator is obtained from the inner product of
    // y^T B x
    // sensitivity is
    //   -dlambda/dp y^T B x = - y^T (dA/dp - lambda dB/dp)
    // or
    //   dlambda/dp = [y^T (dA/dp - lambda dB/dp)]/(y^T B x)
    
    // now calculate the numerator for sensitivity
    // numerator =  ( dA/dp - lambda dB/dp)
    den     = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
    deig_dp = root.eig_vec_left.dot((mat_A_sens.cast<Complex>() -
                                    eig*mat_B_sens.cast<Complex>())*root.eig_vec_right)/den;
    
    // next we need the sensitivity of eigenvalue wrt V
    // identify the sensitivity of solution to be used based on the
    // function arguments
    if (!dXdV) {
        sol_sens = zero_sol_sens.get();
    }
    else
        sol_sens = dXdV;
    
    _initialize_matrix_sensitivity_for_param(*_velocity_param,
                                             *sol_sens,
                                             root.V,
                                             mat_A_sens,
                                             mat_B_sens);

    
    // now calculate the quotient for sensitivity wrt V
    // calculate numerator
    deig_dV = root.eig_vec_left.dot((mat_A_sens.cast<Complex>() -
                                     eig*mat_B_sens.cast<Complex>())*root.eig_vec_right)/den;
    
    // since the constraint that defines flutter speed is that damping = 0,
    // Re(lambda) = 0, then the sensitivity of flutter speed is obtained
    // from the total derivative of this constraint
    //     d Re(lambda)/dp + d Re(lambda)/dV dV/dp = 0
    // or, dV/dp = -[d Re(lambda)/dp] / [d Re(lambda)/dV]
    // finally, the flutter speed sensitivity
    root.V_sens     = - deig_dp.real() / deig_dV.real();
    
    // total sensitivity of the eigenvlaue
    root.root_sens  = deig_dp + deig_dV * root.V_sens;
    
    root.has_sensitivity_data = true;

    libMesh::out
    << "Finished Flutter Sensitivity Solution" << std::endl
    << " ====================================================" << std::endl;
    
}
