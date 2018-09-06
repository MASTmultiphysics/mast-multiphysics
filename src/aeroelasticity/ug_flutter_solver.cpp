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


// MAST includes
#include "aeroelasticity/ug_flutter_solver.h"
#include "aeroelasticity/ug_flutter_solution.h"
#include "aeroelasticity/ug_flutter_root.h"
#include "aeroelasticity/ug_flutter_root_crossover.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "numerics/lapack_zggev_interface.h"
#include "base/parameter.h"
#include "base/nonlinear_system.h"


MAST::UGFlutterSolver::UGFlutterSolver():
MAST::FlutterSolverBase(),
_kr_param(nullptr),
_bref_param(nullptr),
_kr_range(),
_n_kr_divs(0.),
_include_highest_kr_unstable(false) {
    
}



MAST::UGFlutterSolver::~UGFlutterSolver() {
    
    this->clear();
}




void
MAST::UGFlutterSolver::clear() {
    
    this->clear_solutions();
    
    _kr_param          = nullptr;
    _bref_param        = nullptr;
    _kr_range          = std::pair<Real, Real>(0.,0.);
    _n_kr_divs         = 0;
    
    MAST::FlutterSolverBase::clear();
}




void
MAST::UGFlutterSolver::
initialize(MAST::Parameter&                             kr_param,
           MAST::Parameter&                             bref_param,
           Real                                         rho,
           Real                                         kr_lower,
           Real                                         kr_upper,
           unsigned int                                 n_kr_divs,
           std::vector<libMesh::NumericVector<Real> *>& basis) {
    
    
    _kr_param        = &kr_param;
    _bref_param      = &bref_param;
    _rho             = rho;
    _kr_range.first  = kr_lower;
    _kr_range.second = kr_upper;
    _n_kr_divs       = n_kr_divs;

    MAST::FlutterSolverBase::initialize(basis);
}




void
MAST::UGFlutterSolver::clear_solutions() {
    
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
MAST::UGFlutterSolver::n_roots_found() const {
    
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
MAST::UGFlutterSolver::get_root(const unsigned int n) const {
    
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
MAST::UGFlutterSolver::find_next_root(const Real g_tol,
                                              const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
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
MAST::UGFlutterSolver::find_critical_root(const Real g_tol,
                                                  const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
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




void
MAST::UGFlutterSolver::scan_for_roots() {
    
    // if the initial scanning has not been done, then do it now
    if (!_flutter_solutions.size()) {
        // march from the upper limit to the lower to find the roots
        Real
        current_kr  = _kr_range.second,
        delta_kr    = (_kr_range.second - _kr_range.first)/_n_kr_divs;
        
        std::vector<Real> k_vals(_n_kr_divs+1);
        for (unsigned int i=0; i<_n_kr_divs+1; i++) {
            k_vals[i]      = current_kr;
            current_kr    -= delta_kr;
        }
        k_vals[_n_kr_divs] = _kr_range.first; // to get around finite-precision arithmetic
        
        MAST::FlutterSolutionBase* prev_sol = nullptr;
        for (unsigned int i=0; i< _n_kr_divs+1; i++) {
            
            current_kr = k_vals[i];
            std::unique_ptr<MAST::FlutterSolutionBase> sol =
            _analyze(current_kr, prev_sol);
            
            prev_sol = sol.get();
            
            if (_output)
                sol->print(*_output);
            
            // add the solution to this solver
            bool if_success =
            _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                      (current_kr, sol.release())).second;
            
            libmesh_assert(if_success);
        }
        
        _identify_crossover_points();
    }
}





void
MAST::UGFlutterSolver::print_sorted_roots()
{
    
    // write only if the output was set
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
        << std::setw(15) << "kr"
        << std::setw(15) << "g"
        << std::setw(15) << "V" << std::endl;
        
        // update the iterator for this analysis
        sol_it = _flutter_solutions.begin();
        
        // write the data from all solutions
        for ( ; sol_it != sol_end; sol_it++)
        {
            const MAST::FlutterRootBase& root =
            sol_it->second->get_root(i);
            
            *_output
            << std::setw(15) << root.kr
            << std::setw(15) << root.g
            << std::setw(15) << root.V << std::endl;
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
        << "kr     = " << std::setw(15) << std::setprecision(15) << root.kr << std::endl
        << "V      = " << std::setw(15) << std::setprecision(15) << root.V << std::endl
        << "g      = " << std::setw(15) << root.g << std::endl
        << "omega  = " << std::setw(15) << root.omega << std::endl
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
MAST::UGFlutterSolver::print_crossover_points()
{
    // print only if the output was set
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
MAST::UGFlutterSolver::
_bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                  MAST::FlutterSolutionBase*>& ref_sol_range,
                  const unsigned int root_num,
                  const Real g_tol,
                  const unsigned int max_iters) {
    
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real
    lower_kr = ref_sol_range.first->get_root(root_num).kr,
    lower_g  = ref_sol_range.first->get_root(root_num).g,
    upper_kr = ref_sol_range.second->get_root(root_num).kr,
    upper_g  = ref_sol_range.second->get_root(root_num).g,
    new_kr   = 0.;
    unsigned int n_iters = 0;
    
    MAST::FlutterSolutionBase* new_sol = nullptr;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, nullptr);
    
    while (n_iters < max_iters) {
        
        libMesh::out
        << upper_kr << "  "
        << upper_g << "  " << std::endl
        << lower_kr << "  "
        << lower_g << "  " << std::endl;
        
        new_kr    = lower_kr +
        (upper_kr-lower_kr)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        
        new_sol  = _analyze(new_kr, ref_sol_range.first).release();
        
        if (_output)
            new_sol->print(*_output);
        
        // add the solution to this solver
        bool if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (new_kr, new_sol)).second;
        
        libmesh_assert(if_success);
        
        const MAST::UGFlutterRoot& root =
        dynamic_cast<const MAST::UGFlutterRoot&>(new_sol->get_root(root_num));
        
        // check if the new damping value
        if (fabs(root.g) <= g_tol) {
            
            rval.first = true;
            rval.second = new_sol;
            return  rval;
        }
        
        // update the V value
        if (root.g < 0.) {
            
            lower_kr = new_kr;
            lower_g = root.g;
        }
        else {
            
            upper_kr = new_kr;
            upper_g = root.g;
        }
        
        n_iters++;
    }
    
    // return false, along with the latest sol
    rval.first = false;
    rval.second = new_sol;
    
    return rval;
}




std::unique_ptr<MAST::FlutterSolutionBase>
MAST::UGFlutterSolver::_analyze(const Real kr_ref,
                                const MAST::FlutterSolutionBase* prev_sol) {
    
    libMesh::out
    << " ====================================================" << std::endl
    << "Eigensolution" << std::endl
    << "   kr_ref = " << std::setw(10) << kr_ref << std::endl;
    
    ComplexMatrixX
    A,
    B;
    
    // initialize the matrices for the structure.
    _initialize_matrices(kr_ref, A, B);
    
    MAST::LAPACK_ZGGEV ges;
    ges.compute(A, B);
    ges.scale_eigenvectors_to_identity_innerproduct();
    
    MAST::UGFlutterSolution* root = new MAST::UGFlutterSolution;
    root->init(*this, kr_ref, (*_bref_param)(), ges);
    if (prev_sol)
        root->sort(*prev_sol);
    
    libMesh::out
    << "Finished Eigensolution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return std::unique_ptr<MAST::FlutterSolutionBase> (root);
}




void
MAST::UGFlutterSolver::_initialize_matrices(Real kr,
                                            ComplexMatrixX &A,
                                            ComplexMatrixX &B) {
    
    // the UG method equations are
    //
    // ((kr/b)^2 M + rho/2 A(kr))q = lambda K q
    // where M and K are the structural reduced-order mass and stiffness
    // matrices, and A(kr) is the generalized aerodynamic force matrix.
    //
    
    const unsigned int n = (unsigned int)_basis_vectors->size();
    
    RealMatrixX
    m      =  RealMatrixX::Zero(n, n),
    k      =  RealMatrixX::Zero(n, n);
    
    ComplexMatrixX
    a      =  ComplexMatrixX::Zero(n, n);
    
    // now prepare a map of the quantities and ask the assembly object to
    // calculate the quantities of interest.
    std::map<MAST::StructuralQuantityType, RealMatrixX*> qty_map;
    qty_map[MAST::MASS]       = &m;
    qty_map[MAST::STIFFNESS]  = &k;
    
    
    // set the velocity value in the parameter that was provided
    (*_kr_param) = kr;
    
    _assembly->assemble_reduced_order_quantity(*_basis_vectors, qty_map);
    
    dynamic_cast<MAST::FSIGeneralizedAeroForceAssembly*>(_assembly)->
    assemble_generalized_aerodynamic_force_matrix(*_basis_vectors, a);
    
 
    // scale the force vector by -1 since MAST calculates all quantities
    // for a R(X)=0 equation so that matrix/vector quantity is assumed
    // to be on the left of the equality. This is not consistent with
    // the expectation of a flutter solver, which expects the force
    // vector to be defined on the RHS. Hence, we multiply the quantity
    // here to maintain consistency.
    a  *= -1.;

    
    A    = pow(kr/(*_bref_param)(),2) * m.cast<Complex>() + (_rho/2.) * a;
    B    = k.cast<Complex>();
}






void
MAST::UGFlutterSolver::
_initialize_matrix_sensitivity_for_param(const MAST::FunctionBase& f,
                                         const libMesh::NumericVector<Real>& dXdp,
                                         Real kr,
                                         ComplexMatrixX& A,
                                         ComplexMatrixX& B) {
    
    // the UG method equations are
    //
    // ((kr/b)^2 M + rho/2 A(kr))q = lambda K q
    // where M and K are the structural reduced-order mass and stiffness
    // matrices, and A(kr) is the generalized aerodynamic force matrix.
    //
    
    const unsigned int n = (unsigned int)_basis_vectors->size();
    
    RealMatrixX
    m      =  RealMatrixX::Zero(n, n),
    k      =  RealMatrixX::Zero(n, n);

    ComplexMatrixX
    a      =  ComplexMatrixX::Zero(n, n);
    
    // now prepare a map of the quantities and ask the assembly object to
    // calculate the quantities of interest.
    std::map<MAST::StructuralQuantityType, RealMatrixX*> qty_map;
    qty_map[MAST::MASS]       = &m;
    qty_map[MAST::STIFFNESS]  = &k;
    
    
    // set the velocity value in the parameter that was provided
    (*_kr_param) = kr;
    
    _assembly->assemble_reduced_order_quantity_sensitivity(f,
                                                           *_basis_vectors,
                                                           qty_map);

    // currently, sensitivity of generalized aero matrix is available only
    // for freq (ie non-structural parameters). Once the dependence on
    // base sol sensitivity is introduced, this will change.
    //dynamic_cast<MAST::FSIGeneralizedAeroForceAssembly*>(_assembly)->
    //assemble_generalized_aerodynamic_force_matrix(*_basis_vectors, a);
    
    
    // scale the force vector by -1 since MAST calculates all quantities
    // for a R(X)=0 equation so that matrix/vector quantity is assumed
    // to be on the left of the equality. This is not consistent with
    // the expectation of a flutter solver, which expects the force
    // vector to be defined on the RHS. Hence, we multiply the quantity
    // here to maintain consistency.
    a  *= -1.;
    
    A    = pow(kr/(*_bref_param)(),2) * m.cast<Complex>() +  (_rho/2.) * a;
    B    = k.cast<Complex>();
}






void
MAST::UGFlutterSolver::
_initialize_matrix_sensitivity_for_kr(Real kr,
                                      ComplexMatrixX& A,
                                      ComplexMatrixX& B) {
    
    // the UG method equations are
    //
    // ((kr/b)^2 M + rho/2 A(kr))q = lambda K q
    // where M and K are the structural reduced-order mass and stiffness
    // matrices, and A(kr) is the generalized aerodynamic force matrix.
    //
    
    const unsigned int n = (unsigned int)_basis_vectors->size();
    
    RealMatrixX
    m      =  RealMatrixX::Zero(n, n);
    
    ComplexMatrixX
    a      =  ComplexMatrixX::Zero(n, n);
    
    // now prepare a map of the quantities and ask the assembly object to
    // calculate the quantities of interest.
    std::map<MAST::StructuralQuantityType, RealMatrixX*> qty_map;
    qty_map[MAST::MASS]       = &m;
    
    
    // set the velocity value in the parameter that was provided
    (*_kr_param) = kr;
    
    _assembly->assemble_reduced_order_quantity(*_basis_vectors,
                                               qty_map);
    
    dynamic_cast<MAST::FSIGeneralizedAeroForceAssembly*>(_assembly)->
    assemble_generalized_aerodynamic_force_matrix(*_basis_vectors, a, _kr_param);
    
    // scale the force vector by -1 since MAST calculates all quantities
    // for a R(X)=0 equation so that matrix/vector quantity is assumed
    // to be on the left of the equality. This is not consistent with
    // the expectation of a flutter solver, which expects the force
    // vector to be defined on the RHS. Hence, we multiply the quantity
    // here to maintain consistency.
    a  *= -1.;

    
    A    = 2.*kr*pow((*_bref_param)(),-2) * m.cast<Complex>() +  (_rho/2.) * a;
    B    = ComplexMatrixX::Zero(n, n);
}





void
MAST::UGFlutterSolver::_identify_crossover_points() {
    
    // if the initial scanning has not been done, then do it now
    const Real tol = 1.0e-5, max_allowable_g = 0.75;
    
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
        Real max_g_val = 0., val = 0.;
        for ( ; sol_it!=sol_end; sol_it++) {
            val = fabs(sol_it->second->get_root(i).g);
            if (val > max_g_val)
                max_g_val = val;
        }
        // check the maximum damping seen for this mode
        if (max_g_val < tol)
            modes_to_neglect[i] = true;
    }
    
    // identify the flutter cross-overs. For this, move from
    // higher k to lower k, and handle the k=0 cases separately
    {
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
        sol_it    = _flutter_solutions.begin(), // first of the pair
        sol_end   = _flutter_solutions.end();
        if (sol_it == sol_end)
            return;
        
        // check if k=0 exists, and identify
        if (fabs(sol_it->second->ref_val()) < tol) { // k = 0
            
            // k=0 makes sense only for divergence roots. Do not use them
            // for crossover points if a finite damping was seen. Hence,
            // do nothing here, and move to the next iterator
            for (unsigned int i=0; i<nvals; i++) {
                if (!sol_it->second->get_root(i).if_nonphysical_root &&
                    fabs(sol_it->second->get_root(i).g) < tol) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::UGFlutterRootCrossover;
                    cross->crossover_solutions.first = sol_it->second;
                    cross->crossover_solutions.second = sol_it->second;
                    cross->root = &sol_it->second->get_root(i);
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( cross->root->V, cross);
                    _flutter_crossovers.insert(val);
                }
            }
        }
    }
    
    // now look for oscillatory roots crossover points in decreasing
    // order of k_red
    for (unsigned int i=0; i<nvals; i++) {
        std::map<Real, MAST::FlutterSolutionBase*>::const_reverse_iterator
        sol_rit      = _flutter_solutions.rbegin(), // first of the pair
        sol_ritp1    = _flutter_solutions.rbegin(), // first of the pair
        sol_rend     = _flutter_solutions.rend();
        if (sol_rit == sol_rend)
            return;
        
        sol_ritp1++; // increment for the next pair of results
        bool if_process = false;
        
        while (sol_ritp1 != sol_rend) {
            if_process =
            (// both should be valid roots
             (!sol_rit->second->get_root(i).if_nonphysical_root &&
              !sol_ritp1->second->get_root(i).if_nonphysical_root) &&
             // atleast one |g| > tol
             (fabs(sol_rit->second->ref_val()) > tol ||
              fabs(sol_ritp1->second->ref_val()) > tol) &&
             // both |g| < max_g
             (fabs(sol_rit->second->get_root(i).g) < max_allowable_g &&
              fabs(sol_ritp1->second->get_root(i).g) < max_allowable_g) &&
             // if the mode has been identified to be trailing along g =0,
             // neglect it
             !modes_to_neglect[i]);
            
            // do not use k_red = 0, or if the root is invalid
            if (if_process) { // look for the flutter roots
                MAST::FlutterSolutionBase *lower = sol_rit->second,
                *upper = sol_ritp1->second;
                
                if ((lower->get_root(i).g <= 0.) &&
                    (upper->get_root(i).g > 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::UGFlutterRootCrossover;
                    cross->crossover_solutions.first  = lower; // -ve g
                    cross->crossover_solutions.second = upper; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( lower->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
                else if ((lower->get_root(i).g > 0.) &&
                         (upper->get_root(i).g <= 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::UGFlutterRootCrossover;
                    cross->crossover_solutions.first  = upper; // -ve g
                    cross->crossover_solutions.second = lower; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( upper->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
            }
            
            // increment the pointers for next pair of roots
            sol_rit++;
            sol_ritp1++;
        }
        
        // now check to see if the root started out as critical at
        // the highest k value.
        sol_rit    = _flutter_solutions.rbegin();
        sol_ritp1   = _flutter_solutions.rbegin();
        sol_ritp1++;
        
        if_process =
        (// both should be valid roots
         (!sol_rit->second->get_root(i).if_nonphysical_root &&
          !sol_ritp1->second->get_root(i).if_nonphysical_root) &&
         // atleast one |g| > tol
         (fabs(sol_rit->second->ref_val()) > tol ||
          fabs(sol_ritp1->second->ref_val()) > tol) &&
         // both |g| < max_g
         (fabs(sol_rit->second->get_root(i).g) < max_allowable_g &&
          fabs(sol_ritp1->second->get_root(i).g) < max_allowable_g) &&
         // if the mode has been identified to be trailing along g =0,
         // neglect it
         !modes_to_neglect[i]);
        
        Real g_val = sol_rit->second->get_root(i).g;
        if (if_process &&
            g_val > 0 &&
            g_val < max_allowable_g &&
            _include_highest_kr_unstable) {
            MAST::FlutterRootCrossoverBase* cross =
            new MAST::UGFlutterRootCrossover;
            // here, both roots for crossover are set to be the same
            // note that this will not work with bisection search, since
            // it needs a bracket to begin with.
            // However, the Newton search should be able to find the
            // critical root location.
            cross->crossover_solutions.first  = sol_rit->second;
            cross->crossover_solutions.second = sol_rit->second;
            cross->root_num = i;
            std::pair<Real, MAST::FlutterRootCrossoverBase*>
            val( sol_rit->second->get_root(i).V, cross);
            _flutter_crossovers.insert(val);
        }
        
    }
}




void
MAST::UGFlutterSolver::
calculate_sensitivity(MAST::FlutterRootBase& root,
                      const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>* dXdp,
                      libMesh::NumericVector<Real>* dXdkr) {
    
    
    
    libMesh::out
    << " ====================================================" << std::endl
    << "UG Sensitivity Solution" << std::endl
    << "   k_red = " << std::setw(10) << root.kr << std::endl
    << "   V_ref = " << std::setw(10) << root.V << std::endl;
    
    Complex
    eig              = root.root,
    deig_dp          = 0.,
    deig_dkr         = 0.,
    den              = 0.;
    
    Real
    dkr_dp           = 0.,
    dg_dp            = 0.,
    dg_dkr           = 0.;

    
    // get the sensitivity of the matrices
    ComplexMatrixX
    mat_A,
    mat_B,
    mat_A_sens,
    mat_B_sens;
    
    ComplexVectorX v;
    
    // initialize the baseline matrices
    _initialize_matrices(root.kr, mat_A, mat_B);
    
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
                                             root.kr,
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
    
    // next we need the sensitivity of eigenvalue wrt kr
    // identify the sensitivity of solution to be used based on the
    // function arguments
    sol_sens = zero_sol_sens.get();
    
    _initialize_matrix_sensitivity_for_kr(root.kr,
                                          mat_A_sens,
                                          mat_B_sens);
    
    // now calculate the quotient for sensitivity wrt k_red
    // calculate numerator
    deig_dkr = root.eig_vec_left.dot((mat_A_sens.cast<Complex>() -
                                     eig*mat_B_sens.cast<Complex>())*root.eig_vec_right)/den;
    
    // using this, the following quantities are caluclated
    // since eig = lambda =  (1+ig)/V^2.
    // V =  (1/re(eig))^(-1/2)
    // g =  im(eig)/re(eig)
    //
    // Therefore, the sensitivity of g and V are
    // dV/dp =  -1/2 (1/re(eig))^(-3/2) deig_re/dp
    // dg/dp =  deig_im/dp / re(eig) - im(eig)/re(eig)^2 deig_re/dp

    dg_dp  =
    deig_dp.imag()/eig.real()  - eig.imag()/pow(eig.real(),2) * deig_dp.real();
    dg_dkr =
    deig_dkr.imag()/eig.real() - eig.imag()/pow(eig.real(),2) * deig_dkr.real();
    

    // since the constraint that defines flutter speed is that damping = 0,
    //      g(p, kr) = 0,
    // then the total derivative of this constraint is
    //      dg/dp + dg/dkr dkr/dp = 0
    // or,  dkr/dp = -dg/dp / dg/dkr
    dkr_dp       = -dg_dp / dg_dkr;
    root.kr_sens = dkr_dp;
    
    // Using this, the sensitivity of flutter speed if calculated as
    // dV/dp = dV/dp + dV/dkr dkr/dp
    root.root_sens            = deig_dp + deig_dkr * dkr_dp;
    root.V_sens               = -.5*root.root_sens.real()/pow(eig.real(), 1.5);
    root.has_sensitivity_data = true;
    
    libMesh::out
    << "Finished Flutter Sensitivity Solution" << std::endl
    << " ====================================================" << std::endl;
    
}



