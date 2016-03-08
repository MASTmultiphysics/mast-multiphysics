/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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


MAST::UGFlutterSolver::UGFlutterSolver():
MAST::FlutterSolverBase(),
_kr_param(NULL),
_bref_param(NULL),
_kr_range(),
_n_kr_divs(0.)
 {
    
}



MAST::UGFlutterSolver::~UGFlutterSolver() {
    
    this->clear();
}




void
MAST::UGFlutterSolver::clear() {
    
    this->clear_solutions();
    
    _kr_param          = NULL;
    _bref_param        = NULL;
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
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
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
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
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
    std::map<Real, MAST::FlutterRootCrossoverBase*>::iterator
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
    return std::pair<bool, MAST::FlutterRootBase*> (false, NULL);
}



std::pair<bool,  MAST::FlutterRootBase*>
MAST::UGFlutterSolver::find_critical_root(const Real g_tol,
                                                  const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::map<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    if (it == end) // no potential cross-over points were identified
        return std::pair<bool, MAST::FlutterRootBase*> (false, NULL);
    
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
        current_kr = _kr_range.first,
        delta_kr   = (_kr_range.second - _kr_range.first)/_n_kr_divs;
        
        std::vector<Real>
        kr_vals(_n_kr_divs+1);
        
        
        for (unsigned int i=0; i<_n_kr_divs+1; i++) {
            kr_vals[i] = current_kr;
            current_kr += delta_kr;
        }
        
        kr_vals[_n_kr_divs] = _kr_range.second; // to get around finite-precision arithmetic
        
        MAST::FlutterSolutionBase* prev_sol = NULL;
        for (unsigned int i=0; i<_n_kr_divs+1; i++) {
            current_kr    = kr_vals[i];
            std::auto_ptr<MAST::FlutterSolutionBase> sol =
            _analyze(current_kr, prev_sol);
            
            prev_sol = sol.get();
            
            sol->print(_output);
            
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
MAST::UGFlutterSolver::print_sorted_roots(std::ostream* output)
{
    if (!output)
        output = &_output;
    
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
        *output
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
            
            *output
            << std::setw(15) << root.V
            << std::setw(15) << std::real(root.root)
            << std::setw(15) << std::imag(root.root) << std::endl;
        }
        *output << std::endl << std::endl;
    }
    
    
    // write the roots identified using iterative search technique
    std::streamsize prec = output->precision();
    
    unsigned int nroots = this->n_roots_found();
    *output << std::endl
    << "n critical roots identified: " << nroots << std::endl;
    for (unsigned int i=0; i<nroots; i++)
    {
        const MAST::FlutterRootBase& root = this->get_root(i);
        *output
        << "** Root : " << std::setw(5) << i << " **" << std::endl
        << "V      = " << std::setw(15) << std::setprecision(15) << root.V << std::endl
        << "g      = " << std::setw(15) << std::real(root.root) << std::endl
        << "omega  = " << std::setw(15) << std::imag(root.root) << std::endl
        << std::setprecision(prec) // set the precision to the default value
        << "Modal Participation : " << std::endl ;
        for (unsigned int j=0; j<nvals; j++)
            *output
            << "(" << std::setw(5) << j << "): "
            << std::setw(10) << root.modal_participation(j)
            << std::setw(3)  << " ";
        *output << std::endl << std::endl;
    }
}


void
MAST::UGFlutterSolver::print_crossover_points(std::ostream* output)
{
    if (!output)
        output = &_output;
    
    *output << "n crossover points found: "
    << std::setw(5) << _flutter_crossovers.size() << std::endl;
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    unsigned int i=0;
    
    for ( ; it != end; it++) {
        *output << "** Point : " << std::setw(5) << i << " **" << std::endl;
        it->second->print(*output);
        *output << std::endl;
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
    lower_kr = ref_sol_range.first->ref_val(),
    lower_g  = ref_sol_range.first->get_root(root_num).root.real(),
    upper_kr = ref_sol_range.second->ref_val(),
    upper_g  = ref_sol_range.second->get_root(root_num).root.real(),
    new_kr   = 0.;
    unsigned int n_iters = 0;
    
    MAST::FlutterSolutionBase* new_sol = NULL;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, NULL);
    
    while (n_iters < max_iters) {
        
        new_kr    = lower_kr +
        (upper_kr-lower_kr)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        
        new_sol  = _analyze(new_kr, ref_sol_range.first).release();
        
        new_sol->print(_output);
        
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




std::auto_ptr<MAST::FlutterSolutionBase>
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
    
    
    return std::auto_ptr<MAST::FlutterSolutionBase> (root);
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
    
    _assembly->assemble_reduced_order_quantity(*_assembly->system().solution,
                                               *_basis_vectors,
                                               qty_map);
    dynamic_cast<MAST::FSIGeneralizedAeroForceAssembly*>(_assembly)->
    assemble_generalized_aerodynamic_force_matrix(*_assembly->system().solution,
                                                  *_basis_vectors,
                                                  a);
    
    
    A    = pow(kr/(*_bref_param)(),2) * m.cast<Complex>() +  (_rho/2.) * a;
    B    = k.cast<Complex>();
}






void
MAST::UGFlutterSolver::
_initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                         const unsigned int i,
                                         const libMesh::NumericVector<Real>& dXdp,
                                         Real kr,
                                         ComplexMatrixX& A,
                                         ComplexMatrixX& B) {
    
    libmesh_assert(false);
    
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
    
    _assembly->assemble_reduced_order_quantity_sensitivity
    (params,
     i,
     *_assembly->system().solution,
     dXdp,
     *_basis_vectors,
     qty_map);
    
    
    // put the matrices back in the system matrices
    A    = pow(kr/(*_bref_param)(),2) * m.cast<Complex>() +  (_rho/2.) * a;
    B    = k.cast<Complex>();
}





void
MAST::UGFlutterSolver::_identify_crossover_points() {
    
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
                    new MAST::UGFlutterRootCrossover;
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
            new MAST::UGFlutterRootCrossover;
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
MAST::UGFlutterSolver::
calculate_sensitivity(MAST::FlutterRootBase& root,
                      const libMesh::ParameterVector& params,
                      const unsigned int i,
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
    ComplexMatrixX
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
    std::auto_ptr<libMesh::NumericVector<Real> > zero_sol_sens;
    if (!dXdp) {
        zero_sol_sens.reset(_assembly->system().solution->zero_clone().release());
        sol_sens = zero_sol_sens.get();
    }
    else
        sol_sens = dXdp;
    
    // calculate the eigenproblem sensitivity
    _initialize_matrix_sensitivity_for_param(params,
                                             i,
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
    libMesh::ParameterVector param_V;
    param_V.resize(1);
    param_V[0]  =  _kr_param->ptr();
    
    
    // identify the sensitivity of solution to be used based on the
    // function arguments
    if (!dXdV) {
        sol_sens = zero_sol_sens.get();
    }
    else
        sol_sens = dXdV;
    
    _initialize_matrix_sensitivity_for_param(param_V,
                                             0,
                                             *sol_sens,
                                             root.V,
                                             mat_A_sens,
                                             mat_B_sens);
    
    // now calculate the quotient for sensitivity wrt k_red
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



