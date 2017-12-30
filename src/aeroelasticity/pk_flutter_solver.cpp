/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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
#include "aeroelasticity/pk_flutter_solver.h"
#include "aeroelasticity/pk_flutter_solution.h"
#include "aeroelasticity/pk_flutter_root.h"
#include "aeroelasticity/pk_flutter_root_crossover.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "numerics/lapack_zggev_interface.h"
#include "base/parameter.h"


MAST::PKFlutterSolver::PKFlutterSolver():
MAST::FlutterSolverBase(),
_velocity_param(nullptr),
_kred_param(nullptr),
_bref_param(nullptr),
_V_range(std::pair<Real, Real>(0., 0.)),
_n_V_divs(0),
_kr_range(std::pair<Real, Real>(0., 0.)),
_n_k_red_divs(0)
{ }



MAST::PKFlutterSolver::~PKFlutterSolver() {

    this->clear();
}



void
MAST::PKFlutterSolver::clear() {
    
    this->clear_solutions();
    
    _velocity_param   = nullptr;
    _kred_param       = nullptr;
    _bref_param       = nullptr;
    _V_range          = std::pair<Real, Real>(0.,0.);
    _kr_range         = std::pair<Real, Real>(0.,0.);
    _n_V_divs         = 0;
    _n_k_red_divs     = 0;
    
    MAST::FlutterSolverBase::clear();
}



void
MAST::PKFlutterSolver::
initialize(MAST::Parameter&                             V_param,
           MAST::Parameter&                             kr_param,
           MAST::Parameter&                             bref_param,
           Real                                         rho,
           Real                                         V_lower,
           Real                                         V_upper,
           unsigned int                                 n_V_divs,
           Real                                         kr_lower,
           Real                                         kr_upper,
           unsigned int                                 n_kr_divs,
           std::vector<libMesh::NumericVector<Real>*>& basis) {
    
    
    _velocity_param     = &V_param;
    _kred_param         = &kr_param;
    _bref_param         = &bref_param;
    _rho                = rho;
    _V_range.first      = V_lower;
    _V_range.second     = V_upper;
    _n_V_divs           = n_V_divs;
    _kr_range.first     = kr_lower;
    _kr_range.second    = kr_upper;
    _n_k_red_divs       = n_kr_divs;
    
    MAST::FlutterSolverBase::initialize(basis);

}




void
MAST::PKFlutterSolver::clear_solutions() {
    
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
MAST::PKFlutterSolver::n_roots_found() const {
    
    std::multimap<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int n = 0;
    for ( ; it!=end; it++)
        if (it->second->root) // a valid root pointer has been assigned
            n++;
    
    return n;
}



void
MAST::PKFlutterSolver::scan_for_roots() {
    
    // if the initial scanning has not been done, then do it now
    if (!_flutter_solutions.size()) {
        
        // the outer loop consists of the reference velocity at which
        // aerodynamics is calculated
        Real current_k_red = _kr_range.second,
        delta_k_red = (_kr_range.second-_kr_range.first)/_n_k_red_divs;
        
        std::vector<Real> k_red_vals(_n_k_red_divs+1);
        for (unsigned int i=0; i<_n_k_red_divs+1; i++) {
            k_red_vals[i] = current_k_red;
            current_k_red -= delta_k_red;
        }
        k_red_vals[_n_k_red_divs] = _kr_range.first; // to get around finite-precision arithmetic
        
        //
        //  outer loop is on reduced frequency
        //
        for (unsigned int j=0; j<_n_k_red_divs+1; j++) {
            
            current_k_red = k_red_vals[j];
            
            // march from the upper limit to the lower to find the roots
            Real current_v_ref = _V_range.first,
            delta_v_ref = (_V_range.second-_V_range.first)/_n_V_divs;
            
            std::vector<Real> v_ref_vals(_n_V_divs+1);
            for (unsigned int i=0; i<_n_V_divs+1; i++) {
                v_ref_vals[i] = current_v_ref;
                current_v_ref += delta_v_ref;
            }
            v_ref_vals[_n_V_divs] = _V_range.second; // to get around finite-precision arithmetic
            
            MAST::FlutterSolutionBase* prev_sol = nullptr;
            
            //
            // inner loop is on reduced frequencies
            //
            for (unsigned int i=0; i<_n_V_divs+1; i++) {
                current_v_ref = v_ref_vals[i];
                std::unique_ptr<MAST::FlutterSolutionBase> sol =
                _analyze(current_k_red,
                         current_v_ref,
                         prev_sol);
                
                
                if (_output)
                    sol->print(*_output);
                
                // add the solution to this solver
                _insert_new_solution(current_k_red, sol.release());
                
                // now get a pointer to the previous solution
                // get the solution from the database for this reduced frequency
                std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
                _flutter_solutions.find(current_v_ref);
                
                libmesh_assert(it != _flutter_solutions.end());
                prev_sol = it->second;
            }
            
        }
        _identify_crossover_points();
    }
}



void
MAST::PKFlutterSolver::print_sorted_roots()
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
MAST::PKFlutterSolver::print_crossover_points()
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
MAST::PKFlutterSolver::_bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                                        MAST::FlutterSolutionBase*>& ref_sol_range,
                                        const unsigned int root_num,
                                        const Real g_tol,
                                        const unsigned int max_iters) {
    
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real
    lower_k = ref_sol_range.first->get_root(root_num).kr,
    lower_v = ref_sol_range.first->get_root(root_num).V,
    lower_g = ref_sol_range.first->get_root(root_num).g,
    upper_k = ref_sol_range.second->get_root(root_num).kr,
    upper_v = ref_sol_range.second->get_root(root_num).V,
    upper_g = ref_sol_range.second->get_root(root_num).g,
    new_k   = lower_k + (upper_k-lower_k)/(upper_g-lower_g)*(0.-lower_g),
    new_v   = 0.; // linear interpolation
    unsigned int n_iters = 0;
    
    std::unique_ptr<MAST::FlutterSolutionBase> new_sol;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, nullptr);
    
    while (n_iters < max_iters) {
        
        new_v = lower_v + (upper_v-lower_v)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        
        new_sol.reset(_analyze(new_k,
                               new_v,
                               ref_sol_range.first).release());
        
        if (_output)
            new_sol->print(*_output);
        
        // add the solution to this solver
        _insert_new_solution(new_k, new_sol.release());
        
        // get the solution from the database for this reduced frequency
        std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
        _flutter_solutions.find(new_v);
        
        libmesh_assert(it != _flutter_solutions.end());
        rval.second = it->second;
        const MAST::FlutterRootBase& root = rval.second->get_root(root_num);
        
        // use the estimated flutter velocity to get the next
        // V_ref value for aerodynamic matrices.
        new_k = root.kr;
        
        // check if the new damping value
        if (fabs(root.g) <= g_tol) {
            rval.first = true;
            return  rval;
        }
        
        // update the v_val
        if (root.g < 0.) {
            lower_v = new_v;
            lower_g = root.g;
        }
        else {
            upper_v = new_v;
            upper_g = root.g;
        }
        
        n_iters++;
    }
    
    // return false, along with the latest sol
    rval.first = false;
    
    return rval;
}




const MAST::FlutterRootBase&
MAST::PKFlutterSolver::get_root(const unsigned int n) const {
    
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
MAST::PKFlutterSolver::find_next_root(const Real g_tol,
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
MAST::PKFlutterSolver::find_critical_root(const Real g_tol,
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




/*std::pair<bool, MAST::FlutterSolutionBase*>
MAST::PKFlutterSolver::_newton_search(const MAST::FlutterSolutionBase& init_sol,
                                     const unsigned int root_num,
                                     const Real tol,
                                     const unsigned int max_iters) {
    
    libmesh_error(); // to be implemented
    
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, nullptr);
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real k_red, v_ref;
    unsigned int n_iters = 0;
    
    std::unique_ptr<MAST::FlutterSolutionBase> new_sol;
    
    RealVectorX res, sol, dsol;
    RealMatrixX jac, stiff;
    ComplexMatrixX mat_A, mat_B, mat_A_sens, mat_B_sens;
    ComplexVectorX v;
    
    res.resize(2); sol.resize(2); dsol.resize(2);
    jac.resize(2,2);
    
    // initialize the solution to the values of init_sol
    k_red  = init_sol.get_root(root_num).k_red_ref;
    v_ref  = init_sol.get_root(root_num).V_ref;
    sol(0) = k_red;
    sol(1) = v_ref;
    
    const MAST::FlutterSolutionBase* prev_sol = &init_sol;
    
    bool if_continue = true;
    
    while (if_continue) {
        
        // evaluate the residual and Jacobians
        std::unique_ptr<MAST::FlutterSolutionBase> pk_sol =
        this->analyze(k_red, v_ref, prev_sol);
        
        pk_sol->print(_output);
        
        // add the solution to this solver
        _insert_new_solution(v_ref, pk_sol.release());
        
        // now get a pointer to the previous solution
        // get the solution from the database for this reduced frequency
        std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
        _flutter_solutions.find(v_ref);
        
        libmesh_assert(it != _flutter_solutions.end());
        rval.second = it->second;
        prev_sol = it->second;
        
        // solve the Newton update problem
        const MAST::FlutterRootBase& root = prev_sol->get_root(root_num);
        Complex eig = root.root, eig_k_red_sens = 0., den = 0., eig_V_ref_sens = 0.;
        
        // initialize the baseline matrices
        initialize_matrices(k_red, v_ref, mat_A, mat_B, stiff);
        
        // solve the sensitivity problem
        // first with respect to k_red
        // next we need the sensitivity of k_red before we can calculate
        // the sensitivity of flutter eigenvalue
        initialize_matrix_sensitivity_for_reduced_freq(k_red,
                                                       v_ref,
                                                       mat_A_sens,
                                                       mat_B_sens);
        
        // now calculate the quotient for sensitivity wrt k_red
        // calculate numerator
        mat_B_sens *= -eig;
        mat_B_sens += mat_A_sens;
        v = mat_B_sens*root.eig_vec_right;
        den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
        eig_k_red_sens = root.eig_vec_left.dot(v) / den;
        
        //std::unique_ptr<MAST::FlutterSolutionBase> PK_dsol =
        //this->analyze(k_red+.001, v_ref, prev_sol);
        //eig -= PK_dsol->get_root(root_num).root;
        //eig /= -.001;
         
        //std::cout << eig << std::endl;
        //std::cout << eig_k_red_sens << std::endl;
        
        // next, sensitivity wrt V_ref
        initialize_matrix_sensitivity_for_V_ref(k_red,
                                                v_ref,
                                                mat_A_sens,
                                                mat_B_sens);
        
        // now calculate the quotient for sensitivity wrt V_ref
        // calculate numerator
        mat_B_sens *= -eig;
        mat_B_sens += mat_A_sens;
        v = mat_B_sens*root.eig_vec_right;
        den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
        eig_V_ref_sens = root.eig_vec_left.dot(v) / den;
        
 //std::unique_ptr<MAST::FlutterSolutionBase> PK_dsol =
 //this->analyze(k_red, v_ref+.001, prev_sol);
 //eig -= PK_dsol->get_root(root_num).root;
 //eig /= -.001;
 
 //std::cout << eig << std::endl;
 //std::cout << eig_V_ref_sens << std::endl;
 
        // residual
        res(0) = root.root.imag();
        res(1) = v_ref*std::sqrt(root.root.real()) - 1.;
        
        // Jacobian
        jac(0,0) = eig_k_red_sens.imag();
        jac(0,1) = eig_V_ref_sens.imag();
        jac(1,0) = 0.5*v_ref*pow(eig.real(), -0.5)*eig_k_red_sens.real();
        jac(1,1) = std::sqrt(eig.real()) + 0.5*v_ref*pow(eig.real(), -0.5)*eig_V_ref_sens.real();
        
        // now calculate the updates
        //     r0 + J *dx = 0
        // =>  dx = - inv(J) * r0
        // =>  x1 = x0 + dx
        
        dsol = -jac.inverse()*res;
        sol += dsol;
        
        // get the updated parameter values
        k_red = sol(0);
        v_ref = sol(1);
        
        // increment the iteration counter
        n_iters++;
        
        // set the flag
        if (fabs(root.g) < tol) {
            rval.first = true;
            return rval;
        }
        
        if (n_iters >= max_iters)
            if_continue = false;
    }
    
    // return false, along with the latest sol
    rval.first = false;
    
    return rval;
}
*/


void
MAST::PKFlutterSolver::_insert_new_solution(const Real k_red,
                                            MAST::FlutterSolutionBase* sol) {
    
    const Real v_ref = sol->ref_val();
    
    // if the value does not already exist, then insert it in the map
    std::map<Real, MAST::FlutterSolutionBase*>::iterator it =
    _flutter_solutions.find(v_ref);
    
    if ( it == _flutter_solutions.end()) {
        bool if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (v_ref, sol)).second;
        
        libmesh_assert(if_success);
        
        // do not delete the solution since it is stored in the solver
    }
    else {
        // assuming that the roots are sorted, compare them with the existing
        // root and swap the ones that are closer to v_ref
        
        // first make sure that the two have the same number of roots
        libmesh_assert_equal_to(it->second->n_roots(), sol->n_roots());
        
        Real dk_old = 0., k_ref_old = 0., dk_new = 0.;
        MAST::FlutterRootBase *old_root, *new_root;
        
        for (unsigned int i=0; i<sol->n_roots(); i++) {
            old_root  = &(it->second->get_root(i));
            new_root  = &(sol->get_root(i));
            k_ref_old = old_root->kr;
            dk_old    = fabs(fabs(old_root->kr) - k_ref_old);
            dk_new    = fabs(fabs(new_root->kr) - k_red);
            
            if (dk_new < dk_old)
                old_root->copy_root(*new_root);
        }
        
        // delete the solution since its relevants roots have been swapped
        // out
        delete sol;
    }
}





void MAST::PKFlutterSolver::_identify_crossover_points()
{
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
                    new MAST::PKFlutterRootCrossover;
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
        sol_rit    = _flutter_solutions.rbegin(), // first of the pair
        sol_ritp1    = _flutter_solutions.rbegin(), // first of the pair
        sol_rend   = _flutter_solutions.rend();
        if (sol_rit == sol_rend)
            return;
        
        sol_ritp1++; // increment for the next pair of results
        while (sol_ritp1 != sol_rend) {
            // do not use k_red = 0, or if the root is invalid
            if (sol_rit->second->get_root(i).if_nonphysical_root ||
                sol_ritp1->second->get_root(i).if_nonphysical_root ||
                fabs(sol_rit->second->ref_val()) < tol ||
                fabs(sol_ritp1->second->ref_val()) < tol ||
                fabs(sol_rit->second->get_root(i).g) > max_allowable_g ||
                fabs(sol_ritp1->second->get_root(i).g) > max_allowable_g) {
                // do nothing
            }
            else if (!modes_to_neglect[i]) { // look for the flutter roots
                MAST::FlutterSolutionBase *lower = sol_rit->second,
                *upper = sol_ritp1->second;
                
                if ((lower->get_root(i).g <= 0.) &&
                    (upper->get_root(i).g > 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::PKFlutterRootCrossover;
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
                    new MAST::PKFlutterRootCrossover;
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
    }
}



std::unique_ptr<MAST::FlutterSolutionBase>
MAST::PKFlutterSolver::_analyze(const Real k_red,
                                const Real v_ref,
                                const MAST::FlutterSolutionBase* prev_sol) {
    // solve the eigenproblem  L x = lambda R x
    ComplexMatrixX R, L;
    RealMatrixX stiff;
    
    libMesh::out
    << " ====================================================" << std::endl
    << "PK Solution" << std::endl
    << "   k_red = " << std::setw(10) << k_red << std::endl
    << "   V_ref = " << std::setw(10) << v_ref << std::endl;
    
    _initialize_matrices(k_red, v_ref, L, R, stiff);
    LAPACK_ZGGEV ges;
    ges.compute(L, R);
    ges.scale_eigenvectors_to_identity_innerproduct();
    
    MAST::PKFlutterSolution* root = new MAST::PKFlutterSolution;
    root->init(*this,
               k_red, v_ref,
               (*_bref_param)(),
               stiff, ges);
    if (prev_sol)
        root->sort(*prev_sol);
    
    libMesh::out
    << "Finished PK Solution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return std::unique_ptr<MAST::FlutterSolutionBase> (root);
}




void
MAST::PKFlutterSolver::calculate_sensitivity(MAST::FlutterRootBase& root,
                                             const libMesh::ParameterVector& params,
                                             const unsigned int i) {

    /*
    libMesh::out
    << " ====================================================" << std::endl
    << "PK Sensitivity Solution" << std::endl
    << "   k_red = " << std::setw(10) << root.kr << std::endl
    << "   V_ref = " << std::setw(10) << root.V << std::endl;
    
    Complex eig = root.root, sens = 0., k_sens = 0., den = 0.;
    Real par_g_par_alpha = 0., par_g_par_kref = 0., par_k_par_alpha = 0.,
    V_sens=0.;
    
    // get the sensitivity of the matrices
    ComplexMatrixX mat_A, mat_B, mat_A_sens, mat_B_sens;
    ComplexVectorX v;
    RealMatrixX stiff;
    
    // initialize the baseline matrices
    _initialize_matrices(root.kr, root.V, mat_A, mat_B, stiff);
    
    // calculate the eigenproblem sensitivity
    _initialize_matrix_sensitivity_for_param(params, i,
                                             root.kr,
                                             root.V,
                                             mat_A_sens,
                                             mat_B_sens);
    
    // the eigenproblem is     A x - lambda B x = 0
    // therefore, the denominator is obtained from the inner product of
    // x^T B x
    // sensitivity is
    //   -dlambda/dp x^T B x = - x^T (dA/dp - lambda dB/dp)
    // or
    //   dlambda/dp = [x^T (dA/dp - lambda dB/dp)]/(x^T B x)
    
    // now calculate the quotient for sensitivity
    // numerator =  ( dA/dp - lambda dB/dp)
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
    sens = root.eig_vec_left.dot(v)/den;
    
    // now add the correction from sensitivity of g(k) = 0
    par_g_par_alpha =
    sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2) * sens.real();
    
    
    // next we need the sensitivity of k_red before we can calculate
    // the sensitivity of flutter eigenvalue
    _initialize_matrix_sensitivity_for_reduced_freq(root.kr,
                                                    root.V,
                                                    mat_A_sens,
                                                    mat_B_sens);
    
    // now calculate the quotient for sensitivity wrt k_red
    // calculate numerator
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    k_sens = root.eig_vec_left.dot(v) / den;
    
    // use this to calculate the partial derivative of g wrt k_red
    par_g_par_kref =
    k_sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2) * k_sens.real();
    
    // use this to calculate the sensitivity of k_red wrt alpha
    par_k_par_alpha = -par_g_par_alpha / par_g_par_kref;
    
    // finally add the correction to the flutter sensitivity
    sens += k_sens * par_k_par_alpha;
    
    // finally, the flutter speed sensitivity
    V_sens = -.5*sens.real()/pow(eig.real(), 1.5);
    
    // set value in the return root
    root.has_sensitivity_data = true;
    root.root_sens  = sens;
    root.V_sens     = V_sens;
    
    libMesh::out
    << "Finished PK Sensitivity Solution" << std::endl
    << " ====================================================" << std::endl;
 */
}




void
MAST::PKFlutterSolver::_initialize_matrices(const Real k_red,
                                            const Real v_ref,
                                            ComplexMatrixX& A, // stiff, aero, damp
                                            ComplexMatrixX& B, // mass
                                            RealMatrixX& stiff)// stiffness
{
    // the PK method equations are
    //
    //   p [ I  0 ] {  X } =  [ 0      I ] {  X }
    //     [ 0  M ] { pX }    [-K-qA   -C] { pX }
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
    (*_kred_param)      = k_red;
    (*_velocity_param)  = v_ref;
    
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

    A.topRightCorner    (n, n)    =  ComplexMatrixX::Identity(n, n);
    A.bottomLeftCorner  (n, n)    = -k.cast<Complex>() + _rho/2.*v_ref*v_ref*a;
    B.topLeftCorner     (n, n)    = ComplexMatrixX::Identity(n, n);
    B.bottomRightCorner (n, n)    = m.cast<Complex>();
    
}



void
MAST::PKFlutterSolver::
_initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                         unsigned int p,
                                         const Real k_red,
                                         const Real v_ref,
                                         ComplexMatrixX& L,   // stiff, aero, damp
                                         ComplexMatrixX& R) { // mass
    

    libmesh_error(); // to be implemented
}



