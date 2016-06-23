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

/*
// MAST includes
#include "aeroelasticity/pk_flutter_solver.h"
#include "aeroelasticity/pk_flutter_solution.h"
#include "aeroelasticity/pk_flutter_root.h"
#include "aeroelasticity/pk_flutter_root_crossover.h"


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
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
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
            
            MAST::FlutterSolutionBase* prev_sol = NULL;
            
            //
            // inner loop is on reduced frequencies
            //
            for (unsigned int i=0; i<_n_V_divs+1; i++) {
                current_v_ref = v_ref_vals[i];
                std::auto_ptr<MAST::FlutterSolutionBase> sol =
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



std::pair<bool, MAST::FlutterSolutionBase*>
MAST::PKFlutterSolver::bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                                        MAST::FlutterSolutionBase*>& ref_sol_range,
                                        const unsigned int root_num,
                                        const Real g_tol,
                                        const unsigned int max_iters) {
    
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real
    lower_k = ref_sol_range.first->get_root(root_num).k_red,
    lower_v = ref_sol_range.first->get_root(root_num).V,
    lower_g = ref_sol_range.first->get_root(root_num).g,
    upper_k = ref_sol_range.second->get_root(root_num).k_red,
    upper_v = ref_sol_range.second->get_root(root_num).V,
    upper_g = ref_sol_range.second->get_root(root_num).g,
    new_k   = lower_k + (upper_k-lower_k)/(upper_g-lower_g)*(0.-lower_g),
    new_v   = 0.; // linear interpolation
    unsigned int n_iters = 0;
    
    std::auto_ptr<MAST::FlutterSolutionBase> new_sol;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, NULL);
    
    while (n_iters < max_iters) {
        
        new_v = lower_v + (upper_v-lower_v)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        
        new_sol.reset(analyze(new_k,
                              new_v,
                              ref_sol_range.first).release());
        
        new_sol->print(_output);
        
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
        new_k = root.k_red;
        
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



std::pair<bool, MAST::FlutterSolutionBase*>
MAST::PKFlutterSolver::newton_search(const MAST::FlutterSolutionBase& init_sol,
                                     const unsigned int root_num,
                                     const Real tol,
                                     const unsigned int max_iters) {
    
    libmesh_error(); // to be implemented
    
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, NULL);
    // assumes that the upper k_val has +ve g val and lower k_val has -ve
    // k_val
    Real k_red, v_ref;
    unsigned int n_iters = 0;
    
    std::auto_ptr<MAST::FlutterSolutionBase> new_sol;
    
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
        std::auto_ptr<MAST::FlutterSolutionBase> pk_sol =
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
        
        //std::auto_ptr<MAST::FlutterSolutionBase> PK_dsol =
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
        
 //std::auto_ptr<MAST::FlutterSolutionBase> PK_dsol =
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
            k_ref_old = old_root->k_red_ref;
            dk_old    = fabs(fabs(old_root->k_red) - k_ref_old);
            dk_new    = fabs(fabs(new_root->k_red) - k_red);
            
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
                    new MAST::FlutterRootCrossoverBase;
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
                    new MAST::FlutterRootCrossoverBase;
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
                    new MAST::FlutterRootCrossoverBase;
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



std::auto_ptr<MAST::FlutterSolutionBase>
MAST::PKFlutterSolver::analyze(const Real k_red,
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
    
    initialize_matrices(k_red, v_ref, L, R, stiff);
    LAPACK_ZGGEV ges;
    ges.compute(L, R);
    ges.scale_eigenvectors_to_identity_innerproduct();
    
    MAST::PKFlutterSolution* root = new MAST::PKFlutterSolution;
    root->init(*this,
               k_red, v_ref,
               flight_condition->ref_chord,
               stiff, ges);
    if (prev_sol)
        root->sort(*prev_sol);
    
    libMesh::out
    << "Finished PK Solution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return std::auto_ptr<MAST::FlutterSolutionBase> (root);
}




void
MAST::PKFlutterSolver::calculate_sensitivity(MAST::FlutterRootBase& root,
                                             const libMesh::ParameterVector& params,
                                             const unsigned int i) {
    // make sure that the aero_structural_model is a valid pointer
    libmesh_assert(aero_structural_model);
    
    libMesh::out
    << " ====================================================" << std::endl
    << "PK Sensitivity Solution" << std::endl
    << "   k_red = " << std::setw(10) << root.k_red << std::endl
    << "   V_ref = " << std::setw(10) << root.V << std::endl;
    
    Complex eig = root.root, sens = 0., k_sens = 0., den = 0.;
    Real par_g_par_alpha = 0., par_g_par_kref = 0., par_k_par_alpha = 0.,
    V_sens=0.;
    
    // get the sensitivity of the matrices
    ComplexMatrixX mat_A, mat_B, mat_A_sens, mat_B_sens;
    ComplexVectorX v;
    RealMatrixX stiff;
    
    // initialize the baseline matrices
    initialize_matrices(root.k_red_ref, root.V_ref, mat_A, mat_B, stiff);
    
    // calculate the eigenproblem sensitivity
    initialize_matrix_sensitivity_for_param(params, i,
                                            root.k_red_ref,
                                            root.V_ref,
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
    initialize_matrix_sensitivity_for_reduced_freq(root.k_red_ref,
                                                   root.V_ref,
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
    
}




void
MAST::PKFlutterSolver::initialize_matrices(const Real k_red,
                                           const Real v_ref,
                                           ComplexMatrixX& L, // stiff, aero, damp
                                           ComplexMatrixX& R, // mass
                                           RealMatrixX& stiff)// stiffness
{
    
    bool has_matrix = false;
    RealMatrixX mat_r;
    ComplexMatrixX mat_c;
    
    const unsigned int n_dofs = aero_structural_model->n_dofs();
    L.setZero(2*n_dofs, 2*n_dofs);
    R.setZero(2*n_dofs, 2*n_dofs);
    
    
    // stiffness matrix
    has_matrix = aero_structural_model->get_structural_stiffness_matrix(stiff);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++) {
        L(i,n_dofs+i) = 1.;
        R(i,i) = 1.;
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) = -stiff(i,j);
    }
    
    // damping matrix
    has_matrix = aero_structural_model->get_structural_damping_matrix(mat_r);
    
    if (has_matrix) {
        for (unsigned int i=0; i<n_dofs; i++)
            for (unsigned int j=0; j<n_dofs; j++)
                L(n_dofs+i,n_dofs+j) = -mat_r(i,j);
    }
    
    
    // mass matrix
    has_matrix = aero_structural_model->get_structural_mass_matrix(mat_r);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            R(n_dofs+i,n_dofs+j) = mat_r(i,j);
    
    // aerodynamic operator matrix
    has_matrix =
    aero_structural_model->get_aero_operator_matrix(k_red, v_ref, mat_c);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) +=
            0.5 * flight_condition->gas_property.rho * pow(v_ref,2) *
            mat_c(i,j);
}



void
MAST::PKFlutterSolver::
initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                        unsigned int p,
                                        const Real k_red,
                                        const Real v_ref,
                                        ComplexMatrixX& L,   // stiff, aero, damp
                                        ComplexMatrixX& R) { // mass
    
    bool has_matrix = false;
    RealMatrixX mat_r;
    ComplexMatrixX mat_c;
    
    const unsigned int n_dofs = aero_structural_model->n_dofs();
    L.setZero(2*n_dofs, 2*n_dofs);
    R.setZero(2*n_dofs, 2*n_dofs);
    
    
    // mass matrix sensitivity
    has_matrix =
    aero_structural_model->get_structural_stiffness_matrix_sensitivity(params,
                                                                       p,
                                                                       mat_r);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++) {
        L(i,n_dofs+i) = 1.;
        R(i,i) = 1.;
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) = -mat_r(i,j);
    }
    
    
    // mass matrix sensitivity
    has_matrix =
    aero_structural_model->get_structural_damping_matrix_sensitivity(params,
                                                                     p,
                                                                     mat_r);
    
    if (has_matrix) {
        for (unsigned int i=0; i<n_dofs; i++)
            for (unsigned int j=0; j<n_dofs; j++)
                L(n_dofs+i,n_dofs+j) = -mat_r(i,j);
    }
    
    
    // mass matrix sensitivity
    has_matrix =
    aero_structural_model->get_structural_mass_matrix_sensitivity(params,
                                                                  p,
                                                                  mat_r);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            R(n_dofs+i,n_dofs+j) = mat_r(i,j);
    
    
    // aerodynamic operator matrix sensitivity
    has_matrix =
    aero_structural_model->get_aero_operator_matrix_sensitivity(params,
                                                                p,
                                                                k_red, v_ref,
                                                                mat_c);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) +=
            0.5 * flight_condition->gas_property.rho * pow(v_ref,2) *
            mat_c(i,j);
    
}




void
MAST::PKFlutterSolver::
initialize_matrix_sensitivity_for_reduced_freq(const Real k_red,
                                               const Real v_ref,
                                               ComplexMatrixX& L,   // stiff, aero, damp
                                               ComplexMatrixX& R) { // mass
    bool has_matrix = false;
    RealMatrixX mat_r;
    ComplexMatrixX mat_c;
    
    const unsigned int n_dofs = aero_structural_model->n_dofs();
    
    L.setZero(2*n_dofs, 2*n_dofs);
    R.setZero(2*n_dofs, 2*n_dofs);
    
    // aerodynamic operator matrix sensitivity
    has_matrix =
    aero_structural_model->get_aero_operator_matrix_sensitivity_for_reduced_freq(k_red,
                                                                                 v_ref,
                                                                                 mat_c);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) +=
            0.5 * flight_condition->gas_property.rho * pow(v_ref,2) *
            mat_c(i,j);
}



void
MAST::PKFlutterSolver::
initialize_matrix_sensitivity_for_V_ref(const Real k_red,
                                        const Real v_ref,
                                        ComplexMatrixX& L,   // stiff, aero, damp
                                        ComplexMatrixX& R) { // mass
    bool has_matrix = false;
    RealMatrixX mat_r;
    ComplexMatrixX mat_c;
    
    const unsigned int n_dofs = aero_structural_model->n_dofs();
    
    L.setZero(2*n_dofs, 2*n_dofs);
    R.setZero(2*n_dofs, 2*n_dofs);
    
    // aerodynamic operator matrix
    has_matrix =
    aero_structural_model->get_aero_operator_matrix(k_red, v_ref, mat_c);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) +=
            flight_condition->gas_property.rho * v_ref * mat_c(i,j);
    
    
    // aerodynamic operator matrix sensitivity
    has_matrix =
    aero_structural_model->get_aero_operator_matrix_sensitivity_for_V_ref(k_red,
                                                                          v_ref,
                                                                          mat_c);
    libmesh_assert(has_matrix);
    
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<n_dofs; j++)
            L(n_dofs+i,j) +=
            0.5 * flight_condition->gas_property.rho * pow(v_ref,2) *
            mat_c(i,j);
}


*/