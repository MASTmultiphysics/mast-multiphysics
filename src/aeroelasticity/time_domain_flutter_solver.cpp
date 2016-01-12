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
#include "aeroelasticity/time_domain_flutter_solver.h"
#include "aeroelasticity/flutter_solution_base.h"
#include "aeroelasticity/time_domain_flutter_root_base.h"
#include "aeroelasticity/flutter_root_crossover_base.h"
#include "elasticity/structural_fluid_interaction_assembly.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "numerics/lapack_dggev_interface.h"


MAST::TimeDomainFlutterSolver::TimeDomainFlutterSolver() {
    
}



MAST::TimeDomainFlutterSolver::~TimeDomainFlutterSolver() {
    
    this->clear_solutions();
}




void
MAST::TimeDomainFlutterSolver::
initialize(MAST::StructuralFluidInteractionAssembly&   assembly,
           Real                                         V_lower,
           Real                                         V_upper,
           unsigned int                                 n_V_divs,
           std::vector<libMesh::NumericVector<Real> *>& basis) {
    
    _assembly       = &assembly;
    _basis_vectors  = &basis;
    _V_range.first  = V_lower;
    _V_range.second = V_upper;
    _n_V_divs       = n_V_divs;
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
    
    std::map<Real, MAST::FlutterRootCrossoverBase*>::const_iterator
    it = _flutter_crossovers.begin(),
    end = _flutter_crossovers.end();
    
    unsigned int n = 0;
    for ( ; it!=end; it++)
        if (it->second->root) // a valid root pointer has been assigned
            n++;
    
    return n;
}





const MAST::TimeDomainFlutterRootBase&
MAST::TimeDomainFlutterSolver::get_root(const unsigned int n) const {
    
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





std::pair<bool, const MAST::TimeDomainFlutterRootBase*>
MAST::TimeDomainFlutterSolver::find_next_root(const Real g_tol,
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
            return std::pair<bool, const MAST::TimeDomainFlutterRootBase*> (true, cross->root);
        }
        
        it++;
    }
    
    // if it gets here, no new root was found
    return std::pair<bool, MAST::TimeDomainFlutterRootBase*> (false, NULL);
}



std::pair<bool, const MAST::TimeDomainFlutterRootBase*>
MAST::TimeDomainFlutterSolver::find_critical_root(const Real g_tol,
                                                  const unsigned int n_bisection_iters)
{
    // iterate over the cross-over points and calculate the next that has
    // not been evaluated
    std::map<Real, MAST::FlutterRootCrossoverBase*>::iterator
    it = _flutter_crossovers.begin(), end = _flutter_crossovers.end();
    
    if (it == end) // no potential cross-over points were identified
        return std::pair<bool, MAST::TimeDomainFlutterRootBase*> (false, NULL);
    
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
    return std::pair<bool, const MAST::TimeDomainFlutterRootBase*> (true, it->second->root);
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
        
        MAST::FlutterSolutionBase* prev_sol = NULL;
        for (unsigned int i=0; i<_n_V_divs+1; i++) {
            current_V = V_vals[i];
            std::auto_ptr<MAST::FlutterSolutionBase> sol =
            _analyze(current_V, prev_sol);
            
            prev_sol = sol.get();
            
            sol->print(_output);
            
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
MAST::TimeDomainFlutterSolver::print_sorted_roots(std::ostream* output)
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
            const MAST::TimeDomainFlutterRootBase& root =
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
        const MAST::TimeDomainFlutterRootBase& root = this->get_root(i);
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
MAST::TimeDomainFlutterSolver::print_crossover_points(std::ostream* output)
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
    lower_g  = ref_sol_range.first->get_root(root_num).g,
    upper_V  = ref_sol_range.second->ref_val(),
    upper_g  = ref_sol_range.second->get_root(root_num).g,
    new_V    = 0.;
    unsigned int n_iters = 0;
    
    MAST::FlutterSolutionBase* new_sol = NULL;
    std::pair<bool, MAST::FlutterSolutionBase*> rval(false, NULL);
    
    while (n_iters < max_iters) {
        
        new_V    = lower_V +
        (upper_V-lower_V)/(upper_g-lower_g)*(0.-lower_g); // linear interpolation
        
        new_sol  = _analyze(new_V, ref_sol_range.first).release();
        
        new_sol->print(_output);
        
        // add the solution to this solver
        bool if_success =
        _flutter_solutions.insert(std::pair<Real, MAST::FlutterSolutionBase*>
                                  (new_V, new_sol)).second;
        
        libmesh_assert(if_success);
        
        const MAST::TimeDomainFlutterRootBase& root = new_sol->get_root(root_num);
        
        // check if the new damping value
        if (fabs(root.g) <= g_tol) {
            
            rval.first = true;
            rval.second = new_sol;
            return  rval;
        }
        
        // update the V value
        if (root.g < 0.) {
            
            lower_V = new_V;
            lower_g = root.g;
        }
        else {
            
            upper_V = new_V;
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
MAST::TimeDomainFlutterSolver::_analyze(const Real v_ref,
                                       const MAST::FlutterSolutionBase* prev_sol) {
    
    libMesh::out
    << " ====================================================" << std::endl
    << "Eigensolution" << std::endl
    << "   V_ref = " << std::setw(10) << v_ref << std::endl;
    
    
    RealMatrixX m, c, k;
    
    // initialize the matrices for the structure.
    _initialize_matrices(v_ref, m, c, k);
    
    // now use the matrices to create the matrices for first-order model
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
    
    const unsigned int
    nmodes =  (unsigned int)k.rows();

    RealMatrixX
    A      =  RealMatrixX::Zero(2*nmodes, 2*nmodes),
    B      =  RealMatrixX::Zero(2*nmodes, 2*nmodes);

    B.topLeftCorner(nmodes, nmodes)      = RealMatrixX::Identity(nmodes, nmodes);
    B.bottomRightCorner(nmodes, nmodes)  = m;
    
    
    A.topRightCorner(nmodes, nmodes)     = RealMatrixX::Identity(nmodes, nmodes);
    A.bottomLeftCorner(nmodes, nmodes)   = -k;
    A.bottomRightCorner(nmodes, nmodes)  = -c;
    
    
    MAST::LAPACK_DGGEV ges;
    ges.compute(A, B);
    ges.scale_eigenvectors_to_identity_innerproduct();
    for (unsigned int i=0; i<nmodes*2; i++)
        std::cout << ges.alphas()(i)/ges.betas()(i) << std::endl;
    std::cout << std::endl;

    
    MAST::FlutterSolutionBase* root = new MAST::FlutterSolutionBase;
    root->init(*this, v_ref, ges);
    if (prev_sol)
        root->sort(*prev_sol);
    
    libMesh::out
    << "Finished Eigensolution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return std::auto_ptr<MAST::FlutterSolutionBase> (root);
}




void
MAST::TimeDomainFlutterSolver::_initialize_matrices(Real U_inf,
                                                    RealMatrixX &m,
                                                    RealMatrixX &c,
                                                    RealMatrixX &k) {
    
    // get all the side and volume loads for this structure and set the velocity
    // for the piston theory boundary conditions
    
    // firs the side BC
    MAST::SideBCMapType& side_bc = _assembly->discipline().side_loads();
    MAST::SideBCMapType::iterator
    side_it  = side_bc.begin(),
    side_end = side_bc.end();
    
    for ( ; side_it != side_end; side_it++)
        if (side_it->second->type() == MAST::PISTON_THEORY)
            dynamic_cast<MAST::PistonTheoryBoundaryCondition*>(side_it->second)->set_U_inf(U_inf);

    // next, the volume BC
    MAST::VolumeBCMapType& vol_bc = _assembly->discipline().volume_loads();
    MAST::VolumeBCMapType::iterator
    vol_it  = vol_bc.begin(),
    vol_end = vol_bc.end();
    
    for ( ; vol_it != vol_end; vol_it++)
        if (vol_it->second->type() == MAST::PISTON_THEORY)
            dynamic_cast<MAST::PistonTheoryBoundaryCondition*>(vol_it->second)->set_U_inf(U_inf);


    // now prepare a map of the quantities and ask the assembly object to
    // calculate the quantities of interest.
    std::map<MAST::StructuralQuantityType, RealMatrixX*> qty_map;
    qty_map[MAST::MASS]       = &m;
    qty_map[MAST::DAMPING]    = &c;
    qty_map[MAST::STIFFNESS]  = &k;

    // zero the matrices to the size of the number of basis
    const unsigned int n = (unsigned int)_basis_vectors->size();
    m.setZero(n,n);
    c.setZero(n,n);
    k.setZero(n,n);
    
    _assembly->assemble_reduced_order_quantity(*_assembly->system().solution,
                                               *_basis_vectors,
                                               qty_map);
}




void
MAST::TimeDomainFlutterSolver::_identify_crossover_points() {
    
    /////////////////////////////////////////////////////////////////////
    // only non-negative frequencies (>=0) are used, since a time-domain
    // solver provides complex-conjugate roots.
    /////////////////////////////////////////////////////////////////////
    
    // if the initial scanning has not been done, then do it now
    const Real tol = 1.0e-5, max_allowable_g = 100.;
    
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
            
            g_val = fabs(sol_it->second->get_root(i).g);
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
             (fabs(sol_it->second->get_root(i).g) < max_allowable_g &&
              fabs(sol_itp1->second->get_root(i).g) < max_allowable_g) &&
             // if the mode has been identified to be trailing along g =0,
             // neglect it
             !modes_to_neglect[i]);
            

            if (if_process) {
                
                // look for the flutter roots
                MAST::FlutterSolutionBase
                *lower = sol_it->second,
                *upper = sol_itp1->second;
                
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
         (fabs(sol_it->second->get_root(i).g) < max_allowable_g &&
          fabs(sol_itp1->second->get_root(i).g) < max_allowable_g) &&
         // if the mode has been identified to be trailing along g =0,
         // neglect it
         !modes_to_neglect[i]);
        
        Real g_val = sol_it->second->get_root(i).g;
        
        if (if_process &&
            g_val > 0 &&
            g_val < max_allowable_g) {
            
            MAST::FlutterRootCrossoverBase* cross =
            new MAST::FlutterRootCrossoverBase;
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

