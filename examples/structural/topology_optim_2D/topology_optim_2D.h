///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2018  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//#ifndef __mast_topology_optim_2d_optimization_h__
//#define __mast_topology_optim_2d_optimization_h__
//
//// C++ includes
//#include <memory>
//#include <map>
//
//// MAST includes
//#include "base/field_function_base.h"
//#include "base/physics_discipline_base.h"
//#include "elasticity/structural_system_initialization.h"
//#include "property_cards/isotropic_material_property_card.h"
//#include "property_cards/solid_2d_section_element_property_card.h"
//#include "base/parameter.h"
//#include "base/constant_field_function.h"
//#include "optimization/function_evaluation.h"
//#include "boundary_condition/dirichlet_boundary_condition.h"
//
//
//// libMesh includes
//#include "libmesh/libmesh.h"
//#include "libmesh/equation_systems.h"
//#include "libmesh/serial_mesh.h"
//#include "libmesh/mesh_generation.h"
//#include "libmesh/nonlinear_implicit_system.h"
//#include "libmesh/fe_type.h"
//#include "libmesh/dof_map.h"
//#include "libmesh/mesh_function.h"
//#include "libmesh/parameter_vector.h"
//#include "libmesh/getpot.h"
//
//
//// get this from the global namespace
//extern libMesh::LibMeshInit* _init;
//
//
//namespace MAST {
//    
//    
//    // Forward declerations
//    class StructuralSystemInitialization;
//    class Parameter;
//    class PhysicsDisciplineBase;
//    class ConstantFieldFunction;
//    class IsotropicMaterialPropertyCard;
//    class Solid2DSectionElementPropertyCard;
//    class DirichletBoundaryCondition;
//    class BoundaryConditionBase;
//    class NonlinearImplicitAssembly;
//    class StructuralNonlinearAssemblyElemOperations;
//    class RealOutputFunction;
//    
//    
//    /*!
//     *   This class provides the Youngs Modulus that is scaled by the 
//     *   density.
//     */
//    class YoungsModulus: public MAST::FieldFunction<Real> {
//    public:
//        YoungsModulus(const std::string& nm,
//                      MAST::FieldFunction<Real> &rho,
//                      const Real base_modulus,
//                      const Real penalty);
//                
//        virtual ~YoungsModulus();
//        
//    protected:
//        
//        MAST::FieldFunction<Real> &_rho;
//        
//        Real _base_modulus;
//        
//        Real _penalty;
//        
//    public:
//        
//        virtual void operator() (const libMesh::Point& p, Real t, Real& v) const;
//        
//        virtual void derivative(   const MAST::FunctionBase& f,
//                                const libMesh::Point& p,
//                                Real t,
//                                Real& v) const;
//        
//    };
//
//    
//    
//    struct TopologyOptimization2D:
//    public MAST::FunctionEvaluation {
//        
//        
//        TopologyOptimization2D(const libMesh::Parallel::Communicator& comm);
//        
//        
//        ~TopologyOptimization2D();
//        
//        /*!
//         *   initializes the object for specified characteristics
//         */
//        void init(GetPot& infile, libMesh::ElemType e_type, bool if_vk);
//        
//        
//        /*!
//         *   initialize the design variables values and bounds
//         */
//        virtual void init_dvar(std::vector<Real>& x,
//                               std::vector<Real>& xmin,
//                               std::vector<Real>& xmax);
//        
//        
//        /*!
//         *    the core routine that performs the function evaluations
//         */
//        virtual void evaluate(const std::vector<Real>& dvars,
//                              Real& obj,
//                              bool eval_obj_grad,
//                              std::vector<Real>& obj_grad,
//                              std::vector<Real>& fvals,
//                              std::vector<bool>& eval_grads,
//                              std::vector<Real>& grads);
//        
//        /*!
//         *   customized output
//         */
//        virtual void output(unsigned int iter,
//                            const std::vector<Real>& x,
//                            Real obj,
//                            const std::vector<Real>& fval,
//                            bool if_write_to_optim_file) const;
//        
//        /*!
//         *  @returns a pointer to the function that evaluates the objective
//         */
//        virtual MAST::FunctionEvaluation::funobj
//        get_objective_evaluation_function();
//        
//        
//        /*!
//         *  @returns a pointer to the function that evaluates the constraint
//         */
//        virtual MAST::FunctionEvaluation::funcon
//        get_constraint_evaluation_function();
//        
//
//        
//        bool _initialized;
//        
//        // length of domain
//        Real _length;
//        
//        // width of domain
//        Real _width;
//        
//        // penalty value for DVs
//        Real _penalty;
//        
//        // volume fraction
//        Real _volume_fraction;
//        
//        
//        // number of elements and number of stations at which DVs are defined
//        unsigned int
//        _n_divs_x,
//        _n_divs_y,
//        _n_elems;
//        
//        // create the mesh
//        libMesh::SerialMesh*           _mesh;
//        
//        // create the equation system
//        libMesh::EquationSystems*      _eq_sys;
//        
//        // create the libmesh system
//        MAST::NonlinearSystem*  _sys;
//        
//        // system to plot the density variable over the mesh
//        libMesh::ExplicitSystem*           _rho_sys;
//        
//        // initialize the system to the right set of variables
//        MAST::StructuralSystemInitialization* _structural_sys;
//        MAST::PhysicsDisciplineBase*           _discipline;
//        
//        // nonlinear assembly object
//        MAST::NonlinearImplicitAssembly                 *_assembly;
//        MAST::StructuralNonlinearAssemblyElemOperations *_elem_ops;
//        
//        // create the property functions and add them to the
//        MAST::Parameter
//        *_nu,
//        *_kappa,
//        *_press,
//        *_th,
//        *_zero;
//        
//        MAST::ConstantFieldFunction
//        *_nu_f,
//        *_kappa_f,
//        *_th_f,
//        *_hoff_f,
//        *_press_f;
//
//        
//        // vector of elements that is in the same sequence as the DVs. This is
//        // used to map the DV vector to the density parameters.
//        std::vector<const libMesh::Elem*>           _elems;
//
//        // element density parameters
//        std::map<const libMesh::Elem*, MAST::Parameter*>    _elem_rho;
//        
//        // element density functions
//        std::map<const libMesh::Elem*, MAST::ConstantFieldFunction*>  _elem_rho_f;
//
//        // element modulus functions
//        std::map<const libMesh::Elem*, MAST::YoungsModulus*>          _elem_E_f;
//
//        // element material property card, one per element
//        std::map<const libMesh::Elem*, MAST::IsotropicMaterialPropertyCard*> _elem_m_card;
//        
//        // element property card, one per element
//        std::map<const libMesh::Elem*, MAST::Solid2DSectionElementPropertyCard*> _elem_p_card;
//        
//        // create the Dirichlet boundary condition on left edge
//        MAST::DirichletBoundaryCondition*               _dirichlet_left;
//        
//        // create the Dirichlet boundary condition on right edge
//        MAST::DirichletBoundaryCondition*               _dirichlet_right;
//        
//        // create the pressure boundary condition
//        MAST::BoundaryConditionBase*                    _p_load;
//        
//        // output object
//        MAST::RealOutputFunction*                       _output;
//        
//    };
//}
//
//
//#endif /* __mast_topology_optim_2d_optimization_h__ */
//
