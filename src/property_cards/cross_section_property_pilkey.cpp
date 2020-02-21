#include <list>

#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/cross_section_property_pilkey.h"
#include "base/warping_assembly.h"
#include "elasticity/warping_system_initialization.h"
#include "base/nonlinear_system.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "elasticity/structural_nonlinear_assembly.h"

#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/mesh_triangle_holes.h"

// this class allows us to set the solver and preconditioner to be appropriate
// for linear elasticity. This is modeled after libMesh's 
// systems_of_equations_ex6.C example sorce code.

/**
 * A class that is used by libMesh to set options on PETSc's linear solver.
 * Modeled after the code in libMesh's systems_of_equations_ex6.C example 
 * source code.
 */
class PetscLinearSolverConfiguration : public libMesh::SolverConfiguration
{
public:
    
    PetscLinearSolverConfiguration(
        libMesh::PetscLinearSolver<libMesh::Number> & petsc_linear_solver):
        _petsc_linear_solver(petsc_linear_solver)
    {
    }
    
    virtual void configure_solver()
    {
        PetscErrorCode ierr = 0;
        
        /** 
         * Let's PETSc know that it is expecting a zero initial guess. This
         * make sense since direct solvers don't utialize an intial guess.
         */
        ierr = KSPSetInitialGuessNonzero(_petsc_linear_solver.ksp(), PETSC_FALSE);
        CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
        
        /**
         * Set the KSP type to PREONLY which indicates only the preconditoner
         * will be applied, and no iterations performed afterwards. This
         * results in a direct solver when the preconditioner is a complete
         * factorization of the matrix.
         */
        ierr = KSPSetType(_petsc_linear_solver.ksp(), KSPPREONLY);
        CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
        
        /**
         * Set the preconditioner to be a complete LU factorization
         */
        ierr = PCSetType(_petsc_linear_solver.pc(), PCLU);
        CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
        
        /**
         * Have MUMPS perform the preconditioner factorization.
         */
        ierr = PCFactorSetMatSolverType(_petsc_linear_solver.pc(), MATSOLVERMUMPS);
        CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
        
        /**
         * During factorization, to minimize fill-in (an increase in the number
         * of nonzeros (nnz) in a sparse matrix versus the nnz in the 
         * unfactored matrix) which causes increased memory usage and higher
         * computational cost, the matrix is reordered. If multiple matrices
         * are to be factored and they have similar sparsity structures (such 
         * as when performing optimization and doing multiple FEA solves), the
         * reordering that was computed on the first matrix will be reused on
         * all subsequent matrices. This reduces computational cost by avoiding
         * redundant calculation of the reordering permutations.
         */
        ierr = PCFactorSetReuseOrdering(_petsc_linear_solver.pc(), PETSC_TRUE);
        CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
    }
    
    libMesh::PetscLinearSolver<libMesh::Number>& _petsc_linear_solver;
};


/**
 * A class that is used by libMesh to set options on PETSc's nonlinear solver.
 * Modeled after the code in libMesh's systems_of_equations_ex6.C example 
 * source code.
 */
class PetscNonlinearSolverConfiguration : public libMesh::SolverConfiguration
{
public:
    
    PetscNonlinearSolverConfiguration(libMesh::PetscNonlinearSolver<libMesh::Number> & petsc_nonlinear_solver):
    _petsc_nonlinear_solver(petsc_nonlinear_solver)
    {
    }
    
    virtual void configure_solver()
    {
        PetscErrorCode ierr = 0;
        
        /**
         * We can make the nonlinear solver act as a linear solver by specifing
         * to do only one linear solver.   This is useful in MAST where 
         * everything is modeled in a nonlinear system, and we want to force
         * a linear solver if the system is truly linear.
         */
        ierr = SNESSetType(_petsc_nonlinear_solver.snes(), SNESKSPONLY);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        /**
         * We can set the linesearch used by the nonlinear solver.  Here, we
         * used a basic linesearch with a damping value of 1.0 which
         * corresponds to the full Newton-Raphson method (when using a 
         * Newton-Raphson based nonlinear solver).
         */
        SNESLineSearch linesearch;
        ierr = SNESGetLineSearch(_petsc_nonlinear_solver.snes(), &linesearch);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        ierr = SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        ierr = SNESLineSearchSetDamping(linesearch, 1.0);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        /**
         * Recall that a nonlinear solver essentially linearizes a nonlinear
         * system and does iterative linear solves until it converges to the 
         * nonlinear solution.  Below, we configure the options for the linear
         * solver called by the nonlinear solver.
         */
        
        KSP ksp;
        SNESGetKSP(_petsc_nonlinear_solver.snes(), &ksp);
        
        ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        ierr = KSPSetType(ksp, KSPPREONLY);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        PC pc;
        KSPGetPC(ksp, &pc);
        
        ierr = PCSetType(pc, PCLU);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
        
        ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        CHKERRABORT(_petsc_nonlinear_solver.comm().get(), ierr);
    }
    
    libMesh::PetscNonlinearSolver<libMesh::Number>& _petsc_nonlinear_solver;
};


MAST::CrossSection::CrossSection(
    const libMesh::LibMeshInit& init, 
    const uint target_n_elems,
    MAST::Solid1DSectionElementPropertyCard& property,
    const libMesh::ElemType element_type): 
    _init(init), _property(property), _target_n_elems(target_n_elems),
    _element_type(element_type)
{
    _mesh.reset(new libMesh::ReplicatedMesh(_init.comm()));
}


Real MAST::CrossSection::get_area(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    return _A;
}


Real MAST::CrossSection::get_area_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    _update_geometry(f, p, t);
    _calculate_shoelace_geometric_properties_sensitivities(f, p, t);
    return _dA;
}


RealVectorX MAST::CrossSection::get_second_area_moments(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    RealVectorX I(3);
    I << _Ixx, _Iyy, _Ixy;
    return I;
}


RealVectorX MAST::CrossSection::get_second_area_moments_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    _update_geometry(f, p, t);
    RealVectorX dI(3);
    dI << _dIxx, _dIyy, _dIxy;
    return dI;
}


Real MAST::CrossSection::get_first_area_moment_z(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    return _Qx;
}


Real MAST::CrossSection::get_first_area_moment_z_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    _update_geometry(f, p, t);
    return _dQx;
}


Real MAST::CrossSection::get_first_area_moment_y(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    return _Qy;
}


Real MAST::CrossSection::get_first_area_moment_y_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    _update_geometry(f, p, t);
    return _dQy;
}


Real MAST::CrossSection::get_second_polar_area_moment(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    Real Ip = _Ixx + _Iyy;
    return Ip;
}


Real MAST::CrossSection::get_second_polar_area_moment_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    _update_geometry(f, p, t);
    Real dIp = _dIxx + _dIyy;
    return dIp;
}


libMesh::Point MAST::CrossSection::get_centroid(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    return libMesh::Point(_Cx, _Cy, 0.);
}


libMesh::Point MAST::CrossSection::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    _update_geometry(f, p, t);
    return libMesh::Point(_dCx, _dCy, 0.);
}


RealVectorX MAST::CrossSection::get_second_area_moments_about_centroid(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    RealVectorX Ic(3);
    Ic << _Ixxc, _Iyyc, _Ixyc;
    return Ic;
}


RealVectorX MAST::CrossSection::get_radii_of_gyration(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    RealVectorX r(2);
    r << _rx, _ry;
    return r;
}


RealVectorX MAST::CrossSection::get_elastic_section_modulii(const libMesh::Point& p, const Real t)
{
    _update_geometry(p, t);
    RealVectorX Z(4);
    Z << _Z_xx_p, _Z_xx_n, _Z_yy_p, _Z_yy_n;
    return Z;
}


Real MAST::CrossSection::get_torsion_constant(const libMesh::Point& p, const Real t)
{
    _update_warping_geometry(p, t);
    return _wp.J;
}

// TODO: Implement analytical derivative instead of finite difference
Real MAST::CrossSection::get_torsion_constant_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    //Real delta = 1.0e-06;
    //Real delta = 1.4901161193847656e-08; // np.spacing(1)**0.5
    const Real delta = 1.220703125e-04; // np.spacing(1)**0.25
    MAST::Parameter* param = dynamic_cast<MAST::Parameter*>(&f);
    (*param)() += delta;
    _update_warping_geometry(p, t);
    const Real Jh = _wp.J;
    (*param)() -= 2.0*delta;
    _update_warping_geometry(p, t);
    const Real Jn = _wp.J;
    (*param)() += delta;
    return (Jh - Jn)/(2.0*delta);
}


libMesh::Point MAST::CrossSection::get_shear_center(const libMesh::Point& p, const Real t)
{
    _update_warping_geometry(p, t);
    libMesh::Point shear_center(_wp.xs, _wp.ys, 0.);
    return shear_center;
}

// TODO: Implement analytical derivative instead of finite difference
libMesh::Point MAST::CrossSection::get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    //Real delta = 1.0e-06;
    const Real delta = 1.220703125e-04; // np.spacing(1)**0.25
    MAST::Parameter* param = dynamic_cast<MAST::Parameter*>(&f);
    (*param)() += delta;
    _update_warping_geometry(p, t);
    const Real xs_h = _wp.xs;
    const Real ys_h = _wp.ys;
    
//     (*param)() += delta;
//     _update_warping_geometry(p, t);
//     const Real xs_2h = _wp.xs;
//     const Real ys_2h = _wp.ys;
//     (*param)() -= delta;
    
    (*param)() -= 2.0*delta;
    _update_warping_geometry(p, t);
    const Real xs_n = _wp.xs;
    const Real ys_n = _wp.ys;
    
//     (*param)() -= delta;
//     _update_warping_geometry(p, t);
//     const Real xs_2n = _wp.xs;
//     const Real ys_2n = _wp.ys;
//     (*param)() + =delta;
    
    (*param)() += delta;
    
    //return libMesh::Point((xs_h - xs_n)/(2.0*delta), (ys_h - ys_n)/(2.0*delta), 0.);
    Real dxs = (xs_h - xs_n)/(2.0*delta);
    Real dys = (ys_h - ys_n)/(2.0*delta);
//     Real dxs = (xs_2n - 8.*xs_n + 8.*xs_h - xs_2h)/(12.*delta);
//     Real dys = (ys_2n - 8.*ys_n + 8.*ys_h - ys_2h)/(12.*delta);
    return libMesh::Point(dxs, dys, 0.);
}


RealVectorX MAST::CrossSection::get_shear_coefficients(const libMesh::Point& p, const Real t)
{
    _update_warping_geometry(p, t);
    RealVectorX kappa(3);
    kappa << _wp.kappa_x, _wp.kappa_y, _wp.kappa_xy;
    return kappa;
}

// TODO: Implement analytical derivative instead of finite difference
RealVectorX MAST::CrossSection::get_shear_coefficients_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    //Real delta = 1.0e-06;
    const Real delta = 1.220703125e-04; // np.spacing(1)**0.25
    MAST::Parameter* param = dynamic_cast<MAST::Parameter*>(&f);
    
    (*param)() += delta;
    _update_warping_geometry(p, t);
    Real kappa_x_h =  _wp.kappa_x;
    Real kappa_y_h =  _wp.kappa_y;
    Real kappa_xy_h = _wp.kappa_xy;
    
    (*param)() -= 2.0*delta;
    _update_warping_geometry(p, t);
    Real kappa_x_n =  _wp.kappa_x;
    Real kappa_y_n =  _wp.kappa_y;
    Real kappa_xy_n = _wp.kappa_xy;
    
    (*param)() += delta;
    
    Real dkappa_x =  (kappa_x_h -  kappa_x_n)  / (2.0*delta);
    Real dkappa_y =  (kappa_y_h -  kappa_y_n)  / (2.0*delta);
    Real dkappa_xy = (kappa_xy_h - kappa_xy_n) / (2.0*delta);
    
    RealVectorX dkappa(3);
    dkappa << dkappa_x, dkappa_y, dkappa_xy;
    
    return dkappa;
}


Real MAST::CrossSection::get_warping_constant(const libMesh::Point& p, const Real t)
{
    _update_warping_geometry(p, t);
    return _wp.gamma;
}

// TODO: Implement analytical derivative instead of finite difference
Real MAST::CrossSection::get_warping_constant_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    //Real delta = 1.0e-06;
    const Real delta = 1.220703125e-04; // np.spacing(1)**0.25
    MAST::Parameter* param = dynamic_cast<MAST::Parameter*>(&f);
    
    (*param)() += delta;
    _update_warping_geometry(p, t);
    Real gammah = _wp.gamma;
    (*param)() -= 2.0*delta;
    _update_warping_geometry(p, t);
    Real gamman = _wp.gamma;
    (*param)() += delta;
    return (gammah - gamman)/(2.0*delta);
}


void MAST::CrossSection::_update_geometry(const libMesh::Point& p, const Real t)
{
    // Get the points defining the outter geometry from the mesh
    std::vector<libMesh::Point> geometry_points = _property.get_geom_points(p, t);
    std::vector<std::vector<libMesh::Point>> hole_list = _property.get_holes_points(p, t);
    
    // Check if the geometry is the same as the last time it was calculated
    if ((geometry_points == _geometry_points) and (hole_list == _hole_list))
    {
        // Geometry did not change, no need to redo calculations.
        //libMesh::out << "Geometry the same. SKIPPING update." << std::endl;
    }
    else
    {
        // Geometry did change, need to re-run calculations
        _geometry_points = geometry_points;
        _hole_list = hole_list;
        _calculate_shoelace_geometric_properties();
        _build_model(p, t);
        //libMesh::out << "Geometry is different. UPDATING." << std::endl;
    }
}


void MAST::CrossSection::_update_warping_geometry(const libMesh::Point& p, const Real t)
{
    // Get the points defining the outter geometry from the mesh
    std::vector<libMesh::Point> geometry_points = _property.get_geom_points(p, t);
    std::vector<std::vector<libMesh::Point>> hole_list = _property.get_holes_points(p, t);
    
    
    // Check if the geometry is the same as the last time it was calculated
    if ((geometry_points == _geometry_points) and (hole_list == _hole_list) and (_warping_calculated))
    {
        // Geometry did not change, no need to redo calculations.
        libMesh::out << "Geometry the same. SKIPPING warping update." << std::endl;
    }
    else
    {
        libMesh::out << "Geometry is different. UPDATING warping." << std::endl;
        // Geometry did change, need to re-run calculations
        _warping_calculated = false;
        _geometry_points = geometry_points;
        _hole_list = hole_list;
        _calculate_shoelace_geometric_properties();
        _build_model(p, t);
        _solve_warping();
        _wp = _assembly->calculate_warping_properties(*_F_warp, *_Omega, *_Psi, 
                                                *_Phi, _A, _Ixxc, _Iyyc, _Ixyc,
                                                _Cx, _Cy);
        _warping_calculated = true;
        libMesh::out << "Done updating warping." << std::endl;
    }
}


void MAST::CrossSection::_update_geometry(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    // Get the points defining the outter geometry from the mesh
    std::vector<libMesh::Point> geometry_points = _property.get_geom_points(p, t);
    std::vector<std::vector<libMesh::Point>> hole_list = _property.get_holes_points(p, t);
    std::string sensitivity_parameter = f.name();
    
    // Check if the geometry is the same as the last time it was calculated
    if ((geometry_points == _geometry_points) and (hole_list == _hole_list) and (sensitivity_parameter == _sensitivity_parameter))
    {
        // Geometry did not change, no need to redo calculations.
        //libMesh::out << "Geometry the same. SKIPPING update." << std::endl;
    }
    else
    {
        // Geometry did change, need to re-run calculations
        _geometry_points = geometry_points;
        _hole_list = hole_list;
        _sensitivity_parameter = sensitivity_parameter;
        _calculate_shoelace_geometric_properties_sensitivities(f, p, t);
        _build_model(p, t);
        //libMesh::out << "Geometry is different. UPDATING." << std::endl;
    }
}


void MAST::CrossSection::calculate_shoelace_geometric_properties()
{
    _calculate_shoelace_geometric_properties();
    
    libMesh::out << "\nShoelace Geometric Properties" << std::endl;
    libMesh::out << "xc = " << _Cx << std::endl;
    libMesh::out << "yc = " << _Cy << std::endl;
    libMesh::out << "A = " << _A << std::endl;
    libMesh::out << "Ay = " << _Qy << std::endl;
    libMesh::out << "Az = " << _Qx << std::endl;
    libMesh::out << "Iyy = " << _Iyy << std::endl;
    libMesh::out << "Izz = " << _Ixx << std::endl;
    libMesh::out << "Iyz = " << _Ixy << std::endl;
    libMesh::out << "Ip = " << (_Ixx+_Iyy) << std::endl;
}


void MAST::CrossSection::_calculate_shoelace_geometric_properties()
{
    // Shoelace formula
    Real xi, yi, xi1, yi1, yin1, Ibase;
    
    Real A = 0.0;
    Real Cx = 0.0;
    Real Cy = 0.0;
    Real Ixx = 0.0;
    Real Iyy = 0.0;
    Real Ixy = 0.0;
    Real y_max = -1e16;
    Real y_min = 1e16;
    Real x_max = -1e16;
    Real x_min = 1e16;
    
    // Loop through points describing outer geometry boundary
    const uint n = _geometry_points.size();
    for (uint i=0; i<n; i++)
    {
        xi = _geometry_points[i](0);
        yi = _geometry_points[i](1);
        
        // Get max and min x values
        x_max = std::max(x_max, xi);
        x_min = std::min(x_min, xi);
        
        // Get max and min y values
        y_max = std::max(y_max, yi);
        y_min = std::min(y_min, yi);
        
        if (i==0)
        {
            yin1 = _geometry_points[n-1](1);
        }
        else
        {
            yin1 = _geometry_points[i-1](1);
        }
        
        if (i==n-1)
        {
            xi1 = _geometry_points[0](0);
            yi1 = _geometry_points[0](1);
        }
        else
        {
            xi1 = _geometry_points[i+1](0);
            yi1 = _geometry_points[i+1](1);
        }
        
        A += xi*(yi1 - yin1);
        Ibase = (xi*yi1 - xi1*yi);
        Cx += (xi + xi1) * Ibase;
        Cy += (yi + yi1) * Ibase;
        Ixx += Ibase*(yi*yi + yi*yi1 + yi1*yi1);
        Iyy += Ibase*(xi*xi + xi*xi1 + xi1*xi1);
        Ixy += Ibase*(xi*yi1 + 2.0*xi*yi + 2.0*xi1*yi1 + xi1*yi);
    }
    
    // Loop through holes in the geometry
    for (const auto& hole_points : _hole_list)
    {
        // Loop through points describing the hole
        const uint n = hole_points.size();
        for (uint i=0; i<n; i++)
        {
            xi = hole_points[i](0);
            yi = hole_points[i](1);
            
            if (i==0)
            {
                yin1 = hole_points[n-1](1);
            }
            else
            {
                yin1 = hole_points[i-1](1);
            }
            
            if (i==n-1)
            {
                xi1 = hole_points[0](0);
                yi1 = hole_points[0](1);
            }
            else
            {
                xi1 = hole_points[i+1](0);
                yi1 = hole_points[i+1](1);
            }
            
            // NOTE: We are using -= here instead of += above, because this
            // area is considered void so its contribution should be negative.
            A -= xi*(yi1 - yin1);
            Cx -= (xi + xi1) * (xi*yi1 - xi1*yi);
            Cy -= (yi + yi1) * (xi*yi1 - xi1*yi);
            Ibase = (xi*yi1 - xi1*yi);
            Ixx -= Ibase*(yi*yi + yi*yi1 + yi1*yi1);
            Iyy -= Ibase*(xi*xi + xi*xi1 + xi1*xi1);
            Ixy -= Ibase*(xi*yi1 + 2.0*xi*yi + 2.0*xi1*yi1 + xi1*yi);
        } // end loop through hole points
    } // end loop through hole list
    
    _A = abs(A)/2.0;
    _Cx = Cx/(6.0*_A);
    _Cy = Cy/(6.0*_A);
    _Qx = _Cy * _A;
    _Qy = _Cx * _A;
    _Ixx = Ixx / 12.0;
    _Iyy = Iyy / 12.0;
    _Ixy = Ixy / 24.0;
    _Ixxc = _Ixx - _Cy*_Cy * _A;
    _Iyyc = _Iyy - _Cx*_Cx * _A;
    _Ixyc = _Ixy - _Cx*_Cy * _A;
    _rx = sqrt(_Ixx/_A);
    _ry = sqrt(_Iyy/_A);
    _Z_xx_p = _Ixxc / (y_max - _Cy);
    _Z_xx_n = _Ixxc / (_Cy - y_min);
    _Z_yy_p = _Iyyc / (x_max - _Cx);
    _Z_yy_n = _Iyyc / (_Cx - x_max);
    
    _shoelace_calculated = true;
}


void MAST::CrossSection::_calculate_shoelace_geometric_properties_sensitivities(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    // Shoelace formula
    Real xi, yi, xi1, yi1, yin1, Ibase;
    Real dxi, dyi, dxi1, dyi1, dyin1, dIbase;
    
    Real A = 0.0;
    Real Cy = 0.0;
    Real Cx = 0.0;
    Real Ixx = 0.0;
    Real Iyy = 0.0;
    Real Ixy = 0.0;
    Real y_max = -1e16;
    Real y_min = 1e16;
    Real x_max = -1e16;
    Real x_min = 1e16;
    
    Real dA = 0.0;
    Real dCx = 0.0;
    Real dCy = 0.0;
    Real dIxx = 0.0;
    Real dIyy = 0.0;
    Real dIxy = 0.0;
    
    Real dy_max = 0.0;
    Real dy_min = 0.0;
    Real dx_max = 0.0;
    Real dx_min = 0.0;
    
    Real alpha_x  = 1.0;
    Real alpha_y = 1.0;
    
    Real dy_max_temp1;
    Real dy_max_temp2;
    Real dy_min_temp1;
    Real dy_min_temp2;
    Real dx_max_temp1;
    Real dx_max_temp2;
    Real dx_min_temp1;
    Real dx_min_temp2;
    Real pn = 20.0;
    
    // Get the geometry point sensitivity
    const std::vector<libMesh::Point> geometry_point_sensitivity = _property.get_geom_points_derivative(f, p, t); 
    
    // Loop through points describing outer geometry boundary
    const uint n = _geometry_points.size();
    for (uint i=0; i<n; i++)
    {
        xi = _geometry_points[i](0);
        dxi = geometry_point_sensitivity[i](0);
        
        yi = _geometry_points[i](1);
        dyi = geometry_point_sensitivity[i](1);
        
        // Get max and min x values
        x_max = std::max(x_max, xi);
        x_min = std::min(x_min, xi);
        
        // Get max and min y values
        y_max = std::max(y_max, yi);
        y_min = std::min(y_min, yi);
        
        dy_max_temp1 += pow(yi/alpha_y, pn);
        dy_max_temp2 += pow(yi/alpha_y, pn-1) * (dyi/alpha_y);
        
        dy_min_temp1 += pow(-yi/alpha_y, pn);
        dy_min_temp2 += pow(-yi/alpha_y, pn-1) * (-dyi/alpha_y);
        
        dx_max_temp1 += pow(xi/alpha_x, pn);
        dx_max_temp2 += pow(xi/alpha_x, pn-1) * (dxi/alpha_x);
        
        dx_min_temp1 += pow(-xi/alpha_x, pn);
        dx_min_temp2 += pow(-xi/alpha_x, pn-1) * (-dxi/alpha_x);
        
        if (i==0)
        {
            yin1 = _geometry_points[n-1](1);
            dyin1 = geometry_point_sensitivity[n-1](1);
        }
        else
        {
            yin1 = _geometry_points[i-1](1);
            dyin1 = geometry_point_sensitivity[i-1](1);
        }
        
        if (i==n-1)
        {
            xi1 = _geometry_points[0](0);
            dxi1 = geometry_point_sensitivity[0](0);
            
            yi1 = _geometry_points[0](1);
            dyi1 = geometry_point_sensitivity[0](1);
        }
        else
        {
            xi1 = _geometry_points[i+1](0);
            dxi1 = geometry_point_sensitivity[i+1](0);
            
            yi1 = _geometry_points[i+1](1);
            dyi1 = geometry_point_sensitivity[i+1](1);
        }
        
        A += xi*(yi1 - yin1);
        Ibase = (xi*yi1 - xi1*yi);
        Cx += (xi + xi1) * Ibase;
        Cy += (yi + yi1) * Ibase;
        Ixx += Ibase*(yi*yi + yi*yi1 + yi1*yi1);
        Iyy += Ibase*(xi*xi + xi*xi1 + xi1*xi1);
        Ixy += Ibase*(xi*yi1 + 2.0*xi*yi + 2.0*xi1*yi1 + xi1*yi);
        
        dA += dxi*(yi1 - yin1) + xi*(dyi1 - dyin1);
        dIbase = dxi*yi1 + xi*dyi1 - dxi1*yi - xi1*dyi;
        dCx += (dxi + dxi1) * Ibase + (xi + xi1) * dIbase; 
        dCy += (dyi + dyi1) * Ibase + (yi + yi1) * dIbase; 
        dIxx += dIbase*(yi*yi + yi*yi1 + yi1*yi1) + Ibase * ( 2.*yi*dyi + yi*dyi1 + dyi*yi1 + 2.*yi1*dyi1 );
        dIyy += dIbase*(xi*xi + xi*xi1 + xi1*xi1) + Ibase * ( 2.*xi*dxi + xi*dxi1 + dxi*xi1 + 2.*xi1*dxi1 );
        dIxy += dIbase*(xi*yi1 + 2.0*xi*yi + 2.0*xi1*yi1 + xi1*yi) + Ibase * ( xi*dyi1 + dxi*yi1 + 2.*(xi*dyi+dxi*yi) + 2.*(xi1*dyi1+dxi1*yi1) + xi1*dyi + dxi1*yi );
    }
    dy_max =  alpha_y * pow(dy_max_temp1, 1.0/pn - 1.0) * dy_max_temp2;
    dy_min = -alpha_y * pow(dy_min_temp1, 1.0/pn - 1.0) * dy_min_temp2; 
    
    dx_max =  alpha_x * pow(dx_max_temp1, 1.0/pn - 1.0) * dx_max_temp2;
    dx_min = -alpha_x * pow(dx_min_temp1, 1.0/pn - 1.0) * dx_min_temp2; 
    
    const std::vector<std::vector<libMesh::Point>> hole_derivatives_list = _property.get_holes_points_derivative(f, p, t);
    
    // Loop through holes in the geometry
    for (uint j=0; j<_hole_list.size(); j++)
    {
        const std::vector<libMesh::Point> hole_points = _hole_list[j];
        const std::vector<libMesh::Point> hole_point_derivative = hole_derivatives_list[j];
        // Loop through points describing the hole
        const uint n = hole_points.size();
        for (uint i=0; i<n; i++)
        {
            xi = hole_points[i](0);
            dxi = hole_point_derivative[i](0);
            
            yi = hole_points[i](1);
            dyi = hole_point_derivative[i](1);
            
            if (i==0)
            {
                yin1 = hole_points[n-1](1);
                dyin1 = hole_point_derivative[n-1](1);
            }
            else
            {
                yin1 = hole_points[i-1](1);
                dyin1 = hole_point_derivative[i-1](1);
            }
            
            if (i==n-1)
            {
                xi1 = hole_points[0](0);
                dxi1 = hole_point_derivative[0](0);
                
                yi1 = hole_points[0](1);
                dyi1 = hole_point_derivative[0](1);
            }
            else
            {
                xi1 = hole_points[i+1](0);
                dxi1 = hole_point_derivative[i+1](0);
                
                yi1 = hole_points[i+1](1);
                dyi1 = hole_point_derivative[i+1](1);
            }
            
            // NOTE: We are using -= here instead of += above, because this
            // area is considered void so its contribution should be negative.
            A -= xi*(yi1 - yin1);
            Ibase = (xi*yi1 - xi1*yi);
            Cx -= (xi + xi1) * Ibase;
            Cy -= (yi + yi1) * Ibase;
            Ixx -= Ibase*(yi*yi + yi*yi1 + yi1*yi1);
            Iyy -= Ibase*(xi*xi + xi*xi1 + xi1*xi1);
            Ixy -= Ibase*(xi*yi1 + 2.0*xi*yi + 2.0*xi1*yi1 + xi1*yi);
            
            dA -= dxi*(yi1 - yin1) + xi*(dyi1 - dyin1);
            dIbase = dxi*yi1 + xi*dyi1 - dxi1*yi - xi1*dyi;
            dCx -= (dxi + dxi1) * Ibase + (xi + xi1) * dIbase; 
            dCy -= (dyi + dyi1) * Ibase + (yi + yi1) * dIbase; 
            dIxx -= dIbase*(yi*yi + yi*yi1 + yi1*yi1) + Ibase * ( 2.*yi*dyi + yi*dyi1 + dyi*yi1 + 2.*yi1*dyi1 );
            dIyy -= dIbase*(xi*xi + xi*xi1 + xi1*xi1) + Ibase * ( 2.*xi*dxi + xi*dxi1 + dxi*xi1 + 2.*xi1*dxi1 );
            dIxy -= dIbase*(xi*yi1 + 2.0*xi*yi + 2.0*xi1*yi1 + xi1*yi) + Ibase * ( xi*dyi1 + dxi*yi1 + 2.*(xi*dyi+dxi*yi) + 2.*(xi1*dyi1+dxi1*yi1) + xi1*dyi + dxi1*yi );
        } // end loop through hole points
    } // end loop through hole list
    
    Real tCx = Cx;
    Real tCy = Cy;
    A = abs(A)/2.0;
    Cx = tCx/(6.0*A);
    Cy = tCy/(6.0*A);
    Ixx = Ixx / 12.0;
    Iyy = Iyy / 12.0;
    Ixy = Ixy / 24.0;
    Real Ixxc = Ixx - Cy*Cy * A;
    Real Iyyc = Iyy - Cx*Cx * A;
    Real Ixyc = Ixy - Cx*Cy * A;
    
    _dA = 0.5*(abs(A)/A * dA);
    _dCx = (dCx - tCx*_dA/A)/(6.*A);
    _dCy = (dCy - tCy*_dA/A)/(6.*A);
    _dQx = Cy * _dA + _dCy * A;
    _dQy = Cx * _dA + _dCx * A;
    _dIxx = dIxx / 12.0;
    _dIyy = dIyy / 12.0;
    _dIxy = dIxy / 24.0;
    Real _dIxxc = _dIxx - 2.*Cy*_dCy*A - Cy*Cy*_dA;
    Real _dIyyc = _dIyy - 2.*Cx*_dCx*A - Cx*Cx*_dA;
    Real _dIxyc = _dIxy - ((Cx*dCy + dCx*Cy) * A + Cx*Cy * dA);
    Real _drx = ( pow((Ixx*A),-0.5)*_dIxx - pow(Ixx,0.5)*pow(A,-1.5)*_dA )/2.0;
    Real _dry = ( pow((Iyy*A),-0.5)*_dIyy - pow(Iyy,0.5)*pow(A,-1.5)*_dA )/2.0;
    
    Real _dZ_xx_p = ( _dIxxc - Ixxc/(y_max - Cy) * (dy_max - dCy) )  / (y_max - Cy);
    Real _dZ_xx_n = ( _dIxxc - Ixxc/(Cy - y_min) * (dCy - dy_min) )  / (Cy - y_min);
    
    Real _dZ_yy_p = ( _dIyyc - Iyyc/(x_max - Cx) * (dx_max - dCx) )  / (x_max - Cx);
    Real _dZ_yy_n = ( _dIyyc - Iyyc/(Cx - x_min) * (dCx - dx_min) )  / (Cx - x_min);
    
//     libMesh::out << "dA = " << _dA << std::endl;
//     libMesh::out << "dCx = " << _dCx << std::endl;
//     libMesh::out << "dCy = " << _dCy << std::endl;
//     libMesh::out << "dQx = " << _dQx << std::endl;
//     libMesh::out << "dQy = " << _dQy << std::endl;
//     libMesh::out << "dIxx = " << _dIxx << std::endl;
//     libMesh::out << "dIyy = " << _dIyy << std::endl;
//     libMesh::out << "dIxy = " << _dIxy << std::endl;
//     libMesh::out << "dIxxc = " << _dIxxc << std::endl;
//     libMesh::out << "dIyyc = " << _dIyyc << std::endl;
//     libMesh::out << "dIxyc = " << _dIxyc << std::endl;
//     libMesh::out << "drx = " << _drx << std::endl;
//     libMesh::out << "dry = " << _dry << std::endl;
//     libMesh::out << "dZ_xx_p = " << _dZ_xx_p << std::endl;
//     libMesh::out << "dZ_xx_n = " << _dZ_xx_n << std::endl;
//     libMesh::out << "dZ_yy_p = " << _dZ_yy_p << std::endl;
//     libMesh::out << "dZ_yy_n = " << _dZ_yy_n << std::endl;
}



void MAST::CrossSection::calculate_geometric_properties()
{
    // Get the geometric properties which can be explicitly calculated
    _gp = _assembly->calculate_geometric_properties();
    
//     libMesh::out << "\nPilkey Geometric Properties:" << std::endl;
//     libMesh::out << "xc = " << _gp.xc << std::endl;
//     libMesh::out << "yc = " << _gp.yc << std::endl;
//     libMesh::out << "A = " << _gp.A << std::endl;
//     libMesh::out << "Ax = " << _gp.Qx << std::endl;
//     libMesh::out << "Ay = " << _gp.Qy << std::endl;
//     libMesh::out << "Ixx = " << _gp.Ixx << std::endl;
//     libMesh::out << "Iyy = " << _gp.Iyy << std::endl;
//     libMesh::out << "Ixy = " << _gp.Ixy << std::endl;
//     libMesh::out << "Ip = " << _gp.Ip << std::endl;
//     libMesh::out << "Ixxc = " << _gp.Ixxc << std::endl;
//     libMesh::out << "Iyyc = " << _gp.Iyyc << std::endl;
//     libMesh::out << "Ixyc = " << _gp.Ixyc << std::endl;
}

void MAST::CrossSection::calculate_warping_properties()
{
    _wp = _assembly->calculate_warping_properties(*_F_warp, *_Omega, 
                                                  *_Psi, *_Phi, _gp.A, 
                                                  _gp.Ixxc, _gp.Iyyc, _gp.Ixyc,
                                                  _gp.xc, _gp.yc);
    
//     libMesh::out << "xs = " << _wp.xs << std::endl;
//     libMesh::out << "ys = " << _wp.ys << std::endl;
//     libMesh::out << "xs_t = " << _wp.xs_t << std::endl;
//     libMesh::out << "ys_t = " << _wp.ys_t << std::endl;
//     libMesh::out << "J = " << _wp.J << std::endl;
//     libMesh::out << "kappa_x = " << _wp.kappa_x << std::endl;
//     libMesh::out << "kappa_y = " << _wp.kappa_y << std::endl;
//     libMesh::out << "kappa_xy = " << _wp.kappa_xy << std::endl;
//     libMesh::out << "gamma = " << _wp.gamma << std::endl;
}


void MAST::CrossSection::_mesh_cross_section(const libMesh::Point& p, const Real t)
{
    // Clear the mesh
    _mesh->clear();
    
    // Set the Mesh Dimensions
    _mesh->set_mesh_dimension(2);
    _mesh->set_spatial_dimension(2);

    /**
    * Now define the corner points on the boundaries of the mesh. If we need
    * a node to appear in a certain position, for example to apply a boundary
    * condition or a load on that node, we should also create those points as
    * well.  Here we go in a clockwise direction.
    */
    std::vector<libMesh::Point> geom_points = _property.get_geom_points(p, t);
    for (const libMesh::Point point : geom_points)
    {
        _mesh->add_point(point);
    }

    /**
    * Declare the TriangleInterface object. This is where we can set up 
    * parameters of the triangulation and where the actual triangulate 
    * function lives.
    */
    libMesh::TriangleInterface tri(*_mesh);
    
    /**
     *  Create a vector containing holes which should be meshed around during
     *  triangularization.
     */
    std::vector<std::vector<libMesh::Point>> holes_list = _property.get_holes_points(p, t);
    std::vector<libMesh::TriangleInterface::Hole*> holes(holes_list.size());
    uint i=0;
    for (const auto& hole_points : holes_list)
    {
        // Get the centroid of the hole.
        Real c_x = 0.0;
        Real c_y = 0.0;
        for (const libMesh::Point point : hole_points)
        {
            c_x += point(0);
            c_y += point(1);
        }
        c_x /= Real(hole_points.size());
        c_y /= Real(hole_points.size());
        libMesh::Point hole_centroid(c_x, c_y, 0.);
        
        // Create the triangulation hole and add it to the vector of all holes
        holes[i] = new libMesh::TriangleInterface::ArbitraryHole(hole_centroid, hole_points);
        i++;
    }
    
    tri.attach_hole_list(&holes);

    /**
    * Custmoize the variables for the triangulation. desired_area can control
    * the coarseness of the mesh, with smaller values resulting in a finer
    * mesh. elem_type can be used to set the type of element to either constant
    * strain triangles (TRI3) or linear strain triangles (TRI6).
    */
    if (not _shoelace_calculated)
    {
        _calculate_shoelace_geometric_properties();
    }
    tri.desired_area() = 1.5*_A/Real(_target_n_elems);
    tri.elem_type() = _element_type;


    /** A Planar Straight Line Graph (PSLG) is essentially a list
    * of segments which have to exist in the final triangulation.
    * For an L-shaped domain, Triangle will compute the convex
    * hull of boundary points if we do not specify the PSLG.
    * The PSLG algorithm is also required for triangulating domains
    * containing holes.
    */
    tri.triangulation_type() = libMesh::TriangleInterface::PSLG;

    /** Turn on/off Laplacian mesh smoothing after generation.
    * By default this is on.
    */
    tri.smooth_after_generating() = true;

    /**
    * After the triangulation parameters have been defined, we go ahead and
    * triangulate the mesh.
    */
    tri.triangulate();
    
    /**
     * Clean up the holes' raw pointers
     */
    for (uint i=0; i<holes.size(); i++)
    {
        delete holes[i];
    }

    // Print the Information About the Mesh
    //_mesh->print_info();
    //_mesh->boundary_info->print_info();
}


void MAST::CrossSection::_build_model(const libMesh::Point& p, const Real t)
{
    // Mesh the cross section
    _mesh_cross_section(p, t);
    
    // Create EquationSystems object, which is a container for multiple systems of equations that are defined on a given mesh.
    _equation_systems.reset(new libMesh::EquationSystems(*_mesh));
    
    // Add system of type MAST::NonlinearSystem (which wraps libMesh::NonlinearImplicitSystem) to the EquationSystems container.
    // We name the system "structural" and also get a reference to the system so we can easily reference it later.
    _system = &(_equation_systems->add_system<MAST::NonlinearSystem>("warping"));
     
    // Create a finite element type for the system. Here, we use first order Lagrangian-type finite elements.
    _fetype.reset(new libMesh::FEType(_mesh->elem_ptr(0)->default_order(), libMesh::LAGRANGE));
    
    _warping_system.reset(new MAST::WarpingSystemInitialization(*_system, _system->name(), *_fetype));
    
    _discipline.reset(new MAST::PhysicsDisciplineBase(*_equation_systems));
    
    _equation_systems->init();
    
    _parameters["th"].reset(new MAST::Parameter("th", 1.0));
    _parameters["zero"].reset(new MAST::Parameter("zero", 0.0));
    
    // Create field functions to dsitribute these constant parameters throughout the model
    _field_functions["h"].reset(new MAST::ConstantFieldFunction("h", *_parameters["th"]));
    _field_functions["off"].reset(new MAST::ConstantFieldFunction("off", *_parameters["zero"]));
    
    _section.reset(new MAST::Solid2DSectionElementPropertyCard);

    _section->add(*_field_functions["h"]);
    _section->add(*_field_functions["off"]);
    
    // Set the materials for the section
    _section->set_material(_property.get_material());
    _section->set_strain(MAST::LINEAR_STRAIN);
    _section->set_bending_model(MAST::NO_BENDING);
    _section->set_warping_only(true);
    
    // Set the section properties for subdoamin 0 in the discipline
    _discipline->set_property_for_subdomain(0, *_section);
    
    _assembly.reset(new MAST::WarpingAssembly);
    _assembly->set_force_jacobian_symmetry(true);
    _assembly->set_discipline_and_system(*_discipline, *_warping_system);
    
    _elem_ops.reset(new MAST::StructuralNonlinearAssemblyElemOperations);
    _elem_ops->set_discipline_and_system(*_discipline, *_warping_system);
    
    _nonlinear_system = &(_assembly->system());
    
    // Get the geometric properties which can be explicitly calculated
    _gp = _assembly->calculate_geometric_properties();
    
    // Shift the mesh so that the global axis lies on the centroid.
    libMesh::Point centroid(_gp.xc, _gp.yc, 0.0);
    for (auto& node_ptr : _mesh->node_ptr_range())
    {
        node_ptr->subtract(centroid);
    }
    
    // Zero the solution before solving
    _nonlinear_system->solution->zero();
}


void MAST::CrossSection::_solve_warping()
{
    // Get the warping, shear_x, and shear_y loads
    _F_warp = _nonlinear_system->solution->zero_clone();
    _F_shearx = _nonlinear_system->solution->zero_clone();
    _F_sheary = _nonlinear_system->solution->zero_clone();
    _assembly->get_loads(*_F_warp, *_F_shearx, *_F_sheary, _gp.xc, _gp.yc, 
                         _gp.Ixxc, _gp.Iyyc, _gp.Ixyc);
    
    // Configure the nonlinear solver
    libMesh::PetscNonlinearSolver<libMesh::Number>* petsc_nonlinear_solver = dynamic_cast<libMesh::PetscNonlinearSolver<libMesh::Number>*>(_nonlinear_system->nonlinear_solver.get());
    libmesh_assert(petsc_nonlinear_solver);
    PetscNonlinearSolverConfiguration petsc_nonlinear_solver_config(*petsc_nonlinear_solver);
    petsc_nonlinear_solver->set_solver_configuration(petsc_nonlinear_solver_config);
    
    // Solve the system for the warping function
    _nonlinear_system->solve(*_elem_ops, *_assembly);
    
    /** 
     * Get the linear solver object
     * See sensitivity_solve() method in libMesh's implicit_system.C source code
     * for an example of getting the linear solver and using it to solve 
     * additional right hand sides. 
     */
    _nonlinear_system->assemble_before_solve = false;
    libMesh::LinearSolver<libMesh::Number>* linear_solver = _nonlinear_system->get_linear_solver();
        
    /** 
     * Configure the Linear Solver Object
     * See libMesh's systems_of_equations_ex6.C example sorce code as reference.
     */
    libMesh::PetscLinearSolver<libMesh::Number>* petsc_linear_solver = dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number>*>(linear_solver);
    libmesh_assert(petsc_linear_solver);
    PetscLinearSolverConfiguration petsc_solver_config(*petsc_linear_solver);
    petsc_linear_solver->set_solver_configuration(petsc_solver_config);
    
    _Omega = _nonlinear_system->solution->clone();
//     _Omega = _nonlinear_system->solution->zero_clone();
//     linear_solver->solve(*(_nonlinear_system->matrix), *_Omega, *_F_warp, 1e-15, 100000);
    
    // solve the linear system for the shear_x function
    _Psi = _nonlinear_system->solution->zero_clone();
    linear_solver->solve(*(_nonlinear_system->matrix), *_Psi, *_F_shearx, 1e-15, 100000);
    
    
    // solve the linear system for the shear_y function
    _Phi = _nonlinear_system->solution->zero_clone();
    linear_solver->solve(*(_nonlinear_system->matrix), *_Phi, *_F_sheary, 1e-15, 100000);
    
    
    _nonlinear_system->release_linear_solver(linear_solver); // To avoid memory leak
}
