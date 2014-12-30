/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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
#include "property_cards/multilayer_2d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/field_function_base.h"


namespace MAST {
    namespace Multilayer2DSectionProperty {
        
        
        class LayerOffset: public MAST::FieldFunction<Real> {
        public:
            LayerOffset(const Real base,
                        unsigned int layer_num,
                        std::vector<MAST::FieldFunction<Real>*>& layer_h):
            MAST::FieldFunction<Real>("off"),
            _base(base),
            _layer_num(layer_num),
            _layer_h(layer_h) {
                for (unsigned int i=0; i < _layer_h.size(); i++)
                    _functions.insert(_layer_h[i]->master());
            }
            
            LayerOffset(const MAST::Multilayer2DSectionProperty::LayerOffset &f):
            MAST::FieldFunction<Real>(f),
            _base(f._base),
            _layer_num(f._layer_num)
            {
                // initialize the vector
                _layer_h.resize(f._layer_h.size());
                for (unsigned int i=0; i < _layer_h.size(); i++) {
                    _layer_h[i] = f._layer_h[i]->clone().release();
                    _functions.insert(_layer_h[i]->master());
                }
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Multilayer2DSectionProperty::LayerOffset(*this));
            }
            
            virtual ~LayerOffset() {
                // delete all the layer functions
                for (unsigned int i=0; i<_layer_h.size(); i++)
                    delete _layer_h[i];
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real val = 0.;
                m = 0.;
                for (unsigned int i=0; i<_layer_num; i++) {
                    (*_layer_h[i])(p, t, val);
                    m += val; // this calculates offset from h=0, and then it is modified for base
                }
                // finally, add half of the current layer thickness
                (*_layer_h[_layer_num])(p, t, val);
                m += 0.5*val;
                
                // now add the base offset
                for (unsigned int i=0; i<_layer_h.size(); i++) {
                    (*_layer_h[i])(p, t, val);
                    m -= 0.5*(1.+_base)*val; // currently the offset is chosen from h=0;
                }
            }
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                const Real t,
                                Real& m) const {
                Real val = 0.;
                m = 0.;
                for (unsigned int i=0; i<_layer_num; i++) {
                    _layer_h[i]->derivative(d, f, p, t, val);
                    m += val; // this calculates offset from h=0, and then it is modified for base
                }
                // finally, add half of the current layer thickness
                _layer_h[_layer_num]->derivative(d, f, p, t, val);
                m += 0.5*val;
                
                // now add the base offset
                for (unsigned int i=0; i<_layer_h.size(); i++) {
                    _layer_h[i]->derivative(d, f, p, t, val);
                    m -= 0.5*(1.+_base)*val; // currently the offset is chosen from h=0;
                }
            }
            
        protected:
            
            const Real _base;
            const unsigned int _layer_num;
            std::vector<MAST::FieldFunction<Real>*> _layer_h;
        };
        
        
        class Matrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            Matrix(std::vector<MAST::FieldFunction<RealMatrixX>*>& layer_mats):
            MAST::FieldFunction<RealMatrixX>("Matrix2D"),
            _layer_mats(layer_mats) {
                for (unsigned int i=0; i < _layer_mats.size(); i++) {
                    _functions.insert(_layer_mats[i]->master());
                }
            }
            
            Matrix(const MAST::Multilayer2DSectionProperty::Matrix &f):
            MAST::FieldFunction<RealMatrixX>(f) {
                // initialize the vector
                _layer_mats.resize(f._layer_mats.size());
                for (unsigned int i=0; i < _layer_mats.size(); i++) {
                    _layer_mats[i] = f._layer_mats[i]->clone().release();
                    _functions.insert(_layer_mats[i]->master());
                }
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Multilayer2DSectionProperty::Matrix(*this));
            }
            
            virtual ~Matrix() {
                // delete all the layer functions
                for (unsigned int i=0; i<_layer_mats.size(); i++)
                    delete _layer_mats[i];
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                // add the values of each matrix to get the integrated value
                RealMatrixX mi;
                for (unsigned int i=0; i<_layer_mats.size(); i++) {
                    (*_layer_mats[i])(p, t, mi);
                    // use the size of the layer matrix to resize the output
                    // all other layers should return the same sized matrices
                    if (i==0)
                        m.resize(mi.rows(), mi.cols());
                    
                    m += mi;
                }
            }
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                const Real t,
                                RealMatrixX& m) const {
                // add the values of each matrix to get the integrated value
                RealMatrixX mi;
                m.resize(2,2);
                for (unsigned int i=0; i<_layer_mats.size(); i++) {
                    _layer_mats[i]->derivative(d, f, p, t, mi);
                    // use the size of the layer matrix to resize the output
                    // all other layers should return the same sized matrices
                    if (i==0)
                        m.resize(mi.rows(), mi.cols());
                    
                    m += mi;
                }
            }
            
            
        protected:
            
            std::vector<MAST::FieldFunction<RealMatrixX>*> _layer_mats;
        };
        
        
    }
    
}




MAST::Multilayer2DSectionElementPropertyCard::~Multilayer2DSectionElementPropertyCard() {
    // delete the layer offset functions
    for (unsigned int i=0; i<_layer_offsets.size(); i++)
        delete _layer_offsets[i];
}



unsigned int
MAST::Multilayer2DSectionElementPropertyCard::dim() const {
    return 2;
}




void
MAST::Multilayer2DSectionElementPropertyCard::
set_layers(const Real base,
           std::vector<MAST::Solid2DSectionElementPropertyCard*>& layers) {
    
    // make sure that this has not been already set
    libmesh_assert(_layers.size() == 0);
    
    // now create the vector of offsets for each later
    const unsigned n_layers = (unsigned int)layers.size();
    _layers = layers;
    _layer_offsets.resize(n_layers);
    
    
    for (unsigned int i=0; i<n_layers; i++) {
        
        // offsets to be provided as functions to each layer
        std::vector<MAST::FieldFunction<Real>*> layer_h(n_layers);
        for (unsigned int j=0; j<n_layers; j++)
            layer_h[j] = _layers[j]->get<MAST::FieldFunction<Real> >("h").clone().release();
        
        // create the offset function
        _layer_offsets[i] =
        new MAST::Multilayer2DSectionProperty::LayerOffset
        (base, i, layer_h);
        // tell the layer about the offset
        _layers[i]->add(*_layer_offsets[i]);
    }
}



const std::vector<MAST::Solid2DSectionElementPropertyCard*>&
MAST::Multilayer2DSectionElementPropertyCard::get_layers() const {
    // make sure they have been set
    libmesh_assert(_layers.size() > 0);
    return _layers;
}



bool
MAST::Multilayer2DSectionElementPropertyCard::if_isotropic() const {
    return false;
}





bool
MAST::Multilayer2DSectionElementPropertyCard::depends_on(const MAST::FunctionBase& f) const {
    // ask each layer for the dependence
    for (unsigned int i=0; i<_layers.size(); i++)
        if (_layers[i]->depends_on(f))
            return true;
    
    // ask each offset for the dependence
    for (unsigned int i=0; i<_layer_offsets.size(); i++)
        if (_layer_offsets[i]->depends_on(f))
            return true;
    
    // if it gets here, then there is no dependence
    return false;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
stiffness_A_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->stiffness_A_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
stiffness_B_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->stiffness_B_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
stiffness_D_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->stiffness_D_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
damping_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->damping_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
inertia_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->inertia_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
thermal_expansion_A_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->thermal_expansion_A_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
thermal_expansion_B_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->thermal_expansion_B_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->transverse_shear_stiffness_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
prestress_A_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->prestress_A_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Multilayer2DSectionElementPropertyCard::
prestress_B_matrix(const MAST::ElementBase& e) {
    
    // prepare vector of matrix functions from each layer
    std::vector<MAST::FieldFunction<RealMatrixX>*> layer_mats(_layers.size());
    for (unsigned int i=0; i<_layers.size(); i++)
        layer_mats[i] = _layers[i]->prestress_B_matrix(e).release();
    
    // now create the integrated object
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > rval
    (new MAST::Multilayer2DSectionProperty::Matrix(layer_mats));
    
    return rval;
}



