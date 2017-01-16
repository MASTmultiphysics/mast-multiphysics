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
#include "aeroelasticity/flutter_root_base.h"


MAST::FlutterRootBase::FlutterRootBase():
has_sensitivity_data  (false),
if_nonphysical_root   (false),
kr                    (0.),
g                     (0.),
kr_sens               (0.),
V                     (0.),
omega                 (0.),
V_sens                (0.),
root                  (0.),
root_sens             (0.)
{ }



MAST::FlutterRootBase::FlutterRootBase(const FlutterRootBase& f):
has_sensitivity_data   (f.has_sensitivity_data),
if_nonphysical_root    (f.if_nonphysical_root),
kr                     (f.kr),
g                      (f.g),
kr_sens                (f.kr_sens),
V                      (f.V),
omega                  (f.omega),
V_sens                 (f.V_sens),
root                   (f.root),
root_sens              (f.root_sens),
eig_vec_right          (f.eig_vec_right),
eig_vec_left           (f.eig_vec_left),
modal_participation    (f.modal_participation)
{ }



void
MAST::FlutterRootBase::copy_root(const MAST::FlutterRootBase& f) {
    
    has_sensitivity_data   = f.has_sensitivity_data;
    if_nonphysical_root    = f.if_nonphysical_root;
    kr                     = f.kr;
    g                      = f.g;
    kr_sens                = f.kr_sens;
    V                      = f.V;
    omega                  = f.omega;
    V_sens                 = f.V_sens;
    root                   = f.root;
    root_sens              = f.root_sens;
    eig_vec_right          = f.eig_vec_right;
    eig_vec_left           = f.eig_vec_left;
    modal_participation    = f.modal_participation;
}


