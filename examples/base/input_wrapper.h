/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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


#ifndef __mast_input_wrapper_h__
#define __mast_input_wrapper_h__


// libMesh inclues
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

namespace MAST {
    
    namespace Examples {

        class GetPotWrapper {
            
        public:

            /*!
             *   Identifies if an input file name is provided in
             *   \p argc by the parameter name \p input. If found,
             *   the input parameters are extracted from the input file,
             *   otherwise they are obtained from \p argc.
             */
            GetPotWrapper(const int argc,
                          const char* const* argv,
                          const std::string& input):
            _input   (nullptr),
            _print   (false) {
                
                GetPot command_line(argc, argv);
                std::string
                input_name   = command_line(input,    "");

                if (input_name == "")
                    _input = new GetPot(argc, argv);
                else
                    _input = new GetPot(input_name);
                
                _print = (*_input)("print_params", false);
            }


            /*!
             *   creates a wrapper around the input arguments given to
             *   the executable in \p argc and \p argv
             */
            GetPotWrapper(const int argc, const char* const* argv):
            _input   (new GetPot(argc, argv)),
            _print   (false) {}
            
            /*!
             *   creates a wrapper for the input arguments in the input file
             *   named \p file.
             */
            GetPotWrapper(const std::string&  file):
            _input   (new GetPot(file)),
            _print   (false) {
                
            }
            
            /*!
             *   creates a wrapper for the input arguments in the input stream
             *   \p input.
             */
            GetPotWrapper(std::ifstream&     input):
            _input   (new GetPot(input)),
            _print   (false) {
                
            }

            virtual ~GetPotWrapper() {
                delete _input;
            }

            
            GetPot& get() {
                return *_input;
            }
            
            void set_print(bool f) {
                _print = f;
            }
            
            
            const char* operator()(const std::string& nm, const std::string& doc, const char* v) {
                
                if (_print)
                    libMesh::out
                    << nm << " : " << doc
                    << "    [Default: v = " << v << " ]" << std::endl;
                
                return (*_input)(nm, v);
            }

            
            template <typename ValType>
            ValType operator()(const std::string& nm, const std::string& doc, const ValType& v) {
                
                if (_print)
                    libMesh::out
                    << nm << " : " << doc
                    << "    [Default: v = " << v << " ]" << std::endl;
                
                return (*_input)(nm, v);
            }

            template <typename ValType>
            ValType operator()(const std::string& nm, const std::string& doc, const ValType& v, const unsigned int i) {
               
                
                if (_print)
                    libMesh::out
                    << "Vector: " << nm << " : " << doc
                    << "    [ Default: v(i) = " << v << " ] " << std::endl;
                
                return (*_input)(nm, v, i);
            }

            
        protected:
            
            GetPot* _input;
            bool    _print;
        };
        
        
    }
}

#endif // __mast_input_wrapper_h__
