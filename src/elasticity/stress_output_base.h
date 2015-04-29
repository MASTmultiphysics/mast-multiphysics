
#ifndef __mast__stress_output_base__
#define __mast__stress_output_base__

// C++ includes
#include <map>

// MAST includes
#include "base/output_function_base.h"
#include "base/mast_data_types.h"


namespace MAST {
    
    class StressStrainOutputBase:
    public MAST::OutputFunctionBase {
        
    public:
        
        StressStrainOutputBase();
        
        virtual ~StressStrainOutputBase();
        
        
        /*!
         *   sets true/false for evaluation of stress. True by default
         */
        void set_evaluate_stress(bool val);
        
        
        /*!
         *   sets true/false for evaluation of strain. False by default
         */
        void set_evaluate_strain(bool val);

        
        /*!
         *   @returns true/false based on whether or not to evaluate 
         *   the stress
         */
        bool if_evaluate_stress() const;


        /*!
         *   @returns true/false based on whether or not to evaluate
         *   the strain
         */
        bool if_evaluate_strain() const;

        
        /*!
         *   add the stress tensor associated with the qp
         */
        void add_stress_at_qp_location(const libMesh::Point& quadrature_pt,
                                       const libMesh::Point& physical_pt,
                                       const RealVectorX& vec);

        /*!
         *   add the strain tensor associated with the qp
         */
        void add_strain_at_qp_location(const libMesh::Point& quadrature_pt,
                                       const libMesh::Point& physical_pt,
                                       const RealVectorX& vec);

    protected:

        struct Data {
            RealVectorX     stress;
            libMesh::Point  quadrature_pt;
            libMesh::Point  physical_pt;
        };
        

        /*!
         *   if the stress should be evaluated
         */
        bool _if_evaluate_stress;

        
        /*!
         *   if the strain should be evaluated
         */
        bool _if_evaluate_strain;

        /*!
         *    vector of stress with the associated location details
         */
        std::vector<MAST::StressStrainOutputBase::Data> _stress_data;

        /*!
         *    vector of strain with the associated location details
         */
        std::vector<MAST::StressStrainOutputBase::Data> _strain_data;

    };
}

#endif // __mast__stress_output_base__
