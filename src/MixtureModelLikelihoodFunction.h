#ifndef GENEVOL_MIXTUREMODELLIKELIHOODFUNCTION_H
#define GENEVOL_MIXTUREMODELLIKELIHOODFUNCTION_H

// From bpp-core
#include <Bpp/Numeric/Function/Functions.h>

// From bpp-phyl
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>



#include "LikelihoodUtils.h"
#include "ModelParameters.h"

using namespace std;

namespace bpp {
    class AlphaLikelihoodFunction : 
        public Function,
        public AbstractParametrizable
    {
        private:
            ModelParameters* m_;
            PhyloTree* tree_;
            int categories_;
            double fval_;
        private:
            std::map<std::string, std::vector<double>> generateMMValues(int categories, double alphaGain, double betaGain, double shapeLoss) const;
            double calculateFunctionValue() const;
        public:
            AlphaLikelihoodFunction(ModelParameters* m, PhyloTree* tree) : AbstractParametrizable("") {
                m_ = m;
                tree_ = tree;
                categories_ = m->categories_;
                addParameter_(new Parameter("alphaGain0_1", m->alphaGain_, make_shared<IntervalConstraint>(0.5, 3, false, true)));
                addParameter_(new Parameter("alphaLoss0_1", m->alphaLoss_, make_shared<IntervalConstraint>(0.5, 3, false, true)));
                addParameter_(new Parameter("dupl0_1", 3, make_shared<IntervalConstraint>(0.1, 10, false, true)));
                fireParameterChanged(getParameters());
            }
            double getValue() const {return fval_;}
            Clonable* clone() const { return new AlphaLikelihoodFunction(*this); }
            void fireParameterChanged(const ParameterList& parameters) {
                setParameters(getParameters());
                fval_ = calculateFunctionValue();
            }
            void setParameters(const ParameterList& parameters) { matchParametersValues(parameters); }
            double getParameterValueByName(string name) { return getParameterValue(name); }
            std::vector<SingleProcessPhyloLikelihood*> getLikelihoodProcesses() const;

    };
}

#endif // GENEVOL_MIXTUREMODELLIKELIHOODFUNCTION_H
