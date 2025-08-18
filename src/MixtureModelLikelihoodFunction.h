#ifndef GENEVOL_MIXTUREMODELLIKELIHOODFUNCTION_H
#define GENEVOL_MIXTUREMODELLIKELIHOODFUNCTION_H

// From bpp-core
#include <Bpp/Numeric/Function/Functions.h>

// From bpp-phyl
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>


// Local modules
#include "ChromosomeSubstitutionModel.h"
#include "LikelihoodUtils.h"
#include "ModelParameters.h"

using namespace std;

namespace bpp {
    class MixtureModelLikelihoodFunction : 
        public Function,
        public AbstractParametrizable
    {
        private:
            ModelParameters* m_;
            PhyloTree* tree_;
            int categories_;
            double fval_;
            std::map<int, std::vector<std::string>> rateParams;
        private:
            std::map<std::string, std::vector<double>> generateMMValues(int categories, double alphaGain, double betaGain, double alphaLoss, double betaLoss) const;
            double calculateFunctionValue() const;
            std::vector<double> getRateParametersByEvent(int eventType) const;
            void setRateParameters() {
                auto rateChangeType = m_->mixtureRateChangeType_;
                for (int i = 0; i < rateChangeType.size(); i++)
                {
                    if (rateChangeType[i] == ChromosomeNumberDependencyFunction::IGNORE){
                        continue;
                    }
                    int funcType = rateChangeType[i];
                    auto it = m_->expectedNumOfParams.find(funcType);
                    if (it == m_->expectedNumOfParams.end()) {
                        throw std::runtime_error("Unknown function type enum: " + std::to_string(funcType));
                    }
                    int requiredParams = it->second;
                    for (int j = 1; j < requiredParams; j++)
                    {
                        int eventIndex = i+1;
                        string paramName = "rate" + m_->eventTypeToString.at(eventIndex) + std::to_string(j) + "_1";
                        rateParams[eventIndex].push_back(paramName);
                        addParameter_(new Parameter(paramName, 1, make_shared<IntervalConstraint>(0.5, 5, false, true))); //TODO: fix hardcoded interval
                    }
                    
                }
                
            }
        public:
            MixtureModelLikelihoodFunction(ModelParameters* m, PhyloTree* tree) : AbstractParametrizable("") {
                m_ = m;
                tree_ = tree;
                categories_ = m->categories_;
                // Basic parameters
                // TODO: fix hardcoded intervals
                addParameter_(new Parameter("alphaGain0_1", m->alphaGain_, make_shared<IntervalConstraint>(0.5, 45, false, true)));
                addParameter_(new Parameter("betaGain0_1", m->betaGain_, make_shared<IntervalConstraint>(0.5, 5, false, true)));
                addParameter_(new Parameter("alphaLoss0_1", m->alphaLoss_, make_shared<IntervalConstraint>(0.5, 45, false, true)));
                addParameter_(new Parameter("betaLoss0_1", m->betaLoss_, make_shared<IntervalConstraint>(0.5, 5, false, true)));
                
                // Dependency function parameters 
                setRateParameters();
                fireParameterChanged(getParameters());
            }
            double getValue() const {return fval_;}
            Clonable* clone() const { return new MixtureModelLikelihoodFunction(*this); }
            void fireParameterChanged(const ParameterList& parameters) {
                setParameters(getParameters());
                fval_ = calculateFunctionValue();
            }
            void setParameters(const ParameterList& parameters) { matchParametersValues(parameters); }
            double getParameterValueByName(string name) { return getParameterValue(name); }
            int getParametersCount() { return getParameters().size(); }
            std::vector<SingleProcessPhyloLikelihood*> getLikelihoodProcesses() const;

    };
}

#endif // GENEVOL_MIXTUREMODELLIKELIHOODFUNCTION_H
