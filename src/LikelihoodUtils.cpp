#include "LikelihoodUtils.h"

using namespace bpp;
using namespace std;


SingleProcessPhyloLikelihood* LikelihoodUtils::createLikelihoodProcess(ModelParameters* m, PhyloTree* tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType) {

    // Create substitution process
    int baseNum = -999;
    bool demiOnlyForEven = false;
    auto mapOfNodeIds = LikelihoodUtils::getMapOfNodeIds(tree);
    std::shared_ptr<DiscreteDistribution> gammaDist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    auto parTree = std::make_shared<ParametrizablePhyloTree>(tree);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(gammaDist, parTree);
    
    // Create substitution model
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(m->alphabet_, rateParams, baseNum, m->countRange_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, demiOnlyForEven);
    subProcesses->addModel(chrModel, mapOfNodeIds[1]);
    SubstitutionProcess* nsubPro = subProcesses->clone();

    // Create likelihood object
    Context* context = new Context();
    bool weightedRootFreqs = true;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *m->container_->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;
}


std::map<uint, vector<uint>> LikelihoodUtils::getMapOfNodeIds(PhyloTree* tree) {
    std::map<uint, vector<uint>> mapModelNodesIds;
    auto nodes = tree->getAllNodes();
    for (size_t i = 0; i < nodes.size(); i++){
        uint nodeId = tree->getNodeIndex(nodes[i]);
        if (nodeId == tree->getRootIndex()){
            continue;
        }
        mapModelNodesIds[1].push_back(nodeId);
    }
    // mapModelNodesIds[1].push_back(tree->getRootIndex());
    return mapModelNodesIds;
}

void LikelihoodUtils::deleteLikelihoodProcess(SingleProcessPhyloLikelihood* lik) {
    auto sequenceData = lik->getData();
    auto process = &(lik->getSubstitutionProcess());
    auto context = &(lik->getContext());
    delete process;
    delete sequenceData;
    delete context;
    delete lik;
}
