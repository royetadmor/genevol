#include <iostream>
#include <fstream>
#include <string>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSiteSimulator.h>

#include "ModelParameters.h"
#include "LikelihoodUtils.h"
#include "TreeUtils.h"

using namespace std;
using namespace bpp;

int main(int argc, char** argv)
{   
    BppApplication app(argc, argv, "Simulator");
    ModelParameters* m = new ModelParameters(app);
    
    // Sim-specific parameters
    int    numSites    = ApplicationTools::getIntParameter   ("_numSites",   app.getParams(), 1000,   "", true, -1);
    int    seed        = ApplicationTools::getIntParameter   ("_seed",       app.getParams(), 42,     "", true, -1);
    string outputFasta = ApplicationTools::getStringParameter("_outputFasta", app.getParams(), "simulated.fasta", "", true, -1);
    string outputTree  = ApplicationTools::getStringParameter("_outputTree",  app.getParams(), "simulated.nwk",   "", true, -1);

    // Seed the random number generator
    RandomTools::setSeed(static_cast<long>(seed));

    // Read and scale tree
    Newick reader;
    std::shared_ptr<bpp::PhyloTree> tree = std::move(reader.readPhyloTree(m->treeFilePath_));
    m->validateTree(tree);
    if (m->branchMul_ > 0.0) {
        tree->scaleTree(m->branchMul_);
    }

    // Validate WGD: number of zero-length edges must match number of q values
    vector<uint> wgdEdgeIds = TreeUtils::collectWgdEdgesInOrder(tree);
    if (wgdEdgeIds.size() != m->fixedWgdQ_.size()) {
        throw runtime_error(
            "WGD count mismatch: tree has " + to_string(wgdEdgeIds.size())
            + " zero-length edge(s) but _fixedWgdQ has "
            + to_string(m->fixedWgdQ_.size()) + " value(s).");
    }

    // Build WGD q map
    map<uint, double> wgdQMap;
    for (size_t i = 0; i < wgdEdgeIds.size(); ++i) {
        wgdQMap[wgdEdgeIds[i]] = m->fixedWgdQ_[i];
    }

    // Print summary
    cout << "===== Simulator =====" << endl;
    cout << "Tree:         " << m->treeFilePath_ << endl;
    cout << "Min state:    " << m->minState_ << endl;
    cout << "Max state:    " << m->maxState_ << endl;
    cout << "Num sites:    " << numSites << endl;
    cout << "Seed:         " << seed << endl;
    cout << "Root lambda:  " << m->rootLambda_ << endl;
    if (m->branchMul_ > 0.0) {
        cout << "Branch mul:   " << m->branchMul_ << endl;
    }
    cout << "Gain func:    " << ApplicationTools::getStringParameter("_gainFunc", app.getParams(), "CONST", "", true, -1)
         << "  params: " << VectorTools::paste(m->paramMap_.at(1), ",") << endl;
    cout << "Loss func:    " << ApplicationTools::getStringParameter("_lossFunc", app.getParams(), "CONST", "", true, -1)
         << "  params: " << VectorTools::paste(m->paramMap_.at(0), ",") << endl;
    if (!wgdEdgeIds.empty()) {
        cout << "WGD events:   " << wgdEdgeIds.size() << "  q values: ";
        for (size_t i = 0; i < m->fixedWgdQ_.size(); ++i) {
            cout << m->fixedWgdQ_[i];
            if (i + 1 < m->fixedWgdQ_.size()) cout << ",";
        }
        cout << endl;
    }
    cout << "Output FASTA: " << outputFasta << endl;
    cout << "Output tree:  " << outputTree << endl;
    cout << "=====================" << endl;

    // Build substitution process
    auto process = LikelihoodUtils::createSubstitutionProcess(
        m, tree, m->paramMap_, m->rateChangeType_, m->rDist_, wgdQMap, 0.0, m->rootLambda_);

    // Simulate
    SimpleSubstitutionProcessSiteSimulator siteSim(process);
    const vector<string>& seqNames = siteSim.getSequenceNames();

    cout << "Simulating " << numSites << " sites..." << endl;
    vector<unique_ptr<Site>> sites;
    sites.reserve(static_cast<size_t>(numSites));
    for (int i = 0; i < numSites; ++i)
        sites.push_back(siteSim.simulateSite());

    // Write simulated data as TSV (matches genevol input format)
    ofstream out(outputFasta);
    if (!out.is_open())
        throw runtime_error("Could not open output file: " + outputFasta);

    // Header row: Organizem \t Family1 \t Family2 ...
    out << "Organizem";
    for (int i = 1; i <= numSites; ++i)
        out << "\tFamily" << i;
    out << "\n";

    // One row per species
    for (size_t j = 0; j < seqNames.size(); ++j) {
        out << seqNames[j];
        for (int i = 0; i < numSites; ++i)
            out << "\t" << (*sites[i])[static_cast<int>(j)];
        out << "\n";
    }
    cout << "Written: " << outputFasta << endl;

    // Collapse WGD node triplets back into single branches before writing
    TreeUtils::collapseWgdNodes(tree);
    TreeUtils::writeTree(tree, outputTree);

    delete m;
    return 0;
}