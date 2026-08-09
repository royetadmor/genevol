#include <iostream>
#include <fstream>
#include <string>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSiteSimulator.h>
#include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

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
    int    numSites    = ApplicationTools::getIntParameter    ("_numSites",    app.getParams(), 1000,             "", true, -1);
    string outputFasta = ApplicationTools::getStringParameter ("_outputFasta", app.getParams(), "simulated.fasta", "", true, -1);
    string outputTree  = ApplicationTools::getStringParameter ("_outputTree",  app.getParams(), "simulated.nwk",   "", true, -1);
    bool   debugOutput = ApplicationTools::getBooleanParameter("_debugOutput", app.getParams(), false,            "", true, -1);

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

    // Collect non-root node IDs and labels for debug output
    uint rootId = tree->getRootIndex();
    vector<uint>   debugNodeIds;
    vector<string> debugNodeLabels;
    if (debugOutput) {
        for (auto& node : tree->getAllNodes()) {
            uint id = tree->getNodeIndex(node);
            if (id == rootId) continue;
            debugNodeIds.push_back(id);
            string name;
            try { name = node->getName(); } catch (...) {}
            if (name.empty()) name = "node_" + to_string(id);
            debugNodeLabels.push_back(name);
        }
    }

    cout << "Simulating " << numSites << " sites..." << endl;
    vector<unique_ptr<SiteInterface>> sites;
    sites.reserve(static_cast<size_t>(numSites));

    // Per-site ancestral (root) state; per-branch event counts [branch][site]
    // Event matrix: total substitution counts summed across all sites and branches
    vector<size_t>          ancestralStates;
    vector<vector<size_t>>  branchCounts;
    size_t alphabetSize = static_cast<size_t>(m->alphabet_->getSize());
    RowMatrix<double>       eventMatrix(alphabetSize, alphabetSize);
    if (debugOutput) {
        ancestralStates.reserve(static_cast<size_t>(numSites));
        branchCounts.assign(debugNodeIds.size(), vector<size_t>(static_cast<size_t>(numSites), 0));
    }

    for (int i = 0; i < numSites; ++i) {
        if (debugOutput) {
            auto result = siteSim.dSimulateSite();
            ancestralStates.push_back(result->getAncestralState(rootId));
            for (size_t b = 0; b < debugNodeIds.size(); ++b) {
                branchCounts[b][static_cast<size_t>(i)] = result->getSubstitutionCount(debugNodeIds[b]);
                result->getMutationPath(debugNodeIds[b]).getEventCounts(eventMatrix);
            }
            sites.push_back(result->getSite());
        } else {
            sites.push_back(siteSim.simulateSite());
        }
    }

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

    // Write debug files if requested
    if (debugOutput) {
        string outputDir;
        size_t lastSlash = outputFasta.find_last_of("/\\");
        outputDir = (lastSlash != string::npos) ? outputFasta.substr(0, lastSlash) : ".";

        // Ancestral states: one row per site
        string ancestralFile = outputDir + "/ancestral_states.tsv";
        ofstream aout(ancestralFile);
        if (!aout.is_open())
            throw runtime_error("Could not open output file: " + ancestralFile);
        aout << "site_id\troot_state\n";
        for (int i = 0; i < numSites; ++i)
            aout << (i + 1) << "\t" << ancestralStates[static_cast<size_t>(i)] << "\n";
        cout << "Written: " << ancestralFile << endl;

        // Branch event counts: rows = branches, columns = sites
        string branchFile = outputDir + "/branch_event_counts.tsv";
        ofstream bout(branchFile);
        if (!bout.is_open())
            throw runtime_error("Could not open output file: " + branchFile);
        bout << "branch";
        for (int i = 1; i <= numSites; ++i)
            bout << "\tFamily" << i;
        bout << "\n";
        for (size_t b = 0; b < debugNodeIds.size(); ++b) {
            bout << debugNodeLabels[b];
            for (int i = 0; i < numSites; ++i)
                bout << "\t" << branchCounts[b][static_cast<size_t>(i)];
            bout << "\n";
        }
        cout << "Written: " << branchFile << endl;

        // Event matrix: total transitions summed across all sites and branches (sparse, non-zero only)
        string matrixFile = outputDir + "/event_matrix.tsv";
        ofstream mout(matrixFile);
        if (!mout.is_open())
            throw runtime_error("Could not open output file: " + matrixFile);
        mout << "from_state\tto_state\tcount\n";
        for (size_t from = 0; from < alphabetSize; ++from) {
            for (size_t to = 0; to < alphabetSize; ++to) {
                double count = eventMatrix(from, to);
                if (count > 0)
                    mout << (m->minState_ + from) << "\t" << (m->minState_ + to) << "\t" << count << "\n";
            }
        }
        cout << "Written: " << matrixFile << endl;
    }

    // Collapse WGD node triplets back into single branches before writing
    TreeUtils::collapseWgdNodes(tree);
    TreeUtils::writeTree(tree, outputTree);

    delete m;
    return 0;
}