#include "utils.h"

/*
############################################
RDF UTILS
############################################
*/

RNode defineMetadata(RNode df) {
    return df.DefinePerSample("xsec", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("xsec");})
        .DefinePerSample("lumi", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("lumi");})
        .DefinePerSample("nevents", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("nevents");})
        .DefinePerSample("category", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("category");})
        .DefinePerSample("type", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("type");})
        .DefinePerSample("year", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("year");})
        .DefinePerSample("namewithyear", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("namewithyear");})
        .Define("xsec_weight", "1000 * xsec * lumi / nevents")
        .Define("isData", "category == \"data\"")
        .Define("is2016", "year == \"2016preVFP\" || year == \"2016postVFP\"")
        .Define("is2017", "year == \"2017\"")
        .Define("is2018", "year == \"2018\"");
}

RNode removeDuplicates(RNode df){
    return df.Filter(FilterOnePerKind(), {"run", "luminosityBlock", "event"}, "REMOVED DUPLICATES");
}

RNode applyObjectMask(RNode df, const std::string& maskName, const std::string& objectName) {
    auto columnNames = df.GetColumnNames();
    for (const auto& colName : columnNames) {
        if (colName.starts_with(objectName + "_")) {
            df = df.Redefine(colName, colName + "[" + maskName + "]");
        }
    }
    df = df.Redefine("n" + objectName, "Sum(" + maskName + ")");
    return df;
}

/*
############################################
LUMIMASK - GOLDEN JSON
############################################
*/

bool operator< ( const lumiMask::LumiBlockRange& lh, const lumiMask::LumiBlockRange& rh )
{
    return ( lh.run() == rh.run() ) ? ( lh.lastLumi() < rh.firstLumi() ) : lh.run() < rh.run();
}

lumiMask lumiMask::fromJSON(const std::vector<std::string>& files, lumiMask::Run firstRun, lumiMask::Run lastRun)
{
  const bool noRunFilter = ( firstRun == 0 ) && ( lastRun == 0 );

  std::vector<lumiMask::LumiBlockRange> accept;

  for ( const auto& file : files ) {
    boost::property_tree::ptree ptree;
    boost::property_tree::read_json(file, ptree);
    for ( const auto& runEntry : ptree ) {
        const lumiMask::Run run = std::stoul(runEntry.first);
        if ( noRunFilter || ( ( firstRun <= run ) && ( run <= lastRun ) ) ) {
        for ( const auto& lrEntry : runEntry.second ) {
            const auto lrNd = lrEntry.second;
            const lumiMask::LumiBlock firstLumi = std::stoul(lrNd.begin()->second.data());
            const lumiMask::LumiBlock lastLumi  = std::stoul((++lrNd.begin())->second.data());
            accept.emplace_back(run, firstLumi, lastLumi);
        }
        }
    }
  }
  
  return lumiMask(accept);
}

/*
############################################
SELECTION UTILS
############################################
*/

float fdR(float eta1, float phi1, float eta2, float phi2) {
    return ROOT::VecOps::DeltaR(eta1, eta2, phi1, phi2);
}

RVec<float> VdR(const RVec<float>& vec_eta, const RVec<float>& vec_phi, float obj_eta, float obj_phi) {
    RVec<float> out(vec_eta.size());
    if (obj_eta == -999 || obj_phi == -999) {
        std::fill(out.begin(), out.end(), 10.0f);
        return out;
    }
    for (size_t i = 0; i < vec_eta.size(); i++) {
        out[i] = ROOT::VecOps::DeltaR(vec_eta[i], obj_eta, vec_phi[i], obj_phi);
    }
    return out;
}

RVec<float> VVdR(const RVec<float>& vec_eta1, const RVec<float>& vec_phi1, const RVec<float>& vec_eta2, const RVec<float>& vec_phi2) {
    if (vec_eta1.empty()) {
        return RVec<float>();
    }
    if (vec_eta2.empty()) {
        return RVec<float>(vec_eta1.size(), 999.0f);
    }
    RVec<float> out(vec_eta1.size());
    for (size_t i = 0; i < vec_eta1.size(); i++) {
        float mindR = 999.;
        for (size_t j = 0; j < vec_eta2.size(); j++) {
            float dR = ROOT::VecOps::DeltaR(vec_eta1[i], vec_eta2[j], vec_phi1[i], vec_phi2[j]);
            if (dR < mindR) {
                mindR = dR;
            }
        }
        out[i] = mindR;
    }
    return out;
}

float fInvariantMass(float obj1_pt, float obj1_eta, float obj1_phi, float obj1_mass, 
                    float obj2_pt, float obj2_eta, float obj2_phi, float obj2_mass) {
    TLorentzVector obj1, obj2;
    obj1.SetPtEtaPhiM(obj1_pt, obj1_eta, obj1_phi, obj1_mass);
    obj2.SetPtEtaPhiM(obj2_pt, obj2_eta, obj2_phi, obj2_mass);
    return (obj1 + obj2).M();
}

RVec<float> VInvariantMass(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                          const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invMass(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invMass[i] = (obj1 + obj2).M();
    }
    return invMass;
}

RVec<float> VInvariantPt(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                        const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invPt(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invPt[i] = (obj1 + obj2).Pt();
    }
    return invPt;
}

RVec<float> VInvariantPhi(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                         const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invPhi(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invPhi[i] = (obj1 + obj2).Phi();
    }
    return invPhi;
}

RVec<float> VTransverseMass(const RVec<float>& vec_pt, const RVec<float>& vec_phi, float obj_pt, float obj_phi) {
    RVec<float> mt(vec_pt.size());
    for (size_t i = 0; i < vec_pt.size(); i++) {
        mt[i] = std::sqrt(2 * vec_pt[i] * obj_pt * (1 - std::cos(ROOT::VecOps::DeltaPhi(vec_phi[i], obj_phi))));
    }
    return mt;
}

// Return for each ak4 jet, the dR from the closest ak8 jet
RVec<float> dRfromClosestJet(const RVec<float>& ak4_eta, const RVec<float>& ak4_phi, const RVec<float>& ak8_eta, const RVec<float>& ak8_phi) {
    RVec<float> vec_minDR = {};
    for (size_t i = 0; i < ak4_eta.size(); i++)
    {
        float mindR = 999.;
        for (size_t j = 0; j < ak8_eta.size(); j++)
        {
            float dR = ROOT::VecOps::DeltaR(ak4_eta.at(i), ak8_eta.at(j), ak4_phi.at(i), ak8_phi.at(j));
            if (dR < mindR) {
                mindR = dR;
            }
        }
        vec_minDR.push_back(mindR);
    }
    return vec_minDR;
}

RVec<RVec<int>> getJetPairs(const RVec<int>& goodJets) {
    if (Sum(goodJets) >= 2) {
        return ROOT::VecOps::Combinations(goodJets, 2);
    } else {
        RVec<RVec<int>> result;
        result.emplace_back(RVec<int>{999});
        result.emplace_back(RVec<int>{999});
        return result;
    }
}

RVec<int> findJetPairWithMaxDeltaEta(const RVec<float>& jet_pt, const RVec<float>& jet_eta, 
                                     const RVec<float>& jet_phi, const RVec<float>& jet_mass) {
    RVec<int> selected_jet_indices = {};
    int first_jet_idx = -1;
    int second_jet_idx = -1;
    float max_delta_eta = 0;
    
    for (size_t i = 0; i < jet_eta.size(); i++) {
        for (size_t j = i + 1; j < jet_eta.size(); j++) {
            float delta_eta = std::abs(jet_eta[i] - jet_eta[j]);
            if (delta_eta > max_delta_eta) {
                max_delta_eta = delta_eta;
                first_jet_idx = i;
                second_jet_idx = j;
            }
        }
    }
    
    selected_jet_indices.push_back(first_jet_idx);
    selected_jet_indices.push_back(second_jet_idx);
    
    std::sort(selected_jet_indices.begin(), selected_jet_indices.end(), 
              [&jet_pt](int i, int j) { return jet_pt[i] > jet_pt[j]; });
              
    return selected_jet_indices;
}

/*
############################################
SNAPSHOT
############################################
*/

std::string setOutputDirectory(const std::string &ana, const std::string &output_subdir, bool spanet_training) {
    // If on UAF and the USER environment variable is defined, store output on ceph
    const char* userEnv = getenv("USER");
    std::string storage_dir = "/data/userdata/";
    std::string output_dir = "./";
    if (spanet_training) {
        output_dir = storage_dir + std::string(userEnv) + "/vbsvvhAnalysis/spanet_training/";
    }
    else if (userEnv != nullptr && std::filesystem::exists(storage_dir) && std::filesystem::is_directory(storage_dir)) {
        output_dir = storage_dir + std::string(userEnv) + "/vbsvvhAnalysis/preselection/" + ana + "/";
    }

    if (!output_subdir.empty()) {
        output_dir += "/" + output_subdir + "/";
    }

    std::filesystem::path directory_path(output_dir);
    // Check if the directory exists
    if (std::filesystem::exists(directory_path)) {
        std::cerr << "Output directory already exists: " << directory_path << std::endl;
    }
    // Try to create the directory and any missing parent directories
    else if (std::filesystem::create_directories(directory_path)) {
        std::cout << "Created output directory : " << directory_path << std::endl;
    }
    else {
        std::cerr << "Failed to create output directory: " << directory_path << std::endl;
        std::exit(EXIT_FAILURE); 
    }

    return directory_path;
}

void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool isData, bool dumpInput)
{
    auto ColNames = df.GetDefinedColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("event");

    for (auto &&ColName : ColNames) {
        if (ColName.starts_with("_") && !ColName.starts_with("_cut")) {
            continue;
        }
        final_variables.push_back(ColName);
    }

    // add LHE info
    if (!isData) {
        final_variables.push_back("LHEReweightingWeight");
        final_variables.push_back("nLHEReweightingWeight");
    }

    // store all columns from input nanoAOD tree
    if (dumpInput) {
        auto nanoColNames = df.GetColumnNames();
        for (auto &&colName : nanoColNames) {
            if ((std::find(final_variables.begin(), final_variables.end(), colName) == final_variables.end()) &&
                (colName.find("HLT") == std::string::npos) && (colName.find("L1") == std::string::npos)) {
                final_variables.push_back(colName);
            }
        }
    }

    std::string outputFile = outputDir + "/" + outputFileName + ".root";
    df.Snapshot("Events", outputFile, final_variables);
    std::cout << " -> Stored output file: " << outputFile << std::endl;
}


void saveSpanetSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName)
{
    auto ColNames = df.GetColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("event");

    for (auto &&ColName : ColNames) {
        if (ColName.starts_with("Jet_") || 
            ColName.starts_with("FatJet_") ||
            ColName.starts_with("PuppiMET_") ||
            ColName.starts_with("GenPart_") ||  
            ColName.starts_with("Lepton_") ||
            ColName.starts_with("gen_") || 
            ColName.starts_with("truth_")) {
                final_variables.push_back(ColName);
            }
    }

    std::string outputFile = outputDir + "/" + outputFileName + "_spanet_training_data.root";
    df.Snapshot("Events", outputFile, final_variables);
    std::cout << " -> Stored output file: " << outputFile << std::endl;
}
