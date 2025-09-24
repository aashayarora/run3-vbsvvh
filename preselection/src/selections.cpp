#include "selections.h"

RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string>& trigger_map) {
    if (trigger_map.empty()) {
        std::cerr << "Warning: No trigger map provided. Skipping trigger selection." << std::endl;
        return df_;
    }
    if (trigger_map.find(channel) == trigger_map.end()) {
        std::cerr << "Warning: Channel '" << channel << "' not found in trigger map. Skipping trigger selection." << std::endl;
        return df_;
    }

    std::string trigger_condition = trigger_map.at(channel);
    return df_.Filter(trigger_condition, "C1: Trigger Selection");
}

RNode ElectronSelections(RNode df_) {
    auto df = df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
        .Define("_looseElectrons", 
            "Electron_pt > 7 &&"
            "abs(Electron_SC_eta) < 2.5 && "
            "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
            "abs(Electron_sip3d) < 8 && "
            "Electron_cutBased >= 2 && "
            "Electron_pfRelIso03_all < 0.4 && "
            "Electron_lostHits <= 1")
        .Define("_tightElectrons", "_looseElectrons &&" 
            "Electron_pt > 30 && "
            "Electron_cutBased >= 4 && "
            "Electron_pfRelIso03_all < 0.15 && "
            "Electron_hoe < 0.1 && "
            "Electron_eInvMinusPInv > -0.04 && "
            "((abs(Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
            "Electron_convVeto == true && "
            "Electron_tightCharge == 2 && "
            "Electron_lostHits == 0")
        .Define("nElectron_Loose", "nElectron == 0 ? 0 : Sum(_looseElectrons)")
        .Define("nElectron_Tight", "nElectron_Loose == 0 ? 0 : Sum(_tightElectrons)")
        .Define("vvhTightLepMaskElectron", "_tightElectrons");
    return applyObjectMask(df, "_tightElectrons", "Electron");
}

RNode MuonSelections(RNode df_) {
    auto df = df_.Define("_looseMuons", 
            "Muon_pt > 5 && "
            "Muon_pfIsoId >= 2 && "
            "abs(Muon_eta) < 2.4 && "
            "abs(Muon_dxy) < 0.2 && "
            "abs(Muon_dz) < 0.5 && "
            "abs(Muon_sip3d) < 8 && "
            "Muon_looseId")
        .Define("_tightMuons", "_looseMuons && "
            "Muon_pt > 30 && "
            "Muon_pfIsoId > 4 && "
            "Muon_tightCharge == 2 && "
            "Muon_highPurity && "
            "Muon_tightId")
        .Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
        .Define("nMuon_Tight", "nMuon_Loose == 0 ? 0 : Sum(_tightMuons)")
        .Define("vvhTightLepMaskMuon", "_tightMuons");
    return applyObjectMask(df, "_tightMuons", "Muon");
}

RNode LeptonSelections(RNode df_) {
    auto df = ElectronSelections(df_);
    df = MuonSelections(df);
    return df.Define("Lepton_pt", "Concatenate(Electron_pt, Muon_pt)")
        .Define("_LeptonSorted", "Argsort(-Lepton_pt)")
        .Redefine("Lepton_pt", "Take(Lepton_pt, _LeptonSorted)")
        .Define("Lepton_eta", "Take(Concatenate(Electron_eta, Muon_eta), _LeptonSorted)")
        .Define("Lepton_phi", "Take(Concatenate(Electron_phi, Muon_phi), _LeptonSorted)")
        .Define("Lepton_mass", "Take(Concatenate(Electron_mass, Muon_mass), _LeptonSorted)")
        .Define("Lepton_charge", "Take(Concatenate(Electron_charge, Muon_charge), _LeptonSorted)");
}

RNode AK4JetsSelection(RNode df_) {
    auto df = df_.Define("_dR_ak4_lep", VVdR, {"Jet_eta", "Jet_phi", "Lepton_eta", "Lepton_phi"})
        .Define("_good_ak4jets", " _dR_ak4_lep > 0.4 && "
            "Jet_pt > 20 && "
            "abs(Jet_eta) < 5.0 && "
            "Jet_jetId >= 2")
        .Define("Jet_isTightBTag", "Jet_btagUParTAK4B > 0.4648")
        .Define("Jet_isMediumBTag", "Jet_btagUParTAK4B > 0.1272")
        .Define("Jet_isLooseBTag", "Jet_btagUParTAK4B > 0.0246");
    df = applyObjectMask(df, "_good_ak4jets", "Jet");
    return df;
}

RNode AK8JetsSelection(RNode df_) {
    auto df = df_.Define("_dR_ak8_lep", VVdR, {"FatJet_eta", "FatJet_phi", "Lepton_eta", "Lepton_phi"})
        .Define("_good_ak8jets", "_dR_ak8_lep > 0.8 && "
            "FatJet_pt > 250 && "
            "abs(FatJet_eta) <= 2.5 && "
            "FatJet_msoftdrop > 40 && "
            "FatJet_jetId > 0")
        .Define("FatJet_HvsQCD", "FatJet_globalParT3_Xbb / (FatJet_globalParT3_Xbb + FatJet_globalParT3_QCD)")
        .Define("FatJet_WvsQCD", "(FatJet_globalParT3_Xqq + FatJet_globalParT3_Xcs) / (FatJet_globalParT3_Xqq + FatJet_globalParT3_Xcs + FatJet_globalParT3_QCD)")
        .Define("FatJet_ZvsQCD", "(FatJet_globalParT3_Xqq + FatJet_globalParT3_Xcc + FatJet_globalParT3_Xbb) / (FatJet_globalParT3_Xqq + FatJet_globalParT3_Xcc + FatJet_globalParT3_Xbb + FatJet_globalParT3_QCD)");
    df = applyObjectMask(df, "_good_ak8jets", "FatJet");
    return df;
}

RNode OneLepBoostedAnalysis(RNode df_) {
    auto df = df_.Filter("Jet_pt.size() >= 2 && FatJet_pt.size() >= 2", "C3: At least 2 AK4 jets and 2 AK8 jet"); 

    df = df.Define("_vbs_candidate_jets", findJetPairWithMaxDeltaEta, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
            .Define("vbs_jet1_pt", "Jet_pt[_vbs_candidate_jets[0]]")
            .Define("vbs_jet1_eta", "Jet_eta[_vbs_candidate_jets[0]]")
            .Define("vbs_jet1_phi", "Jet_phi[_vbs_candidate_jets[0]]")
            .Define("vbs_jet1_mass", "Jet_mass[_vbs_candidate_jets[0]]")
            .Define("vbs_jet2_pt", "Jet_pt[_vbs_candidate_jets[1]]")
            .Define("vbs_jet2_eta", "Jet_eta[_vbs_candidate_jets[1]]")
            .Define("vbs_jet2_phi", "Jet_phi[_vbs_candidate_jets[1]]")
            .Define("vbs_jet2_mass", "Jet_mass[_vbs_candidate_jets[1]]")
            .Define("vbs_mjj", "(ROOT::Math::PtEtaPhiMVector(vbs_jet1_pt, vbs_jet1_eta, vbs_jet1_phi, vbs_jet1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(vbs_jet2_pt, vbs_jet2_eta, vbs_jet2_phi, vbs_jet2_mass)).M()")
            .Define("vbs_detajj", "abs(vbs_jet1_eta - vbs_jet2_eta)")
            .Define("vbs_candidate_found", "_vbs_candidate_jets[0] != -1 && _vbs_candidate_jets[1] != -1");

    df = df.Define("_fatjet_vbs1_dR", VdR, {"FatJet_eta", "FatJet_phi", "vbs_jet1_eta", "vbs_jet1_phi"})
            .Define("_fatjet_vbs2_dR", VdR, {"FatJet_eta", "FatJet_phi", "vbs_jet2_eta", "vbs_jet2_phi"})
            .Define("_boosted_h_candidate_jets", 
                "_fatjet_vbs1_dR >= 0.8 && "
                "_fatjet_vbs2_dR >= 0.8")
            .Define("_best_h_idx", "FatJet_HvsQCD.size() != 0 ? ArgMax(FatJet_HvsQCD[_boosted_h_candidate_jets]) : 999.0")
            .Define("boosted_h_candidate_score", "_best_h_idx != 999.0 ? FatJet_HvsQCD[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_found", "boosted_h_candidate_score > 0")
            .Define("boosted_h_candidate_eta", "boosted_h_candidate_found ? FatJet_eta[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_phi", "boosted_h_candidate_found ? FatJet_phi[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_mass", "boosted_h_candidate_found ? FatJet_mass[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_pt", "boosted_h_candidate_found ? FatJet_pt[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_tau21", "boosted_h_candidate_found ? FatJet_tau2[_boosted_h_candidate_jets][_best_h_idx] / FatJet_tau1[_boosted_h_candidate_jets][_best_h_idx] : -999.0f");

    df = df.Define("_fatjet_h_dR", VdR, {"FatJet_eta", "FatJet_phi", "boosted_h_candidate_eta", "boosted_h_candidate_phi"})
            .Define("_boosted_v_candidate_jets", 
                "_fatjet_h_dR >= 0.8 && "
                "_fatjet_vbs1_dR >= 0.8 && "
                "_fatjet_vbs2_dR >= 0.8")
            .Define("_best_w_idx", "FatJet_WvsQCD.size() != 0 ? ArgMax(FatJet_WvsQCD[_boosted_v_candidate_jets]) : -1")
            .Define("_best_z_idx", "FatJet_ZvsQCD.size() != 0 ? ArgMax(FatJet_ZvsQCD[_boosted_v_candidate_jets]) : -1")
            .Define("boosted_w_candidate_score", "_best_w_idx != -1 ? FatJet_WvsQCD[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_z_candidate_score", "_best_z_idx != -1 ? FatJet_ZvsQCD[_boosted_v_candidate_jets][_best_z_idx] : -999.0f")
            .Define("boosted_w_candidate_found", "boosted_w_candidate_score > 0 && boosted_w_candidate_score > boosted_z_candidate_score")
            .Define("boosted_z_candidate_found", "boosted_z_candidate_score > 0 && boosted_z_candidate_score > boosted_w_candidate_score")
            .Define("boosted_w_candidate_eta", "boosted_w_candidate_found ? FatJet_eta[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_w_candidate_phi", "boosted_w_candidate_found ? FatJet_phi[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_w_candidate_mass", "boosted_w_candidate_found ? FatJet_mass[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_w_candidate_pt", "boosted_w_candidate_found ? FatJet_pt[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_w_candidate_tau21", "boosted_w_candidate_found ? FatJet_tau2[_boosted_v_candidate_jets][_best_w_idx] / FatJet_tau1[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_z_candidate_eta", "boosted_z_candidate_found ? FatJet_eta[_boosted_v_candidate_jets][_best_z_idx] : -999.0f")
            .Define("boosted_z_candidate_phi", "boosted_z_candidate_found ? FatJet_phi[_boosted_v_candidate_jets][_best_z_idx] : -999.0f")
            .Define("boosted_z_candidate_mass", "boosted_z_candidate_found ? FatJet_mass[_boosted_v_candidate_jets][_best_z_idx] : -999.0f")
            .Define("boosted_z_candidate_pt", "boosted_z_candidate_found ? FatJet_pt[_boosted_v_candidate_jets][_best_z_idx] : -999.0f")
            .Define("boosted_z_candidate_tau21", "boosted_z_candidate_found ? FatJet_tau2[_boosted_v_candidate_jets][_best_z_idx] / FatJet_tau1[_boosted_v_candidate_jets][_best_z_idx] : -999.0f");

    df = df.Define("fully_reconstructed", "vbs_candidate_found && boosted_h_candidate_found && (boosted_w_candidate_found || boosted_z_candidate_found)");
    return df;
}

RNode OneLepResolvedAnalysis(RNode df_) {
    auto df = df_.Filter("Jet_pt.size() >= 4 && FatJet_pt.size() == 1", "C3: At least 4 AK4 jets and 1 AK8 jet");

    df = df.Define("_vbs_candidate_jets", findJetPairWithMaxDeltaEta, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
            .Define("vbs_jet1_pt", "Jet_pt[_vbs_candidate_jets[0]]")
            .Define("vbs_jet1_eta", "Jet_eta[_vbs_candidate_jets[0]]")
            .Define("vbs_jet1_phi", "Jet_phi[_vbs_candidate_jets[0]]")
            .Define("vbs_jet1_mass", "Jet_mass[_vbs_candidate_jets[0]]")
            .Define("vbs_jet2_pt", "Jet_pt[_vbs_candidate_jets[1]]")
            .Define("vbs_jet2_eta", "Jet_eta[_vbs_candidate_jets[1]]")
            .Define("vbs_jet2_phi", "Jet_phi[_vbs_candidate_jets[1]]")
            .Define("vbs_jet2_mass", "Jet_mass[_vbs_candidate_jets[1]]")
            .Define("vbs_mjj", "(ROOT::Math::PtEtaPhiMVector(vbs_jet1_pt, vbs_jet1_eta, vbs_jet1_phi, vbs_jet1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(vbs_jet2_pt, vbs_jet2_eta, vbs_jet2_phi, vbs_jet2_mass)).M()")
            .Define("vbs_detajj", "abs(vbs_jet1_eta - vbs_jet2_eta)")
            .Define("vbs_candidate_found", "_vbs_candidate_jets[0] != -1 && _vbs_candidate_jets[1] != -1");

    auto fatJetIndices = [](const RVec<float>& fatjet_hvsqcd, const RVec<float>& fatjet_wvsqcd, const RVec<float>& fatjet_zvsqcd) {
        if (fatjet_hvsqcd.empty()) {
            return -1;
        }
        auto fatjet_scores = {fatjet_hvsqcd[0], fatjet_wvsqcd[0], fatjet_zvsqcd[0]};
        auto max_it = std::max_element(fatjet_scores.begin(), fatjet_scores.end());
        return static_cast<int>(std::distance(fatjet_scores.begin(), max_it));
    };

    df = df.Define("fatjet_indices", fatJetIndices, {"FatJet_HvsQCD", "FatJet_WvsQCD", "FatJet_ZvsQCD"})
            .Define("boosted_h_candidate_found", "fatjet_indices == 0")
            .Define("boosted_h_candidate_eta", "boosted_h_candidate_found ? FatJet_eta[0] : -999.0f")
            .Define("boosted_h_candidate_phi", "boosted_h_candidate_found ? FatJet_phi[0] : -999.0f")
            .Define("boosted_h_candidate_mass", "boosted_h_candidate_found ? FatJet_mass[0] : -999.0f")
            .Define("boosted_h_candidate_pt", "boosted_h_candidate_found ? FatJet_pt[0] : -999.0f")
            .Define("boosted_h_candidate_tau21", "boosted_h_candidate_found ? FatJet_tau2[0] / FatJet_tau1[0] : -999.0f")
            .Define("boosted_h_candidate_score", "boosted_h_candidate_found ? FatJet_HvsQCD[0] : -999.0f")
            .Define("boosted_w_candidate_found", "fatjet_indices == 1")
            .Define("boosted_w_candidate_eta", "boosted_w_candidate_found ? FatJet_eta[0] : -999.0f")
            .Define("boosted_w_candidate_phi", "boosted_w_candidate_found ? FatJet_phi[0] : -999.0f")
            .Define("boosted_w_candidate_mass", "boosted_w_candidate_found ? FatJet_mass[0] : -999.0f")
            .Define("boosted_w_candidate_pt", "boosted_w_candidate_found ? FatJet_pt[0] : -999.0f")
            .Define("boosted_w_candidate_tau21", "boosted_w_candidate_found ? FatJet_tau2[0] / FatJet_tau1[0] : -999.0f")
            .Define("boosted_w_candidate_score", "boosted_w_candidate_found ? FatJet_WvsQCD[0] : -999.0f")
            .Define("boosted_z_candidate_found", "fatjet_indices == 2")
            .Define("boosted_z_candidate_eta", "boosted_z_candidate_found ? FatJet_eta[0] : -999.0f")
            .Define("boosted_z_candidate_phi", "boosted_z_candidate_found ? FatJet_phi[0] : -999.0f")
            .Define("boosted_z_candidate_mass", "boosted_z_candidate_found ? FatJet_mass[0] : -999.0f")
            .Define("boosted_z_candidate_pt", "boosted_z_candidate_found ? FatJet_pt[0] : -999.0f")
            .Define("boosted_z_candidate_tau21", "boosted_z_candidate_found ? FatJet_tau2[0] / FatJet_tau1[0] : -999.0f")
            .Define("boosted_z_candidate_score", "boosted_z_candidate_found ? FatJet_ZvsQCD[0] : -999.0f");

    df = df.Define("_jet_w_dR", VdR, {"Jet_eta", "Jet_phi", "boosted_w_candidate_eta", "boosted_w_candidate_phi"})
            .Define("_jet_z_dR", VdR, {"Jet_eta", "Jet_phi", "boosted_z_candidate_eta", "boosted_z_candidate_phi"})
            .Define("_jet_vbs1_dR", VdR, {"Jet_eta", "Jet_phi", "vbs_jet1_eta", "vbs_jet1_phi"})
            .Define("_jet_vbs2_dR", VdR, {"Jet_eta", "Jet_phi", "vbs_jet2_eta", "vbs_jet2_phi"})
            .Define("_resolved_h_candidate_jets", 
                "_jet_w_dR >= 0.8 && "
                "_jet_z_dR >= 0.8 && "
                "_jet_vbs1_dR >= 0.4 && "
                "_jet_vbs2_dR >= 0.4 && "
                "Jet_isTightBTag")
            .Define("_hbb_candidate_pairs", getJetPairs, {"_resolved_h_candidate_jets"})
            .Define("hbb_candidate_pairs1_pt", "Take(Jet_pt, _hbb_candidate_pairs[0], -999.0f)")
            .Define("hbb_candidate_pairs2_pt", "Take(Jet_pt, _hbb_candidate_pairs[1], -999.0f)")
            .Define("hbb_candidate_pairs1_eta", "Take(Jet_eta, _hbb_candidate_pairs[0], -999.0f)")
            .Define("hbb_candidate_pairs2_eta", "Take(Jet_eta, _hbb_candidate_pairs[1], -999.0f)")
            .Define("hbb_candidate_pairs1_phi", "Take(Jet_phi, _hbb_candidate_pairs[0], -999.0f)")
            .Define("hbb_candidate_pairs2_phi", "Take(Jet_phi, _hbb_candidate_pairs[1], -999.0f)")
            .Define("hbb_candidate_pairs1_mass", "Take(Jet_mass, _hbb_candidate_pairs[0], -999.0f)")
            .Define("hbb_candidate_pairs2_mass", "Take(Jet_mass, _hbb_candidate_pairs[1], -999.0f)")
            .Define("hbb_candidate_pairs1_btag", "Take(Jet_btagUParTAK4B, _hbb_candidate_pairs[0], -999.0f)")
            .Define("hbb_candidate_pairs2_btag", "Take(Jet_btagUParTAK4B, _hbb_candidate_pairs[1], -999.0f)")
            .Define("_hbb_candidate_mjj", "InvariantMasses(hbb_candidate_pairs1_pt, hbb_candidate_pairs1_eta, hbb_candidate_pairs1_phi, hbb_candidate_pairs1_mass, hbb_candidate_pairs2_pt, hbb_candidate_pairs2_eta, hbb_candidate_pairs2_phi, hbb_candidate_pairs2_mass)")
            .Define("_hbb_best_pair_idx", "_hbb_candidate_mjj.size() != 0 ? ArgMin(abs(_hbb_candidate_mjj - 125.0)) : -1")
            .Define("resolved_h_candidate_found", "_hbb_best_pair_idx != -1 && _hbb_candidate_mjj[_hbb_best_pair_idx] > 50")
            .Define("resolved_h_candidate_pt1", "resolved_h_candidate_found ? hbb_candidate_pairs1_pt[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_pt2", "resolved_h_candidate_found ? hbb_candidate_pairs2_pt[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_eta1", "resolved_h_candidate_found ? hbb_candidate_pairs1_eta[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_eta2", "resolved_h_candidate_found ? hbb_candidate_pairs2_eta[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_phi1", "resolved_h_candidate_found ? hbb_candidate_pairs1_phi[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_phi2", "resolved_h_candidate_found ? hbb_candidate_pairs2_phi[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_mass1", "resolved_h_candidate_found ? hbb_candidate_pairs1_mass[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_mass2", "resolved_h_candidate_found ? hbb_candidate_pairs2_mass[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_score1", "resolved_h_candidate_found ? hbb_candidate_pairs1_btag[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_score2", "resolved_h_candidate_found ? hbb_candidate_pairs2_btag[_hbb_best_pair_idx] : -999.0f")
            .Define("resolved_h_candidate_mjj", "resolved_h_candidate_found ? _hbb_candidate_mjj[_hbb_best_pair_idx] : -999.0f");

    df = df.Define("_jet_h_dR", VdR, {"Jet_eta", "Jet_phi", "boosted_h_candidate_eta", "boosted_h_candidate_phi"})
            .Define("_resolved_v_candidate_jets", 
                "_jet_h_dR >= 0.8 && "
                "_jet_vbs1_dR >= 0.4 && "
                "_jet_vbs2_dR >= 0.4")
            .Define("_v_candidate_pairs", getJetPairs, {"_resolved_v_candidate_jets"})
            .Define("v_candidate_pairs1_pt", "Take(Jet_pt, _v_candidate_pairs[0], -999.0f)")
            .Define("v_candidate_pairs2_pt", "Take(Jet_pt, _v_candidate_pairs[1], -999.0f)")
            .Define("v_candidate_pairs1_eta", "Take(Jet_eta, _v_candidate_pairs[0], -999.0f)")
            .Define("v_candidate_pairs2_eta", "Take(Jet_eta, _v_candidate_pairs[1], -999.0f)")
            .Define("v_candidate_pairs1_phi", "Take(Jet_phi, _v_candidate_pairs[0], -999.0f)")
            .Define("v_candidate_pairs2_phi", "Take(Jet_phi, _v_candidate_pairs[1], -999.0f)")
            .Define("v_candidate_pairs1_mass", "Take(Jet_mass, _v_candidate_pairs[0], -999.0f)")
            .Define("v_candidate_pairs2_mass", "Take(Jet_mass, _v_candidate_pairs[1], -999.0f)")
            .Define("_v_candidate_mjj", "InvariantMasses(v_candidate_pairs1_pt, v_candidate_pairs1_eta, v_candidate_pairs1_phi, v_candidate_pairs1_mass, v_candidate_pairs2_pt, v_candidate_pairs2_eta, v_candidate_pairs2_phi, v_candidate_pairs2_mass)")
            .Define("_w_best_pair_idx", "_v_candidate_mjj.size() != 0 ? ArgMin(abs(_v_candidate_mjj - 80.0)) : -1")
            .Define("_z_best_pair_idx", "_v_candidate_mjj.size() != 0 ? ArgMin(abs(_v_candidate_mjj - 91.0)) : -1")
            .Define("resolved_w_candidate_found", "_w_best_pair_idx != -1 && abs(_v_candidate_mjj[_w_best_pair_idx] - 80.0f) < abs(_v_candidate_mjj[_z_best_pair_idx] - 91.0f) && abs(_v_candidate_mjj[_w_best_pair_idx]) > 25.0f")
            .Define("resolved_z_candidate_found", "_z_best_pair_idx != -1 && abs(_v_candidate_mjj[_z_best_pair_idx] - 91.0f) < abs(_v_candidate_mjj[_w_best_pair_idx] - 80.0f) && abs(_v_candidate_mjj[_z_best_pair_idx]) > 25.0f")
            .Define("resolved_w_candidate_pt1", "resolved_w_candidate_found ? v_candidate_pairs1_pt[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_pt2", "resolved_w_candidate_found ? v_candidate_pairs2_pt[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_eta1", "resolved_w_candidate_found ? v_candidate_pairs1_eta[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_eta2", "resolved_w_candidate_found ? v_candidate_pairs2_eta[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_phi1", "resolved_w_candidate_found ? v_candidate_pairs1_phi[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_phi2", "resolved_w_candidate_found ? v_candidate_pairs2_phi[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_mass1", "resolved_w_candidate_found ? v_candidate_pairs1_mass[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_mass2", "resolved_w_candidate_found ? v_candidate_pairs2_mass[_w_best_pair_idx] : -999.0f")
            .Define("resolved_w_candidate_mjj", "resolved_w_candidate_found ? _v_candidate_mjj[_w_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_pt1", "resolved_z_candidate_found ? v_candidate_pairs1_pt[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_pt2", "resolved_z_candidate_found ? v_candidate_pairs2_pt[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_eta1", "resolved_z_candidate_found ? v_candidate_pairs1_eta[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_eta2", "resolved_z_candidate_found ? v_candidate_pairs2_eta[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_phi1", "resolved_z_candidate_found ? v_candidate_pairs1_phi[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_phi2", "resolved_z_candidate_found ? v_candidate_pairs2_phi[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_mass1", "resolved_z_candidate_found ? v_candidate_pairs1_mass[_z_best_pair_idx] : -999.0f")
            .Define("resolved_z_candidate_mass2", "resolved_z_candidate_found ? v_candidate_pairs2_mass[_z_best_pair_idx] : -999.0f") 
            .Define("resolved_z_candidate_mjj", "resolved_z_candidate_found ? _v_candidate_mjj[_z_best_pair_idx] : -999.0f");  

    df = df.Define("fully_reconstructed", "vbs_candidate_found && ((boosted_h_candidate_found && (resolved_w_candidate_found || resolved_z_candidate_found)) ^ ((boosted_w_candidate_found || boosted_z_candidate_found) && resolved_h_candidate_found))");

    return df;
}

RNode runPreselection(RNode df_, std::string channel, bool noCut) {
    auto df = LeptonSelections(df_);
    df = AK4JetsSelection(df);
    df = AK8JetsSelection(df);
    
    if (noCut) return df; // for spanet training data
    
    df = TriggerSelections(df, channel, TriggerMap);
    // channel-specific selections
    if (channel == "1Lep2FJ") {
        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
            "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
            "(Lepton_pt[0] > 40)", "C2: 1-lepton selection");
            df = OneLepBoostedAnalysis(df);
    }
    else if (channel == "1Lep1FJ") {
        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
            "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
            "(Lepton_pt[0] > 40)", "C2: 1-lepton selection");
            df = OneLepResolvedAnalysis(df);
    }
    else if (channel == "0Lep3FJ") {
        df = df.Filter("nMuon_Loose == 0 && nElectron_Loose == 0", "C2: 0-lepton selection");
    }

    return df;
}