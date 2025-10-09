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
        .Define("FatJet_VvsQCD", "(FatJet_globalParT3_Xqq + FatJet_globalParT3_Xcs) / (FatJet_globalParT3_Xqq + FatJet_globalParT3_Xcs + FatJet_globalParT3_QCD)");
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
            .Define("_best_w_idx", "FatJet_VvsQCD.size() != 0 ? ArgMax(FatJet_VvsQCD[_boosted_v_candidate_jets]) : -1")
            .Define("boosted_v_candidate_score", "_best_w_idx != -1 ? FatJet_VvsQCD[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_found", "boosted_v_candidate_score > 0")
            .Define("boosted_v_candidate_eta", "boosted_v_candidate_found ? FatJet_eta[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_phi", "boosted_v_candidate_found ? FatJet_phi[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_mass", "boosted_v_candidate_found ? FatJet_mass[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_pt", "boosted_v_candidate_found ? FatJet_pt[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_tau21", "boosted_v_candidate_found ? FatJet_tau2[_boosted_v_candidate_jets][_best_w_idx] / FatJet_tau1[_boosted_v_candidate_jets][_best_w_idx] : -999.0f");

    df = df.Define("fully_reconstructed", "vbs_candidate_found && boosted_h_candidate_found && boosted_v_candidate_found");
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

    df = df.Define("boosted_h_candidate_found", "FatJet_HvsQCD[0] > FatJet_VvsQCD[0]")
            .Define("boosted_v_candidate_found", "FatJet_VvsQCD[0] > FatJet_HvsQCD[0]")
            .Define("boosted_h_candidate_eta", "boosted_h_candidate_found ? FatJet_eta[0] : -999.0f")
            .Define("boosted_h_candidate_phi", "boosted_h_candidate_found ? FatJet_phi[0] : -999.0f")
            .Define("boosted_h_candidate_mass", "boosted_h_candidate_found ? FatJet_mass[0] : -999.0f")
            .Define("boosted_h_candidate_pt", "boosted_h_candidate_found ? FatJet_pt[0] : -999.0f")
            .Define("boosted_h_candidate_tau21", "boosted_h_candidate_found ? FatJet_tau2[0] / FatJet_tau1[0] : -999.0f")
            .Define("boosted_h_candidate_score", "boosted_h_candidate_found ? FatJet_HvsQCD[0] : -999.0f")
            .Define("boosted_v_candidate_eta", "boosted_v_candidate_found ? FatJet_eta[0] : -999.0f")
            .Define("boosted_v_candidate_phi", "boosted_v_candidate_found ? FatJet_phi[0] : -999.0f")
            .Define("boosted_v_candidate_mass", "boosted_v_candidate_found ? FatJet_mass[0] : -999.0f")
            .Define("boosted_v_candidate_pt", "boosted_v_candidate_found ? FatJet_pt[0] : -999.0f")
            .Define("boosted_v_candidate_tau21", "boosted_v_candidate_found ? FatJet_tau2[0] / FatJet_tau1[0] : -999.0f")
            .Define("boosted_v_candidate_score", "boosted_v_candidate_found ? FatJet_VvsQCD[0] : -999.0f");

    df = df.Define("jet_v_dR", VdR, {"Jet_eta", "Jet_phi", "boosted_v_candidate_eta", "boosted_v_candidate_phi"})
            .Define("jet_h_dR", VdR, {"Jet_eta", "Jet_phi", "boosted_h_candidate_eta", "boosted_h_candidate_phi"})
            .Define("jet_vbs1_dR", VdR, {"Jet_eta", "Jet_phi", "vbs_jet1_eta", "vbs_jet1_phi"})
            .Define("jet_vbs2_dR", VdR, {"Jet_eta", "Jet_phi", "vbs_jet2_eta", "vbs_jet2_phi"})
            .Define("_resolved_candidate_jets",
                "Jet_eta <= 2.5 && "
                "jet_v_dR >= 0.8 && "
                "jet_h_dR >= 0.8 && "
                "jet_vbs1_dR >= 0.4 && "
                "jet_vbs2_dR >= 0.4 && "
                "Jet_isTightBTag")
            .Define("_resolved_candidate_pt", "Jet_pt[_resolved_candidate_jets]")
            .Define("_resolved_candidate_eta", "Jet_eta[_resolved_candidate_jets]")
            .Define("_resolved_candidate_phi", "Jet_phi[_resolved_candidate_jets]")
            .Define("_resolved_candidate_mass", "Jet_mass[_resolved_candidate_jets]")
            .Define("_resolved_candidate_btag", "Jet_btagUParTAK4B[_resolved_candidate_jets]")
            .Define("_resolved_candidate_pairs", getJetPairs, {"_resolved_candidate_pt"})
            .Define("_resolved_candidate_pairs1_pt", "Take(_resolved_candidate_pt, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_pt", "Take(_resolved_candidate_pt, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_eta", "Take(_resolved_candidate_eta, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_eta", "Take(_resolved_candidate_eta, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_phi", "Take(_resolved_candidate_phi, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_phi", "Take(_resolved_candidate_phi, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_mass", "Take(_resolved_candidate_mass, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_mass", "Take(_resolved_candidate_mass, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_btag", "Take(_resolved_candidate_btag, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_btag", "Take(_resolved_candidate_btag, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_mjj", "InvariantMasses(_resolved_candidate_pairs1_pt, _resolved_candidate_pairs1_eta, _resolved_candidate_pairs1_phi, _resolved_candidate_pairs1_mass, _resolved_candidate_pairs2_pt, _resolved_candidate_pairs2_eta, _resolved_candidate_pairs2_phi, _resolved_candidate_pairs2_mass)")
            .Define("_resolved_candidate_dR", "DeltaR(_resolved_candidate_pairs1_eta, _resolved_candidate_pairs2_eta, _resolved_candidate_pairs1_phi, _resolved_candidate_pairs2_phi)")
            .Define("_sorted_resolved_dR", "Argsort(-_resolved_candidate_dR)")
            .Define("resolved_candidate_mjj", "Take(_resolved_candidate_mjj, _sorted_resolved_dR)")
            .Define("resolved_candidate_btag1", "Take(_resolved_candidate_pairs1_btag, _sorted_resolved_dR)")
            .Define("resolved_candidate_btag2", "Take(_resolved_candidate_pairs2_btag, _sorted_resolved_dR)");

    df = df.Define("fully_reconstructed", "vbs_candidate_found && (boosted_h_candidate_found || boosted_v_candidate_found) && (resolved_candidate_mjj.size() > 0 && resolved_candidate_mjj[0] > 0)");

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