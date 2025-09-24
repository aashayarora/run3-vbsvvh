#ifndef SELECTIONS_H
#define SELECTIONS_H

#pragma once

#include <unordered_map>

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

const std::unordered_map<std::string, std::string> TriggerMap = {
    {"1Lep2FJ", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24"},
    {"1Lep1FJ", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24"},
    {"0Lep3FJ", "HLT_PFHT1050"},
    // {"0Lep2FJ", "HLT_BLAH"},
    // {"0Lep2FJMET", "HLT_BLAH"},
};

RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string>& trigger_map);
RNode ElectronSelections(RNode df);
RNode MuonSelections(RNode df);
RNode AK8JetsSelection(RNode df);
RNode AK4JetsSelection(RNode df);

RNode runPreselection(RNode df_, std::string channel, bool noCut);
RNode OneLepBoostedAnalysis(RNode df);
RNode OneLepResolvedAnalysis(RNode df);

#endif // SELECTIONS_H