#!/usr/bin/env python3

import logging
import sys
from datetime import datetime

from argparse import ArgumentParser
from plotter import Hist1D, Plotter

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

NUM_BINS = 30
BASE_PATH = "/data/userdata/aaarora/vbsvvhAnalysis/preselection/1Lep1FJ-bdt/"

SIGNAL_SCALE = 10000.0

hists = [
    Hist1D("Muon_pt", r"Muon $p_{T}$ [GeV]", (NUM_BINS, 0, 500), scale=SIGNAL_SCALE),
    Hist1D("Muon_eta", r"Muon $\eta$", (NUM_BINS, -2.5, 2.5), scale=SIGNAL_SCALE),
    Hist1D("Muon_phi", r"Muon $\phi$", (NUM_BINS, -3.5, 3.5), scale=SIGNAL_SCALE),

    Hist1D("Electron_pt", r"Electron $p_{T}$ [GeV]", (NUM_BINS, 0, 500), scale=SIGNAL_SCALE),
    Hist1D("Electron_eta", r"Electron $\eta$", (NUM_BINS, -2.5, 2.5), scale=SIGNAL_SCALE),
    Hist1D("Electron_phi", r"Electron $\phi$", (NUM_BINS, -3.5, 3.5), scale=SIGNAL_SCALE),
    
    Hist1D("boosted_v_candidate_score", "ParT V tagger score", (NUM_BINS, 0, 1), scale=SIGNAL_SCALE),
    Hist1D("boosted_h_candidate_score", "ParT H tagger score", (NUM_BINS, 0, 1), scale=SIGNAL_SCALE),

    Hist1D("resolved_mjj_1", r"Resolved V candidate $m_{jj}$ [GeV]", (NUM_BINS, 0, 200), scale=SIGNAL_SCALE),
    Hist1D("resolved_ptjj_1", r"Resolved V candidate $p_{T,jj}$ [GeV]", (NUM_BINS, 0, 500), scale=SIGNAL_SCALE),
    Hist1D("resolved_dR_1", r"Resolved V candidate $\Delta R_{jj}$", (NUM_BINS, 0, 2), scale=SIGNAL_SCALE),
    
    Hist1D("boosted_v_candidate_mass", "V candidate mass [GeV]", (NUM_BINS, 0, 200), scale=SIGNAL_SCALE),
    Hist1D("boosted_h_candidate_mass", "H candidate mass [GeV]", (NUM_BINS, 0, 200), scale=SIGNAL_SCALE),

    Hist1D("boosted_h_candidate_eta", "H candidate $\eta$", (NUM_BINS, -2.5, 2.5), scale=SIGNAL_SCALE),
    Hist1D("boosted_v_candidate_eta", "V candidate $\eta$", (NUM_BINS, -2.5, 2.5), scale=SIGNAL_SCALE),

    Hist1D("boosted_h_candidate_pt", r"H candidate $p_{T}$ [GeV]", (NUM_BINS, 250, 500), scale=SIGNAL_SCALE),
    Hist1D("boosted_v_candidate_pt", r"V candidate $p_{T}$ [GeV]", (NUM_BINS, 250, 500), scale=SIGNAL_SCALE),

    Hist1D("vbs_jet1_pt", r"VBS Jet 1 $p_{T}$ [GeV]", (NUM_BINS, 0, 500), scale=SIGNAL_SCALE),
    Hist1D("vbs_jet2_pt", r"VBS Jet 2 $p_{T}$ [GeV]", (NUM_BINS, 0, 200), scale=SIGNAL_SCALE),
    Hist1D("vbs_jet1_eta", r"VBS Jet 1 $\eta$", (NUM_BINS, -5, 5), scale=SIGNAL_SCALE),
    Hist1D("vbs_jet2_eta", r"VBS Jet 2 $\eta$", (NUM_BINS, -5, 5), scale=SIGNAL_SCALE),

    Hist1D("vbs_jet1_phi", r"VBS Jet 1 $\phi$", (NUM_BINS, -3.5, 3.5), scale=SIGNAL_SCALE),
    Hist1D("vbs_jet2_phi", r"VBS Jet 2  $\phi$", (NUM_BINS, -3.5, 3.5), scale=SIGNAL_SCALE),
    
    Hist1D("vbs_mjj", r"$m_{jj}$ [GeV]", (NUM_BINS, 0, 2000), scale=SIGNAL_SCALE),
    Hist1D("vbs_detajj", r"$\Delta \eta_{jj}$", (NUM_BINS, 0, 10), scale=SIGNAL_SCALE),
    Hist1D("vbs_candidate_score", "VBS tagger score", (NUM_BINS, 0, 1), scale=SIGNAL_SCALE),
]

bkg_samples_labels = {
    "TTbar": r"$t\bar{t}$",
    "WJets": "W + Jets",
    "DY": "Drell-Yan",
    "Boson": "Diboson",
    "Other": "Other",
}

sig_samples_labels = [
    r"VBS VVH ($\kappa_V$ = 2.0)",
]

SELECTION_CUT = "fully_reconstructed"

def main():
    sig_files = [
        BASE_PATH + "sig_infer_2p0.root"
    ]
    bkg_files = [
        BASE_PATH + "bkg-ttbar.root",
        BASE_PATH + "bkg-wjets.root",
        BASE_PATH + "bkg-dy.root",
        BASE_PATH + "bkg-boson.root"
    ]   
    data_files = [
        BASE_PATH + "data-2024c.root",
        BASE_PATH + "data-2024d.root",
        BASE_PATH + "data-2024e.root",
        BASE_PATH + "data-2024f.root",
        BASE_PATH + "data-2024g.root",
        BASE_PATH + "data-2024h.root",
        BASE_PATH + "data-2024i.root",
    ]
    
    try:
        logger.info("Initializing plotter...")
        plotter = Plotter(
            sig=sig_files,
            bkg=bkg_files,
            data=data_files,
            bkg_samples_labels=bkg_samples_labels,
            sig_samples_labels=sig_samples_labels,
            cut=SELECTION_CUT,
            year=2024
        )
        logger.info("Plotter initialized successfully")
    except Exception as e:
        logger.error(f"Failed to initialize plotter: {e}")
        sys.exit(1)
    
    output_configs = [
        {
            "name": "Standard plots",
            "density": False,
            "savePath": OUTPUT_DIR
        }
    ]
    
    for config in output_configs:
        logger.info("-" * 80)
        logger.info(f"Creating {config['name']}...")
        logger.info(f"Output directory: {config['savePath']}")
        
        try:
            plotter.make_plots(
                hists,
                save=True,
                density=config["density"],
                savePath=config["savePath"]
            )
            logger.info(f"✓ {config['name']} completed successfully")
        except Exception as e:
            logger.error(f"✗ Failed to create {config['name']}: {e}")
            continue
    

if __name__ == "__main__":
    parser = ArgumentParser(description="Plotter for physics analysis")
    parser.add_argument("--output", type=str, default=f"plots_{datetime.now().strftime('%Y%m%d_%H%M%S')}", help="Output directory for plots")
    args = parser.parse_args()
    OUTPUT_DIR = args.output

    try:
        main()
    except KeyboardInterrupt:
        logger.warning("Plot generation interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)
