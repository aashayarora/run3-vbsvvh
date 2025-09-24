#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights.h"
#include "corrections.h"
#include "utils.h"
#include "selections.h"
#include "spanet.h"
#include "genSelections.h"

#include "argparser.hpp"

struct MyArgs : public argparse::Args {
    std::string &spec = kwarg("i,input", "spec.json path");
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection").set_default("");
    std::string &output = kwarg("o,output", "output root file").set_default("");
    std::string &output_subdir = kwarg("outdir", "output project subdirectory").set_default("");
    
    int &batch_size = kwarg("b,batch_size", "batch size for spanet inference").set_default(64);
    int &nthread = kwarg("n,nthread", "number of threads for ROOT").set_default(0);

    bool &debug = flag("debug", "enable debug mode").set_default(false);
    bool &dumpInput = flag("dump_input", "Dump all input branches to output ROOT file").set_default(false);
    bool &makeSpanetTrainingdata = flag("spanet_training", "Only make training data for SPANet").set_default(false);
};

RNode runAnalysis(RNode df, std::string ana, SPANet::SPANetInference &spanet_inference, bool makeSpanetTrainingdata = false) {
    std::cout << " -> Run " << ana << "::runAnalysis()" << std::endl;
    df = runPreselection(df, ana, makeSpanetTrainingdata);
    
    if (!makeSpanetTrainingdata) {
	    // df = spanet_inference.RunSPANetInference(df);
        // df = spanet_inference.ParseSpanetInference(df);
    }
    return df;
}

int main(int argc, char** argv) {
    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;

    std::vector<std::string> channels = {"0Lep3FJ", "0Lep2FJ", "0Lep2FJMET", "1Lep2FJ", "1Lep1FJ"};
    if (std::find(channels.begin(), channels.end(), args.ana) == channels.end()) {
        if (!args.makeSpanetTrainingdata) {
            std::cerr << "Did not recognize analysis tag: " << args.ana << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // Create output directory
    std::string output_dir = setOutputDirectory(args.ana, args.output_subdir, args.makeSpanetTrainingdata);

    SPANet::SPANetInference spanet_inference("spanet/v2/model.onnx", args.batch_size);

    // add debugging
    if (args.debug) {
        std::cout << " -> Debug mode enabled" << std::endl;
        auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    }

    if (args.nthread > 1) {
        ROOT::EnableImplicitMT(args.nthread);
        ROOT::EnableThreadSafety();
    }

    // Load df
    ROOT::RDataFrame df_ = ROOT::RDF::Experimental::FromSpec(input_spec);
    ROOT::RDF::Experimental::AddProgressBar(df_);

    // Define metadata
    auto df = defineMetadata(df_);

    // Set output file name and input type
    bool isData = false;
    bool isSignal = false;
    if (input_spec.find("data") != std::string::npos) {
        isData = true;
        if (output_file.empty()) {
            output_file = "data";
        }
    }
    else if (input_spec.find("bkg") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "bkg";
        }
    }
    else if (input_spec.find("sig") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "sig";
        }
        isSignal = true;
    }
    else {
        std::cerr << "Could not guess output name from spec name, file must contain sig, bkg or data." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    bool makeSpanetTrainingdata = args.makeSpanetTrainingdata;
    if (!isSignal) {
        makeSpanetTrainingdata = false; // do not make training data for non-signal samples
    }

    // Run analysis
    if (isData) {
        std::cout << " -> Running data analysis" << std::endl;
        df = runAnalysis(df, args.ana, spanet_inference);
        df = applyDataWeights(df);
    } else {
        std::cout << " -> Running MC analysis" << std::endl;
        df = runAnalysis(df, args.ana, spanet_inference, makeSpanetTrainingdata);
        df = applyMCWeights(df);
    }

    if (isSignal) {
        df = GenSelections(df);
        if (makeSpanetTrainingdata) {
            std::cout << " -> Saving SPANet training data" << std::endl;
            saveSpanetSnapshot(df, output_dir, output_file);
            return 0; // Exit after saving training data
        }
    }

    saveSnapshot(df, output_dir, output_file, isData, args.dumpInput);

    return 0;
}
