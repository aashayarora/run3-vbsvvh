from glob import glob
import os
import json
import re
import argparse
import uproot
import concurrent.futures

class Config:
    def __init__(self, sample_category : str, channel : str, samples: str, xsecs: dict, nthreads: int = 8):
        self.sample_category = sample_category
        self.samples = sorted(glob(samples))
        self.config = {"samples": {}}

        self.process_samples(xsecs, nthreads)
        self.write_config(f"{channel}-{sample_category}.json")

    def write_config(self, output_file):
        with open(output_file, "w") as f:
            json.dump(self.config, f, indent=4)

    @staticmethod
    def extract_sample_year(sample):
        # Data
        if any(run in sample for run in ["Run2022C", "Run2022D"]):
            return "2022Re-recoBCD"
        elif any(run in sample for run in ["Run2022E", "Run2022F", "Run2022G"]):
            return "2022Re-recoE+PromptFG"
        # MC
        elif "22EE" in sample:
            return "2022Re-recoE+PromptFG"
        elif "22" in sample:
            return "2022Re-recoBCD"
        elif "24" in sample:
            return "2024Prompt"
        elif "UL16" in sample and "APV" in sample:
            return "2016preVFP"
        elif "UL16" in sample and not "APV" in sample:
            return "2016postVFP"
        elif "UL17" in sample or "UL2017" in sample or "Run2017" in sample:
            return "2017"
        elif "UL18" in sample or "UL2018" in sample or "Run2018" in sample:
            return "2018"
        elif ("Run2016B" in sample or "Run2016C" in sample or "Run2016D" in sample or "Run2016E" in sample or "Run2016F" in sample) and "HIPM" in sample:
            return "2016preVFP"
        elif ("Run2016F" in sample or "Run2016G" in sample or "Run2016H" in sample) and not "HIPM" in sample:
            return "2016postVFP"
        else:
            raise ValueError(f"Error: year not found for {sample}")


    # From a dataset name, get the short version (as defined in the xsec dict)
    @staticmethod
    def get_sample_name_and_xsec(dataset_name,xsec_dict,is_data):
        # Data is data
        if is_data:
            return("data",1)

        # The short name for this dataset, as defined in the xsec dict
        dataset_name_short = ""

        # Loop through the xsec dict and get the name that matches this dataset
        # Raise error if no matches, or if more than one matches
        match_xsec_name = 0
        for xsec_name in xsec_dict:
            if dataset_name.startswith(xsec_name):
                match_xsec_name += 1
                if match_xsec_name > 1:
                    raise Exception(f"More than one xsec name matches the dataset \"{dataset_name}\"")
                else:
                    dataset_name_short = xsec_name
                    dataset_xsec = xsec_dict[xsec_name]
        if match_xsec_name < 1:
            raise Exception(f"Failed to find xsec name match for the dataset \"{dataset_name}\"")

        xsec = xsec_dict[dataset_name_short]

        return(dataset_name_short,xsec)


    @staticmethod
    def get_xsec_weight(xsecs, sample):
        if sample in xsecs:
            return xsecs[sample]
        else:
            raise ValueError(f"xsec not found for {sample}")

    @staticmethod
    def get_lumi(year):
        lumi = {
            "2016preVFP": 19.52,
            "2016postVFP": 16.81,
            "2017": 41.48,
            "2018": 59.83,
            "2022Re-recoBCD": 7.9804,
            "2022Re-recoE+PromptFG": 26.6717,
            "2024Prompt": 109.08
        }
        if year in lumi:
            return lumi[year]
        else:
            raise ValueError(f"lumi not found for {year}")

    @staticmethod
    def extract_mc_sample_type(sample_name):
        sample_type_mapping = {
            "DY": "DY",
            "TTto": "TTbar",
            "Wto": "WJets",
            "VBS": "VBS", 
            "QCD": "QCD",
            "WW": "Boson",
            "WZ": "Boson"
        }
        for key, value in sample_type_mapping.items():
            if key in sample_name:
                return value
        return "Other"

    def get_sample_name(self, sample):
        if self.sample_category == "data":
            match = re.search(r"Run(\d{4}[A-Z])/(.+)", sample)
            return f"Run{match.group(1)}-{match.group(2)}"
        else:
            return re.search(r"/([^/]+)_TuneCP5", sample).group(1)

    def process_samples(self, xsecs, nthreads):
        for sample in self.samples:
            # Get the dataset name
            #dataset_name = os.path.basename(os.path.dirname(sample))
            dataset_name = os.path.basename(sample)
            # Get the info about the sample
            sample_name = self.get_sample_name(sample)
            try:
                process_name_sync_with_xsec_name, xsec = self.get_sample_name_and_xsec(dataset_name,xsecs,is_data=self.sample_category=="data")
            except Exception as e:
                print(f"    -> Skipping {dataset_name} as {e}.")
                continue
            sample_year = self.extract_sample_year(sample)
            num_events = 0
            files_path = f"{sample}/*.root"
            if self.sample_category != "data":
                files = glob(files_path)
                with concurrent.futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
                    results = list(executor.map(self._process_file, files))
                num_events = sum(results)
            else:
                num_events = 1.0
            self.config["samples"].update(
                {
                    f"{sample_name}_{sample_year}": {
                        "trees": ["Events"],
                        "files": [files_path],
                        "metadata": {
                            "namewithyear": f"{sample_year}_{process_name_sync_with_xsec_name}",
                            "category": self.sample_category,
                            "year": sample_year,
                            "type": self.extract_mc_sample_type(sample_name) if self.sample_category != "data" else "Muon" if "Muon" in sample_name else "Electron",
                            "xsec": float(xsec),
                            "lumi": self.get_lumi(sample_year) if self.sample_category != "data" else 1.0,
                            "nevents": num_events
                        }
                    }
                }
            )
    
    @staticmethod
    def _process_file(file):
        try:
            with uproot.open(file) as upf:
                # check if file is corrupted
                assert upf["Events"]
                return sum(upf["Runs"]["genEventSumw"].array())
        except Exception as e:
            with open("corrupt_files", "a") as f:
                f.write(file + "\n")
                return 0
            
if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--channel", type=str, help="channel", required=True)
    argparser.add_argument("--category", help="categories: bkg, sig or data", required=True)
    argparser.add_argument("--config", help="paths file", default="paths.json")
    argparser.add_argument("--xsecs", help="xsec file", default="xsecs.json")
    argparser.add_argument("-n", "--nthreads", type=int, default=8, help="number of threads to use for processing files")
    args = argparser.parse_args()

    if not args.category or not any(cat in args.category for cat in ["bkg", "sig", "data"]):
        raise ValueError("Please provide a valid category")

    with open(args.xsecs, "r") as f_xsecs:
        xsecs = json.load(f_xsecs)

    with open(args.config, "r") as f:
        skim_paths = json.load(f)

    if args.channel not in skim_paths:
        raise ValueError(f"No skim paths found for channel {args.channel}")
    if args.category not in skim_paths[args.channel]:
        raise ValueError(f"No skim paths found for category {args.category}")
    
    config = Config(args.category, args.channel, skim_paths[args.channel][args.category], xsecs, args.nthreads)
