
###################################################################
# Description:                                                    #
#     draw comparision profile                                    #
#     draw value of mean of y per x bin                           #
#     if you use CONDITIONS dictionary, draw these in same canvas #
###################################################################
import ROOT
from ROOT import RDataFrame
import sys
import os
import gc
import threading
from root_tool import load_file
from root_tool import Profile1D_filterx_def

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True) # don't show pop-up (can cause seg fault with multi-thread)
ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(48) # activate multi-thread and set a number of threads

data_entry_numbers = []
log_lock = threading.Lock()
profile_lock = threading.Lock()
canvas_lock = threading.Lock()

##############################
# Set user define parameters #
##############################

# directory, tree, branch name and histogram parameters setting
DIR_NAMES = ["AnalysisMC"]   # set directory of the root file, if no directory : []
TREE_NAMES = ["Analysis"]    # set tree path of the root file, you can use multiple trees.
CONDITIONS = {
    "muon_isRecoZ": {"name": "muon_", "ymean" : "zMass_Reco", "condition": lambda i: "muon_isRecoZ[i] && zMass_Reco > 86.5 && zMass_Reco <96.5", "color": ROOT.kBlack, "marker": 1},
    "muon_isGenZ": {"name": "genMuon_", "ymean" : "zMass_Gen", "condition": lambda i: "genMuon_isGenZ[i] && zMass_Gen > 86.5 && zMass_Gen < 96.5", "color": ROOT.kRed, "marker": 1},
}

BRANCHES = [ # see details at Profile1D_def in the module interface/root_tool.py
    {"output": "muoneta_matchingEff", "name": "eta", "ymean" : None, "bins": [-3, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, 0, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 3], "xmin": -3.0, "xmax": 3.0, "xtitle": "muonEta", "ytitle": "Multiplicity", "condition": None},
    {"output": "muonphi_matchingEff", "name": "phi", "ymean" : None, "bins": 20, "xmin": -3.2, "xmax": 3.2, "xtitle": "muonPhi", "ytitle": "Multiplicity", "condition": None},
]

################################################################
# Usage of this code                                           #
# python3 script.py <path/to/data.root> [path/to/output/folder]#
################################################################

if len(sys.argv) > 1:
    root_path = sys.argv[1]
    try:
        output_path = sys.argv[2]
        if not os.path.exists(output_path): 
            os.makedirs(output_path)
    except IndexError:
        output_path = "."
else:
    print("Usage: script.py <root_path> [output_path]")
    sys.exit(1)

######################################################################################
# helper functions                                                                   #
# log_corrupted_entry(event_number)                                                  #
#     Logs the corrupted event number to a file                                      #
#     parameters:                                                                    #
#         event_number (int) : Corrupted event, which can't be read a entry          #
#                                                                                    #
# draw_histogram                                                                     #
#    draw histogram with data and mc                                                 #
#    parameters:                                                                     #
#        root_rdf (RDataFrame): RDataFrame of data                                   #
#        tree_name (str): current tree                                               #
#        branch (dict): branch info with (bins, xmin, xmax, xtitle, ytitle)          #
#        output_dir (str): histogram output directory                                #
# ####################################################################################
def log_corrupted_entry(event_number): 
    with log_lock: 
        with open("corrupted_entries.log", "a") as log_file: 
            log_file.write(f"Corrupted event: {event_number}\n")



def draw_profile(root_rdf, tree_name, branch, output_dir="."):
    profiles = {}
    origial_branch = branch["name"]
    # make Profile for all CONDITIONS
    for name, style in CONDITIONS.items():
        condition = style["condition"]
        branch["condition"] = condition # set condition
        ymean = style["ymean"]
        branch["ymean"] = ymean
        branch_name = style["name"]
        branch["name"] = branch_name + origial_branch
        profiles[name] = Profile1D_filterx_def(root_rdf, branch)
        profile = profiles[name].GetPtr()

        # Bin uncertainty calculate
        with profile_lock:
            print(type(profile))
            for bin_idx in range(1, profile.GetNbinsX() + 1):
                bin_stat_unc = profile.GetBinError(bin_idx)
                total_unc = bin_stat_unc
                profile.SetBinError(bin_idx, total_unc)

            # set Profile styles
            profile.SetMarkerStyle(style["marker"])
            profile.SetMarkerColor(style["color"])
            profile.SetLineColor(style["color"])

    with canvas_lock:
        canvas = ROOT.TCanvas("canvas", "", 800, 600)

        # pad1
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0, 1, 1)
        pad1.Draw()
        pad1.cd()

        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

        for idx, (name, style) in enumerate(CONDITIONS.items()):
            profile = profiles[name].GetPtr()
            draw_option = "E same" if idx > 0 else "E"
            profile.SetMinimum(86)
            profile.SetMaximum(97)
            profile.Draw(draw_option)
            if isinstance(profile, ROOT.TObject):
                legend.AddEntry(profile, name, "lep")
            else:
                print("Error: Profile is not a valid ROOT TObject")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.035)
        latex.SetTextFont(62)
        latex.DrawLatex(0.15, 0.91, "CMS")
        latex.SetTextFont(42)
        latex.DrawLatex(0.21, 0.91, "#it{In Progress}")
        latex.DrawLatex(0.70, 0.91, "#sqrt{s} = 13.6 TeV, L = 7.6 /fb")

        line = ROOT.TLine(branch['xmin'], 91.1876, branch['xmax'], 91.1876)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.Draw()

        legend.Draw()

        canvas.SaveAs(f"{output_dir}/{branch['output']}_combined.png")
        del canvas



##########################################################################
# main function : Open Data and MC file and call draw_histogram function #
#     Parameters:                                                        #
#         filename (str): data.root file directory                       #
#         output_dir (str): histogram output directory                   #
##########################################################################

def main(filename,output_dir="."):

    root_file = None
    gc.collect()

    root_file, root_dir = load_file(filename, DIR_NAMES[0])

    for tree_name in TREE_NAMES:
        root_tree = root_dir.Get(tree_name)
        if not root_tree:
            print(f"Tree '{tree_name}' not found in data directory '{DIR_NAMES}'!")
            continue

        root_rdf = ROOT.RDataFrame(root_tree)

        for branch in BRANCHES:
            print(f"Draw '{branch}' in '{tree_name}'")
            try:
                draw_profile(root_rdf, tree_name, branch, output_dir)
                        
            except Exception as e:
                # Check for specific error and log the problematic entry
                print(f"Error initializing RDataFrame for tree '{tree_name}': {str(e)}")
                for entry in root_tree:
                    try:
                        event_number = getattr(entry, "eventNumber", None)
                        if event_number:
                            data_entry_numbers.append(event_number)
                    except:
                        log_corrupted_entry(event_number)
                        continue
                continue

    root_file.Close()

main(root_path, output_path)
