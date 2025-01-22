
#############################################
# Description:                              #
#     draw simple histogram with RdataFrame #
#############################################
import ROOT
from ROOT import RDataFrame
import sys
import os
import gc
import threading
from root_tool import load_file
from root_tool import Histo1D_def

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True) # don't show pop-up (can cause seg fault with multi-thread)
ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(48) # activate multi-thread and set a number of threads

data_entry_numbers = []
log_lock = threading.Lock()
canvas_lock = threading.Lock()

##############################
# Set user define parameters #
##############################

# directory, tree, branch name and histogram parameters setting
DIR_NAMES = ["AnalysisMC"]   # set directory of the root file ["data","mc" ], if no directory : []
TREE_NAMES = ["Analysis"]               # set tree path of the root file, you must use same Treename for data and mc.
BRANCHES = [                            # see details at Histo1D_def in the module root_tool.py
#    {"output": "muonpt", "name": "muon_pt", "bins": 200, "xmin": 0.0, "xmax": 200.0, "xtitle": "muonPt [GeV]", "ytitle": "Multiplicity", "condition": lambda i: "muon_pt[i] && fabs(muon_eta[i]) < 2.4"},
#    {"output": "muoneta", "name": "muon_eta", "bins": 60, "xmin": -3.0, "xmax": 3.0, "xtitle": "muonEta", "ytitle": "Multiplicity", "condition": lambda i: "muon_pt[i] && fabs(muon_eta[i]) < 2.4"},
#    {"output": "muonphi", "name": "muon_phi", "bins": 70, "xmin": -3.5, "xmax": 3.5, "xtitle": "muonPhi", "ytitle": "Multiplicity", "condition": lambda i: "muon_pt[i] && fabs(muon_eta[i]) < 2.4"},
#    {"output": "muonpt_nofilter", "name": "muon_pt", "bins": 200, "xmin": 0.0, "xmax": 200.0, "xtitle": "muonPt [GeV]", "ytitle": "Multiplicity"},
#    {"output": "muoneta_nofilter", "name": "muon_eta", "bins": 60, "xmin": -3.0, "xmax": 3.0, "xtitle": "muonEta", "ytitle": "Multiplicity"},
#    {"output": "muonphi_nofilter", "name": "muon_phi", "bins": 70, "xmin": -3.5, "xmax": 3.5, "xtitle": "muonPhi", "ytitle": "Multiplicity"},
    {"output": "runNumber", "name": "runNumber", "bins": 5, "xmin": 0, "xmax": 5, "xtitle": "runNumber", "ytitle": "Multiplicity"},
    {"output": "muonphi_muonMCMatchPass", "name": "muonMCMatchPass", "bins": 5, "xmin": 0, "xmax": 5, "xtitle": "muonMCMatchPass", "ytitle": "Multiplicity"}
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
#        data_rdf (RDataFrame): RDataFrame of data                                   #
#        mc_rdf (RDataFrame): RDataFrame of MC                                       #
#        tree_name (str): current tree                                               #
#        branch (dict): branch info with (bins, xmin, xmax, xtitle, ytitle)          #
#        output_dir (str): histogram output directory                                #
# ####################################################################################
def log_corrupted_entry(event_number): 
    with log_lock: 
        with open("corrupted_entries.log", "a") as log_file: 
            log_file.write(f"Corrupted event: {event_number}\n")


def draw_histogram(root_rdf, tree_name, branch, output_dir="."):

    hist_root = Histo1D_def(root_rdf, branch)

    # draw and decorate canvas
    with canvas_lock: # deactivate multi-thread
        canvas = ROOT.TCanvas("canvas", "", 800, 600)
        hist_root.SetMarkerStyle(20)
        hist_root.SetMarkerColor(ROOT.kBlack)
        hist_root.SetLineColor(ROOT.kBlack)
        
        hist_root.Draw()
        
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.035)
        latex.SetTextFont(62)  # bold font for "CMS"
        latex.DrawLatex(0.15, 0.91, "CMS")
        latex.SetTextFont(42)  # normal font
        latex.DrawLatex(0.21, 0.91, "#it{In Progress}")
        latex.DrawLatex(0.70, 0.91, "#sqrt{s} = 13.6 TeV, L = 7.6 /fb")

        entries = hist_root.Integral()
        print (entries)
        
        
        canvas.SaveAs(f"{output_dir}/{branch['output']}.png")
        del canvas


##########################################################################
# main function : Open Data and MC file and call draw_histogram function #
#     Parameters:                                                        #
#         data_filename (str): data.root file directory                  #
#         mc_filename (str): mc.root file directory                      #
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
                draw_histogram(root_rdf, tree_name, branch, output_dir)
                        
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