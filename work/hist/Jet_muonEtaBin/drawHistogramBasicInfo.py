
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
DIR_NAMES = []   # set directory of the root file ["data","mc" ], if no directory : []
TREE_NAMES = ["Events"]               # set tree path of the root file, you must use same Treename for data and mc.
BRANCHES = [                            # see details at Histo1D_def in the module root_tool.py
    {"output": "Jet_nMuons_09", "name": "Jet_nMuons", "bins": 10, "xmin": 0.0, "xmax": 10.0, "xtitle": "Jet_nMuons (fabs(Jet_eta[i]) < 0.9)", "ytitle": "Multiplicity", "condition": lambda i: "fabs(Jet_eta[i]) < 0.9"},
    {"output": "Jet_muEF_09", "name": "Jet_muEF", "bins": 100, "xmin": 0.0, "xmax": 1.0, "xtitle": "Jet_muEF (Jet_nMuons[i] > 0 && fabs(Jet_eta[i]) < 0.9)", "ytitle": "Multiplicity", "condition":  lambda i: "Jet_nMuons[i] > 0 && fabs(Jet_eta[i]) < 0.9"},
    {"output": "Jet_nMuons_12", "name": "Jet_nMuons", "bins": 10, "xmin": 0.0, "xmax": 10.0, "xtitle": "Jet_nMuons (0.9 < fabs(Jet_eta[i]) < 1.2)", "ytitle": "Multiplicity", "condition": lambda i: "fabs(Jet_eta[i]) >= 0.9 && fabs(Jet_eta[i]) < 1.2"},
    {"output": "Jet_muEF_12", "name": "Jet_muEF", "bins": 100, "xmin": 0.0, "xmax": 1.0, "xtitle": "Jet_muEF (Jet_nMuons[i] > 0  && 0.9 < fabs(Jet_eta[i]) < 1.2)", "ytitle": "Multiplicity", "condition":  lambda i: "Jet_nMuons[i] > 0 && fabs(Jet_eta[i]) >= 0.9 && fabs(Jet_eta[i]) < 1.2"},
    {"output": "Jet_nMuons_21", "name": "Jet_nMuons", "bins": 10, "xmin": 0.0, "xmax": 10.0, "xtitle": "Jet_nMuons (1.2 < fabs(Jet_eta[i]) < 2.1)", "ytitle": "Multiplicity", "condition": lambda i: "fabs(Jet_eta[i]) >= 1.2 && fabs(Jet_eta[i]) < 2.1"},
    {"output": "Jet_muEF_21", "name": "Jet_muEF", "bins": 100, "xmin": 0.0, "xmax": 1.0, "xtitle": "Jet_muEF (Jet_nMuons[i] > 0  && 1.2 < fabs(Jet_eta[i]) < 2.1)", "ytitle": "Multiplicity", "condition":  lambda i: "Jet_nMuons[i] > 0 && fabs(Jet_eta[i]) >= 1.2 && fabs(Jet_eta[i]) < 2.1"},
    {"output": "Jet_nMuons_24", "name": "Jet_nMuons", "bins": 10, "xmin": 0.0, "xmax": 10.0, "xtitle": "Jet_nMuons (2.1 < fabs(Jet_eta[i]) < 2.4)", "ytitle": "Multiplicity", "condition": lambda i: "fabs(Jet_eta[i]) >= 2.1 && fabs(Jet_eta[i]) < 2.4"},
    {"output": "Jet_muEF_24", "name": "Jet_muEF", "bins": 100, "xmin": 0.0, "xmax": 1.0, "xtitle": "Jet_muEF (Jet_nMuons[i] > 0  && 2.1 < fabs(Jet_eta[i]) < 2.4)", "ytitle": "Multiplicity", "condition":  lambda i: "Jet_nMuons[i] > 0 && fabs(Jet_eta[i]) >= 2.1 && fabs(Jet_eta[i]) < 2.4"},
    {"output": "Jet_nMuons_inf", "name": "Jet_nMuons", "bins": 10, "xmin": 0.0, "xmax": 10.0, "xtitle": "Jet_nMuons (2.4 < fabs(Jet_eta[i]))", "ytitle": "Multiplicity", "condition": lambda i: "fabs(Jet_eta[i]) >= 2.4"},
    {"output": "Jet_muEF_inf", "name": "Jet_muEF", "bins": 100, "xmin": 0.0, "xmax": 1.0, "xtitle": "Jet_muEF (Jet_nMuons[i] > 0  && 2.4 < fabs(Jet_eta[i]))", "ytitle": "Multiplicity", "condition":  lambda i: "Jet_nMuons[i] > 0 && fabs(Jet_eta[i]) >= 2.4"},
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
#        latex.DrawLatex(0.70, 0.91, "#sqrt{s} = 13.6 TeV, L = 7.6 /fb")
        latex.DrawLatex(0.70, 0.91, "#sqrt{s} = 13.6 TeV")

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

    root_file, root_dir = load_file(filename, DIR_NAMES)

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