
#####################################################
# Description:                                      #
#     draw data and mc comparision histogram.       #
#     draw both data and mc histogram in same pad 1 #
#     draw data/mc ratio in pad 2                   #
#####################################################
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
DIR_NAMES = ["Analysis","AnalysisMC"] # set directory of the root file ["data","mc" ], if no directory : []
TREE_NAMES = ["Analysis"] # set tree path of the root file, you  
BRANCHES = [
    {"output": "muonpt_isTightRecoZ", "name": "muon_pt", "bins": 200, "xmin": 20.0, "xmax": 120.0, "xtitle": "muonPt [GeV]", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRecoZ[i]"},
    {"output": "muoneta_isTightRecoZ", "name": "muon_eta", "bins": 60, "xmin": -3.0, "xmax": 3.0, "xtitle": "muonEta", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRecoZ[i]"},
    {"output": "muonphi_isTightRecoZ", "name": "muon_phi", "bins": 70, "xmin": -3.5, "xmax": 3.5, "xtitle": "muonPhi", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRecoZ[i]"},
    {"output": "muoniso_isTightRecoZ", "name": "muon_iso", "bins": 150, "xmin": 0.0, "xmax": 0.15, "xtitle": "muonIso", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRecoZ[i]"},
    {"output": "zMass_isTightReco", "name": "zMass_TightReco", "bins": 120, "xmin": 60.0, "xmax": 120.0, "xtitle": "dimuon mass", "ytitle": "Multiplicity"},
    {"output": "muonpt_isTightRPCZ", "name": "muon_pt", "bins": 2000, "xmin": 0.0, "xmax": 200.0, "xtitle": "muonPt [GeV]", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRPCZ[i]"},
    {"output": "muoneta_isTightRPCZ", "name": "muon_eta", "bins": 60, "xmin": -3.0, "xmax": 3.0, "xtitle": "muonEta", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRPCZ[i]"},
    {"output": "muonphi_isTightRPCZ", "name": "muon_phi", "bins": 70, "xmin": -3.5, "xmax": 3.5, "xtitle": "muonPhi", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRPCZ[i]"},
    {"output": "muoniso_isTightRPCZ", "name": "muon_iso", "bins": 150, "xmin": 0.0, "xmax": 0.15, "xtitle": "muonIso", "ytitle": "Multiplicity", "condition": lambda i: "muon_isTightRPCZ[i]"},
    {"output": "zMass_isTightRPC", "name": "zMass_TightRPC", "bins": 120, "xmin": 60.0, "xmax": 120.0, "xtitle": "dimuon mass", "ytitle": "Multiplicity"},
]

# parameters for Normalization factor
# Run 22C,D
mc_cross_section = 2219.0  
mc_cross_section_unc = 0.2327  
data_lumi = 7.9804
gen_weight = 2219
n_mc = 2924957
equi_lumi = 0.4327

# Run 22E,F,G
#   mc_cross_section = 2219.0  
#   mc_cross_section_unc = 0.2327  
#   data_lumi = 5.8070 + 17.7819 + 3.0828
#   gen_weight = 2209
#   n_mc = 10148870
#   equi_lumi = 0.4327

###################################################################################
# Usage of this code                                                              #
# python3 script.py <path/to/data.root> <path/to/mc.root> [path/to/output/folder] #
###################################################################################

if len(sys.argv) > 2:
    data_path = sys.argv[1]
    mc_path = sys.argv[2]
    try:
        output_path = sys.argv[3]
        if not os.path.exists(output_path): 
            os.makedirs(output_path)
    except IndexError:
        output_path = "."
else:
    print("Usage: script.py <data_path> <mc_path> [output_path]")
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



def draw_histogram(data_rdf, mc_rdf, tree_name, branch, output_dir="."):

    hist_data = Histo1D_def(data_rdf, branch)
    hist_mc = Histo1D_def(mc_rdf, branch)

    mc_events_total = n_mc if n_mc > 0 else 1.0
    nom_factor = equi_lumi * gen_weight * data_lumi * mc_cross_section / mc_events_total
    hist_mc.Scale(nom_factor)

    # Bin uncertainty calculate
    for bin_idx in range(1, hist_mc.GetNbinsX() + 1):
        bin_content = hist_mc.GetBinContent(bin_idx)
        bin_stat_unc = hist_mc.GetBinError(bin_idx)
        bin_scale_unc = bin_content * (mc_cross_section_unc / mc_cross_section)
        total_unc = ROOT.TMath.Sqrt(bin_stat_unc**2 + bin_scale_unc**2)
        hist_mc.SetBinError(bin_idx, total_unc)

    # draw and decorate canvas
    with canvas_lock: # deactivate multi-thread
        canvas = ROOT.TCanvas("canvas", "", 800, 800)
        hist_data.SetMarkerStyle(20)
        hist_data.SetMarkerColor(ROOT.kBlack)
        hist_data.SetLineColor(ROOT.kBlack)
        hist_mc.SetFillColor(ROOT.kBlue)
        hist_mc.SetLineColor(ROOT.kBlue)
        
        # pad1
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.Draw()
        pad1.cd()
        
        hist_mc.Draw("HIST")
        hist_data.Draw("E same")
        
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.035)
        latex.SetTextFont(62)  # bold font for "CMS"
        latex.DrawLatex(0.15, 0.91, "CMS")
        latex.SetTextFont(42)  # normal font
        latex.DrawLatex(0.21, 0.91, "#it{In Progress}")
        latex.DrawLatex(0.70, 0.91, "#sqrt{s} = 13.6 TeV, L = 7.98/fb")
        
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        legend.AddEntry(hist_data.GetValue(), "Data", "lep")
        legend.AddEntry(hist_mc.GetValue(), "MC", "f")
        legend.Draw()
        
        # pad2
        canvas.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.2)
        pad2.Draw()
        pad2.cd()
        
        ratioHist = hist_data.Clone("ratioHist")
        ratioHist.Divide(hist_mc.GetValue())
        ratioHist.SetTitle("")
        ratioHist.GetYaxis().SetTitle("Data / MC")
        ratioHist.GetYaxis().SetNdivisions(505)
        ratioHist.GetYaxis().SetTitleSize(20)
        ratioHist.GetYaxis().SetTitleFont(43)
        ratioHist.GetYaxis().SetTitleOffset(1.55)
        ratioHist.GetYaxis().SetLabelFont(43)
        ratioHist.GetYaxis().SetLabelSize(15)
        ratioHist.GetXaxis().SetTitleSize(20)
        ratioHist.GetXaxis().SetTitleFont(43)
        ratioHist.GetXaxis().SetTitleOffset(4.0)
        ratioHist.GetXaxis().SetLabelFont(43)
        ratioHist.GetXaxis().SetLabelSize(15)
        ratioHist.Draw("ep")
        
        line = ROOT.TLine(branch['xmin'], 1, branch['xmax'], 1)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.Draw()
        
        canvas.SaveAs(f"{output_dir}/{branch['output']}.png")
        del canvas


##########################################################################
# main function : Open Data and MC file and call draw_histogram function #
#     Parameters:                                                        #
#         data_filename (str): data.root file directory                  #
#         mc_filename (str): mc.root file directory                      #
#         output_dir (str): histogram output directory                   #
##########################################################################

def main(data_filename, mc_filename, output_dir="."):

    data_file = None
    mc_file = None
    gc.collect()

    data_dir, mc_dir, data_file, mc_file = load_file(data_filename, mc_filename, DIR_NAMES)

    for tree_name in TREE_NAMES:
        data_tree = data_dir.Get(tree_name)
        mc_tree = mc_dir.Get(tree_name)
        if not data_tree:
            print(f"Tree '{tree_name}' not found in data directory '{DIR_NAMES}'!")
            continue

        if not mc_tree:
            print(f"Tree '{tree_name}' not found in MC directory '{DIR_NAMES}'!")
            continue

        data_rdf = ROOT.RDataFrame(data_tree)
        mc_rdf = ROOT.RDataFrame(mc_tree)

        for branch in BRANCHES:
            print(f"Draw '{branch}' in '{tree_name}'")
            try:
                draw_histogram(data_rdf, mc_rdf, tree_name, branch, output_dir)
                        
            except Exception as e:
                # Check for specific error and log the problematic entry
                print(f"Error initializing RDataFrame for tree '{tree_name}': {str(e)}")
                for entry in data_tree:
                    try:
                        event_number = getattr(entry, "eventNumber", None)
                        if event_number:
                            data_entry_numbers.append(event_number)
                    except:
                        log_corrupted_entry(event_number)
                        continue
                continue

    data_file.Close()
    mc_file.Close()

main(data_path, mc_path, output_path)
