import ROOT
from ROOT import RDataFrame
import sys
import gc
import threading

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True) #don't show pop-up
ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(48)


# directory, tree, branch name and histogram parameters setting
DIR_NAME = "Analysis"
TREE_NAMES = ["RecoMuons", "RPCMuons", "TightRPCMuons", "GEMMuons", "TightGEMMuons"] #["RecoMuons", "RPCMuons", "TightRPCMuons", "GEMMuons", "TightGEMMuons"]
BRANCHES = [
    {"name": "zBosonMass", "bins": 60, "xmin": 60.0, "xmax": 120.0, "xtitle": "diMuon Mass [GeV]", "ytitle": "Events"},
    {"name": "vertexdz", "bins": 100, "xmin": 0.0, "xmax": 0.1, "xtitle": "vertexdz [cm]", "ytitle": "Events"},
    {"name": "muonSize", "bins": 5, "xmin": 0, "xmax": 5, "xtitle": "muonSize", "ytitle": "Multiplicity"},
    {"name": "muonPt", "bins": 2000, "xmin": 0.0, "xmax": 200.0, "xtitle": "muonPt [GeV]", "ytitle": "Multiplicity"},
    {"name": "muonEta", "bins": 60, "xmin": -3.0, "xmax" : 3.0, "xtitle": "muonEta", "ytitle": "Multiplicity"},
    {"name": "muonPhi", "bins": 70, "xmin": -3.5, "xmax": 3.5, "xtitle": "muonPhi", "ytitle": "Multiplicity"},
    {"name": "muonIso", "bins": 150, "xmin": 0.0, "xmax": 0.15, "xtitle": "muonIso", "ytitle": "Multiplicity"},
]
data_entry_numbers = []
log_lock = threading.Lock()
canvas_lock = threading.Lock()


if len(sys.argv) > 2:
    data_path = sys.argv[1]
    mc_path = sys.argv[2]
    try:
        output_path = sys.argv[3]
    except IndexError:
        output_path = "."
else:
    print("Usage: script.py <data_path> <mc_path> [output_path]")
    sys.exit(1)

def log_corrupted_entry(event_number): 
    """Logs the corrupted event number to a file."""
    with log_lock: 
        with open("corrupted_entries.log", "a") as log_file: 
            log_file.write(f"Corrupted event: {event_number}\n")

def draw_histogram(data_rdf, mc_rdf, tree_name, branch, output_dir="."):
    """
    draw histogram with data and mc
    
    Parameters:
        data_rdf (RDataFrame): RDataFrame of data
        mc_rdf (RDataFrame): RDataFrame of MC
        tree_name (str): current tree
        branch (dict): branch info with (bins, xmin, xmax, xtitle, ytitle)
        output_dir (str): histogram output directory
    """

    # Scale factor    
    mc_cross_section = 6688.0  
    mc_cross_section_unc = 83.99  
    data_lumi = 34652.1
    n_mc = 72909628
    mc_events_total = n_mc if n_mc > 0 else 1.0
    scale_factor = mc_cross_section * data_lumi / mc_events_total
    
    print("flag1")
    
    hist_data = data_rdf.Histo1D((
        f"{tree_name}_data_{branch['name']}", 
        f"{branch['xtitle']}", 
        branch['bins'], branch['xmin'], branch['xmax']
    ), branch['name'])

    hist_mc = mc_rdf.Histo1D((
        f"{tree_name}_mc_{branch['name']}", 
        f"{branch['xtitle']}", 
        branch['bins'], branch['xmin'], branch['xmax']
    ), branch['name'])
    print("flag2")
    hist_mc.Scale(scale_factor)
    print("flag3")
    # Bin uncertainty calculate
    for bin_idx in range(1, hist_mc.GetNbinsX() + 1):
        bin_content = hist_mc.GetBinContent(bin_idx)
        bin_stat_unc = hist_mc.GetBinError(bin_idx)
        bin_scale_unc = bin_content * (mc_cross_section_unc / mc_cross_section)
        total_unc = ROOT.TMath.Sqrt(bin_stat_unc**2 + bin_scale_unc**2)
        hist_mc.SetBinError(bin_idx, total_unc)
    print("flag4")
    with canvas_lock:
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
        latex.SetTextFont(42)  # normal font for "Preliminary"
        latex.DrawLatex(0.21, 0.91, "#it{In Progress}")
        latex.DrawLatex(0.70, 0.91, "#sqrt{s} = 13.6 TeV, L = 11.06/fb")
        
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
        
        canvas.SaveAs(f"{output_dir}/{tree_name}_{branch['name']}.png")
        del canvas


def main(data_filename, mc_filename, output_dir="."):
    """
    main function : Open Data and MC file and call draw_histogram function

    Parameters:
        data_filename (str): data.root file directory
        mc_filename (str): mc.root file directory
        output_dir (str): histogram output directory
    """
    data_file = None
    mc_file = None
    gc.collect()
    
    data_file = ROOT.TFile.Open(data_filename)
    mc_file = ROOT.TFile.Open(mc_filename)
    print("files have been opened")
    
    if not data_file or data_file.IsZombie():
        print(f"Error: Unable to open data file {data_filename}")
        return

    if not mc_file or mc_file.IsZombie():
        print(f"Error: Unable to open MC file {mc_filename}")
        return

    data_dir = data_file.Get(DIR_NAME)
    mc_dir = mc_file.Get(DIR_NAME)

    if not data_dir:
        print(f"Error: Directory '{DIR_NAME}' not found in data file!")
        return

    if not mc_dir:
        print(f"Error: Directory '{DIR_NAME}' not found in MC file!")
        return

    for tree_name in TREE_NAMES:
        data_tree = data_dir.Get(tree_name)
        mc_tree = mc_dir.Get(tree_name)

        if not data_tree:
            print(f"Tree '{tree_name}' not found in data directory '{DIR_NAME}'!")
            continue

        if not mc_tree:
            print(f"Tree '{tree_name}' not found in MC directory '{DIR_NAME}'!")
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
