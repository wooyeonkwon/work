import ROOT
from ROOT import RDataFrame

data_filename = "/data1/users/dndus0107/AnalysisResults/processed_data/Analysis_Data_22EFG_2.root"
data_file = ROOT.TFile.Open(data_filename)
data_dir = data_file.Get("Analysis")
data_tree = data_dir.Get("Analysis")
data_rdf = ROOT.RDataFrame(data_tree)
hist = data_rdf.Histo1D((
        f"title",  #
        "xname",               #
        60, 60, 120  # bin, xmin, xmax
    ), "zMass_TightReco")

canvas = ROOT.TCanvas("c", "c", 400, 400)
hist.Draw("hist")
canvas.Print("./testHist_0107.pdf")
canvas.Clear()
