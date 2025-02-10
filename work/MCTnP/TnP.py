###################################################################
# Description:                                                    #
#     Tag and Probe with Z to mumu channel                        #
#     Tag Selection : TightID & TightIso & pt >24 & eta < 2.4     #
#     Probe Selection : Tracker & pt > 15 & eta < 2.4             #     
#     event Selection : Z mass window == [60,120] , dvz < 0.5     #
#     passing : Global & RPC ,Global & GEM                        #
###################################################################
import ROOT
from ROOT import RDataFrame
import sys
import os
import gc
import threading
from root_tool import load_file


ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True) # don't show pop-up (can cause seg fault with multi-thread)
ROOT.EnableImplicitMT(16)

data_entry_numbers = []
log_lock = threading.Lock()

##############################
# Set user define parameters #
##############################
N_THREAD = 1
# directory, tree, branch name and histogram parameters setting
DIR_NAMES = ["AnalysisMC"]   # set directory of the root file, if no directory : []
TREE_NAMES = ["Analysis"]    # set tree path of the root file, you can use multiple trees.
IS_MC = True                 # if you are using MC, set True for MCMathcing
BRANCHES = [ # set branch to read
    "eventNumber",
    "runNumber", 
    "lumiSection",
    "muon_charge",
    "muon_pt",
    "muon_eta",
    "muon_phi",
    "muon_iso",
    "muon_vz",
    "muon_isGlobal",
    "muon_isTracker",
    "muon_isRPC",
    "muon_isGEM",
    "muon_isTight"
]
if IS_MC :
    BRANCHES = BRANCHES + ["muonMCMatched", "muonMCMatchPass","genWeight"]

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
# make_root_file(root_rdf, tree_name, filename, output_dir)                          #                                         #
#    make empty root file with brach structure                                       #
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

def z_mass_reconstruction(root_rdf, tree_name, output_dir):
    """ Reconstruct Z mass and filter events for RPC and GEM probes """
    output_file = ROOT.TFile(f"{output_dir}/{tree_name}_TnP.root", "RECREATE")
    output_tree = ROOT.TTree("TnPTree", "Tag and Probe Tree")

    # Define variables for the output tree
    zMassRPC = ROOT.std.vector('double')()
    zDvzRPC = ROOT.std.vector('double')()
    zMassGEM = ROOT.std.vector('double')()
    zDvzGEM = ROOT.std.vector('double')()
    muon_isRPCPassProbe = ROOT.std.vector('bool')()
    muon_isGEMPassProbe = ROOT.std.vector('bool')()
    print("flag1")
    # Associate variables with branches
    output_tree.Branch("zMassRPC", zMassRPC)
    output_tree.Branch("zDvzRPC", zDvzRPC)
    output_tree.Branch("zMassGEM", zMassGEM)
    output_tree.Branch("zDvzGEM", zDvzGEM)
    output_tree.Branch("muon_isRPCPassProbe", muon_isRPCPassProbe)
    output_tree.Branch("muon_isGEMPassProbe", muon_isGEMPassProbe)
    print("flag2")
    # Iterate over each event
    for event in root_rdf.AsNumpy().all():
        print("flag3")
        tag_muon_indices = [
            i for i, (pt, eta, isTight, isGlobal) in enumerate(zip(
                event["muon_pt"], event["muon_eta"], event["muon_isTight"], event["muon_isGlobal"]
            )) if pt > 24 and abs(eta) < 2.4 and isTight and isGlobal
        ]

        for tag_idx in tag_muon_indices:
            print(f"{tag_idx}")
            for probe_idx in range(len(event["muon_pt"])):
                if tag_idx == probe_idx:  # Same muon
                    continue
                if event["muon_charge"][tag_idx] == event["muon_charge"][probe_idx]:  # Same charge
                    continue

                # Calculate Z mass and dvz
                z_mass = calculate_z_mass(event["muon_pt"][tag_idx], event["muon_eta"][tag_idx], event["muon_phi"][tag_idx],
                                          event["muon_pt"][probe_idx], event["muon_eta"][probe_idx], event["muon_phi"][probe_idx])
                dvz = abs(event["muon_vz"][tag_idx] - event["muon_vz"][probe_idx])

                if 60 <= z_mass <= 120 and dvz < 0.5:
                    zMassRPC.clear()
                    zDvzRPC.clear()
                    zMassGEM.clear()
                    zDvzGEM.clear()
                    muon_isRPCPassProbe.clear()
                    muon_isGEMPassProbe.clear()

                    zMassRPC.push_back(z_mass)
                    zDvzRPC.push_back(dvz)
                    zMassGEM.push_back(z_mass)
                    zDvzGEM.push_back(dvz)
                    muon_isRPCPassProbe.push_back(event["muon_isRPC"][probe_idx])
                    muon_isGEMPassProbe.push_back(event["muon_isGEM"][probe_idx])

                    output_tree.Fill()
    
    output_file.Write()
    output_file.Close()


def calculate_z_mass(pt1, eta1, phi1, pt2, eta2, phi2):
    """ Calculate invariant mass of two muons """
    vec1 = ROOT.TLorentzVector()
    vec2 = ROOT.TLorentzVector()
    vec1.SetPtEtaPhiM(pt1, eta1, phi1, 0.105)  # Muon mass = 0.105 GeV
    vec2.SetPtEtaPhiM(pt2, eta2, phi2, 0.105)
    return (vec1 + vec2).M()




##########################################################################
# main function : Open Data and MC file and call draw_histogram function #
#     Parameters:                                                        #
#         filename (str): data.root file directory                       #
#         output_dir (str): histogram output directory                   #
##########################################################################

def main(filename, output_dir="."):

    root_file = None
    gc.collect()

    root_file = ROOT.TFile.Open(filename)
    root_dir = root_file.Get(DIR_NAMES[0])
    for tree_name in TREE_NAMES:
        root_tree = root_dir.Get(tree_name)
        if not root_tree:
            print(f"Tree '{tree_name}' not found in data directory '{DIR_NAMES}'!")
            continue
        root_rdf = ROOT.RDataFrame(root_tree, BRANCHES)
        print("flag0")
        if IS_MC :
            root_rdf = root_rdf.Filter("muonMCMatchPass == 1")
        try:
            z_mass_reconstruction(root_rdf, tree_name, output_dir)
                    
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
