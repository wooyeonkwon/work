import ROOT
from ROOT import RDataFrame
import sys

# Load and check directory in root files and return direcrtory
def load_file(data_filename, mc_filename, DIR_NAMES):
        data_file = ROOT.TFile.Open(data_filename)
        mc_file = ROOT.TFile.Open(mc_filename)
        print("files have been opened")
        
        if not data_file or data_file.IsZombie():
            print(f"Error: Unable to open data file {data_filename}")
            return sys.exit(1)

        if not mc_file or mc_file.IsZombie():
            print(f"Error: Unable to open MC file {mc_filename}")
            return sys.exit(1)
        
        if DIR_NAMES != []:
            data_dir = data_file.Get(DIR_NAMES[0])
            mc_dir = mc_file.Get(DIR_NAMES[1])
            
            ROOT.SetOwnership(data_file, False)
            ROOT.SetOwnership(mc_file, False)
            ROOT.SetOwnership(data_dir, False)
            ROOT.SetOwnership(mc_dir, False)

            if not data_dir:
                print(f"Error: Directory '{DIR_NAMES}' not found in data file!")
                # List all keys in the data file
                for key in data_file.GetListOfKeys():
                    print(key.GetName())
                return sys.exit(1)

            if not mc_dir:
                print(f"Error: Directory '{DIR_NAMES}' not found in MC file!")
                # List all keys in the data file
                for key in mc_file.GetListOfKeys():
                    print(key.GetName())
                return sys.exit(1)
            return data_dir, mc_dir, data_file, mc_file
        else :
            return data_file, mc_file, data_file, mc_file


#define Rdataframe with conditions and make it a histogram
#use with branch define format like example below
#"condition": lambda i: <"condition">
#condition must returns boolean type and must be wrapped with double quote
#BRANCHES = [
#    {"name": "zMass_TightReco", "bins": 120, "xmin": 60.0, "xmax": 120.0, "xtitle": "dimuon mass", "ytitle": "Multiplicity"}, #if no condition delete condition element
#    {"name": "zMass_TightRPC", "bins": 120, "xmin": 60.0, "xmax": 120.0, "xtitle": "dimuon mass", "ytitle": "Multiplicity", "condition": lambda i: "(muon_isTightReco[i] == 1) && (muon_pt[i] > 40)"}, 
#    ...
#]
def Histo1D_def(rdf, branch):
    variable_name = branch["name"]
    condition = branch.get("condition", None)
    new_var_name = variable_name

    if condition:
        condition_str = condition("i")
        new_var_name = f"filtered_{variable_name}"
        rdf = rdf.Define(new_var_name,
                         f"""
                         std::vector<double> filtered;
                         for (size_t i = 0; i < {variable_name}.size(); ++i) {{
                             if ({condition_str}) {{
                                 filtered.push_back({variable_name}[i]);
                             }}
                         }}
                         return filtered;
                         """)

    hist = rdf.Histo1D((
        f"{variable_name}_hist",      
        branch["xtitle"],             
        branch["bins"],               
        branch["xmin"], branch["xmax"]  
    ), new_var_name)

    ROOT.SetOwnership(hist, False)

    return hist
