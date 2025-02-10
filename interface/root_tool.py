import ROOT
from ROOT import RDataFrame
import sys
import numpy as np

# Load and check directory in root files and return direcrtory
def load_file(root_filename, DIR_NAMES):
        root_file = ROOT.TFile.Open(root_filename)
        print(f"{root_filename} file has been opened")
        
        if not root_file or root_file.IsZombie():
            print(f"Error: Unable to open data file {root_filename}")
            return sys.exit(1)
        
        if DIR_NAMES != []:
            root_dir = root_file.Get(DIR_NAMES)            
            ROOT.SetOwnership(root_file, False)
            ROOT.SetOwnership(root_dir, False)

            if not root_dir:
                print(f"Error: Directory '{DIR_NAMES}' not found in data file!")
                # List all keys in the data file
                for key in root_file.GetListOfKeys():
                    print(key.GetName())
                return sys.exit(1)
            
            return root_file, root_dir
        else :
            return root_file, root_dir


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
    condition = branch.get("condition", None)  # if there is no condition, None
    new_var_name = variable_name

    # check if scalar or vector
    # nevermind error messages like the one below
    # error: member reference base type 'const Double_t' (aka 'const double') is not a structure or union for (size_t i = 0; i < var1.size(); ++i) {

    is_vector = True
    try:
        scalar_check_code = f"return {variable_name}.size();"
        rdf.Define("check", scalar_check_code)
    except:
        is_vector = False

    if condition:
        if is_vector:
            condition_str = condition("i") if callable(condition) else condition
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
        else:
            condition_str = condition("") if callable(condition) else condition
            new_var_name = f"filtered_{variable_name}"
            rdf = rdf.Define(new_var_name,
                             f"""
                             return {condition_str} ? {variable_name} : std::numeric_limits<double>::quiet_NaN();
                             """)
    else:
        if is_vector:
            new_var_name = variable_name
        else:
            new_var_name = variable_name

    # create histogram
    hist = rdf.Histo1D((
        f"{variable_name}_hist",
        branch["xtitle"],        
        branch["bins"], branch["xmin"], branch["xmax"]
    ), new_var_name)

    ROOT.SetOwnership(hist, False)

    return hist

# make profile with filter muon
def Profile1D_filterx_def(rdf, branch, weight = False):
    variable_name = branch["name"]
    condition = branch.get("condition", None)  # if there is no condition, None
    new_var_name = variable_name

    # check if scalar or vector
    # nevermind error messages like the one below
    # error: member reference base type 'const Double_t' (aka 'const double') is not a structure or union for (size_t i = 0; i < var1.size(); ++i) {

    is_vector = True
    try:
        scalar_check_code = f"return {variable_name}.size();"
        rdf.Define("check", scalar_check_code)
    except:
        is_vector = False

    if condition:
        if is_vector:
            condition_str = condition("i") if callable(condition) else condition
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
        else:
            condition_str = condition("") if callable(condition) else condition
            new_var_name = f"filtered_{variable_name}"
            rdf = rdf.Define(new_var_name,
                             f"""
                             return {condition_str} ? {variable_name} : std::numeric_limits<double>::quiet_NaN();
                             """)
    else:
        if is_vector:
            new_var_name = variable_name
        else:
            new_var_name = variable_name
    # create profile
    if isinstance(branch["bins"], list):
        bins = np.array(branch["bins"], dtype=np.float64)
        profile = rdf.Profile1D((
            f"{variable_name}_profile",
            branch["xtitle"],
            len(bins) - 1, 
            bins
        ), new_var_name, branch["ymean"])
    else:
        profile = rdf.Profile1D((
            f"{variable_name}_profile",
            branch["xtitle"],        
            branch["bins"], branch["xmin"], branch["xmax"]
        ),  new_var_name, branch["ymean"])

    ROOT.SetOwnership(profile, False)

    return profile

# make profile with filter for change y value (result) (use for efficiency)
def Profile1D_filtery_def(rdf, branch, weight = False):
    variable_ymean = branch["ymean"]
    condition = branch.get("condition", None)  # if there is no condition, None
    new_var_ymean = variable_ymean

    # check if scalar or vector
    # nevermind error messages like the one below
    # error: member reference base type 'const Double_t' (aka 'const double') is not a structure or union for (size_t i = 0; i < var1.size(); ++i) {

    is_vector = True
    try:
        scalar_check_code = f"return {variable_ymean}.size();"
        rdf.Define("check", scalar_check_code)
    except:
        is_vector = False

    if condition:
        if is_vector:
            condition_str = condition("i") if callable(condition) else condition
            new_var_ymean = f"filtered_{variable_ymean}"
            rdf = rdf.Define(new_var_ymean,
                             f"""
                             std::vector<double> filtered;
                             for (size_t i = 0; i < {variable_ymean}.size(); ++i) {{
                                 if ({condition_str}) {{
                                     filtered.push_back({variable_ymean}[i]);
                                 }}
                                 else {{
                                     filtered.push_back(false);
                                 }}
                             }}
                             return filtered;
                             """)
        else:
            condition_str = condition("") if callable(condition) else condition
            new_var_ymean = f"filtered_{variable_ymean}"
            rdf = rdf.Define(new_var_ymean,
                    f"""
                    return {condition_str} ? {variable_ymean} : std::numeric_limits<double>::quiet_NaN();
                    """)

    else:
        if is_vector:
            new_var_ymean = variable_ymean
        else:
            new_var_ymean = variable_ymean


    if weight:
        variable_weight = variable_ymean
        # create profile
        if isinstance(branch["bins"], list):
            bins = np.array(branch["bins"], dtype=np.float64)
            profile = rdf.Profile1D((
                f"{variable_ymean}_profile",
                branch["xtitle"],
                len(bins) - 1, 
                bins
            ), branch["name"], new_var_ymean, variable_weight)
        else:
            profile = rdf.Profile1D((
                f"{variable_ymean}_profile",
                branch["xtitle"],        
                branch["bins"], branch["xmin"], branch["xmax"]
            ),  branch["name"], new_var_ymean, variable_weight)

    else :
        # create profile
        if isinstance(branch["bins"], list):
            bins = np.array(branch["bins"], dtype=np.float64)
            profile = rdf.Profile1D((
                f"{variable_ymean}_profile",
                branch["xtitle"],
                len(bins) - 1, 
                bins
            ), branch["name"], new_var_ymean)
        else:
            profile = rdf.Profile1D((
                f"{variable_ymean}_profile",
                branch["xtitle"],        
                branch["bins"], branch["xmin"], branch["xmax"]
            ),  branch["name"], new_var_ymean)
        


    ROOT.SetOwnership(profile, False)

    return profile
