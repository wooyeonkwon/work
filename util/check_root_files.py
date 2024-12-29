import ROOT
import sys

ROOT.EnableImplicitMT()

ROOT.gErrorIgnoreLevel = ROOT.kError

def check_root_files(file_list):
    """
    Function to check ROOT files for corrupted branches or leaves that may cause
    `R__unzip: error -5 in inflate (zlib)` errors.

    Specifically checks if eventNumber 1034034509 exists, and if so, verifies the branches of that entry.

    Parameters:
        file_list (list): List of paths to ROOT files.
    """
    TARGET_EVENT_NUMBER = 90742544

    for file_path in file_list:
        if not file_path.endswith(".root"):
            print(f"Skipping non-ROOT file: {file_path}")
            continue

        root_file = ROOT.TFile.Open(file_path)

        if not root_file or root_file.IsZombie():
            print(f"Error: Unable to open file {file_path}")
            continue

        print(f"Checking file: {file_path}")
        corrupted = False

        for key in root_file.GetListOfKeys():
            obj = key.ReadObj()

            # Check if the object is a directory
            if obj.IsA().InheritsFrom(ROOT.TDirectory.Class()):
                dir_name = obj.GetName()
                directory = root_file.Get(dir_name)

                for sub_key in directory.GetListOfKeys():
                    tree = sub_key.ReadObj()
                    if not tree.InheritsFrom(ROOT.TTree.Class()):
                        continue

                    tree_name = tree.GetName()
                    print(f"  Checking Tree: {tree_name}")

                    try:
                        # Use Filter to find TARGET_EVENT_NUMBER
                        filtered_tree = tree.CopyTree(f"eventNumber == {TARGET_EVENT_NUMBER}")
                        if filtered_tree.GetEntries() > 0:
                            print(f"    Found eventNumber {TARGET_EVENT_NUMBER} in Tree: {tree_name}")
                            
                            for branch in tree.GetListOfBranches():
                                branch_name = branch.GetName()
                                try:
                                    for entry in filtered_tree:
                                        value = getattr(entry, branch_name, None)
                                except Exception as e:
                                    corrupted = True
                                    print(f"    Corrupted Branch: {branch_name} in Tree: {tree_name}")
                                    print(f"      Error: {e}")

                    except Exception as e:
                        print(f"Error while scanning tree {tree_name}: {e}")

        if not corrupted:
            print(f"File {file_path} has no corrupted branches for eventNumber {TARGET_EVENT_NUMBER}.")
        else:
            print(f"File {file_path} contains corrupted branches for eventNumber {TARGET_EVENT_NUMBER}.")

        root_file.Close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_root_files.py <file1.root> <file2.root> ...")
        sys.exit(1)

    file_list = sys.argv[1:]
    check_root_files(file_list)

"""
TTree *tree = (TTree*)_file0->GetDirectory("Analysis")->Get("GEMMuons")
Long64_t entry = 2900294;
TBranch *branch = tree->GetBranch("eventNumber");
int eventNumber;
branch->SetAddress(&eventNumber);
branch->GetEntry(entry);
eventNumber
"""
