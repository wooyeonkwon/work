import ROOT
import sys
import os

ROOT.EnableImplicitMT()

ROOT.gErrorIgnoreLevel = ROOT.kError

def check_root_files(file_list):
    """
    Function to check ROOT files for errors, including `R__unzip` issues, and validate tree entries and leaves.
    
    Parameters:
        file_list (list): List of paths to ROOT files.
    """
    for file_path in file_list:
        if not file_path.endswith(".root"):
            print(f"Skipping non-ROOT file: {file_path}")
            continue

        root_file = ROOT.TFile.Open(file_path)

        if not root_file or root_file.IsZombie():
            print(f"Error: Unable to open file {file_path}")
            continue

        print(f"Checking file: {file_path}")

        for key in root_file.GetListOfKeys():
            obj = key.ReadObj()
            if obj.IsA().InheritsFrom(ROOT.TDirectory.Class()):
                dir_name = obj.GetName()
                directory = root_file.Get(dir_name)

                for sub_key in directory.GetListOfKeys():
                    tree = sub_key.ReadObj()
                    if not tree.InheritsFrom(ROOT.TTree.Class()):
                        continue

                    tree_name = tree.GetName()
                    print(f"  Tree: {tree_name}, Entries: {tree.GetEntries()}")

                    for entry in range(tree.GetEntries()):
                        try:
                            tree.GetEntry(entry)
                            for branch in tree.GetListOfBranches():
                                branch_name = branch.GetName()
                                #print(f"  branch: {branch_name}")
                                for leaf in branch.GetListOfLeaves():
                                    leaf_name = leaf.GetName()
                                    #print(f"  leaf: {leaf_name}")
                                    leaf_value = leaf.GetValue()
                        except Exception as e:
                            print(f"    Error in tree '{tree_name}' at entry {entry}: {e}")

        root_file.Close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_root_files.py <file1.root> <file2.root> ...")
        sys.exit(1)

    file_list = sys.argv[1:]
    check_root_files(file_list)
