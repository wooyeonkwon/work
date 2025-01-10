[folder discription]

skim : 
    using EDFilter, drop useless branches and use HLT trigger. usually running with CRAB in T3
Analysis :
    using EDAnalyzer, MuonSelection and Z reconstruction for "data".
mc :
    using EDAnalyzer, MuonSelection and Z reconstruction for "mc".
        
        output file path : /data1/users/dndus0107/skimmed_data/ (see details at README file in the path)

hist :
    various histogram drawing codes for

MuonCorrection : 
    muon pt correction for data (made by Analysis) and mc (made by mc)

util :
    utility files for handle root files

scratch, Analysis_tnp : 
    ignore these