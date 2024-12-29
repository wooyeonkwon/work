import FWCore.ParameterSet.Config as cms
import os
import sys

process = cms.Process("Analysis")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))


if len(sys.argv) > 1:
    sublist_file = sys.argv[1] 
    if os.path.exists(sublist_file):
        with open(sublist_file, 'r') as f:
            file_paths = [f'file://{line.strip()}' for line in f.readlines()]
    else:
        raise FileNotFoundError(f"Input list file '{sublist_file}' not found.")

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(*file_paths)
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)


if len(sys.argv) > 2:
    output_file_name = sys.argv[2]
else:
    output_file_name = "data_2022C.root"

process.TFileService = cms.Service("TFileService",
    fileName=cms.string(output_file_name)
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons=cms.InputTag("muons"),
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
    vertices=cms.InputTag("offlinePrimaryVertices"),
    hltPath=cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.Analysis)

process.options.wantSummary = True
