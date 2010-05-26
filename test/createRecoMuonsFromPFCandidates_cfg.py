import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)


## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/relval/CMSSW_3_6_0_pre6/RelValTTbar/GEN-SIM-RECO/START36_V4-v1/0010/323F69C9-9F44-DF11-8BA0-002618943943.root'
    )
)


# load the module which reads the collection of PFCandidates, and
# creates a reco::Muon for each PFCandidate of type muon
process.load("RecoParticleFlow.PFProducer.recoMuonFromPFProducer_cfi")

# put it in the path
process.p = cms.Path(process.recoMuonFromPFProducer)

process.load("FastSimulation.Configuration.EventContent_cff")
process.out = cms.OutputModule("PoolOutputModule",
# don't forget to keep all collections of recoMuons in the event.
# no need to keep the recoPFCandidates for a muon analysis
    outputCommands = cms.untracked.vstring('keep recoMuons_*_*_*',
                                           'keep recoPFCandidates_particleFlow_*_*'),
    fileName = cms.untracked.string('muonsFromPF.root')
)

process.outpath = cms.EndPath(process.out )




