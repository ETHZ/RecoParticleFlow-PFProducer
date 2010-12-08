import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")


# load your files
process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring(
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_10_0_pre5/RelValSingleGammaPt35/GEN-SIM-RECO/MC_39Y_V5-v1/0085/F0B2257F-0EF5-DF11-A259-0018F3D096A2.root',
    'rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_10_0_pre5/RelValSingleGammaPt35/GEN-SIM-RECO/MC_39Y_V5-v1/0085/96A8B56F-F7F4-DF11-BBC8-001A92811718.root'
    )
                             )                         


# Event output
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )


# get the proper global tag automatically
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['mc']



# Local re-reco
process.localReReco = cms.Sequence(process.siPixelRecHits+
                                   process.siStripMatchedRecHits+
                                   process.particleFlowCluster)

# Track re-reco
process.globalReReco =  cms.Sequence(process.offlineBeamSpot+
                                     process.recopixelvertexing+
                                     process.ckftracks+
                                     process.caloTowersRec+
                                     process.vertexreco+
                                     process.recoJets+
                                     process.muonrecoComplete+
                                     process.electronGsfTracking+
#                                     process.trackerOnlyConversionSequence+
                                     process.metreco)

# Particle Flow re-processing
process.pfReReco = cms.Sequence(process.particleFlowReco)

process.p = cms.Path(process.localReReco+
                     process.globalReReco+
                     process.pfReReco)

## IMPORTANT SWITCH ON THE NANCY CONVERSIONS IN THE BLOCK!
process.particleFlowBlock.useConversions = True



# define your analyzer
process.load("RecoParticleFlow.PFProducer.pfBlockAnalyzer_cff")

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("histo.root")
    )
process.p5 = cms.Path(process.pfBlockAnalyzer)
    

process.schedule = cms.Schedule(process.p,process.p5)






