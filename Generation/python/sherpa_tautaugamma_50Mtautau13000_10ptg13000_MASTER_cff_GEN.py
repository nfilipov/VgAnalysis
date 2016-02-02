# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/sherpa_tautaugamma_50Mtautau13000_10ptg13000_MASTER_cff -s GEN --conditions auto:mc --eventcontent RAWSIM -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/Generator/python/sherpa_tautaugamma_50Mtautau13000_10ptg13000_MASTER_cff nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('sherpa_tautaugamma_50Mtautau13000_10ptg13000_MASTER_cff_GEN.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("SherpaGeneratorFilter",
    FetchSherpack = cms.bool(False),
    SherpaDefaultWeight = cms.double(1.0),
    SherpaParameters = cms.PSet(
        MPI_Cross_Sections = cms.vstring(' MPIs in Sherpa, Model = Amisic:', 
            ' semihard xsec = 73.8879 mb,', 
            ' non-diffractive xsec = 18.1593 mb with nd factor = 0.335.'),
        Run = cms.vstring('(run){', 
            ' EVENTS = 100;', 
            ' EVENT_MODE = HepMC;', 
            ' HEPMC2_GENEVENT_OUTPUT = hepmc;', 
            ' # avoid comix re-init after runcard modification', 
            ' WRITE_MAPPING_FILE 3;', 
            ' ME_SIGNAL_GENERATOR Comix Amegic LOOPGEN;', 
            ' EVENT_GENERATION_MODE Weighted;', 
            ' LOOPGEN:=BlackHat;', 
            '}(run)', 
            '(beam){', 
            ' BEAM_1 = 2212; BEAM_ENERGY_1 = 6500.;', 
            ' BEAM_2 = 2212; BEAM_ENERGY_2 = 6500.;', 
            '}(beam)', 
            '(processes){', 
            ' Process 93 93 -> 15 -15 22;', 
            ' Order_EW 3;', 
            ' CKKW sqr(20./E_CMS);', 
            ' Print_Graphs MyGraphs;', 
            ' #Integration_Error 0.01;', 
            ' ME_Generator Amegic;', 
            ' End process;', 
            '}(processes)', 
            '(selector){', 
            ' Mass 15 -15 50 E_CMS;', 
            ' PT 22 10. E_CMS;', 
            ' PseudoRapidity 22 -3.0 3.0;', 
            ' PseudoRapidity 15 -3.0 3.0;', 
            ' PseudoRapidity -15 -3.0 3.0;', 
            ' DeltaR 22 15 0.6 100000;', 
            ' DeltaR 22 -15 0.6 100000;', 
            '}(selector)', 
            '(shower){', 
            ' CSS_EW_MODE = 1', 
            '}(shower)', 
            '(integration){', 
            ' FINISH_OPTIMIZATION = Off', 
            '}(integration)', 
            '(isr){', 
            ' PDF_LIBRARY     = LHAPDFSherpa', 
            ' PDF_SET         = CT10', 
            ' PDF_SET_VERSION = 0', 
            '}(isr)', 
            '(mi){', 
            ' MI_HANDLER = Amisic  # None or Amisic', 
            '}(mi)'),
        parameterSets = cms.vstring('MPI_Cross_Sections', 
            'Run')
    ),
    SherpaPath = cms.string('./'),
    SherpaPathPiece = cms.string('./'),
    SherpaProcess = cms.string('tautaugamma_50Mtautau13000_10ptg13000'),
    SherpaResultDir = cms.string('Result'),
    SherpackChecksum = cms.string('6cf02ad3dbf9206e54e69c2860cc249c'),
    SherpackLocation = cms.string(''),
    crossSection = cms.untracked.double(-1),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.int32(0)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


