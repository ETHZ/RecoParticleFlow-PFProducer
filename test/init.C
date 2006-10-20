{

gSystem->Load("libFWCoreFWLite.so");
gSystem->Load("libRecoParticleFlowPFProducer.so");
AutoLibraryLoader::enable();
gSystem->Load("libCintex.so");
ROOT::Cintex::Cintex::Enable();

}
