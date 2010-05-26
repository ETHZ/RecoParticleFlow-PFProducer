{
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gSystem->Load("libCintex.so");
  ROOT::Cintex::Cintex::Enable();
  TFile f("muonsFromPF.root");


  Events.Draw("recoPFCandidates_particleFlow__RECO.obj.charge()","recoPFCandidates_particleFlow__RECO.obj.particleId()==3","");
  Events.SetLineColor(4);
  Events.Draw("recoMuons_recoMuonFromPFProducer__ANALYSIS.obj.charge()","","same");
  
}
