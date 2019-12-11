R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(/data4/Users/snuintern1/SKFlatAnalyzer/lib/libDataFormats.so)
R__LOAD_LIBRARY(/data4/Users/snuintern1/SKFlatAnalyzer/lib/libAnalyzerTools.so)
R__LOAD_LIBRARY(/data4/Users/snuintern1/SKFlatAnalyzer/lib/libAnalyzers.so)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/lib/libLHAPDF.so)


void run_1(){

  SSdimuon m;

  m.SetTreeName("recoTree/SKFlat");

  m.LogEvery = 1000;
  m.IsDATA = true;
  m.DataStream = "SingleMuon";
  m.DataYear = 2016;
  m.AddFile("/data7//DATA/SKFlat/Run2Legacy_v3/2016/DATA/SingleMuon/periodC/190421_220953/0000/SKFlatNtuple_2016_DATA_1.root");
  m.AddFile("/data7//DATA/SKFlat/Run2Legacy_v3/2016/DATA/SingleMuon/periodC/190421_220953/0000/SKFlatNtuple_2016_DATA_10.root");
  m.AddFile("/data7//DATA/SKFlat/Run2Legacy_v3/2016/DATA/SingleMuon/periodC/190421_220953/0000/SKFlatNtuple_2016_DATA_100.root");
  m.AddFile("/data7//DATA/SKFlat/Run2Legacy_v3/2016/DATA/SingleMuon/periodC/190421_220953/0000/SKFlatNtuple_2016_DATA_101.root");
  m.SetOutfilePath("hists.root");
  m.Init();
  m.initializeAnalyzerTools();
  m.initializeAnalyzer();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}
