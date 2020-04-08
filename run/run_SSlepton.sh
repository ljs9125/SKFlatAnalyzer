#DATA
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i SingleMuon -y 2016 -n 20 &

#MC
#QCD
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-1000toInf_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-800to1000_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-600to800_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-470to600_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-300to470_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-170to300_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-120to170_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-80to120_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-50to80_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-30to50_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-20to30_MuEnrichedPt5 -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i QCD_Pt-15to20_MuEnrichedPt5 -y 2016 -n 5 &

#VV
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i WW_pythia -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i WZ_pythia -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i ZZ_pythia -y 2016 -n 5 &

#W,Z
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i WJets_MG -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i DYJets_MG -y 2016 -n 5 &

#Top
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i SingleTop_tW_antitop_NoFullyHad -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i SingleTop_tW_top_NoFullyHad -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i TTLJ_powheg -y 2016 -n 5 &
python python/SKFlat.py -a SSlepton --skim SkimTree_NIsoMuon -i TTLL_powheg -y 2016 -n 5 &







