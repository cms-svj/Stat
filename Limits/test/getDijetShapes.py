import os
import ROOT as r

# from https://github.com/CMSDIJET/DijetRootTreeAnalyzer/blob/master/python/WriteDataCard.py
def getDownFromUpNom(hUp,hNom):
    hDown = hUp.Clone(hUp.GetName()+'Down')
    for i in range(1,hDown.GetNbinsX()+1):
        nom = hNom.GetBinContent(i)
        up = hUp.GetBinContent(i)
        if up > 0:
            down = nom*nom / up 
            hDown.SetBinContent(i, down)
        else:
            hDown.SetBinContent(i, 0)
    return hDown

systs = [
    "",
    "_jesUp",
    "_jesDown",
    "_jerUp",
    "_jerDown",
]

# check for dijet files
jnames = ["ResonanceShapes_qq_13TeV_Spring16{}.root".format(x) for x in [y.upper() for y in systs]]
for jname in jnames:
    if not os.path.isfile(jname):
        # make JERDOWN from JERUP w/ getDownFromUpNom()
        if "JERDOWN" in jname:
            fNom = r.TFile.Open(jnames[0])
            fUp = r.TFile.Open(next(name for name in jnames if "JERUP" in name))
            fDown = r.TFile.Open(jname,"RECREATE")
            # loop over all histos
            hnames = [k.GetName() for k in fNom.GetListOfKeys()]
            for hname in hnames:
                hNom = fNom.Get(hname)
                hUp = fUp.Get(hname)
                hDown = getDownFromUpNom(hUp,hNom)
                fDown.cd()
                hDown.Write(hname)
            fDown.Close()
        else:
            os.system("wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/{}".format(jname))
