import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import getParamNames, makeSigDict, getSignameLong, paramVal

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--dir", dest="dir", type=str, default="root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/Run2ProductionV17/scan", help="ntuple directory")
parser.add_argument("-s", "--signals", dest="signals", type=str, required=True, help="file w/ list of signals")
parser.add_argument("-c", "--cuts", dest="cuts", type=str, default="", help="cuts to apply")
parser.add_argument("-x", "--suffix", dest="suffix", type=str, required=True, help="suffix for output filename")
parser.add_argument("-r", "--reco", dest="reco", default=False, action="store_true", help="use reco instead of gen quantities")
args = parser.parse_args()

import ROOT as r
r.gInterpreter.ProcessLine('#include "{}/src/Analysis/KCode/KMath.h"'.format(os.environ["CMSSW_BASE"]))
drawname = "GenMT>>h(80,0,8000)"
if args.reco: drawname = drawname.replace("GenMT","MT_AK8")

param_names = getParamNames()
param_values = []
with open(args.signals,'r') as sfile:
    for line in sfile:
        line = line.rstrip()
        if len(line)==0: continue
        param_values.append(line.split())

signals = [makeSigDict(param_values[i],param_names) for i in range(len(param_values))]

# make a Combine-esque tree in order to reuse plotLimit code
base_qtys = ["quantileExpected","limit"]
qtys = base_qtys + ["trackedParam_{}".format(q) for q in param_names+["xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")
qobj = r.quantile_t()
qobj.quantileExpected = 0.5
qobj.trackedParam_xsec = 1.0

file = r.TFile.Open("sigWidth_{}.root".format(args.suffix),"RECREATE")
tree = r.TTree("limit","limit")
for qty in qtys:
    tree.Branch(qty, r.AddressOf(qobj,qty), '{}/D'.format(qty))

for signal in signals:
    sfile = r.TFile.Open("{}/{}_MC2018_scan.root".format(args.dir,getSignameLong(signal)))
    stree = sfile.Get("tree")
    stree.SetAlias("GenMT","KMath::TransverseMass(GenJetsAK8[0].Px()+GenJetsAK8[1].Px(),GenJetsAK8[0].Py()+GenJetsAK8[1].Py(),sqrt((GenJetsAK8[0].E()+GenJetsAK8[1].E())^2-(GenJetsAK8[0].Px()+GenJetsAK8[1].Px())^2-(GenJetsAK8[0].Py()+GenJetsAK8[1].Py())^2-(GenJetsAK8[0].Pz()+GenJetsAK8[1].Pz())^2),GenMET*cos(GenMETPhi),GenMET*sin(GenMETPhi),0)")
    stree.Draw(drawname,args.cuts,"goff")
    shist = r.gDirectory.Get("h")
    qobj.limit = shist.GetRMS()/shist.GetMean()
    for q in getParamNames():
        setattr(qobj,"trackedParam_{}".format(q),paramVal(q,signal[q]))
    # skip xsec
    tree.Fill()

file.Write()
file.Close()
