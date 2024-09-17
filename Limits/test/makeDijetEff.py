import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import getParamNames, makeSigDict, getSignameShort, paramVal

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", dest="input", type=str, default="DijetShapeInterpolator/inputs/svj_dijet.root", help="input file w/ histograms")
parser.add_argument("-s", "--signals", dest="signals", type=str, required=True, help="file w/ list of signals")
parser.add_argument("-x", "--suffix", dest="suffix", type=str, default="", help="suffix for output filename")
args = parser.parse_args()

import ROOT as r

param_names = getParamNames()
param_values = []
with open(args.signals,'r') as sfile:
    for line in sfile:
        line = line.rstrip()
        if len(line)==0: continue
        param_values.append(line.split())

signals = [makeSigDict(param_values[i],param_names) for i in range(len(param_values))]

# make a Combine-esque tree in order to reuse plotLimit code
base_qtys = ["quantileExpected"]
keys = ["cmsdijet","cmsdijetCRmiddle","cmsdijetCRhigh"]
qtys = base_qtys + keys + ["trackedParam_{}".format(q) for q in param_names+["xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")
qobj = r.quantile_t()
qobj.quantileExpected = 0.5
qobj.trackedParam_xsec = 1.0

file = r.TFile.Open("dijetEff{}.root".format("_"+args.suffix if len(args.suffix)>0 else ""),"RECREATE")
tree = r.TTree("limit","limit")
for qty in qtys:
    tree.Branch(qty, r.AddressOf(qobj,qty), '{}/D'.format(qty))

# open input file
input_file = r.TFile.Open(args.input)

for signal in signals:
    signame = getSignameShort(signal)
    for key in keys:
        # get efficiency averaged over all three years
        numer = 0
        denom = 0
        for year in ["MC2016","MC2017","MC2018"]:
            hsuff = "{}_{}_{}".format(signame,year,key)
            cutflow = input_file.Get("cutflow_{}".format(hsuff))
            nEventProc = input_file.Get("nEventProc_{}".format(hsuff))
            numer += cutflow.GetBinContent(cutflow.GetNbinsX())
            denom += nEventProc.GetBinContent(1)
        eff = float(numer)/float(denom)
        setattr(qobj,key,eff)
    for q in getParamNames():
        setattr(qobj,"trackedParam_{}".format(q),paramVal(q,signal[q]))
    tree.Fill()

file.Write()
file.Close()
