import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import getParamNames, paramVal
from runLimitsPool import getXsecs

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", dest="input", type=str, default="monojet.txt", help="input file w/ columns: mass rinv obs exp")
parser.add_argument("-x", "--suffix", dest="suffix", type=str, default="", help="suffix for output filename")
args = parser.parse_args()

import ROOT as r

param_names = getParamNames()
xsecs = getXsecs('dict_xsec_Zprime.txt')

# make a Combine-esque tree in order to reuse plotLimit code
base_qtys = ["quantileExpected"]
keys = ["limit"]
qtys = base_qtys + keys + ["trackedParam_{}".format(q) for q in param_names+["xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")
qobj = r.quantile_t()
# default values
qobj.trackedParam_mDark = 20
qobj.trackedParam_alpha = -2

file = r.TFile.Open("limit_monojet{}.root".format("_"+args.suffix if len(args.suffix)>0 else ""),"RECREATE")
tree = r.TTree("limit","limit")
for qty in qtys:
    tree.Branch(qty, r.AddressOf(qobj,qty), '{}/D'.format(qty))

# read txt file
with open(args.input,'r') as input_file:
    skipped_header = False
    for line in input_file:
        if not skipped_header:
            skipped_header = True
            continue
        # parse columns
        mZprime, rinv, obs, exp = line.rstrip().split()
        # no yield
        if rinv=="0": continue
        xsec = xsecs[mZprime]

        # fill tree
        qobj.trackedParam_mZprime = float(mZprime)
        qobj.trackedParam_rinv = float(rinv)
        qobj.trackedParam_xsec = float(xsec)

        for lim,qtl in zip([exp,obs],[0.5,-1]):
            qobj.limit = float(lim)
            qobj.quantileExpected = qtl
            tree.Fill()

file.Write()
file.Close()
