import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import getParamNames, paramVal

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", dest="input", type=str, required=True, help="input file w/ columns: mZprime mDark rinv alpha mean mean_err stdev stdev_err")
parser.add_argument("-x", "--suffix", dest="suffix", type=str, required=True, help="suffix for output filename")
args = parser.parse_args()

import ROOT as r

param_names = getParamNames()

# make a Combine-esque tree in order to reuse plotLimit code
base_qtys = ["quantileExpected"]
keys = ["mean","meanErr","std","stdErr"]
qtys = base_qtys + keys + ["trackedParam_{}".format(q) for q in param_names+["xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")
qobj = r.quantile_t()
# default values
qobj.quantileExpected = 0.5
qobj.trackedParam_xsec = 1.0

file = r.TFile.Open("bias_dijet{}.root".format("_"+args.suffix if len(args.suffix)>0 else ""),"RECREATE")
tree = r.TTree("limit","limit")
for qty in qtys:
    tree.Branch(qty, r.AddressOf(qobj,qty), '{}/D'.format(qty))

# read txt file
with open(args.input,'r') as input_file:
    for line in input_file:
        # parse columns
        linesplit = line.rstrip().split()
        sig = {j:0 for j in param_names}
        sig["mZprime"], sig["mDark"], sig["rinv"], sig["alpha"], mean, meanErr, std, stdErr = linesplit

        # fill tree
        for q in param_names:
            setattr(qobj,"trackedParam_{}".format(q),paramVal(q,sig[q]))
        qobj.mean = float(mean)
        qobj.meanErr = float(meanErr)
        qobj.std = float(std)
        qobj.stdErr = float(stdErr)
        tree.Fill()

file.Write()
file.Close()
