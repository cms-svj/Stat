import math
from collections import OrderedDict
from copy import deepcopy

# sources:
# SM Z: https://pdg.lbl.gov/2021/reviews/rpp2021-rev-standard-model.pdf
# Zprime: https://pdg.lbl.gov/2021/reviews/rpp2021-rev-zprime-searches.pdf, https://arxiv.org/abs/1010.6058, https://arxiv.org/abs/0801.1345

# constants
sin2thW = 0.2229
vev = 246.22
mZ = 91.1876 # PDG Table 10.5
# move other factors into couplings: (mZ/vev)^2 (from eq. 10.61)
width0 = mZ/vev

# baseline cross section values (delta eta acceptance, branching fractions divided out)
mZprime = 3100
sigma_SSM = 0.1338E-01/.41/0.66
sigma_universal = 0.2397E-01/.41/(5./6.)
g_universal = 0.25

# PDG eq. 10.6
def gV(t3, Q):
    return t3-2*Q*sin2thW

def gA(t3):
    return t3

# PDG eq. 10.61
# equiv. to PDG eq. 87.6 w/ "other factors" removed
# tree-level: assume R_V,A = Nc
# generalize coupling
widthD = 1/(12*math.pi)
def partial(mBoson, Nc, g, units="GeV"):
    width = widthD*Nc*mBoson*(g**2)
    if units=="MeV": width *= 1000
    return width

def total(mZ, partial, models, units="GeV", verbose=False):
    total_width = 0
    branchings = []
    for model in models:
        model_width = 0
        for part,prop in model.items():
            width = partial(mZ, prop["Nc"], prop["g"], units)
            if verbose: print("width ({}) = {:.0f}".format(part,width))
            model_width += width
        branchings.append(model_width)
        total_width += model_width
    branchings = [x/total_width for x in branchings]
    return total_width, branchings

# compute from SM functions
def get_couplings(model):
    for part,props in model.items():
        model[part]["gV"] = gV(props["t3"],props["Q"])
        model[part]["gA"] = gA(props["t3"])
        model[part]["g0"] = math.sqrt(model[part]["gV"]**2 + model[part]["gA"]**2)
        model[part]["g"] = model[part]["g0"]*width0
    return model

# to set an arbitrary value
def set_couplings(model,val):
    for part,props in model.items():
        model[part]["g"] = val
    return model

SM_leptons = OrderedDict()
SM_leptons["e"] = dict(Nc = 1, t3 = -0.5, Q = -1)
SM_leptons["mu"] = deepcopy(SM_leptons["e"])
SM_leptons["tau"] = deepcopy(SM_leptons["e"])
SM_leptons["nu_e"] = dict(Nc = 1, t3 = 0.5, Q = 0)
SM_leptons["nu_mu"] = deepcopy(SM_leptons["nu_e"])
SM_leptons["nu_tau"] = deepcopy(SM_leptons["nu_e"])
SM_leptons = get_couplings(SM_leptons)

SM_quarks = OrderedDict()
SM_quarks["u"] = dict(Nc = 3, t3 = 0.5, Q = 2./3.)
SM_quarks["c"] = deepcopy(SM_quarks["u"])
SM_quarks["t"] = deepcopy(SM_quarks["u"])
SM_quarks["d"] = dict(Nc = 3, t3 = -0.5, Q = -1./3.)
SM_quarks["s"] = deepcopy(SM_quarks["d"])
SM_quarks["b"] = deepcopy(SM_quarks["d"])
SM_quarks = get_couplings(SM_quarks)

# universal-style (leptophobic) SM
SM_universal = deepcopy(SM_quarks)

# nominal dark sector
HV = OrderedDict()
HV["x1"] = dict(Nc = 2, g = 1.0)
HV["x2"] = dict(Nc = 2, g = 1.0)

# dark sector w/ minimum flavor coupling
HV1 = OrderedDict()
HV1["x1"] = dict(Nc = 2, g = 1.0)

# ROOT stuff
import ROOT as r
base_qtys = ["quantileExpected","rate","width"]
param_names = ["g_q","g_dark","xsec","maxwidth"]
qtys = base_qtys + ["trackedParam_{}".format(q) for q in param_names]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")
qobj = r.quantile_t()
qobj.quantileExpected = 0.5
qobj.trackedParam_xsec = 1.0
qobj.trackedParam_maxwidth = 10.0

# function to make a limit-style tree with production rate info
def make_tree(suff, looper):
    file = r.TFile.Open("sigRate_{}.root".format(suff),"RECREATE")
    tree = r.TTree("limit","limit")
    for qty in qtys:
        tree.Branch(qty, r.AddressOf(qobj,qty), '{}/D'.format(qty))
    looper(tree)
    file.Write()
    file.Close()

# parameter ranges
g_q_values = [x/10. for x in range(1,11)]
g_q_values.insert(2, 0.25) # benchmark
g_dark_values = [x/10. for x in range(1,21)]

def looper_SSM(tree,HV_model):
    for g_dark in g_dark_values:
        # reset HV couplings
        HV_model = set_couplings(HV_model, g_dark)

        # get widths
        total_width, BRs = total(mZprime, partial, [SM_leptons,SM_quarks,HV_model])

        # fill quantities
        qobj.trackedParam_g_q = SM_quarks["u"]["g"] # arbitrary choice of value
        qobj.trackedParam_g_dark = g_dark
        qobj.rate = sigma_SSM*BRs[-1]
        qobj.width = total_width/mZprime*100

        # fill tree
        tree.Fill()

def looper_SSM_HV(tree):
    looper_SSM(tree,HV)

def looper_SSM_HV1(tree):
    looper_SSM(tree,HV1)

def looper_universal(tree,HV_model):
    # universal-style (leptophobic) SM
    SM_universal = deepcopy(SM_quarks)

    for g_q in g_q_values:
        # reset SM couplings
        SM_universal = set_couplings(SM_universal, g_q)

        # update cross section
        sigma_universal_now = sigma_universal*(g_q/g_universal)**2

        for g_dark in g_dark_values:
            # reset HV couplings
            HV_model = set_couplings(HV_model, g_dark)

            # get widths
            total_width, BRs = total(mZprime, partial, [SM_universal,HV_model])

            # fill quantities
            qobj.trackedParam_g_q = g_q
            qobj.trackedParam_g_dark = g_dark
            qobj.rate = sigma_universal_now*BRs[-1]
            # divide out maxwidth since it is used as xsec to multiply limit and find intersection
            qobj.width = total_width/mZprime*100/qobj.trackedParam_maxwidth

            # fill tree
            tree.Fill()

def looper_universal_HV(tree):
    looper_universal(tree,HV)

def looper_universal_HV1(tree):
    looper_universal(tree,HV1)

make_tree("SSM", looper_SSM_HV)
make_tree("SSMHV1", looper_SSM_HV1)
make_tree("universal", looper_universal_HV)
make_tree("universalHV1", looper_universal_HV1)

