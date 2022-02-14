import os,sys,subprocess,shlex
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from Stat.Limits.bruteForce import makeVarInfoList
from copy import deepcopy

# make status messages useful
def fprint(msg):
    import sys
    print(msg)
    sys.stdout.flush()

# system interaction
def runCmd(args,realtime=False):
    output = ""
    if realtime:
        p = subprocess.Popen(shlex.split(args), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, bufsize=1)
        for line in iter(p.stdout.readline, b''):
            fprint(line.rstrip())
        p.stdout.close()
        p.wait()
    else:
        try:
            output += subprocess.check_output(shlex.split(args),stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            output += e.output
    return output

def runCmds(commands):
    for command in commands:
        fprint(command)
        output = runCmd(command)
        fprint(output)

# basic (analysis-specific) helpers
def alphaVal(val):
    result = 0
    if val=="peak": result = -2
    elif val=="high": result = -1
    elif val=="low": result = -3
    else: result = float(val)
    return result

def rinvVal(val):
    result = None
    if len(val)>1 and val[0]=="0" and val[1]!=".": result = float("0."+val[1:])
    elif len(val)==1: result = int(val)
    else: result = float(val)
    return result

def paramVal(key,val):
    funcs = {
        "rinv": rinvVal,
        "alpha": alphaVal,
    }
    if key in funcs: return funcs[key](val)
    else: return float(val)

def getCombos():
    combos = {
        "cut": ["highCut","lowCut"],
        "bdt": ["highSVJ2","lowSVJ2"],
    }
    return combos

def getChannel(region):
    return "ch1" if "high" in region else "ch2"

def PoissonErrorUp(N):
    alpha = 1 - 0.6827 #1 sigma interval
    import ROOT as r
    U = r.Math.gamma_quantile_c(alpha/2,N+1,1.)
    return U-N

def OpenFile(fname):
    import ROOT as r
    f = r.TFile.Open(fname)
    if f==None: raise RuntimeError("Could not open ROOT file: {}".format(fname))
    return f

# various signal name handling
def getParamNames():
    params = ["mZprime", "mDark", "rinv", "alpha"]
    return params

def makeSigDict(values,names=None):
    if names is None: names=getParamNames()
    sigdict = {names[j]:values[j] for j in range(len(values))}
    return sigdict

def getSigname(sig):
    params = getParamNames()
    return "SVJ_"+"_".join("{}{}".format(key,sig[key]) for key in params)

def getSignameCheck(sig):
    return sig if isinstance(sig,str) else getSigname(sig)

def getSignameShort(sig):
    params = getParamNames()
    sig2 = deepcopy(sig)
    sig2["rinv"] = str(rinvVal(sig2["rinv"]))
    return "SVJ_"+"_".join(sig2[key] for key in params)

def getSignameLong(sig):
    params = getParamNames()
    sig2 = deepcopy(sig)
    sig2["rinv"] = str(rinvVal(sig2["rinv"]))
    return "SVJ_"+"_".join("{}-{}".format(key,sig2[key]) for key in params)

# generic check for signal directory
def getXname(sig,xname):
    signame = getSignameCheck(sig)
    if not os.path.basename(os.getcwd())==signame: xname = signame+"/"+xname
    return xname

# workspace filename
def getWname(sig,region):
    signame = getSignameCheck(sig)
    wname = "ws_{}_{}_2018_template.root".format(signame,region)
    return getXname(signame,wname)

# individual datacard
def getDname(sig,region):
    signame = getSignameCheck(sig)
    dname = "{}_{}_2018_template_bias.txt".format(signame,region)
    return getXname(signame,dname)

# combined datacard
def getDCname(sig,combo):
    dcname = "datacard_{}_{}.txt".format(getSignameShort(sig), combo)
    return getXname(sig,dcname)

# parameter names
def getPname(region, alt=False):
    pname = "Bkg{}_{}_2018".format("_Alt" if alt else "", region)
    return pname

# combine filename
def getFname(name, method, combo, sig=None, prefix="higgsCombine", seed=None):
    isCombine = prefix=="higgsCombine"
    fname = "{}{}{}{}{}{}.root".format(
        prefix,
        name,
        ".{}".format(method) if isCombine else "",
        ".mH120",
        ".ana{}".format(combo) if len(combo)>0 else "",
        ".{}".format(seed) if seed is not None else ""
    )
    if sig is not None:
        signame = getSignameCheck(sig)
        return getXname(signame,fname)
    else:
        return fname

def getRname(region, alt, npars):
    rname = "fitresult_{}{}_data_obs".format(getPname(region,alt),npars)
    return rname

# extract tracked params from trees
def getBranches(tree, matches=None, exact=False):
    import ROOT as r
    if not exact and not isinstance(matches,list): matches = [matches]
    elif exact and isinstance(matches,list): matches = matches[0]
    branches = []
    for b in tree.GetListOfBranches():
        bname = b.GetName()
        leaf = b.GetLeaf(bname)
        if matches is None or (exact and matches==b.GetName()) or (not exact and all(m in bname for m in matches)):
            branches.append(bname)
    return branches

def getParamsText(params):
    return ["{}={}".format(p.replace('trackedParam_','').replace('trackedError_',''),v) for p,v in sorted(params.iteritems()) if any(x in p for x in ['high','low','shapeBkg'])]

def getParamsTracked(fname, quantile, includeParam=True, includeErr=False, extraCondition=""):
    import ROOT as r

    condition = "abs(quantileExpected-{})<0.001".format(quantile)
    if len(extraCondition)>0: condition += "&&{}".format(extraCondition)
    results = {}
    file = OpenFile(fname)
    tree = file.Get("limit")

    matches = []
    if includeParam and not includeErr: matches = ["trackedParam"]
    elif not includeParam and includeErr: matches = ["trackedError"]
    params = []
    # background normalization factors
    params.extend(getBranches(tree, matches+["shapeBkg"]))
    # background & signal final normalizations
    params.extend(getBranches(tree, matches+["n_exp_final"]))
    # signal strength (from MDF)
    params.extend(getBranches(tree, "r", exact=True))

    for region in ["high","low"]:
        # background fit parameters
        params.extend(sorted(getBranches(tree, matches+[region])))

    # deliver dict of params:values
    n = tree.Draw(':'.join(params),condition,"goff")
    if n<=0:
        raise RuntimeError("Could not get:\n\tparameters: {}\n\tcondition: {}\n\tfrom limit tree in file {}".format(params,condition,fname))

    results = {p:tree.GetVal(i)[0] for i,p in enumerate(params)}
    return results

# extract brute force initial values
def getInitFromBF(fname, wsname, pdfname, region=None):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    file = OpenFile(fname)
    ws = file.Get(wsname)
    pdf = ws.pdf(pdfname)
    pars = makeVarInfoList(pdf.getPars())
    npars = len(pars)

    if region is not None:
        # need to use fit results because later saved workspaces don't track the errors
        fname2 = "fitResults_{}.root".format(region)
        if not os.path.isfile(fname2): fname2 = "../"+fname2
        file2 = OpenFile(fname2)
        result = file2.Get(getRname(region,"Alt" in pdfname, npars))

        setargs = {p.name:result.floatParsFinal().find(p.name).getValV() for p in pars}
        errargs = {p.name:result.floatParsFinal().find(p.name).getError() for p in pars}
        return setargs, errargs
    else:
        setargs = ["{}={}".format(p.name,p.val) for p in pars]
        return setargs

