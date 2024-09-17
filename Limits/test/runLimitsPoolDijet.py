import os,sys,shlex,traceback
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool
from StringIO import StringIO
from array import array
from string import Template
from paramUtils import fprint, runCmd, alphaVal, makeSigDict, getParamNames, getSignameShort, OpenFile, paramVal
from Stat.Limits.bruteForce import silence

def DtoF(hD):
    from ROOT import TH1F
    hF = TH1F()
    hD.Copy(hF)
    return hF

def getMethods():
    return ["fit","ratio"]

def getXsecs(fname):
    with open(fname,'r') as xfile:
        xsecs = {xline.split('\t')[0]: float(xline.split('\t')[1]) for xline in xfile}
    return xsecs

def fill_template(inname, outname, **kwargs):
    if inname==outname:
        raise ValueError("Attempt to overwrite template file: "+inname)
    with open(inname,'r') as temp:
        old_lines = Template(temp.read())
        new_lines = old_lines.substitute(**kwargs)
    with open(outname,'w') as temp:
        temp.write(new_lines)

def transform(iname,wname,w2name,template,sig,fullname,extargs,region,acc,args):
    lumi = 36330+41530+59740
    acc_svj_fname = "dijetEff_svj.root"
    dmWeight = 1 if "DM" in args.finalState else 0
    smWeight = 1 if "SM" in args.finalState else 0

    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    # get SVJ eff for this signal in region: use graph to interpolate for intermediate mass values
    acc_svj_file = OpenFile(acc_svj_fname)
    acc_svj_tree = acc_svj_file.Get("limit")
    drawname = "trackedParam_mZprime:{}".format(region)
    cutname = "&&".join("abs(trackedParam_{}-{})<0.0001".format(key,paramVal(key,sig[key])) for key in ["mDark","rinv","alpha"])
    n = acc_svj_tree.Draw(drawname,cutname,"goff")
    acc_graph = r.TGraph(n, acc_svj_tree.GetV1(), acc_svj_tree.GetV2())
    acc_svj = acc_graph.Eval(float(sig["mZprime"]))

    signame = getSignameShort(sig)
    tname = template.format("template")
    ofname = template.format(fullname).replace(".txt",".root")
    hofname = "hist_"+ofname
    dcfname = template.format(fullname)

    if args.dry_run: return dcfname

    systs = [
        "",
        "_jesUp",
        "_jesDown",
        "_jerUp",
        "_jerDown",
    ]
    # dijet search used SR shapes for CRs
    jnames = ["ResonanceShapes_qq_13TeV_Spring16{}.root".format(x) for x in [y.upper() for y in systs]]
    svjnames = ["ResonanceShapes_svj_13TeV_{}{}.root".format(region,x) for x in systs]

    ifile = OpenFile(iname)
    jfiles = [OpenFile(jname) for jname in jnames]
    svjfiles = [OpenFile(svjname) for svjname in svjnames]
    w = ifile.Get(wname)

    xlist = [2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5455,5663,5877,6099,6328,6564,6808,7060,7320,7589,7866,8152,8752]
    if args.method=="fit":
        xlist = [1530,1607,1687,1770,1856,1945,2037,2132,2231,2332]+xlist
    xbins = array('d',xlist)
    th1x = w.var("th1x")
    bins = th1x.getBinning().array()
    bins.SetSize(th1x.getBinning().numBins()+1)
    bins = array('d',bins)

    # make new workspace, import contents from prev (but take new data)
    w2 = r.RooWorkspace(w2name)
    getattr(w2,"import")(w.allVars())
    getattr(w2,"import")(w.allPdfs())
    getattr(w2,"import")(w.allFunctions())
    getattr(w2,"import")(w.data("data_obs"))

    # get sig histos
    hyield = 0
    for isyst,syst in enumerate(systs):
        if syst=="":
            hofile = r.TFile.Open(hofname,"RECREATE")

        # SVJ histo
        histname = "h_svj_{}".format(signame.replace("SVJ_","").replace(".",""))
        hist = svjfiles[isyst].Get(histname)
        if not isinstance(hist,r.TH1): raise RuntimeError("Could not get {} from {}".format(histname,svjnames[isyst]))
        hist.SetDirectory(0)
        hist = DtoF(hist)
        hist.Scale(dmWeight*sig["xsecSVJ"]*acc_svj*lumi)
        hist_rebinned = hist.Rebin(len(xbins)-1,histname+"_rebin",xbins)

        # add dijet sig histo (scale by xsec * BR * A * lumi)
        dhistname = "h_qq_{}".format(sig["mZprime"])
        dhist = jfiles[isyst].Get(dhistname)
        if not isinstance(dhist,r.TH1): raise RuntimeError("Could not get {} from {}".format(dhistname,jnames[isyst]))
        dhist.SetDirectory(0)
        dhist = DtoF(dhist)
        dhist.Scale(smWeight*sig["xsecDijet"]*acc*lumi)
        dhist_rebinned = dhist.Rebin(len(xbins)-1,dhistname+"_rebin",xbins)

        shist = hist.Clone(histname.replace("svj","sum"))
        shist.Add(dhist)
        shist_rebinned = shist.Rebin(len(xbins)-1,histname.replace("svj","sum")+"_rebin",xbins)

        # write out before resetting bins
        if syst=="":
            hyield = shist_rebinned.Integral()
            hist.Write()
            hist_rebinned.Write()
            dhist.Write()
            dhist_rebinned.Write()
            shist.Write()
            shist_rebinned.Write()

        # reset sig histo bins
        shist_rebinned.GetXaxis().Set(len(bins)-1,bins)

        # import new pdf
        hname2 = signame+syst
        pdf = r.RooDataHist(hname2, hname2, r.RooArgList(th1x), shist_rebinned, 1.);
        getattr(w2,"import")(pdf)

        if syst=="":
            hofile.Close()

    # make datacard
    fill_template(
        tname,
        dcfname,
        signame = signame,
        sigrate = hyield,
        wsfile = ofname,
        wsname = w2name,
        extargs = extargs,
    )

    # write out modified workspace
    if os.path.exists(ofname): os.remove(ofname)
    w2.writeToFile(ofname)

    return dcfname

def getFullname(sig,method,finalState):
    if isinstance(sig,str): signame = sig
    else: signame = getSignameShort(sig)
    return signame+"_"+method+"_"+"_".join(sorted(finalState))

def doLimit(info):
    outputs = []
    args = info["args"]
    sig = info["sig"]
    fullname = getFullname(sig,args.method,args.finalState)

    BR_qq = 5./6.
    BR_dark = 0.47
    sig["xsecSVJ"] = sig["xsec"]
    sig["xsecDijet"] = sig["xsec"]/BR_dark*(1-BR_dark)*BR_qq
    if args.xsec is not None: sig["xsecDijet"] = args.xsec
    # update w/ total xsec
    sig["xsec"] = 0
    if "DM" in args.finalState: sig["xsec"] += sig["xsecSVJ"]
    if "SM" in args.finalState: sig["xsec"] += sig["xsecDijet"]

    params = {key:paramVal(key,val) for key,val in sig.iteritems()}
    params["method"] = int(args.method=="ratio")+1
    setargs = []
    trkargs = []
    extargs = ""
    for p,v in params.iteritems():
        setargs.append(p+"="+str(v))
        trkargs.append(p)
        extargs = extargs+p+" extArg "+str(v)+"\n"
    frzargs = trkargs[:]

    dijet_acc = args.acc
    if args.method=="fit":
        fitparams = ["p1_PFDijet2017","p2_PFDijet2017","p3_PFDijet2017","shapeBkg_PFDijet2017_bkg_PFDijet2017__norm"]
        trkargs.extend(fitparams)
        if args.freezeNorm: frzargs.append("shapeBkg_PFDijet2017_bkg_PFDijet2017__norm")

        # datacard setup
        dcfname = transform("dijet_combine_qq_1900_lumi-137.500_PFDijet2017.root","wPFDijet2017","wPFDijet2018","dijet_combine_fit_{}.txt",sig,fullname,extargs,"cmsdijet",dijet_acc,args)

    elif args.method=="ratio":
        dcfnameSR = transform("dijet_combine_qq_5000_lumi-136.600_PFDijet2017MC.root","wPFDijet2017MC","wPFDijet2018MC","dijet_combine_SR_{}.txt",sig,fullname,extargs,"cmsdijet",dijet_acc,args)
        dcfnameCR = transform("dijet_combine_qq_5000_lumi-136.600_PFDijet2017MCCR.root","wPFDijet2017MCCR","wPFDijet2018MCCR","dijet_combine_CRhigh_{}.txt",sig,fullname,extargs,"cmsdijetCRhigh",dijet_acc*0.45,args)
        dcfnameCRmid = transform("dijet_combine_qq_5000_lumi-136.600_PFDijet2017MCCRmid.root","wPFDijet2017MCCRmid","wPFDijet2018MCCRmid","dijet_combine_CRmid_{}.txt",sig,fullname,extargs,"cmsdijetCRmiddle",dijet_acc*0.27,args)

        # use combined card
        dcfname = dcfnameSR.replace("SR","ratio")
        command = "combineCards.py "+" ".join([dcfnameSR,dcfnameCR,dcfnameCRmid])+" > "+dcfname
        outputs.append(command)
        if not args.dry_run: os.system(command)
        # add ext args
        with open(dcfname,'a') as dcfile:
            dcfile.write(extargs)

    cargs = [
        "--setParameters {}".format(','.join(setargs)),
        "--freezeParameters {}".format(','.join(frzargs)),
        "--trackParameters {}".format(','.join(trkargs)),
        "--keyword-value sig={}".format(fullname),
        "-n {}".format(args.cname),
        "-d {}".format(dcfname),
    ]
    if len(args.args)>0: cargs.append(args.args)

    # run combine
    if args.bias:
        bname = "DijetSig{}".format(args.injectSignal)
        bargs = [
            "--expectSignal {}".format(args.injectSignal),
            "-t {}".format(args.ntoys),
            "--toysFrequentist",
            "--saveToys",
            "--bypassFrequentistFit",
        ]

        gname = args.cname.replace(args.name,bname+"Gen")
        gargs = [v for i,v in enumerate(cargs) if i not in [2,4]]
        gargs.extend(bargs)
        gargs.extend([
            "-n {}".format(gname),
            "-s {}".format(args.seed),
            "--saveWorkspace",
        ])

        fname = args.cname.replace(args.name,bname+"Fit")
        fargs = [v for i,v in enumerate(cargs) if i not in [4]]
        fargs.extend(bargs)
        fargs.extend([
            "-n {}".format(fname),
            "--toysFile higgsCombine{}.GenerateOnly.mH120.sig{}.{}.root".format(gname,fullname,args.seed),
            "--rMin {}".format(-args.rmax),
            "--rMax {}".format(args.rmax),
            "--savePredictionsPerToy",
        ])

        gargs = ' '.join(gargs)
        fprint("Generating toys for {}...".format(fullname))
        gcommand = "combine -M GenerateOnly {}".format(gargs)
        outputs.append(gcommand)
        if not args.dry_run: outputs.append(runCmd(gcommand))

        fargs = ' '.join(fargs)
        fprint("Running fits for {}...".format(fullname))
        fcommand = "combine -M FitDiagnostics {}".format(fargs)
        outputs.append(fcommand)
        if not args.dry_run: outputs.append(runCmd(fcommand))
    else:
        cargs = ' '.join(cargs)
        fprint("Calculating limit for {}...".format(fullname))
        command = "combine -M AsymptoticLimits {}".format(cargs)
        outputs.append(command)
        if not args.dry_run: outputs.append(runCmd(command))

    return outputs

def main(args):
    cname = args.name[:]
    if args.freezeNorm: cname += "Frz"
    args.cname = cname

    # assumes getDijetShapes has already been called

    if not args.just_hadd:
        if args.npool==0:
            for outputs in [doLimit({"args":args,"sig":sig}) for sig in args.signals]:
                fprint('\n'.join(outputs))
        else:
            p = Pool(args.npool if not args.dry_run else 1)
            for outputs in p.imap_unordered(doLimit, [{"args":args,"sig":sig} for sig in args.signals]):
                fprint('\n'.join(outputs))
            p.close()
            p.join()

    import ROOT as r
    if len(args.updateXsec)>0:
        xsecs = getXsecs(args.updateXsec[0])
        try:
            r.xsec_t
        except:
            r.gROOT.ProcessLine("struct xsec_t { Float_t mZprime; Float_t xsec; Double_t limit; };")
    outtrees = []
    outtreesroot = r.TList()
    if args.bias: fprint("\nBIAS FIT RESULTS")
    for sig in args.signals:
        if args.bias:
            bname = "DijetSig{}".format(args.injectSignal)
            fname = args.cname.replace(args.name,bname+"Fit")
            fullname = getFullname(sig,args.method,args.finalState)
            fitfname = "fitDiagnostics{}.mH120.sig{}.{}.root".format(fname,fullname,args.seed)
            fitf = r.TFile.Open(fitfname)
            if fitf!=None:
                fittree = fitf.Get("tree_fit_sb")
                if fittree!=None:
                    minPull = -5
                    maxPull = 5
                    nBins = int(args.ntoys/5.0)
                    selection = "fit_status==0"
                    hname = "pull_{}".format(fullname)
                    fittree.Draw("(r-{})/rErr>>{}({},{},{})".format(args.injectSignal,hname,nBins,minPull,maxPull),selection,"goff")
                    gaus = r.TF1("gaus","gaus(0)",minPull,maxPull)
                    hpull = r.gDirectory.Get(hname)
                    hpull.Fit(gaus,"RQN")
                    fprint("{} {} {} {} {}".format(
                        ' '.join(["{}".format(sig[key]) for key in getParamNames()]),
                        gaus.GetParameter(1),
                        gaus.GetParError(1),
                        gaus.GetParameter(2),
                        gaus.GetParError(2),
                    ))
                    continue
            fprint("ERROR: bias fit for {} did not converge".format(fullname))
            continue
        if args.method.startswith("best"):
            # decide which method to use based on criterion
            fullname = getFullname(sig,"ratio",args.finalState)
            histfname = "hist_dijet_combine_SR_{}.root".format(fullname)
            if len(args.hadd_dir)>0: histfname = args.hadd_dir+"/"+histfname
            histf = r.TFile.Open(histfname)
            signame = getSignameShort(sig)
            if args.method=="best":
                shistname = "h_sum_{}".format(signame.replace("SVJ_","").replace(".",""))
                shist = histf.Get(shistname)
            elif args.method=="best2":
                shistname = "h_qq_{}".format(sig["mZprime"])
                shist = histf.Get(shistname)
            # best3 uses a different criterion
            if args.method=="best3":
                eff = 1
                if int(sig["mZprime"])<3000:
                    method = "fit"
                else:
                    method = "ratio"
            else:
                # criterion: >1 sigma, one-sided
                min_eff = 0.84135
                min_mjj = 2438
                eff = shist.Integral(shist.FindBin(min_mjj),-1)/shist.Integral(-1,-1)
                if eff > min_eff:
                    method = "ratio"
                else:
                    method = "fit"
            print("{}: chose method {} (eff = {:.4f})".format(signame, method, eff))
        else:
            method = args.method
        fullname = getFullname(sig,method,args.finalState)
        fname = "higgsCombine{}.AsymptoticLimits.mH120.sig{}.root".format(cname,fullname)
        if len(args.hadd_dir)>0: fname = args.hadd_dir+"/"+fname
        if args.dry_run:
            outtrees.append(fname)
        else:
            # check if limit converged
            f = r.TFile.Open(fname)
            if f!=None:
                t = f.Get("limit")
                if t!=None:
                    # 5 expected + 1 observed (+ prefit sometimes)
                    if t.GetEntries() >= 6:
                        t.SetDirectory(0)
                        # can't update hadded tree because of ROOT bug, so update first
                        if len(args.updateXsec)>0:
                            xobj = r.xsec_t()
                            t.SetBranchAddress("trackedParam_mZprime",r.AddressOf(xobj,'mZprime'))
                            t.SetBranchAddress("trackedParam_xsec",r.AddressOf(xobj,'xsec'))
                            t.SetBranchAddress("limit",r.AddressOf(xobj,'limit'))
                            nt = t.CloneTree(0)
                            nt.SetDirectory(0)
                            for i in range(t.GetEntries()):
                                t.GetEntry(i)
                                xsec_orig = xobj.xsec
                                xobj.xsec = xsecs[str(int(xobj.mZprime))]
                                xobj.limit = xobj.limit*xsec_orig/xobj.xsec
                                nt.Fill()
                            outtrees.append(nt)
                            outtreesroot.Add(nt)
                        else:
                            outtrees.append(t)
                            outtreesroot.Add(t)
                    else:
                        fprint("Warning: limit for {} did not converge".format(fullname))

    # combine outfiles
    if not args.no_hadd:
        if args.dry_run: fprint(outtrees)
        else:
            outname = "limit_dijet{}{}{}.root".format(
                "_"+cname[4:] if len(cname[4:])>0 else "",
                args.updateXsec[1] if len(args.updateXsec)>0 else "",
                getFullname("",args.method,args.finalState),
            )
            outfile = r.TFile.Open(outname,"RECREATE")
            outtree = r.TTree.MergeTrees(outtreesroot)
            outtree.Write()
            outfile.Close()

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--npool", dest="npool", type=int, default=6, help="number of processes")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-f", "--freezeNorm", dest="freezeNorm", default=False, action="store_true", help="freeze bkg normalization to data")
    hadd_group = parser.add_mutually_exclusive_group()
    hadd_group.add_argument("-j", "--just-hadd", dest="just_hadd", default=False, action="store_true", help="don't run any combine commands, just hadd")
    hadd_group.add_argument("--no-hadd", dest="no_hadd", default=False, action="store_true", help="don't hadd")
    parser.add_argument("--hadd-dir", dest="hadd_dir", type=str, default="", help="directory for files to be hadded if not local")
    sig_group = parser.add_mutually_exclusive_group()
    sig_group.add_argument("--signal", dest="signals", metavar=tuple(getParamNames()), type=str, default=[], nargs=4, help="signal parameters")
    sig_group.add_argument("--signals", dest="signals", type=str, default="", help="text file w/ list of signal parameters")
    parser.set_defaults(signals="default_signals.txt")
    parser.add_argument("-N", "--name", dest="name", type=str, default="Test", help="name for combine files")
    parser.add_argument("-a", "--args", dest="args", type=str, default="", help="extra args for combine")
    parser.add_argument("-u", "--update-xsec", dest="updateXsec", type=str, metavar=('filename','suffix'), default=[], nargs=2, help="info to update cross sections when hadding")
    parser.add_argument("-m", "--method", dest="method", type=str, required=True, choices=["fit","ratio","best","best2","best3"], help="dijet background prediction method")
    parser.add_argument("-s", "--final-state", dest="finalState", type=str, nargs='+', choices=["DM","SM"], help="signal final state(s)")
    parser.add_argument("-A", "--acc", dest="acc", type=float, default=0.41, help="dijet SR acceptance")
    parser.add_argument("-X", "--xsec", dest="xsec", type=float, default=None, help="manual (dijet) cross section")
    parser.add_argument("-b", "--bias", dest="bias", default=False, action="store_true", help="perform bias study")
    parser.add_argument("--inject-signal", dest="injectSignal", type=int, default=0, help="signal strength to inject for bias study")
    parser.add_argument("--ntoys", dest="ntoys", type=int, default=100, help="number of toys for bias study")
    parser.add_argument("--rmax", dest="rmax", type=float, default=80, help="max signal strength allowed for bias study")
    parser.add_argument("--seed", dest="seed", type=int, default=123456, help="seed for toy generation in bias study")
    args = parser.parse_args()

    if args.bias:
        args.just_hadd = False
        args.no_hadd = True

    if not args.just_hadd and args.method.startswith("best"):
        parser.error('Can only use "best" method with --just-hadd')

    # parse signal info
    xsecs = getXsecs('dict_xsec_Zprime.txt')
    param_names = getParamNames()+["xsec"]
    param_values = []
    if isinstance(args.signals,list):
        param_values.append(args.signals)
        param_values[-1].append(xsecs[param_values[-1][0]])
    else:
        with open(args.signals,'r') as sfile:
            for line in sfile:
                line = line.rstrip()
                if len(line)==0: continue
                param_values.append(line.split())
                param_values[-1].append(xsecs[param_values[-1][0]])
    args.signals = [makeSigDict(param_values[i],param_names) for i in range(len(param_values))]

    main(args)

