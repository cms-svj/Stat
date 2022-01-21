import os,sys,shlex,traceback
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool
from StringIO import StringIO
from array import array
from string import Template
from paramUtils import fprint, runCmd, alphaVal, makeSigDict, getParamNames, getSignameShort, OpenFile, paramVal
from Stat.Limits.bruteForce import silence

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

def transform(hname,iname,wname,w2name,template,signame,extargs,dryrun):
    tname = template.format("template")
    ofname = template.format(signame).replace(".txt",".root")
    dcfname = template.format(signame)

    if dryrun: return dcfname
    
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    hfile = OpenFile(hname)
    ifile = OpenFile(iname)
    w = ifile.Get(wname)

    th1x = w.var("th1x")
    bins = th1x.getBinning().array()
    bins.SetSize(th1x.getBinning().numBins()+1)
    bins = array('d',bins)

    # make new workspace, import contents from prev (but take new data)
    w2 = r.RooWorkspace(w2name)
    getattr(w2,"import")(w.allVars())
    getattr(w2,"import")(w.allPdfs())
    getattr(w2,"import")(w.data("data_obs"))
    
    # get sig histo
    systs = [
        "",
        "_jesUp",
        "_jesDown",
        "_jerUp",
        "_jerDown",
    ]
    histname = "widejetmass_"+signame+"_cmsdijet"
    hyield = 0
    for syst in systs:
        hist = hfile.Get(histname+syst)
        if not isinstance(hist,r.TH1): raise RuntimeError("Could not get {} from {}".format(histname+syst,hname))
        hist.SetDirectory(0)
        if syst=="": hyield = hist.Integral()
    
        # reset sig histo bins
        hist.GetXaxis().Set(len(bins)-1,bins)
    
        # import new pdf
        hname2 = signame+syst
        pdf = r.RooDataHist(hname2, hname2, r.RooArgList(th1x), hist, 1.);
        getattr(w2,"import")(pdf)

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

def doLimit(info):
    args = info["args"]
    sig = info["sig"]
    signame = getSignameShort(sig)

    params = {key:paramVal(key,val) for key,val in sig.iteritems()}
    setargs = []
    trkargs = []
    extargs = ""
    for p,v in params.iteritems():
        setargs.append(p+"="+str(v))
        trkargs.append(p)
        extargs = extargs+p+" extArg "+str(v)+"\n"
    frzargs = trkargs[:]

    fitparams = ["p1_PFDijet2017","p2_PFDijet2017","p3_PFDijet2017","shapeBkg_PFDijet2017_bkg_PFDijet2017__norm"]
    trkargs.extend(fitparams)
    if args.freezeNorm: frzargs.append("shapeBkg_PFDijet2017_bkg_PFDijet2017__norm")

    # datacard setup
    dcfname = transform("svj_dijet.root","dijet_combine_qq_1900_lumi-137.500_PFDijet2017.root","wPFDijet2017","wPFDijet2018","dijet_combine_{}.txt",signame,extargs,args.dry_run)

    cargs = args.args
    if len(cargs)>0: cargs += " "
    cargs += "--setParameters {} --freezeParameters {} --trackParameters {} --keyword-value sig={} -n {} -d {}".format(
        ','.join(setargs), ','.join(frzargs), ','.join(trkargs), signame, args.cname, dcfname
    )

    # run combine
    outputs = []
    fprint("Calculating limit for {}...".format(signame))
    command = "combine -M AsymptoticLimits {}".format(cargs)
    outputs.append(command)
    if not args.dry_run: outputs.append(runCmd(command))

    return outputs

def main(args):
    cname = args.name[:]
    if args.freezeNorm: cname += "Frz"
    args.cname = cname

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
    for sig in args.signals:
        signame = getSignameShort(sig)
        fname = "higgsCombine{}.AsymptoticLimits.mH120.sig{}.root".format(cname,signame)
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
                        fprint("Warning: limit for {} did not converge".format(signame))

    # combine outfiles
    if not args.no_hadd:
        if args.dry_run: fprint(outtrees)
        else:
            outname = "limit_dijet{}{}.root".format("_"+cname[4:] if len(cname[4:])>0 else "",args.updateXsec[1] if len(args.updateXsec)>0 else "")
            outfile = r.TFile.Open(outname,"RECREATE")
            outtree = r.TTree.MergeTrees(outtreesroot)
            outtree.Write()
            outfile.Close()

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--npool", dest="npool", type=int, default=6, help="number of processes")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-f", "--freezeNorm", dest="freezeNorm", default=False, action="store_true", help="freeze bkg normalization to data")
    parser.add_argument("-j", "--just-hadd", dest="just_hadd", default=False, action="store_true", help="don't run any combine commands, just hadd")
    parser.add_argument("--no-hadd", dest="no_hadd", default=False, action="store_true", help="don't hadd")
    parser.add_argument("--hadd-dir", dest="hadd_dir", type=str, default="", help="directory for files to be hadded if not local")
    sig_group = parser.add_mutually_exclusive_group()
    sig_group.add_argument("--signal", dest="signals", metavar=tuple(getParamNames()), type=str, default=[], nargs=4, help="signal parameters")
    sig_group.add_argument("--signals", dest="signals", type=str, default="", help="text file w/ list of signal parameters")
    parser.set_defaults(signals="default_signals.txt")
    parser.add_argument("-N", "--name", dest="name", type=str, default="Test", help="name for combine files")
    parser.add_argument("-a", "--args", dest="args", type=str, default="", help="extra args for combine")
    parser.add_argument("-u", "--update-xsec", dest="updateXsec", type=str, metavar=('filename','suffix'), default=[], nargs=2, help="info to update cross sections when hadding")
    args = parser.parse_args()

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

