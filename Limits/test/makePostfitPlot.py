import os,sys,shutil
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from glob import glob
from collections import OrderedDict
from array import array
from paramUtils import getParamsTracked, getFname, makeSigDict, getSigname, getSignameShort, getCombos, runCmds, getChannel, fprint, getPname, getRname, getInitFromBF, OpenFile
from Stat.Limits.bruteForce import silence

input_template = """INPUT
input/input_svj_stack_dijetmtdetahad_2017.txt
input/input_svj_mt_postfit_options.txt
input/input_svj_mt_hist_full{}.txt
"""

ofile_prefix = "test/fit"
bandname = "errorband"

options_template = """OPTION
string+:printsuffix[{psuff}]
vstring:extra_text[{etxt}]
vstring:fits[{fitlist}]
vstring+:chosensets[{signames}]
vstring+:numers[{signames}]
string:rootfile[{ofile}]
bool:treesuffix[0]
string:prelim_text[{prelim}]
"""

fit_template = "{fitname}\ts:fn[{fnname}]\tvd:pars[1,{pvals}]\td:yield[{yieldval}]\ts:legname[{legname}]\tin:input[input/input_svj_mt_fit_opt.txt]\tb:legpars[0]\tc:linecolor[{fitcol}]"
band_template = "\ts:bandfile[test/{signame}/{bandfile}]\ts:bandname[{bandname}]\tc:glinecolor[{fitcol}]\tc:gfillcolor[kGray + 1]"

set_template = """hist\tmc\t{signamefull}\ts:legname[{legname}]\tc:color[{sigcol}]\ti:linestyle[7]\ti:panel[0]\tvs:fits[]\t{signorm}
\tbase\text\t{signamefull}\ts:extfilename[{sigfile}]\ts:exthisto_dir[{hdir}]\tvs:exthisto_in[{signamesafe}]\tvs:exthisto_out[MTAK8]"""

data_template = """hist\tdata\tdata\ti:panel[1]\ts:legname[{dtype} data{inj}]\tb:yieldnorm[0]\tb:poisson[1]
\tbase\text\tdata\ts:extfilename[{dfile}]\ts:exthisto_dir[{hdir}]\tvs:exthisto_in[data_{dtype}]\tvs:exthisto_out[MTAK8]"""

quantile_info = {
    -4: {"legname": "asimov", "name": "asimov", "color": "kYellow + 3", "sigcolor": "kCyan + 2"},
    -3: {"legname": "bestfit", "name": "bestfit", "color": "kOrange + 2", "sigcolor": "kCyan + 2"},
    -2: {"legname": "b-only", "name": "bonly", "color": "kRed", "sigcolor": "kCyan + 2"},
    -1: {"legname": "postfit (obs)", "name": "postfitobs", "color": "kBlue", "sigcolor": "kCyan + 2"},
    0.5: {"legname": "median exp.", "name": "medexp", "color": "kMagenta + 2", "sigcolor": "kCyan + 2"},
}

function_info = {
    ("alt",2):  {"formula": "[0]*exp([1]*x/13000)*(x/13000)^([2])", "legname": "g_{2}(x)", "name": "g"},
    ("alt",3):  {"formula": "[0]*exp([1]*x/13000)*(x/13000)^([2]*(1+[3]*log(x/13000)))", "legname": "g_{3}(x)", "name": "g"},
    ("main",2): {"formula": "([0]*(1-x/13000)^[1])*((x/13000)^(-[2]))", "legname": "f_{1,1}(x)", "name": "f"},
    ("main",3): {"formula": "([0]*(1-x/13000)^[1])*((x/13000)^(-[2]-[3]*log(x/13000)))", "legname": "f_{1,2}(x)", "name": "f"},
}

region_info = {
    "highCut": {"alt": 3, "main": 3, "legname": "high-^{}R_{T}"},
    "lowCut": {"alt": 2, "main": 2, "legname": "low-^{}R_{T}"},
    "highSVJ2": {"alt": 2, "main": 2, "legname": "high-SVJ2"},
    "lowSVJ2": {"alt": 2, "main": 2, "legname": "low-SVJ2"},
}

def makeBandFileName(iname):
    return iname.replace("input_","errorband_",1).replace(".txt",".root")

def makeErrorBand(iname,ws,fitres,region,ftype,norm):
    import ROOT as r
    silence()

    pdfname = "Bkg{}_{}_2018".format("_Alt" if ftype=="alt" else "",region)
    pdf = ws.pdf(pdfname)
    var = ws.var("mH{}_2018".format(region))
    data = ws.data("data_obs")
    # combine output is always RooDataSet w/ both regions (not RooDataHist)
    if type(data)==r.RooDataSet:
        ch = getChannel(region)
        data = data.reduce("CMS_channel==CMS_channel::{}".format(ch))

    # borrowed from ftest.py
    frame = var.frame(r.RooFit.Title(""))
    data.plotOn(frame)
    pdf.plotOn(frame, r.RooFit.VisualizeError(fitres, 1, False), r.RooFit.Normalization(norm, r.RooAbsReal.NumEvent))
    # find the error band
    band = None
    for i in range(int(frame.numItems())):
        bname = frame.nameOf(i)
        if "errorband" in bname:
            band = frame.getCurve(bname)
            break

    # get rid of unwanted end points
    pruned = [(band.GetX()[j],band.GetY()[j]) for j in range(band.GetN()) if band.GetY()[j]>0]
    bx, by = zip(*pruned)
    bx = array('d',bx)
    by = array('d',by)
    gband = r.TGraph(len(bx),bx,by)

    oname = makeBandFileName(iname)
    ofile = r.TFile.Open(oname,"RECREATE")
    ofile.cd()
    gband.Write(bandname)
    ofile.Close()
    return oname

def actuallyPlot(signame,input,postfname,data_file):
    bandfname = makeBandFileName(input)
    current_dir = os.getcwd()
    analysis_dir = os.path.expandvars("$CMSSW_BASE/src/Analysis/")
    # check for data file in main dir
    data_xrd = data_file.startswith("root://")
    if not data_xrd and not os.path.isfile(data_file):
        data_file = "../"+data_file
    # get paths before changing directories
    input_path = os.path.abspath(input)
    bandfname_path = os.path.abspath(bandfname)
    postfname_path = os.path.abspath(postfname) if postfname is not None else ""
    data_file_path = os.path.abspath(data_file) if not data_xrd else data_file
    # link files to Analysis folders, run root command
    os.chdir(analysis_dir)
    commands = [
        "ln -sf {} {}/input".format(input_path,analysis_dir),
        "mkdir -p {}/test/{}".format(analysis_dir,signame),
        "ln -sf {} {}/test/{}".format(bandfname_path,analysis_dir,signame)
    ]
    if postfname is not None:
        commands.append(
            "ln -sf {} {}/test/{}".format(postfname_path,analysis_dir,signame)
        )
    if not data_xrd:
        commands.append(
            "ln -sf {} {}/test".format(data_file_path,analysis_dir)
        )
    commands.append(
        """root -b -l -q 'KPlotDriver.C+(".",{{"{}"}},{{}},1)'""".format("input/{}".format(input))
    )
    runCmds(commands)
    outputs = glob("*.png") + glob("*.pdf") + glob("{}*.root".format(ofile_prefix))
    for output in outputs:
        shutil.move(output,os.path.join(current_dir,os.path.basename(output)))
    fprint(' '.join(outputs))
    os.chdir(current_dir)

def makePostfitPlot(sig, name, method, quantile, data_file, datacard_dir, obs, injected, combo, region, init_dir=None, init_ftype=None, do_plot=False):
    ch = getChannel(region)
    seedname = None
    if not obs: seedname = data_file.split('.')[-2]

    signame = getSignameShort(sig)
    signamesafe = getSigname(sig)
    dtype = "obs" if obs else "toy"
    iname = "input_svj_mt_fit_{dtype}_{region}_{name}_{qname}_{sig}.txt".format(dtype=dtype,region=region,name=name,qname=quantile_info[quantile]["name"],sig=signamesafe)
    rinfo = region_info[region]
    ftype = ""
    finfo = None

    fits = OrderedDict()
    sigs = OrderedDict()
    for q in [quantile]:
        qinfo = quantile_info[q]
        import ROOT as r
        r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

        # using initial fits from outside combine
        if init_dir is not None:
            # input sources
            postfname = None
            fname = "{}/ws_{}.root".format(init_dir,region)
            rname = "{}/fitResults_{}.root".format(init_dir,region)
            wsname = "BackgroundWS"
            ftype = init_ftype

            # items to obtain: params, norm, workspace, RooFitResult
            params = getInitFromBF(fname, wsname, getPname(region, ftype=="alt"))
            params = {item.split("=")[0]:float(item.split("=")[1]) for item in params}

            ftmp = OpenFile(fname)
            ws = ftmp.Get(wsname)
            norm = ws.data("data_obs").sumEntries()
            rtmp = OpenFile(rname)
            fitres = rtmp.Get(getRname(region, ftype=="alt", len(params)))
        else:
            # input sources
            qfitname = "Postfit{:.3f}".format(q)
            fname = getFname(name, method, combo, sig=sig, seed=seedname)
            indfname = getFname(qfitname+name, "MultiDimFit", combo, sig=sig, seed=seedname)
            postfname = getFname(qfitname+name, "", combo, prefix="multidimfit", sig=sig, seed=seedname)

            # items to obtain: params, norm, workspace, RooFitResult
            params = getParamsTracked(fname, quantile)
            if len(params)==0: raise RuntimeError("Could not get tracked parameters for quantile {} from file: {}".format(quantile, fname))
            ftype = "alt" if any(region in p and "_alt" in p for p in params) else "main"
            norm = params["trackedParam_n_exp_final_bin{}_proc_roomultipdf".format(ch)]

            iftmp = OpenFile(indfname)
            ws = iftmp.Get("w")
            pftmp = OpenFile(postfname)
            fitres = pftmp.Get("fit_mdf")

        # get error band
        bandfname = None
        try:
            bandfname = makeErrorBand(iname,ws,fitres,region,ftype,norm)
        except TypeError:
            # in case of: "[#0] ERROR:Eval -- RooFitResult::createHessePdf() ERROR: covariance matrix is not positive definite (|V|=0) cannot construct p.d.f"
            pass

        # common stuff
        pfit = {p:v for p,v in params.iteritems() if region in p}
        pvals = [str(v) for p,v in sorted(pfit.iteritems())]
        finfo = function_info[(ftype,rinfo[ftype])]
        fitname = "{}_{}".format(finfo["name"],qinfo["name"])

        fits[fitname] = fit_template.format(
            fitname = fitname,
            fnname = finfo["formula"],
            pvals = ','.join(pvals),
            # yieldval must be multiplied by bin width
            yieldval = str(norm*100),
            legname = "{}, {}".format(finfo["legname"],qinfo["legname"]),
            fitcol = qinfo["color"],
        )
        if bandfname is not None:
            fits[fitname] += band_template.format(
                fitcol = qinfo["color"],
                bandfile = bandfname,
                bandname = bandname,
                signame = signamesafe,
            )

        # no need to show signal for b-only
#        if qinfo["name"]=="bonly": continue

        signamefull = "{}_{}".format(signame,qinfo["name"])
        sigfileorig = "{}/datacard_final_{}.root".format(datacard_dir,signame)
        hdir = "{}_2018".format(region)
        if qinfo["name"]=="bonly":
            legname = "{} (r = {:.2g})".format(signame,1)
            sigfile = sigfileorig
            hdirsig = hdir
            signorm = "b:yieldnorm[0]"
        else:
            legname = "{} (r = {:.2g})".format(signame,params["r"])
            sigfile = "test/{}".format(postfname if signamesafe in postfname else signamesafe+"/"+postfname)
            hdirsig = "shapes_fit/{}".format(ch)
            signorm = "d:yieldnormval[{}]".format(params["trackedParam_n_exp_final_bin{}_proc_{}".format(ch,signamesafe)])

        sigs[signamefull] = set_template.format(
            signamefull = signamefull,
            legname = legname,
            sigcol = qinfo["sigcolor"],
            signorm = signorm,
            sigfile = sigfile,
            signame = signame,
            signamesafe = signamesafe,
            hdir = hdirsig,
        )

    data = data_template.format(
        dtype = dtype,
        dfile = sigfileorig if obs else "test/"+data_file if do_plot else data_file,
        hdir = hdir,
        inj = "" if obs else " (no signal)" if injected==0 else " (^{{}}m_{{Z'}} = {:.1f} TeV)".format(float(injected)/1000.)
    )

    options = options_template.format(
        # should quantiles be included here?
        psuff = "_{}_fit_{}_{}_{}_{}".format(dtype,region,name,quantile_info[quantile]["name"],signamesafe),
        etxt = rinfo["legname"],
        fitlist = ','.join(fits),
        signames = ','.join(sigs),
        ofile = "{}_{}_{}_{}_{}_{}".format(ofile_prefix,dtype,signamesafe,region,name,quantile_info[quantile]["name"]),
        prelim = "" if obs else "Simulation",
    )

    with open(iname,'w') as ifile:
        lines = [
            input_template.format("_bdt" if combo=="bdt" else ""),
            options,
            "FIT",
            '\n'.join(fits.values()),
            '',
            "SET",
            '\n'.join(sigs.values()),
            data,
        ]
        ifile.write('\n'.join(lines))

    if do_plot: actuallyPlot(signamesafe,iname,postfname,sigfileorig if obs else data_file)

    return iname

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-S","--signal", dest="signal", metavar=("mZprime","mDark","rinv","alpha"), type=str, required=True, nargs=4, help="signal parameters")
    parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[name])")
    parser.add_argument("-M", "--method", dest="method", type=str, default="AsymptoticLimits", help="method name (higgsCombineTest.[method])")
    parser.add_argument("-d", "--data", dest="data", type=str, default="", help="data file name (taken from datacards if obs)")
    parser.add_argument("-o", "--obs", dest="obs", default=False, action="store_true", help="using observed rather than toy data")
    parser.add_argument("-i", "--injected", dest="injected", type=int, default=0, help="injected Zprime mass")
    parser.add_argument("-D", "--datacards", dest="datacards", type=str, default="root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig8/sigfull/", help="datacard histogram location (for prefit)")
    parser.add_argument("-q", "--quantile", dest="quantile", type=float, default=-1, choices=list(quantile_info), help="quantile to plot fits")
    parser.add_argument("-c", "--combo", dest="combo", type=str, required=True, choices=sorted(list(getCombos())), help="combo to plot")
    parser.add_argument("-I", "--use-init", dest="init", type=str, metavar=("dir","fit"), default=[], nargs=2, help="directory from which to use initial fits (outside combine) and fit type (alt or main)")
    parser.add_argument("-p", "--plot", dest="plot", default=False, action="store_true", help="actually make plot(s)")
    args = parser.parse_args()

    if not args.obs and len(args.data)==0:
        parser.error("Data file must be specified if using toy data")
    if len(args.init)>0:
        if not args.obs:
            parser.error("Initial fits can only be used for observed data")
        elif args.quantile!=-2:
            parser.error("Initial fits can only be used for prefit plot")
        elif args.init[1] not in ["alt","main"]:
            parser.error("Unknown fit type {}".format(args.init[1]))
    else: args.init = [None,None]

    args.signal = makeSigDict(args.signal)

    combos = getCombos()

    input_files = []
    for region in combos[args.combo]:
        tmp = makePostfitPlot(args.signal,args.name,args.method,args.quantile,args.data,args.datacards,args.obs,args.injected,args.combo,region,args.init[0],args.init[1],args.plot)
        input_files.append(tmp)
    print ' '.join(input_files)

