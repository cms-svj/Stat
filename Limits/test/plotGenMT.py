import os
import sys

ifname = ""
if len(sys.argv)>1:
    ifname = sys.argv[1]

import ROOT as r
r.gROOT.SetBatch(True)

def signame(width):
    template = "step1_GEN_s-channel_mMed-3100_mDark-20_rinv-0.3_alpha-peak_{}13TeV-pythia8_n-50000_part-1.root"
    return template.format("" if width is None else "width-{}_".format(width))

dir = "root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/width2/GEN/"
mZprime = 3100
widths = [None, 0.01, 0.05, 0.10, 0.15, 0.20, 0.30]
colors = [r.kBlack, r.kBlue, r.kMagenta, r.kOrange+2, r.kCyan+1, r.kYellow+3, r.kRed]

r.gInterpreter.ProcessLine('#include "{}/src/Analysis/KCode/KMath.h"'.format(os.environ["CMSSW_BASE"]))

hists = {
    "mt": [],
    "mzp": [],
}
if len(ifname)==0:
    ofile = r.TFile.Open("hists_genmt.root","RECREATE")
    for iw,width in enumerate(widths):
        sig = signame(width)
        file = r.TFile.Open(dir+sig)
        tree = file.Get("Events")
        tree.SetLineColor(colors[iw])

        tree.SetAlias("GenJetsAK8","recoGenJets_ak8GenJetsNoNu__GEN.obj")
        tree.SetAlias("GenMET","recoGenMETs_genMetTrue__GEN.obj")
        tree.SetAlias("GenMT","KMath::TransverseMass(GenJetsAK8[0].px()+GenJetsAK8[1].px(),GenJetsAK8[0].py()+GenJetsAK8[1].py(),sqrt((GenJetsAK8[0].energy()+GenJetsAK8[1].energy())^2-(GenJetsAK8[0].px()+GenJetsAK8[1].px())^2-(GenJetsAK8[0].py()+GenJetsAK8[1].py())^2-(GenJetsAK8[0].pz()+GenJetsAK8[1].pz())^2),GenMET.px(),GenMET.py(),0)")
        histname = "h{}".format(iw)
        tree.Draw("GenMT>>{}(40,0,4000)".format(histname),"","goff")
        hist = r.gDirectory.Get(histname)
        hist.SetDirectory(0)
        hist.Scale(1./hist.Integral(-1,-1))
        hist.SetTitle("")
        hist.GetXaxis().SetTitle("m_{T}^{gen} [GeV]")
        hist.GetYaxis().SetTitle("Arbitrary units")
        hists["mt"].append(hist)

        tree.SetAlias("GenParticles","recoGenParticles_genParticles__GEN.obj")
        histname2 = "p{}".format(iw)
        tree.Draw("GenParticles.mass()>>{}(200,2000,4000)".format(histname2),"GenParticles.pdgId()==4900023","goff")
        hist2 = r.gDirectory.Get(histname2)
        hist2.SetDirectory(0)
        hist2.Scale(1./hist2.Integral(-1,-1))
        hist2.SetTitle("")
        hist2.GetXaxis().SetTitle("m_{Z'} [GeV]")
        hist2.GetYaxis().SetTitle("Arbitrary units")
        if iw>0: hists["mzp"].append(hist2)

        ofile.cd()
        hist.Write()
        if iw>0: hist2.Write()
else:
    ifile = r.TFile.Open(ifname)
    for iw,width in enumerate(widths):
        hist = ifile.Get("h{}".format(iw))
        hist.SetDirectory(0)
        hists["mt"].append(hist)

        if iw>0:
            hist2 = ifile.Get("p{}".format(iw))
            hist2.SetDirectory(0)
            hists["mzp"].append(hist2)

for htype in hists:
    offset = 1 if htype=="mzp" else 0
    hlist = hists[htype]
    can = r.TCanvas()
    can.SetLogy()
    leg = r.TLegend(0.65,0.5,0.95,0.9)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6);
    leg.SetBorderSize(0)
    leg.SetTextSize(0.037)
    leg.SetTextFont(42)
    haxis = hlist[0].Clone("haxis_{}".format(htype))
    haxis.Reset()
    haxis.GetXaxis().SetLimits(hlist[0].GetXaxis().GetXmin(),hlist[0].GetXaxis().GetXmax()+1000)
    ymin = min([hist.GetMinimum(0) for hist in hlist])
    ymax = max([hist.GetMaximum() for hist in hlist])*1.1
    haxis.GetYaxis().SetRangeUser(ymin,ymax)
    haxis.Draw("hist")
    for ih,hist in enumerate(hlist):
        iho = ih+offset
        leg.AddEntry(hist,"Z' width = {} GeV".format(widths[iho]*mZprime if iho>0 else "0.01"),"l")
        #hist.Draw("hist"+(" same" if ih>0 else ""))
        hist.Draw("hist same")
    leg.Draw()
    can.Print("{}_width.png".format(htype))

