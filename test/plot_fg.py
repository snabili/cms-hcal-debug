import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f = r.TFile(infile)
tps = f.Get("compareL1T/tps")

c = r.TCanvas()
c.SetRightMargin(c.GetRightMargin() * 1.5)
c.SaveAs(outfile + '[')
tps.Draw("iphi:ieta>>histdet(85,-42,42,73,0,72)", "(version==1) * (fg1-fg1_emul)", "COLZ")
r.gDirectory.Get("histdet").SetTitle(";ieta;iphi;#sum L1T FG - FG reemulated")
c.SaveAs(outfile)
tps.Draw("fg1-fg1_emul:soi>>histsoi", "version==1", "COLZ")
r.gDirectory.Get("histsoi").SetTitle(";SOI;L1T FG - FG reemulated;")
c.SetLogz()
c.SaveAs(outfile)
c.SaveAs(outfile + ']')
