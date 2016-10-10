import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f1 = r.TF1("f1", "x", 0, 200)
f1.SetLineColor(r.kRed)

f = r.TFile(infile)
tps = f.Get("chainplotter/matches")

r.gStyle.SetTitleFontSize(.06)
r.gStyle.SetTitleXSize(.05)
r.gStyle.SetTitleYSize(.05)
r.gStyle.SetPadBottomMargin(.13)
r.gStyle.SetPadLeftMargin(.12)
r.gStyle.SetPadRightMargin(.09)
r.gStyle.SetPadTopMargin(.1)

c = r.TCanvas()
c.SetLogz(True)
c.SetRightMargin(c.GetRightMargin() * 1.5)
c.SaveAs(outfile + '[')
tps.Draw("RH_energy:TP_energy>>cmphb(100, 0, 200, 100, 0, 200)", "abs(ieta) <= 16", "COLZ")
f1.Draw("same")
r.gDirectory.Get("cmphb").SetTitle("HB energy comparison;TrigPrim E_{T};RecHit E_{T}")
c.SaveAs(outfile)
tps.Draw("RH_energy:TP_energy>>cmphe(100, 0, 100, 100, 0, 100)", "abs(ieta) > 16 && abs(ieta) < 29", "COLZ")
f1.Draw("same")
r.gDirectory.Get("cmphe").SetTitle("HE energy comparison;TrigPrim E_{T};RecHit E_{T}")
c.SaveAs(outfile)
tps.Draw("RH_energy:TP_energy/0.7>>cmphf(80, 0, 80, 80, 0, 80)", "abs(ieta) > 29", "COLZ")
f1.Draw("same")
r.gDirectory.Get("cmphf").SetTitle("HF energy comparison;TrigPrim E_{T} / 0.7;RecHit E_{T}")
c.SaveAs(outfile)
c.SaveAs(outfile + ']')
