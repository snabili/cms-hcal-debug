import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f = r.TFile(infile)
tps = f.Get("chainplotter/matches")

c = r.TCanvas()
c.SetLogz(True)
c.SetRightMargin(c.GetRightMargin() * 1.5)
c.SaveAs(outfile + '[')
tps.Draw("RH_energy:TP_energy>>cmphb(100, 0, 200, 100, 0, 200)", "abs(ieta) <= 16", "COLZ")
r.gDirectory.Get("cmphb").SetTitle("HB energy comparison;TP E_{T};RH E_{T}")
c.SaveAs(outfile)
tps.Draw("RH_energy:TP_energy>>cmphe(100, 0, 100, 100, 0, 100)", "abs(ieta) > 16 && abs(ieta) < 29", "COLZ")
r.gDirectory.Get("cmphe").SetTitle("HE energy comparison;TP E_{T};RH E_{T}")
c.SaveAs(outfile)
tps.Draw("RH_energy:TP_energy/0.7>>cmphf(80, 0, 80, 80, 0, 80)", "abs(ieta) > 29", "COLZ")
r.gDirectory.Get("cmphf").SetTitle("HF energy comparison;TP E_{T};RH E_{T} / 0.7")
c.SaveAs(outfile)
c.SaveAs(outfile + ']')
