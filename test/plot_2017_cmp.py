import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f = r.TFile(infile)
tps = f.Get("chainplotter/matches")

c = r.TCanvas()
c.SetRightMargin(c.GetRightMargin() * 1.5)
c.SaveAs(outfile + '[')
tps.Draw("RH_energy:TP_energy>>cmphb", "abs(ieta) <= 16", "COLZ")
r.gDirectory.Get("cmphb").SetTitle("HB energy comparison;TP E_{T};RH E_{T}")
c.SaveAs(outfile)
tps.Draw("RH_energy:TP_energy>>cmphe", "abs(ieta) > 16 && abs(ieta) < 29", "COLZ")
r.gDirectory.Get("cmphe").SetTitle("HE energy comparison;TP E_{T};RH E_{T}")
c.SaveAs(outfile)
tps.Draw("RH_energy:TP_energy>>cmphf", "abs(ieta) > 29", "COLZ")
r.gDirectory.Get("cmphf").SetTitle("HF energy comparison;TP E_{T};RH E_{T}")
c.SaveAs(outfile)
c.SaveAs(outfile + ']')
