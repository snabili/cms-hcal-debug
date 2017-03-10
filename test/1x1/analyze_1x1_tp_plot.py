import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f = r.TFile(infile)

e1 = f.Get("analyze/evs")
e2 = f.Get("analyzeOld/evs")
e1.AddFriend(e2, "old")

t1 = f.Get("analyze/tps")
t2 = f.Get("analyzeOld/tps")

c = r.TCanvas()
c.SaveAs(outfile + '[')
e1.Draw("tp_v1_et:old.tp_v0_et>>hist", "", "COLZ")
r.gDirectory.Get("hist").SetTitle("HF 1x1 TP vs 2x3 TP (from RAW);#sum E_{T} HF, v0;#sum E_{T} HF, v1")
c.SaveAs(outfile)
e1.Draw("tp_v1_et:tp_v0_et>>hist", "", "COLZ")
r.gDirectory.Get("hist").SetTitle("HF 1x1 TP vs 2x3 TP (reemulated);#sum E_{T} HF, v0;#sum E_{T} HF, v1")
c.SaveAs(outfile)
t1.Draw("ieta>>hist2", "et")
r.gDirectory.Get("hist2").SetTitle("HF 1x1 TP;ieta;#sum E_{T}")
c.SaveAs(outfile)
t2.Draw("ieta>>hist3", "et")
r.gDirectory.Get("hist3").SetTitle("HF 2x3 TP;ieta;#sum E_{T}")
c.SaveAs(outfile)
c.SaveAs(outfile + ']')
