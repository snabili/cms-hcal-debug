import ROOT as r

r.gROOT.SetBatch()

f = r.TFile("legacy.root")
tree = f.Get("chainplotter/matches")

hist = r.TH2F("hist", "", 40, 0, 40, 40, 0, 40)
tree.Draw("RH_energy:TP_energy>>hist")

diag = r.TF1("diag", "x", 0, 40)
diag.SetLineColor(r.kBlack)

c = r.TCanvas()
hist.Draw("COLZ")
diag.Draw("same")
c.SetLogz()
c.SaveAs("legacy.pdf")
