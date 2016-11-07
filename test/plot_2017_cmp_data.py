import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f1 = r.TF1("f1", "x", 0, 200)
f1.SetLineColor(r.kRed)

f2 = r.TF1("f2", "2*x", 0, 200)
f2.SetLineColor(r.kGreen)

f3 = r.TF1("f3", "3*x", 0, 200)
f3.SetLineColor(r.kBlue)

f = r.TFile(infile)
tps = f.Get("compare/tps")

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
for ieta in range(30, 42):
    if ieta == 30:
        binning = "48,0,24,96,0,48"
    elif ieta < 33:
        binning = "32,0,16,64,0,32"
    else:
        binning = "16,0,8,32,0,16"
    # When taking out the division by valid fiber count
    binning = ','.join(binning.split(',')[3:] * 2)
    tps.Draw("et:et_emul>>hist({})".format(binning), "iphi == 39 && ieta > 0 && abs(ieta) == {}".format(ieta), "COLZ")
    f1.Draw("same")
    f2.Draw("same")
    f3.Draw("same")
    r.gDirectory.Get("hist").SetTitle("Energy comparison for i#eta {};Emul TrigPrim E_{{T}};RAW TrigPrim E_{{T}}".format(ieta))
    c.SaveAs(outfile)
c.SaveAs(outfile + ']')
