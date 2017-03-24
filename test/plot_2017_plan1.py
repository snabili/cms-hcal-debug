import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outstub = sys.argv[1:]

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

hep17 = "(ieta<=29 && ieta>=17 && iphi>=63 && iphi<=66)"
nhep17 = "!(ieta<=29 && ieta>=17 && iphi>=63 && iphi<=66) * (abs(ieta)<=29 && abs(ieta)>=17)"
hb = "(abs(ieta)<=16)"
hf = "(abs(ieta)>29)"

for m in "M0 M2 M3".split():
    rh_en = "RH_energy" + m
    outfile = "{}{}.pdf".format(outstub, m)

    c = r.TCanvas()
    c.SetLogz(True)
    c.SetRightMargin(c.GetRightMargin() * 1.5)
    c.SaveAs(outfile + '[')

    r.TH1.SetDefaultSumw2()
    tps.Draw(rh_en + ">>cmphep17r(10,0,100)", rh_en + " * (TP_energy>0) * " + hep17, "E")
    tps.Draw(rh_en + ">>cmphep17t(10,0,100)", "TP_energy * (TP_energy>0) * " + hep17, "E")
    h17 = r.gDirectory.Get("cmphep17r")
    h17.Divide(r.gDirectory.Get("cmphep17t"))
    h17.SetTitle("HEP17 energy ratio;RecHit E_{T};TrigPrim E_{T} / RecHit E_{T}")
    h17.Draw("E")
    c.SaveAs(outfile)

    tps.Draw(rh_en + ">>cmpnhep17r(10,0,100)", rh_en + " * (TP_energy>0) * " + nhep17, "E")
    tps.Draw(rh_en + ">>cmpnhep17t(10,0,100)", "TP_energy * (TP_energy>0) * " + nhep17, "E")
    h17 = r.gDirectory.Get("cmpnhep17r")
    h17.Divide(r.gDirectory.Get("cmpnhep17t"))
    h17.SetTitle("HE w/o HEP17 energy ratio;RecHit E_{T};TrigPrim E_{T} / RecHit E_{T}")
    h17.Draw("E")
    c.SaveAs(outfile)

    tps.Draw(rh_en + ">>cmphbr(10,0,100)", rh_en + " * (TP_energy>0) * " + hb, "E")
    tps.Draw(rh_en + ">>cmphbt(10,0,100)", "TP_energy * (TP_energy>0) * " + hb, "E")
    h17 = r.gDirectory.Get("cmphbr")
    h17.Divide(r.gDirectory.Get("cmphbt"))
    h17.SetTitle("HB energy ratio;RecHit E_{T};TrigPrim E_{T} / RecHit E_{T}")
    h17.Draw("E")
    c.SaveAs(outfile)

    tps.Draw(rh_en + ">>cmphfr(10,0,100)", rh_en + " * (TP_energy>0) * " + hf, "E")
    tps.Draw(rh_en + ">>cmphft(10,0,100)", "TP_energy * (TP_energy>0) * " + hf, "E")
    h17 = r.gDirectory.Get("cmphfr")
    h17.Divide(r.gDirectory.Get("cmphft"))
    h17.SetTitle("HF energy ratio;RecHit E_{T};TrigPrim E_{T} / RecHit E_{T}")
    h17.Draw("E")
    c.SaveAs(outfile)

    tps.Draw("(TP_energy-{}):ieta>>diffieta".format(rh_en), "TP_energy>0", "COLZ")
    r.gDirectory.Get("diffieta").SetTitle("Energy Difference RH vs TP;ieta;TrigPrim E_{T} - RecHit E_{T}")
    c.SaveAs(outfile)

    tps.Draw("(TP_energy-{}):ieta>>diffieta2".format(rh_en), "(TP_energy>0) * (abs(TP_energy-{})<10)".format(rh_en), "COLZ")
    r.gDirectory.Get("diffieta2").SetTitle("Energy Difference RH vs TP;ieta;TrigPrim E_{T} - RecHit E_{T}")
    c.SaveAs(outfile)

    for cut in "1 2 5".split():
        tps.Draw(rh_en + "/TP_energy:ieta>>diffieta3", "TP_energy>0 && {}>{}".format(rh_en, cut), "COLZ")
        r.gDirectory.Get("diffieta3").SetTitle("Energy Difference RH vs TP (RH E_{T} > " + cut + ");ieta;RecHit E_{T} / TrigPrim E_{T}")
        c.SaveAs(outfile)

        tps.Draw(rh_en + "/TP_energy:ieta>>diffieta4", "(TP_energy>0 && {}>{}) * ".format(rh_en, cut) + hep17, "COLZ")
        r.gDirectory.Get("diffieta4").SetTitle("Energy Difference RH vs TP (HEP17) (RH E_{T} > " + cut + ");ieta;RecHit E_{T} / TrigPrim E_{T}")
        c.SaveAs(outfile)

        tps.Draw(rh_en + "/TP_energy:ieta>>diffieta5", "(TP_energy>0 && {}>{}) * !".format(rh_en, cut) + hep17, "COLZ")
        r.gDirectory.Get("diffieta5").SetTitle("Energy Difference RH vs TP (not HEP17) (RH E_{T} > " + cut + ");ieta;RecHit E_{T} / TrigPrim E_{T}")
        c.SaveAs(outfile)

    tps.Draw(rh_en + ":TP_energy>>cmphb(100, 0, 200, 100, 0, 200)", "abs(ieta) <= 16", "COLZ")
    f1.Draw("same")
    r.gDirectory.Get("cmphb").SetTitle("HB energy comparison;TrigPrim E_{T};RecHit E_{T}")
    c.SaveAs(outfile)

    tps.Draw(rh_en + ":TP_energy>>cmphe(100, 0, 100, 100, 0, 100)", "(abs(ieta)>16 && abs(ieta)<29) * " + nhep17, "COLZ")
    f1.Draw("same")
    r.gDirectory.Get("cmphe").SetTitle("HE (no HEP17) energy comparison;TrigPrim E_{T};RecHit E_{T}")
    c.SaveAs(outfile)

    tps.Draw(rh_en + ":TP_energy>>cmphep17(100, 0, 100, 100, 0, 100)", hep17, "COLZ")
    f1.Draw("same")
    r.gDirectory.Get("cmphep17").SetTitle("HEP17 energy comparison;TrigPrim E_{T};RecHit E_{T}")
    c.SaveAs(outfile)

    tps.Draw(rh_en + ":TP_energy/0.7>>cmphf(80, 0, 80, 80, 0, 80)", "abs(ieta) > 29", "COLZ")
    f1.Draw("same")
    r.gDirectory.Get("cmphf").SetTitle("HF energy comparison;TrigPrim E_{T} / 0.7;RecHit E_{T}")
    c.SaveAs(outfile)

# for i in range(1, 42):
#     tps.Draw(rh_en + ":TP_energy>>cmp(100, 0, 200, 100, 0, 200)", "ieta == {}".format(i), "COLZ")
#     f1.Draw("same")
#     r.gDirectory.Get("cmp").SetTitle("Energy comparison (i#eta {});TrigPrim E_{{T}};RecHit E_{{T}}".format(i))
#     c.SaveAs(outfile)
#     tps.Draw(rh_en + ":TP_energy>>cmp(100, 0, 200, 100, 0, 200)", "ieta == -{}".format(i), "COLZ")
#     f1.Draw("same")
#     r.gDirectory.Get("cmp").SetTitle("Energy comparison (i#eta -{});TrigPrim E_{{T}};RecHit E_{{T}}".format(i))
#     c.SaveAs(outfile)

    c.SaveAs(outfile + ']')
