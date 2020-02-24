from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLine
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad, TLatex
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kSpring, kBlack
import numpy as np

from array import array
import math
import sys,os

x_elas={}
x_radf={}
x_radc={}
x_radf, mu_radf, sig_radf = array('d'), array('d'), array('d')
x_radc, mu_radc, sig_radc = array('d'), array('d'), array('d')


trash_can=[]
g_elas=[]
for ii in range(1,7):
    fin1 = open('p_electron_delta_p_electron_{}_sim.txt'.format(ii),'r')
    x_elas, mu_elas, sig_elas, zeros = array('d'), array('d'), array('d'), array('d')
    # look at ISR events for now

    for line in fin1:
        row = line.split()        
        x=float(row[0])
        m=float(row[1])
        s=float(row[2])
        print(' {} {} {} '.format(x,m,s))
        x_elas.append(x)
        mu_elas.append(-m)
        sig_elas.append(s)
        zeros.append(0.0)

    graph = TGraphErrors(len(x_elas), x_elas, mu_elas, zeros, sig_elas)
    graph.SetName('g_elastic_elfd_prctof_sim_s{}'.format(ii))
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(kBlack)
    trash_can.append(graph)
    g_elas.append(graph)

#rf is rad elastic with proton in FTOF
g_rf=[]
for ii in range(1,7):
    fin2 = open('p_ele_dp_ele_from_angles_FTOF_{}_isr_bcsim.txt'.format(ii),'r')
    x_rf, mu_rf, sig_rf, zeros = array('d'), array('d'), array('d'), array('d')
    # look at ISR events for now

    for line in fin2:
        row = line.split()        
        x=float(row[0])
        m=float(row[1])
        s=float(row[2])
        print(' {} {} {} '.format(x,m,s))
        x_rf.append(x)
        mu_rf.append(m)
        sig_rf.append(s)
        zeros.append(0.0)

    graph = TGraphErrors(len(x_rf), x_rf, mu_rf, zeros, sig_rf)
    graph.SetName('g_rf_elfd_prftof_sim_s{}'.format(ii))
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(kRed)    
    trash_can.append(graph)
    g_rf.append(graph)

g_rc=[]
for ii in range(1,7):
    fin3 = open('p_ele_dp_ele_from_angles_CTOF_{}_isr_bcsim.txt'.format(ii),'r')
    x_rc, mu_rc, sig_rc, zeros = array('d'), array('d'), array('d'), array('d')
    # look at ISR events for now
    x_temp=[] 
    for line in fin3:
        row = line.split()        
        x=float(row[0])
        m=float(row[1])
        s=float(row[2])
        print(' {} {} {} '.format(x,m,s))
        #x_rc.append(x)
        mu_rc.append(m)
        sig_rc.append(s)
        zeros.append(0.0)
        x_temp.append(x)
    np_x_temp = np.array(x_temp)

    np_x_temp += 0.5 *( np_x_temp[1]-np_x_temp[0])
    for xx in np_x_temp:
        x_rc.append(xx)

    graph = TGraphErrors(len(x_rc), x_rc, mu_rc, zeros, sig_rc)
    graph.SetName('g_rc_elfd_prctof_sim_s{}'.format(ii))    
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(kRed+3)
    trash_can.append(graph)
    g_rc.append(graph)

label = TLatex()
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.05)
label.SetTextColor(1)


out_name='comp_simulation_elastic_to_rad_elastic.pdf'
can = TCanvas('can','can',1800,600)
can.Print('{}['.format(out_name))
can.Divide(3,2)
for ii in range(1,7):
    can.cd(ii)
    ss = ii-1
    mg = TMultiGraph()
    mg.Add(g_elas[ss])
    mg.Add(g_rf[ss])
    mg.Add(g_rc[ss])
    mg.Draw('AP')
    mg.GetHistogram().GetXaxis().SetRangeUser(2.5,10)
    l1 = TLine( 2.5, 0.0, 10.0, 0.0 )
    l1.SetLineColor(kRed)    
    can.Update()
    mg.GetHistogram().SetMinimum(-0.35)
    mg.GetHistogram().SetMaximum(0.35)
    mg.Draw('AP')
    l1.Draw('same')
    label.DrawLatex(0.1, 0.925, '#Delta P_{e} vs P_{e} Rad Elas. and Elas. S' + str(ii))
    label.DrawLatex(0.5, 0.015, 'P_{e} (GeV)')
    label.SetTextAngle(90)
    label.DrawLatex(0.035, 0.5, '#Delta P_{e}')
    label.SetTextAngle(0)
    trash_can.append(label)
    trash_can.append(mg)
    trash_can.append(l1)



can.Print(out_name)
can.Print('{}]'.format(out_name))

    
'''
    
for line in f:
    # Split on any whitespace (including tab characters)
    row = line.split()
    # Convert strings to numeric values:
    row[1] = float(row[1])
    row[2] = int(row[2])
    # Append to our list of lists:
    rows.append(row)
'''
