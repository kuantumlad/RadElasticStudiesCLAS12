from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLine
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad, TLatex
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kSpring
import numpy as np

from array import array
import math
import sys,os

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)


def getFitLimits(h_temp, sig):
    temp_mean = h_temp.GetMean()
    temp_rms = h_temp.GetRMS()
    temp_min_fit_range = temp_mean - sig*temp_rms
    temp_max_fit_range = temp_mean + sig*temp_rms
    return [ temp_min_fit_range, temp_max_fit_range, temp_mean, temp_rms ]


def getFitFractionalHeight( h_temp, h_name, percent_max, save_fits):
    # fitting routine adopted from Andrew Puckett.
    print(' special fitting routine ')
    xlow, xhigh, histmax = 0, 0, 0
    binlow, binhigh, binmax = 0, 0, 0

    
    binmax = h_temp.GetMaximumBin()
    histmax = h_temp.GetMaximum()#GetBinContent(h_temp.GetMaximumBin())
    #print('histmax ' + str(histmax))
    
    binlow = binmax
    binhigh = binmax

    while h_temp.GetBinContent(binhigh) >= percent_max*histmax and binhigh <= h_temp.GetNbinsX() : binhigh+=1
    while h_temp.GetBinContent(binlow) >= percent_max*histmax and binlow > 1 : binlow-=1
    
    xlow = h_temp.GetBinLowEdge(binlow)
    xhigh = h_temp.GetBinLowEdge(binhigh+1)
    
    print(h_name)
    print(' >> bin high values ' + str(binhigh) + ' bin low ' + str(binlow) )
    print(" >> values used " + str(xlow) + " " + str(xhigh) + " " + str(histmax) )
    
    fit_temp = TF1(h_name,'gaus', xlow, xhigh )
    fit_temp.SetParameter(0, histmax)
    fit_temp.SetParameter(1, h_temp.GetBinCenter(h_temp.GetMaximumBin()))
    
                          
    
    h_temp.Fit(h_name,"R") # add Q for quiet
    save_fits.append(fit_temp)
    
    temp_mean = fit_temp.GetParameter(1)
    temp_rms = fit_temp.GetParameter(2)
    
    return fit_temp#[temp_mean, temp_rms ]


def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False,
                     y_fit_range=None, landscape=False, x_range=None, vline=None,
                     hline=None):

    root_garbage_can = []
    
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    for i in range(1,7):
        canvas.cd(i)
        
        #if isinstance(histos[title_formatter.format(i)], TH1F):
        if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
            if y_fit_range:
                fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])
                
            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                
            histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
            histos.get(title_formatter.format(i), default_histo).Draw()

            if y_fit_range:
                label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
        #elif isinstance(histos[title_formatter.format(i)], TH2F):
        elif isinstance(histos.get(title_formatter.format(i), default_histo), TH2F):
            # histos[title_formatter.format(i)].Draw('colz')
            histos.get(title_formatter.format(i), default_histo).Draw('colz')
            if log:
                gPad.SetLogz() 
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass

        if vline:
            line = TLine(vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmin(),
                         vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmax())
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)

        if hline:
            xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 
            line = TLine(xmin, hline, xmax, hline)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
            
        if title:
            histos.get(title_formatter.format(i), default_histo).SetTitle('')            
            label.DrawLatex(0.1, 0.925, title + ' S' + str(i))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

    

# script to plot, view, correct electron and proton momentum
# in each theta and phi bin

# load data and file
datatype = sys.argv[2]
hhs={}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj


base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_elastic_kin_corr_'+datatype+'.pdf'
can = TCanvas('can','can')
can.SetCanvasSize(1200,1200)
can.Print('{}['.format(mon_out_file_name))

lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)


gStyle.SetOptStat(00000)

#load theta dependence of deltap/p
dp_ele_from_angles_ctof_isr={}
dp_ele_from_angles_ftof_isr={}

dp_ele_from_angles_ctof_fsr={}
dp_ele_from_angles_ftof_fsr={}

#load proton info
dp_pro_from_angles_ctof_isr={}
dp_pro_from_angles_ftof_isr={}

dp_pro_from_angles_ctof_fsr={}
dp_pro_from_angles_ftof_fsr={}

#load electron info per theta phi bin
dp_ele_ctof_isr_theta_phi = {}
dp_ele_ftof_isr_theta_phi = {}
dp_ele_ctof_fsr_theta_phi = {}
dp_ele_ftof_fsr_theta_phi = {}

dp_pro_ctof_isr_theta_phi = {}
dp_pro_ftof_isr_theta_phi = {}
dp_pro_ctof_fsr_theta_phi = {}
dp_pro_ftof_fsr_theta_phi = {}



for ss in range(0,7):
    temp_h = []
    temp_h2 = []
    temp_h3 = []
    temp_h4 = []

    temp_pro_h = []
    temp_pro_h2 = []
    temp_pro_h3 = []
    temp_pro_h4 = []
    for tb in range(0,10):
        #dp_ele_from_angles_FTOF_1_thetabin4_isr
        if 'dp_ele_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr' in hhs:
            print('hi')
            temp_h.append(hhs['dp_ele_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr'])
        if 'dp_ele_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr' in hhs:
            temp_h2.append(hhs['dp_ele_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr'])
        if 'dp_ele_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr' in hhs:
            temp_h3.append(hhs['dp_ele_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr'])
        if 'dp_ele_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr' in hhs:
            temp_h4.append(hhs['dp_ele_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr'])

    for tb in range(0,10):
        if 'dp_pro_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr' in hhs:
            temp_pro_h.append(hhs['dp_pro_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr'])
        if 'dp_pro_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr' in hhs:
            temp_pro_h2.append(hhs['dp_pro_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_isr'])
        if 'dp_pro_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr' in hhs:
            temp_pro_h3.append(hhs['dp_pro_from_angles_CTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr'])
        if 'dp_pro_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr' in hhs:
            temp_pro_h4.append(hhs['dp_pro_from_angles_FTOF_'+ str(ss) + '_thetabin' + str(tb) + '_fsr'])



    temp_tb_pb_ctof_isr=[]
    temp_tb_pb_ftof_isr=[]
    temp_tb_pb_ctof_fsr=[]
    temp_tb_pb_ftof_fsr=[]
    for tb in range(0,10):
        temp_pb_ctof_isr = []
        temp_pb_ftof_isr = []
        temp_pb_ctof_fsr = []
        temp_pb_ftof_fsr = []

        for pb in range(0,13):
            if 'dp_ele_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_isr' in hhs:
                temp_pb_ctof_isr.append(hhs['dp_ele_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_isr'])            
            if 'dp_ele_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_isr' in hhs:
                temp_pb_ftof_isr.append(hhs['dp_ele_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_isr'])
            if 'dp_ele_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_fsr' in hhs:
                temp_pb_ctof_fsr.append(hhs['dp_ele_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_fsr'])
            if 'dp_ele_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_fsr' in hhs:
                temp_pb_ftof_fsr.append(hhs['dp_ele_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + '_fsr'])
        
        temp_tb_pb_ctof_isr.append(temp_pb_ctof_isr)
        temp_tb_pb_ftof_isr.append(temp_pb_ftof_isr)
        temp_tb_pb_ctof_fsr.append(temp_pb_ctof_fsr)
        temp_tb_pb_ftof_fsr.append(temp_pb_ftof_fsr)


    dp_ele_ctof_isr_theta_phi[ss]=temp_tb_pb_ctof_isr
    dp_ele_ftof_isr_theta_phi[ss]=temp_tb_pb_ftof_isr
    dp_ele_ctof_fsr_theta_phi[ss]=temp_tb_pb_ctof_fsr
    dp_ele_ftof_fsr_theta_phi[ss]=temp_tb_pb_ftof_fsr

    ##############
    ## for protons 
    temp_pro_tb_pb_ctof_isr=[]
    temp_pro_tb_pb_ftof_isr=[]
    temp_pro_tb_pb_ctof_fsr=[]
    temp_pro_tb_pb_ftof_fsr=[]
    for tb in range(0,5):
        temp_pro_pb_ctof_isr = []
        temp_pro_pb_ftof_isr = []
        temp_pro_pb_ctof_fsr = []
        temp_pro_pb_ftof_fsr = []

        for pb in range(0,10):
            if 'dp_pro_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'isr' in hhs:
                temp_pro_pb_ctof_isr.append(hhs['dp_pro_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'isr'])
            if 'dp_pro_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'isr' in hhs:
                temp_pro_pb_ftof_isr.append(hhs['dp_pro_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'isr'])
            if 'dp_pro_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'fsr' in hhs:
                temp_pro_pb_ctof_fsr.append(hhs['dp_pro_from_angles_CTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'fsr'])
            if 'dp_pro_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'fsr' in hhs:
                temp_pro_pb_ftof_fsr.append(hhs['dp_pro_from_angles_FTOF_' + str(ss) + '_thetabin' + str(tb) + '_phibin' + str(pb) + 'fsr'])
        temp_pro_tb_pb_ctof_isr.append(temp_pro_pb_ctof_isr)
        temp_pro_tb_pb_ftof_isr.append(temp_pro_pb_ftof_isr)
        temp_pro_tb_pb_ctof_fsr.append(temp_pro_pb_ctof_fsr)
        temp_pro_tb_pb_ftof_fsr.append(temp_pro_pb_ftof_fsr)
    dp_pro_ctof_isr_theta_phi[ss]=temp_tb_pb_ctof_isr
    dp_pro_ftof_isr_theta_phi[ss]=temp_tb_pb_ftof_isr
    dp_pro_ctof_fsr_theta_phi[ss]=temp_tb_pb_ctof_fsr
    dp_pro_ftof_fsr_theta_phi[ss]=temp_tb_pb_ftof_fsr

        
        
    dp_ele_from_angles_ctof_isr[ss]=temp_h
    dp_ele_from_angles_ftof_isr[ss]=temp_h2
    dp_ele_from_angles_ctof_fsr[ss]=temp_h3
    dp_ele_from_angles_ftof_fsr[ss]=temp_h4

    dp_pro_from_angles_ctof_isr[ss]=temp_pro_h
    dp_pro_from_angles_ftof_isr[ss]=temp_pro_h2
    dp_pro_from_angles_ctof_fsr[ss]=temp_pro_h3
    dp_pro_from_angles_ftof_fsr[ss]=temp_pro_h4
#########################################################################
#########################################################################


print(dp_ele_from_angles_ctof_isr)


#########################################################################
 # first plot event selection related items
 
plot_sector_page(can, hhs, 'w_CTOF_{}', lab, save_name=mon_out_file_name,
                 title='Electron (Forward), Proton (CTOF) W', xtitle='W (GeV)')
 
plot_sector_page(can, hhs, 'w_FTOF_{}', lab, save_name=mon_out_file_name,
                 title='Electron (Forward), Proton (FTOF) W', xtitle='W (GeV)')


can.Clear()
can.Divide(2,2)
can.cd(1)
h_w_eep_ftof = hhs['w_eep_FTOF_elastic']
h_w_eep_ftof.SetTitle('Elastic: Electron (Forward), Proton (FTOF), W Pass All But W Cut; W (GeV); Counts')
h_w_eep_ftof.Draw()
can.Update()
ymax = gPad.GetUymax()
lw_min = TLine(1.15, 0.0, 1.15, ymax)
lw_min.SetLineColor(kRed)
lw_min.Draw('same')

can.cd(2)
h_angleep_eep_ftof = hhs['angle_ep_eep_FTOF_elastic']
h_angleep_eep_ftof.SetTitle('Elastic: Electron (Forward), Proton (FTOF), #Phi_{ep} Pass All But #Phi_{ep} Cut; #Phi_{ep} (deg); Counts')
h_angleep_eep_ftof.Draw()
can.Update()
ymax = gPad.GetUymax()
laep_min = TLine(178, 0.0, 178, ymax)
laep_min.SetLineColor(kRed)
laep_min.Draw('same')
can.cd(3)
h_theta_gamma_eep_ftof = hhs['theta_gamma_eep_FTOF_elastic']
h_theta_gamma_eep_ftof.SetTitle('Elastic: Electron (Forward), Proton (FTOF), #Theta_{beam#gamma} Pass All But #Theta_{beam#gamma} Cut; #Theta_{beam#gamma} (deg) ; Counts')
h_theta_gamma_eep_ftof.Draw()
can.Update()
ymax = gPad.GetUymax()
lg_min = TLine(3, 0.0, 3, ymax)
lg_min.SetLineColor(kRed)
lg_min.Draw('same')
can.cd(4)
h_miss_mass_eep_ftof = hhs['missing_mass_eep_FTOF_elastic']
h_miss_mass_eep_ftof.SetTitle('Elastic: Electron (Forward), Proton (FTOF), MM2 Pass All But MM2 Cut; MM2 (GeV^2); Counts')
h_miss_mass_eep_ftof.Draw()
can.Update()
ymax = gPad.GetUymax()
lmm_min = TLine(-0.1, 0.0, -0.1, ymax)
lmm_min.SetLineColor(kRed)
lmm_min.Draw('same')
lmm_max = TLine(0.1, 0.0, 0.1, ymax)
lmm_max.SetLineColor(kRed)
lmm_max.Draw('same')

can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,2)
can.cd(1)
h_w_eep_ctof = hhs['w_eep_CTOF_elastic']
h_w_eep_ctof.SetTitle('Elastic: Electron (Forward), Proton (CTOF), W Pass All But W Cut; W (GeV); Counts')
h_w_eep_ctof.Draw()
can.Update()
ymax = gPad.GetUymax()
lw_min = TLine(1.15, 0.0, 1.15, ymax)
lw_min.SetLineColor(kRed)
lw_min.Draw('same')
can.cd(2)
h_angleep_eep_ctof = hhs['angle_ep_eep_CTOF_elastic']
h_angleep_eep_ctof.SetTitle('Elastic: Electron (Forward), Proton (CTOF), #Phi_{ep} Pass All But #Phi_{ep} Cut; #Phi_{ep} (deg); Counts')
h_angleep_eep_ctof.Draw()
can.Update()
ymax = gPad.GetUymax()
laep_min = TLine(178, 0.0, 178, ymax)
laep_min.SetLineColor(kRed)
laep_min.Draw('same')
can.cd(3)
h_theta_gamma_eep_ctof = hhs['theta_gamma_eep_CTOF_elastic']
h_theta_gamma_eep_ctof.SetTitle('Elastic: Electron (Forward), Proton (CTOF), #Theta_{#gamma} Pass All But #Theta_{#gamma} Cut; #Theta_{#gamma} (deg) ; Counts')
h_theta_gamma_eep_ctof.Draw()
can.Update()
ymax = gPad.GetUymax()
lg_min = TLine(3, 0.0, 3, ymax)
lg_min.SetLineColor(kRed)
lg_min.Draw('same')
can.cd(4)
h_miss_mass_eep_ctof = hhs['missing_mass_eep_CTOF_elastic']
h_miss_mass_eep_ctof.SetTitle('Elastic: Electron (Forward), Proton (CTOF), MM2 Pass All But MM2 Cut; MM2 (GeV^2); Counts')
h_miss_mass_eep_ctof.Draw()
can.Update()
ymax = gPad.GetUymax()
lmm_min = TLine(-0.4, 0.0, -0.4, ymax)
lmm_min.SetLineColor(kRed)
lmm_min.Draw('same')
lmm_max = TLine(0.4, 0.0, 0.4, ymax)
lmm_max.SetLineColor(kRed)
lmm_max.Draw('same')

can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,1)
can.cd(1)
h_theta_egamma_eep = hhs['theta_egamma_eep_FTOF_elastic']
h_theta_egamma_eep.SetTitle('Elastic: Electron (Forward), Proton (FTOF), Pass All But #theta_{e #gamma}; #theta_{e #gamma}; counts')
h_theta_egamma_eep.Draw()
can.Update()
ymax = gPad.GetUymax()
leg_min = TLine(3, 0.0, 3, ymax)
leg_min.SetLineColor(kRed)
leg_min.Draw('same')
can.cd(2)
h_theta_egamma_eep_ctof = hhs['theta_egamma_eep_CTOF_elastic']
h_theta_egamma_eep_ctof.SetTitle('Elastic: Electron (Forward), Proton (CTOF), Pass All But #theta_{e #gamma}; #theta_{e #gamma}; counts')
h_theta_egamma_eep_ctof.Draw()
can.Update()
ymax = gPad.GetUymax()
leg_min = TLine(3, 0.0, 3, ymax)
leg_min.SetLineColor(kRed)
leg_min.Draw('same')
can.Print(mon_out_file_name)

# draw the W vs sum of angles for all events with electron and at least one proton
can.Clear()
can.Divide(2,1)
can.cd(1)
h_w_theta_sum_ftof = hhs['w_theta_sum_FTOF']
h_w_theta_sum_ftof.SetTitle('Electron (Forward), Proton (FTOF), All, W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ftof.Draw('colz')
can.cd(2)
h_w_theta_sum_ctof = hhs['w_theta_sum_CTOF']
h_w_theta_sum_ctof.SetTitle('Electron (Forward), Proton (CTOF), All, W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ctof.Draw('colz')
can.Print(mon_out_file_name)


can.Clear()
can.Divide(2,1)
can.cd(1)
h_w_theta_sum_ftof_isr = hhs['w_theta_sum_pass_w_elastic_FTOF']
h_w_theta_sum_ftof_isr.SetTitle('Elastic: El. (Forward), Pr. (FTOF), Pass W, #gamma (IRS), W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ftof_isr.Draw('colz')
can.cd(2)
h_w_theta_sum_ctof_isr = hhs['w_theta_sum_pass_w_elastic_CTOF']
h_w_theta_sum_ctof_isr.SetTitle('Elastic: El. (Forward), Pr. (CTOF), Pass W, #gamma (ISR), W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ctof_isr.Draw('colz')
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,2)
can.cd(1)
h_w_theta_sum_ftof_isr_final = hhs['w_theta_sum_pass_all_FTOF_isr']
h_w_theta_sum_ftof_isr_final.SetTitle('El (Forward), Pr (FTOF), #gamma (IRS), Pass All, W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ftof_isr_final.Draw('colz')
can.cd(2)
h_w_theta_sum_ctof_isr_final = hhs['w_theta_sum_pass_all_CTOF_isr']
h_w_theta_sum_ctof_isr_final.SetTitle('El (Forward), Pr (CTOF), #gamma (ISR), Pass All, W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ctof_isr_final.Draw('colz')
can.cd(3)
h_w_theta_sum_ftof_fsr_final = hhs['w_theta_sum_pass_all_FTOF_fsr']
h_w_theta_sum_ftof_fsr_final.SetTitle('El (Forward), Pr (FTOF), #gamma (FRS), Pass All, W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ftof_fsr_final.Draw('colz')
can.cd(4)
h_w_theta_sum_ctof_fsr_final = hhs['w_theta_sum_pass_all_CTOF_fsr']
h_w_theta_sum_ctof_fsr_final.SetTitle('El (Forward), Pr. (CTOF), #gamma (FSR), Pass All, W vs #theta_{el} + #theta_{pr}; W (GeV); #theta_{el} + #theta_{pr} (deg)')
h_w_theta_sum_ctof_fsr_final.Draw('colz')
can.Print(mon_out_file_name)


####################################################################################################
## check the final kinematic distributions of the proton and electron in each of the four cases
plot_sector_page(can, hhs, 'p_ele_theta_ele_CTOF_{}_isr', lab, save_name=mon_out_file_name,
                 title='Proton (CTOF), #gamma (ISR)', xtitle='p_{e}', ytitle='#theta_{e}', log=True)

plot_sector_page(can, hhs, 'p_ele_theta_ele_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 title='Proton (FTOF), #gamma (ISR)', xtitle='p_{e}', ytitle='#theta_{e}', log=True)

plot_sector_page(can, hhs, 'p_pro_theta_pro_CTOF_{}_isr', lab, save_name=mon_out_file_name,
                 title='Proton (CTOF) #gamma (ISR)', xtitle='p_{p}', ytitle='#theta_{p}', log=True)

plot_sector_page(can, hhs, 'p_pro_theta_pro_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 title='Proton (FTOF), #gamma (ISR)', xtitle='p_{p}', ytitle='#theta_{p}', log=True)

##################################################################################################
## same but for FSR 

####################################################################################################
## check the final kinematic distributions of the proton and electron in each of the four cases
plot_sector_page(can, hhs, 'p_ele_theta_ele_CTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 title='Proton (CTOF), #gamma (FSR)', xtitle='p_{e}', ytitle='#theta_{e}', log=True)

plot_sector_page(can, hhs, 'p_ele_theta_ele_FTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 title='Proton (FTOF), #gamma (FSR)', xtitle='p_{e}', ytitle='#theta_{e}', log=True)

plot_sector_page(can, hhs, 'p_pro_theta_pro_CTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 title='Proton (CTOF), #gamma (FSR)', xtitle='p_{p}', ytitle='#theta_{p}', log=True)

plot_sector_page(can, hhs, 'p_pro_theta_pro_FTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 title='Proton (FTOF), #gamma (FSR)', xtitle='p_{p}', ytitle='#theta_{p}', log=True)

can.Clear()
can.Divide(7,6)

ii=0

dp_ele_theta_ctof_isr = {}
dp_ele_theta_ftof_isr = {}
dp_ele_theta_ctof_fsr = {}
dp_ele_theta_ftof_fsr = {}

# fit the dp for electron per theta bin in each sector
can_skip=7
store_fit =[]
for ss in range(1,7):
    delta_theta = 3
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_from_angles_ctof_isr[ss])

    dp_ele_theta_x, dp_ele_theta_errx = array('d'), array('d')
    dp_ele_fit_y, dp_ele_fit_erry = array('d'), array('d')


    for tb in range(0,max_tb):
        can.cd(ii+1)
        temp_dp_ele_fit = TF1('dp_ele_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb),'gaus',-0.075,0.075)
        dp_ele_from_angles_ctof_isr[ss][tb].SetTitle('Electron (Forward), Proton (CTOF), #gamma (ISR) #Delta P/P_{rec} Sect.' + str(ss) + ' Theta Bin min' + str(theta_min) + ' max ' + str(theta_max) + ' ; #Delta P/P_{rec}; counts')
        dp_ele_from_angles_ctof_isr[ss][tb].RebinX(4)
        dp_ele_from_angles_ctof_isr[ss][tb].Draw()
        #dp_ele_from_angles_ctof_isr[ss][tb].Fit('dp_ele_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_ele_fit = getFitFractionalHeight(dp_ele_from_angles_ctof_isr[ss][tb],'dp_ele_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb), 0.4, store_fit)

        print( " TESTING MY NEW FIT " )

        dp_ele_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_ele_theta_errx.append(0.0)
        print(' check fit par ' )
        print(dp_ele_fit.GetParameter(1))
        dp_ele_fit_y.append(dp_ele_fit.GetParameter(1))#[0])
        dp_ele_fit_erry.append(dp_ele_fit.GetParameter(2))#[1])

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    ii=can_skip
    can_skip+=7
    temp_dp_ele = []
    temp_dp_ele.append(dp_ele_theta_x)
    temp_dp_ele.append(dp_ele_theta_errx)
    temp_dp_ele.append(dp_ele_fit_y)
    temp_dp_ele.append(dp_ele_fit_erry)
    
    dp_ele_theta_ctof_isr[ss] = temp_dp_ele


can.Print(mon_out_file_name)

can.Clear()
can.Divide(7,6)
#######
## proton is in FTOF and ISR 
# fit the dp for electron per theta bin in each sector
ii=0
can_skip=7
for ss in range(1,7):
    delta_theta = 5
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_from_angles_ftof_isr[ss])

    dp_ele_theta_x, dp_ele_theta_errx = array('d'), array('d')
    dp_ele_fit_y, dp_ele_fit_erry = array('d'), array('d')


    for tb in range(0,max_tb):
        can.cd(ii+1)
        temp_fit_limits = getFitLimits(dp_ele_from_angles_ftof_isr[ss][tb],3.0)
        #test_fit_results = getFitFractionalHeight(dp_ele_from_angles_ftof_isr[ss][tb],'dp_ele_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb), 0.4, store_fit)
        print( " TESTING MY NEW FIT " )
        

        dp_ele_fit = TF1('dp_ele_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb),'gaus',temp_fit_limits[0],temp_fit_limits[1])
        dp_ele_fit.SetParameter(1,temp_fit_limits[2])
        dp_ele_fit.SetParameter(2,temp_fit_limits[3])

        dp_ele_from_angles_ftof_isr[ss][tb].SetTitle('Electron #Delta P/P_{rec} FTOF ISR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_ele_from_angles_ftof_isr[ss][tb].RebinX(4)
        dp_ele_from_angles_ftof_isr[ss][tb].Draw()
        dp_ele_from_angles_ftof_isr[ss][tb].Fit('dp_ele_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb),'R')        

        dp_ele_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_ele_theta_errx.append(0.0)
        dp_ele_fit_y.append(dp_ele_fit.GetParameter(1))#test_fit_results.GetParameter(1))#
        dp_ele_fit_erry.append(dp_ele_fit.GetParameter(2))#test_fit_results.GetParameter(2))#dp_ele_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta
    
    #can.cd(can_skip)
    ii=can_skip
    can_skip+=7

    temp_dp_ele = []
    temp_dp_ele.append(dp_ele_theta_x)
    temp_dp_ele.append(dp_ele_theta_errx)
    temp_dp_ele.append(dp_ele_fit_y)
    temp_dp_ele.append(dp_ele_fit_erry)
    
    dp_ele_theta_ftof_isr[ss] = temp_dp_ele
    
can.Print(mon_out_file_name)


#########################################
## proton is in the CTOF with FSR
can.Clear()
can.Divide(6,6)

# fit the dp for electron per theta bin in each sector
can_skip=6
ii=0
for ss in range(1,7):
    delta_theta = 3
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_from_angles_ctof_fsr[ss])

    dp_ele_theta_x, dp_ele_theta_errx = array('d'), array('d')
    dp_ele_fit_y, dp_ele_fit_erry = array('d'), array('d')


    for tb in range(0,max_tb):
        can.cd(ii+1)
        temp_fit_limits = getFitLimits(dp_ele_from_angles_ctof_fsr[ss][tb],2.5)
        dp_ele_fit = TF1('dp_ele_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'landau',temp_fit_limits[0], temp_fit_limits[1])
        dp_ele_from_angles_ctof_fsr[ss][tb].SetTitle('Electron #Delta P/P_{rec} CTOF FSR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_ele_from_angles_ctof_fsr[ss][tb].Draw()
        dp_ele_from_angles_ctof_fsr[ss][tb].Fit('dp_ele_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_ele_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_ele_theta_errx.append(0.0)
        dp_ele_fit_y.append(dp_ele_fit.GetParameter(1))
        dp_ele_fit_erry.append(dp_ele_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    ii=can_skip
    can_skip+=6
    temp_dp_ele = []
    temp_dp_ele.append(dp_ele_theta_x)
    temp_dp_ele.append(dp_ele_theta_errx)
    temp_dp_ele.append(dp_ele_fit_y)
    temp_dp_ele.append(dp_ele_fit_erry)
    
    dp_ele_theta_ctof_fsr[ss] = temp_dp_ele
    
can.Print(mon_out_file_name)


#########################################
## proton is in the FTOF with FSR
can.Clear()
can.Divide(6,6)
ii=0
can_skip=6
# fit the dp for electron per theta bin in each sector
for ss in range(1,7):
    delta_theta = 5
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_from_angles_ftof_fsr[ss])

    dp_ele_theta_x, dp_ele_theta_errx = array('d'), array('d')
    dp_ele_fit_y, dp_ele_fit_erry = array('d'), array('d')


    for tb in range(0,max_tb):
        can.cd(ii+1)
        dp_ele_fit = TF1('dp_ele_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'landau',-0.075,0.075) # Change back?
        dp_ele_from_angles_ftof_fsr[ss][tb].SetTitle('Electron #Delta P/P_{rec} FTOF FSR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_ele_from_angles_ftof_fsr[ss][tb].RebinX(2)
        dp_ele_from_angles_ftof_fsr[ss][tb].Draw()
        dp_ele_from_angles_ftof_fsr[ss][tb].Fit('dp_ele_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_ele_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_ele_theta_errx.append(0.0)
        dp_ele_fit_y.append(dp_ele_fit.GetParameter(1))
        dp_ele_fit_erry.append(dp_ele_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    ii=can_skip
    can_skip+=6
    temp_dp_ele = []
    temp_dp_ele.append(dp_ele_theta_x)
    temp_dp_ele.append(dp_ele_theta_errx)
    temp_dp_ele.append(dp_ele_fit_y)
    temp_dp_ele.append(dp_ele_fit_erry)
    
    dp_ele_theta_ftof_fsr[ss] = temp_dp_ele
    
can.Print(mon_out_file_name)
        

can.Clear()
can.Divide(7,6)

ii=0
can_skip=7
dp_pro_theta_ctof_isr = {}
dp_pro_theta_ftof_isr = {}
dp_pro_theta_ctof_fsr = {}
dp_pro_theta_ftof_fsr = {}

#delta p for proton with proton in CTOF and ISR
# fit the dp for electron per theta bin in each sector
for ss in range(1,7):
    delta_theta = 5
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_from_angles_ctof_isr[ss])

    dp_pro_theta_x, dp_pro_theta_errx = array('d'), array('d')
    dp_pro_fit_y, dp_pro_fit_erry = array('d'), array('d')
    
    
    for tb in range(0,max_tb):
        can.cd(ii+1)

        temp_fit_limits = getFitLimits(dp_pro_from_angles_ctof_isr[ss][tb], 2.5)

        dp_pro_fit = TF1('dp_pro_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb),'gaus',temp_fit_limits[0], temp_fit_limits[1])
        dp_pro_from_angles_ctof_isr[ss][tb].SetTitle('Proton #Delta P/P_{rec} CTOF ISR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_pro_from_angles_ctof_isr[ss][tb].Draw()
        dp_pro_from_angles_ctof_isr[ss][tb].Fit('dp_pro_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_pro_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_pro_theta_errx.append(0.0)
        dp_pro_fit_y.append(dp_pro_fit.GetParameter(1))
        dp_pro_fit_erry.append(dp_pro_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    ii=can_skip
    can_skip+=7
    temp_dp_pro = []
    temp_dp_pro.append(dp_pro_theta_x)
    temp_dp_pro.append(dp_pro_theta_errx)
    temp_dp_pro.append(dp_pro_fit_y)
    temp_dp_pro.append(dp_pro_fit_erry)
    
    dp_pro_theta_ctof_isr[ss] = temp_dp_pro

can.Print(mon_out_file_name)

can.Clear()
can.Divide(7,6)

#delta p for proton with proton in FTOF and ISR
# fit the dp for electron per theta bin in each sector
ii=0
can_skip=7
for ss in range(1,7):
    delta_theta = 5
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_from_angles_ftof_isr[ss])

    dp_pro_theta_x, dp_pro_theta_errx = array('d'), array('d')
    dp_pro_fit_y, dp_pro_fit_erry = array('d'), array('d')

    for tb in range(0,max_tb):
        can.cd(ii+1)
        
        temp_fit_limits = getFitLimits(dp_pro_from_angles_ftof_isr[ss][tb],2)
        dp_pro_fit = TF1('dp_pro_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb),'gaus',temp_fit_limits[0], temp_fit_limits[1])
        dp_pro_from_angles_ftof_isr[ss][tb].SetTitle('Proton #Delta P/P_{rec} FTOF ISR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_pro_from_angles_ftof_isr[ss][tb].RebinX(2)
        dp_pro_from_angles_ftof_isr[ss][tb].Draw()
        dp_pro_from_angles_ftof_isr[ss][tb].Fit('dp_pro_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_pro_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_pro_theta_errx.append(0.0)
        dp_pro_fit_y.append(dp_pro_fit.GetParameter(1))
        dp_pro_fit_erry.append(dp_pro_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    ii=can_skip
    can_skip+=7
    temp_dp_pro = []
    temp_dp_pro.append(dp_pro_theta_x)
    temp_dp_pro.append(dp_pro_theta_errx)
    temp_dp_pro.append(dp_pro_fit_y)
    temp_dp_pro.append(dp_pro_fit_erry)
    
    dp_pro_theta_ftof_isr[ss] = temp_dp_pro

can.Print(mon_out_file_name)

can.Clear()
can.Divide(6,6)

ii=0
can_skip=6
#delta p for proton with proton in CTOF and FSR
# fit the dp for electron per theta bin in each sector
for ss in range(1,6):
    delta_theta = 3
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_from_angles_ctof_fsr[ss])

    dp_pro_theta_x, dp_pro_theta_errx = array('d'), array('d')
    dp_pro_fit_y, dp_pro_fit_erry = array('d'), array('d')


    for tb in range(0,max_tb):
        can.cd(ii+1)
        
        temp_fit_limits = getFitLimits(dp_pro_from_angles_ctof_fsr[ss][tb],2)
        dp_pro_fit = TF1('dp_pro_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'landau', temp_fit_limits[0], temp_fit_limits[1])
        dp_pro_from_angles_ctof_fsr[ss][tb].SetTitle('Proton #Delta P/P_{rec} CTOF FSR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_pro_from_angles_ctof_fsr[ss][tb].RebinX(2)
        dp_pro_from_angles_ctof_fsr[ss][tb].Draw()
        dp_pro_from_angles_ctof_fsr[ss][tb].Fit('dp_pro_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_pro_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_pro_theta_errx.append(0.0)
        dp_pro_fit_y.append(dp_pro_fit.GetParameter(1))
        dp_pro_fit_erry.append(dp_pro_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    ii=can_skip
    can_skip+=6
    temp_dp_pro = []
    temp_dp_pro.append(dp_pro_theta_x)
    temp_dp_pro.append(dp_pro_theta_errx)
    temp_dp_pro.append(dp_pro_fit_y)
    temp_dp_pro.append(dp_pro_fit_erry)
    
    dp_pro_theta_ctof_fsr[ss] = temp_dp_pro

can.Print(mon_out_file_name)

can.Clear()
can.Divide(7,6)
ii=0

# delta p for proton with proton in FTOF and FSR
# fit the dp for electron per theta bin in each sector
for ss in range(1,7):
    delta_theta = 5
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_from_angles_ftof_fsr[ss])

    dp_pro_theta_x, dp_pro_theta_errx = array('d'), array('d')
    dp_pro_fit_y, dp_pro_fit_erry = array('d'), array('d')

    for tb in range(0,max_tb):
        can.cd(ii+1)
        
        temp_fit_limits = getFitLimits(dp_pro_from_angles_ftof_fsr[ss][tb],3.0)
        dp_pro_fit = TF1('dp_pro_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'landau',temp_fit_limits[0], temp_fit_limits[1])
        dp_pro_from_angles_ftof_fsr[ss][tb].SetTitle('Proton #Delta P/P_{rec} FTOF FSR Sector '+str(ss)+' Theta ['+str(theta_min)+','+str(theta_max)+']; #Delta P/P_{rec}; counts')
        dp_pro_from_angles_ftof_fsr[ss][tb].RebinX(2)
        dp_pro_from_angles_ftof_fsr[ss][tb].Draw()
        dp_pro_from_angles_ftof_fsr[ss][tb].Fit('dp_pro_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb),'R')        
        dp_pro_theta_x.append( theta_min + (theta_max - theta_min)/2.0)
        dp_pro_theta_errx.append(0.0)
        dp_pro_fit_y.append(dp_pro_fit.GetParameter(1))
        dp_pro_fit_erry.append(dp_pro_fit.GetParameter(2))

        ii+=1
        theta_min=theta_max
        theta_max+=delta_theta

    temp_dp_pro = []
    temp_dp_pro.append(dp_pro_theta_x)
    temp_dp_pro.append(dp_pro_theta_errx)
    temp_dp_pro.append(dp_pro_fit_y)
    temp_dp_pro.append(dp_pro_fit_erry)
    
    dp_pro_theta_ftof_fsr[ss] = temp_dp_pro

can.Print(mon_out_file_name)


store_graphs=[]
store_line=[]
can.Clear()
can.Divide(2,3)
#plot delta p of electron when proton is in ctof with ISR
for ss in range(1, len(dp_ele_theta_ctof_isr)+1):
    can.cd(ss)
    #print(dp_ele_theta_ctof_isr[ss][2])
    dp_res_theta_ctof_isr_dependence = TGraphErrors(len(dp_ele_theta_ctof_isr[ss][0]), dp_ele_theta_ctof_isr[ss][0], dp_ele_theta_ctof_isr[ss][2], dp_ele_theta_ctof_isr[ss][1], dp_ele_theta_ctof_isr[ss][3])
    dp_res_theta_ctof_isr_dependence.SetTitle("Electron #Delta p/p Sector " + str(ss) + " vs Theta CTOF ISR; #theta (deg); #Delta p/p")
    dp_res_theta_ctof_isr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ctof_isr_dependence) # save graphs in memory or root complains
    dp_res_theta_ctof_isr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ctof_isr_dependence.GetHistogram().SetMinimum(-0.25)
    l0 = TLine( dp_ele_theta_ctof_isr[ss][0][0], 0.0, dp_ele_theta_ctof_isr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    dp_res_theta_ctof_isr_dependence.Draw("AP")
    l0.Draw('same')

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,3)

#plot delta p of electron when proton is in ftof with ISR
for ss in range(1, len(dp_ele_theta_ftof_isr)+1):
    can.cd(ss)
    if len(dp_ele_theta_ftof_isr[ss][2]) == 0 : continue
    dp_res_theta_ftof_isr_dependence = TGraphErrors(len(dp_ele_theta_ftof_isr[ss][0]), dp_ele_theta_ftof_isr[ss][0], dp_ele_theta_ftof_isr[ss][2], dp_ele_theta_ftof_isr[ss][1], dp_ele_theta_ftof_isr[ss][3])
    dp_res_theta_ftof_isr_dependence.SetTitle("Electron #Delta p/p Sector " + str(ss) + " vs Theta FTOF ISR; #theta (deg); #Delta p/p")
    dp_res_theta_ftof_isr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ftof_isr_dependence) # save graphs in memory or root complains
    dp_res_theta_ftof_isr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ftof_isr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ftof_isr_dependence.Draw("AP")
    l0 = TLine( dp_ele_theta_ftof_isr[ss][0][0], 0.0, dp_ele_theta_ftof_isr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,3)

#plot delta p of electron when proton is in ctof with FSR
for ss in range(1, len(dp_ele_theta_ctof_fsr)+1):
    can.cd(ss)
    if len(dp_ele_theta_ctof_fsr[ss][2]) == 0 : continue
    dp_res_theta_ctof_fsr_dependence = TGraphErrors(len(dp_ele_theta_ctof_fsr[ss][0]), dp_ele_theta_ctof_fsr[ss][0], dp_ele_theta_ctof_fsr[ss][2], dp_ele_theta_ctof_fsr[ss][1], dp_ele_theta_ctof_fsr[ss][3])
    dp_res_theta_ctof_fsr_dependence.SetTitle("Electron #Delta p/p Sector " + str(ss) + " vs Theta CTOF FSR; #theta (deg); #Delta p/p")
    dp_res_theta_ctof_fsr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ctof_fsr_dependence) # save graphs in memory or root complains
    dp_res_theta_ctof_fsr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ctof_fsr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ctof_fsr_dependence.Draw("AP")
    l0 = TLine( dp_ele_theta_ctof_fsr[ss][0][0], 0.0, dp_ele_theta_ctof_fsr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,3)

#plot delta p of electron when proton is in ftof with FSR
for ss in range(1, len(dp_ele_theta_ftof_fsr)+1):
    can.cd(ss)
    if len(dp_ele_theta_ftof_fsr[ss][2]) == 0: continue
    dp_res_theta_ftof_fsr_dependence = TGraphErrors(len(dp_ele_theta_ftof_fsr[ss][0]), dp_ele_theta_ftof_fsr[ss][0], dp_ele_theta_ftof_fsr[ss][2], dp_ele_theta_ftof_fsr[ss][1], dp_ele_theta_ftof_fsr[ss][3])
    dp_res_theta_ftof_fsr_dependence.SetTitle("Electron #Delta p/p Sector " + str(ss) + " vs Theta FTOF FSR; #theta (deg); #Delta p/p")
    dp_res_theta_ftof_fsr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ftof_fsr_dependence) # save graphs in memory or root complains
    dp_res_theta_ftof_fsr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ftof_fsr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ftof_fsr_dependence.Draw("AP")
    l0 = TLine( dp_ele_theta_ftof_fsr[ss][0][0], 0.0, dp_ele_theta_ftof_fsr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')

can.Print(mon_out_file_name)

####################################################
## check the proton dependence
can.Clear()
can.Divide(2,3)
#plot delta p of proton when proton is in ctof with ISR
for ss in range(1, len(dp_pro_theta_ctof_isr)+1):
    can.cd(ss)
    #print(dp_ele_theta_ctof_isr[ss][2])
    if len(dp_pro_theta_ctof_isr[ss][2]) == 0 : continue
    dp_res_theta_ctof_isr_dependence = TGraphErrors(len(dp_pro_theta_ctof_isr[ss][0]), dp_pro_theta_ctof_isr[ss][0], dp_pro_theta_ctof_isr[ss][2], dp_pro_theta_ctof_isr[ss][1], dp_pro_theta_ctof_isr[ss][3])
    dp_res_theta_ctof_isr_dependence.SetTitle("Proton #Delta p/p Sector " + str(ss) + " vs Theta CTOF ISR; #theta (deg); #Delta p/p")
    dp_res_theta_ctof_isr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ctof_isr_dependence) # save graphs in memory or root complains
    dp_res_theta_ctof_isr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ctof_isr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ctof_isr_dependence.Draw("AP")
    l0 = TLine( dp_pro_theta_ctof_isr[ss][0][0], 0.0, dp_pro_theta_ctof_isr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,3)

#plot delta p of proton when proton is in ftof with ISR
for ss in range(1, len(dp_pro_theta_ftof_isr)+1):
    can.cd(ss)
    if len(dp_pro_theta_ftof_isr[ss][2]) == 0 : continue
    dp_res_theta_ftof_isr_dependence = TGraphErrors(len(dp_pro_theta_ftof_isr[ss][0]), dp_pro_theta_ftof_isr[ss][0], dp_pro_theta_ftof_isr[ss][2], dp_pro_theta_ftof_isr[ss][1], dp_pro_theta_ftof_isr[ss][3])
    dp_res_theta_ftof_isr_dependence.SetTitle("Proton #Delta p/p Sector " + str(ss) + " vs Theta FTOF ISR; #theta (deg); #Delta p/p")
    dp_res_theta_ftof_isr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ftof_isr_dependence) # save graphs in memory or root complains
    dp_res_theta_ftof_isr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ftof_isr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ftof_isr_dependence.Draw("AP")
    l0 = TLine( dp_pro_theta_ftof_isr[ss][0][0], 0.0, dp_pro_theta_ftof_isr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')

    

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,3)

#plot delta p of proton when proton is in ctof with FSR
# bugged need to check
for ss in range(1, len(dp_pro_theta_ctof_fsr)+1):
    can.cd(ss)
    if len(dp_pro_theta_ctof_fsr[ss][2]) == 0 : continue
    dp_res_theta_ctof_fsr_dependence = TGraphErrors(len(dp_pro_theta_ctof_fsr[ss][0]), dp_pro_theta_ctof_fsr[ss][0], dp_pro_theta_ctof_fsr[ss][2], dp_pro_theta_ctof_fsr[ss][1], dp_pro_theta_ctof_fsr[ss][3])
    dp_res_theta_ctof_fsr_dependence.SetTitle("Proton #Delta p/p Sector " + str(ss) + " vs Theta CTOF FSR; #theta (deg); #Delta p/p")
    dp_res_theta_ctof_fsr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ctof_fsr_dependence) # save graphs in memory or root complains
    dp_res_theta_ctof_fsr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ctof_fsr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ctof_fsr_dependence.Draw("AP")
    l0 = TLine( dp_pro_theta_ctof_fsr[ss][0][0], 0.0, dp_pro_theta_ctof_fsr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')


can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,3)

#plot delta p of proton when proton is in ftof with FSR
for ss in range(1, len(dp_pro_theta_ftof_fsr)+1):
    can.cd(ss)
    if len(dp_pro_theta_ftof_fsr[ss][2]) == 0: continue
    dp_res_theta_ftof_fsr_dependence = TGraphErrors(len(dp_pro_theta_ftof_fsr[ss][0]), dp_pro_theta_ftof_fsr[ss][0], dp_pro_theta_ftof_fsr[ss][2], dp_pro_theta_ftof_fsr[ss][1], dp_pro_theta_ftof_fsr[ss][3])
    dp_res_theta_ftof_fsr_dependence.SetTitle("Proton #Delta p/p Sector " + str(ss) + " vs Theta FTOF FSR; #theta (deg); #Delta p/p")
    dp_res_theta_ftof_fsr_dependence.SetMarkerStyle(21)
    store_graphs.append(dp_res_theta_ftof_fsr_dependence) # save graphs in memory or root complains
    dp_res_theta_ftof_fsr_dependence.GetHistogram().SetMaximum(0.25)
    dp_res_theta_ftof_fsr_dependence.GetHistogram().SetMinimum(-0.25)
    dp_res_theta_ftof_fsr_dependence.Draw("AP")
    l0 = TLine( dp_pro_theta_ftof_fsr[ss][0][0], 0.0, dp_pro_theta_ftof_fsr[ss][0][-1], 0.0)
    l0.SetLineColor(kRed)
    store_line.append(l0)
    l0.Draw('same')

can.Print(mon_out_file_name)



####################################################
####################################################
####################################################
## After checking the dependence of delta p/p vs theta,
## it is important to check it per theta phi bin

# delta p for electron with proton in CTOF and ISR
# fit the dp for electron per theta bin in each sector
g_dp_ele_theta_phi_ctof_isr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_ctof_isr_theta_phi[ss])

    can.Clear()
    can.Divide(8,6)
    ii=0

    temp_dp_phi_theta=[]
    
    for tb in range(0,max_tb):
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')

        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range
        print(' number of phi bins for ele ctof isr ' + str( len(dp_ele_ctof_isr_theta_phi[ss][tb]) ) )
        for pb in range(0, len(dp_ele_ctof_isr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits = getFitLimits(dp_ele_ctof_isr_theta_phi[ss][tb][pb],2.5)
            dp_ele_ctof_isr_theta_phi[ss][tb][pb].RebinX(6)
            dp_fit = getFitFractionalHeight(dp_ele_ctof_isr_theta_phi[ss][tb][pb],'dp_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb), 0.4, store_fit)
            #TF1('dp_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'gaus',temp_fit_limits[0], temp_fit_limits[1])
            #dp_fit.SetParameter(1,temp_fit_limits[2])
            #dp_fit.SetParameter(2,temp_fit_limits[3])
            #print(dp_ele_ctof_isr_theta_phi[ss][tb])
            
            if len(dp_ele_ctof_isr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_ele_ctof_isr_theta_phi[ss][tb][pb].SetTitle('Electron #Delta P/P_{rec} CTOF ISR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_ele_ctof_isr_theta_phi[ss][tb][pb].Draw()
            #dp_ele_ctof_isr_theta_phi[ss][tb][pb].Fit('dp_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
    
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        
        temp_dp_phi_theta.append(temp_dp_phi)
        #print(temp_dp_phi_theta)
    g_dp_ele_theta_phi_ctof_isr[ss] = temp_dp_phi_theta
    
    can.Print(mon_out_file_name)


# delta p for electron with proton in CTOF and FSR
# fit the dp for electron per theta bin in each sector
g_dp_ele_theta_phi_ctof_fsr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_ctof_fsr_theta_phi[ss])

    can.Clear()
    can.Divide(8,6)
    ii=0

    temp_dp_phi_theta=[]

    for tb in range(0,max_tb):        
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')


        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range
        for pb in range(0, len(dp_ele_ctof_fsr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits=getFitLimits(dp_ele_ctof_fsr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'landau',temp_fit_limits[0], temp_fit_limits[1])  #gaus(0)*pol2(3)',temp_fit_limits[0], temp_fit_limits[1])
                      
            if len(dp_ele_ctof_fsr_theta_phi[ss][tb][pb]) == 0 : 
                phi_max=phi_low
                phi_low-=phi_range
                continue

            dp_ele_ctof_fsr_theta_phi[ss][tb][pb].SetTitle('Electron #Delta P/P_{rec} CTOF FSR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_ele_ctof_fsr_theta_phi[ss][tb][pb].RebinX(6)
            dp_ele_ctof_fsr_theta_phi[ss][tb][pb].Draw()
            dp_ele_ctof_fsr_theta_phi[ss][tb][pb].Fit('dp_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))
            phi_max=phi_low
            phi_low-=phi_range

            ii+=1
        
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        temp_dp_phi_theta.append(temp_dp_phi)
    
    g_dp_ele_theta_phi_ctof_fsr[ss] = temp_dp_phi_theta

    can.Print(mon_out_file_name)

# delta p for electron with proton in FTOF and ISR
# fit the dp for electron per theta bin in each sector
g_dp_ele_theta_phi_ftof_isr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_ftof_isr_theta_phi[ss])
    
    can.Clear()
    can.Divide(8,6)
    ii=0
    
    temp_dp_phi_theta=[]
    
    for tb in range(0,max_tb):
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')

        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range
        for pb in range(0, len(dp_ele_ftof_isr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits = getFitLimits(dp_ele_ftof_isr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'gaus',temp_fit_limits[0], temp_fit_limits[1])
            dp_fit.SetParameter(1,temp_fit_limits[2])
            dp_fit.SetParameter(2,temp_fit_limits[3])
            
            if len(dp_ele_ftof_isr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_ele_ftof_isr_theta_phi[ss][tb][pb].SetTitle('Electron #Delta P/P_{rec} FTOF ISR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_ele_ftof_isr_theta_phi[ss][tb][pb].RebinX(4)
            dp_ele_ftof_isr_theta_phi[ss][tb][pb].Draw()
            dp_ele_ftof_isr_theta_phi[ss][tb][pb].Fit('dp_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
        
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        temp_dp_phi_theta.append(temp_dp_phi)
        
    g_dp_ele_theta_phi_ftof_isr[ss] = temp_dp_phi_theta
    can.Print(mon_out_file_name)

# delta p for electron with proton in FTOF and FSR
# fit the dp for electron per theta bin in each sector
g_dp_ele_theta_phi_ftof_fsr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_ele_ftof_fsr_theta_phi[ss])

    can.Clear()
    can.Divide(8,6)
    ii=0

    temp_dp_phi_theta=[]
    
    for tb in range(0,max_tb):
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')

        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range

        for pb in range(0, len(dp_ele_ftof_fsr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits=getFitLimits(dp_ele_ftof_fsr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'landau',temp_fit_limits[0], temp_fit_limits[1])
            #dp_fit.SetParameter(1,temp_fit_limits[2])
            #dp_fit.SetParameter(2,temp_fit_limits[3])

            if len(dp_ele_ftof_fsr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_ele_ftof_fsr_theta_phi[ss][tb][pb].SetTitle('Electron #Delta P/P_{rec} FTOF FSR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_ele_ftof_fsr_theta_phi[ss][tb][pb].RebinX(4)
            dp_ele_ftof_fsr_theta_phi[ss][tb][pb].Draw()
            dp_ele_ftof_fsr_theta_phi[ss][tb][pb].Fit('dp_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
        
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        temp_dp_phi_theta.append(temp_dp_phi)

    g_dp_ele_theta_phi_ftof_fsr[ss] = temp_dp_phi_theta
    
    can.Print(mon_out_file_name)


##########################################################################################
##########################################################################################
# plot the graph delta p/p for electron for each case per sector
# i.e. use the lists in g_dp_ele_theta_phi_*
l_mg_dp_res_phi_ele_ctof_isr = []
le0 = []
for ss in range(1,len(g_dp_ele_theta_phi_ctof_isr)+1):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,9)
    ii=1

    mg_dp_res_phi_ele_ctof_isr = TMultiGraph()

    for tb in range(1, len(g_dp_ele_theta_phi_ctof_isr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_ele_theta_phi_ctof_isr[ss][tb]
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_ele_ctof_isr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_ele_ctof_isr_dep.SetTitle("Electron #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " CTOF ISR ; #phi (deg); #Delta p/p ")
        dp_res_phi_ele_ctof_isr_dep.SetMarkerStyle(21)
        dp_res_phi_ele_ctof_isr_dep.SetMarkerColor(kSpring+tb)
        store_graphs.append(dp_res_phi_ele_ctof_isr_dep)
        dp_res_phi_ele_ctof_isr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_ele_ctof_isr_dep.GetHistogram().SetMinimum(-0.15)
        mg_dp_res_phi_ele_ctof_isr.Add(dp_res_phi_ele_ctof_isr_dep)        
        dp_res_phi_ele_ctof_isr_dep.Draw("AP")
        l0 = TLine( dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        le0.append(l0)
        l0.Draw('same')
        store_line.append(l0)

        ii+=1
    
    store_graphs.append(mg_dp_res_phi_ele_ctof_isr)
    can.cd(5)
    mg_dp_res_phi_ele_ctof_isr.SetTitle("Electron Sector " + str(ss) + "#Delta P/P per Theta Bin vs #phi CTOF ISR; #phi (deg); #Delta p/p" )
    l_mg_dp_res_phi_ele_ctof_isr.append(mg_dp_res_phi_ele_ctof_isr)

    can.Print(mon_out_file_name)

# delta p / p for proton in ftof with isr 
l_mg_dp_res_phi_ele_ftof_isr = []
le1=[]
for ss in range(1,len(g_dp_ele_theta_phi_ftof_isr)+1 ):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,7)
    ii=1
    mg_dp_res_phi_ele_ftof_isr = TMultiGraph()

    for tb in range(1, len(g_dp_ele_theta_phi_ftof_isr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_ele_theta_phi_ftof_isr[ss][tb]
        print(dp_phi_dep[0])
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_ele_ftof_isr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_ele_ftof_isr_dep.SetTitle("Electron #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " FTOF ISR ; #phi (deg); #Delta p/p ")
        dp_res_phi_ele_ftof_isr_dep.SetMarkerStyle(21)
        dp_res_phi_ele_ftof_isr_dep.SetMarkerColor(kSpring+tb)
        
        store_graphs.append(dp_res_phi_ele_ftof_isr_dep)
        mg_dp_res_phi_ele_ftof_isr.Add(dp_res_phi_ele_ftof_isr_dep)
        dp_res_phi_ele_ftof_isr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_ele_ftof_isr_dep.GetHistogram().SetMinimum(-0.15)
        dp_res_phi_ele_ftof_isr_dep.Draw("AP")
        l0 = TLine( dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        le1.append(l0)
        l0.Draw('same')
        store_line.append(l0)

        ii+=1
    store_graphs.append(mg_dp_res_phi_ele_ftof_isr)
    can.cd(5)
    mg_dp_res_phi_ele_ftof_isr.SetTitle("Sector " + str(ss) + "#Delta P/P per Theta Bin vs #phi FTOF ISR; #phi (deg); #Delta p/p" )
    l_mg_dp_res_phi_ele_ftof_isr.append(mg_dp_res_phi_ele_ftof_isr)
    
    can.Print(mon_out_file_name)

# delta p / p for proton in ctof with fsr 
l_mg_dp_res_phi_ele_ctof_fsr = []
le2=[]
for ss in range(1,len(g_dp_ele_theta_phi_ctof_fsr)+1 ):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,4)
    ii=1
    mg_dp_res_phi_ele_ctof_fsr = TMultiGraph()
    for tb in range(1, len(g_dp_ele_theta_phi_ctof_fsr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_ele_theta_phi_ctof_fsr[ss][tb]
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_ele_ctof_fsr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_ele_ctof_fsr_dep.SetTitle("Electron #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " CTOF FSR ; #phi (deg); #Delta p/p ")
        dp_res_phi_ele_ctof_fsr_dep.SetMarkerStyle(21)
        dp_res_phi_ele_ctof_fsr_dep.SetMarkerColor(kSpring+tb)

        store_graphs.append(dp_res_phi_ele_ctof_fsr_dep)
        dp_res_phi_ele_ctof_fsr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_ele_ctof_fsr_dep.GetHistogram().SetMinimum(-0.15)
        dp_res_phi_ele_ctof_fsr_dep.Draw("AP")
        l0 = TLine( dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        le2.append(l0)
        l0.Draw('same')
        store_line.append(l0)


        ii+=1
    store_graphs.append(mg_dp_res_phi_ele_ctof_fsr)
    can.cd(5)
    mg_dp_res_phi_ele_ctof_fsr.SetTitle("Electron S " + str(ss) + "#Delta P/P per Theta Bin vs #phi CTOF FSR; #phi (deg); #Delta p/p" )
    l_mg_dp_res_phi_ele_ctof_fsr.append(mg_dp_res_phi_ele_ctof_fsr)

    can.Print(mon_out_file_name)


# delta p / p for proton in ftof with fsr 

l_mg_dp_res_phi_ele_ftof_fsr = []
le3=[]
for ss in range(1,len(g_dp_ele_theta_phi_ftof_fsr)+1 ):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,6)
    ii=1
    mg_dp_res_phi_ele_ftof_fsr = TMultiGraph()
    for tb in range(1, len(g_dp_ele_theta_phi_ftof_fsr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_ele_theta_phi_ftof_fsr[ss][tb]
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_ele_ftof_fsr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_ele_ftof_fsr_dep.SetTitle("Electron #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " FTOF FSR ; #phi (deg); #Delta p/p ")
        dp_res_phi_ele_ftof_fsr_dep.SetMarkerStyle(21)
        dp_res_phi_ele_ftof_fsr_dep.SetMarkerColor(kSpring+tb)
        store_graphs.append(dp_res_phi_ele_ftof_fsr_dep)
        mg_dp_res_phi_ele_ftof_fsr.Add(dp_res_phi_ele_ftof_fsr_dep)

        dp_res_phi_ele_ftof_isr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_ele_ftof_isr_dep.GetHistogram().SetMinimum(-0.15)
        dp_res_phi_ele_ftof_fsr_dep.Draw("AP")
        l0 = TLine( dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        l0.Draw('same')
        le3.append(l0)
        store_line.append(l0)


        ii+=1
    store_graphs.append(mg_dp_res_phi_ele_ftof_fsr)
    can.cd(5)
    mg_dp_res_phi_ele_ftof_fsr.SetTitle("Sector " + str(ss) + "#Delta P/P per Theta Bin vs #phi FTOF FSR; #phi (deg) ; #Delta p/p" )
    l_mg_dp_res_phi_ele_ftof_fsr.append(mg_dp_res_phi_ele_ftof_fsr)
    can.Print(mon_out_file_name)


##########################################################################33
## check the proton delta p/ p theta and phi dependence here
g_dp_pro_theta_phi_ctof_isr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_ctof_isr_theta_phi[ss])

    can.Clear()
    can.Divide(8,6)
    ii=0

    temp_dp_phi_theta=[]
    
    for tb in range(0,max_tb):
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')

        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range
        print(' number of phi bins for pro ctof isr ' + str( len(dp_pro_ctof_isr_theta_phi[ss][tb]) ) )
        for pb in range(0, len(dp_pro_ctof_isr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits = getFitLimits(dp_pro_ctof_isr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'gaus',temp_fit_limits[0], temp_fit_limits[1])
            dp_fit.SetParameter(1,temp_fit_limits[2])
            dp_fit.SetParameter(2,temp_fit_limits[3])
            #print(dp_pro_ctof_isr_theta_phi[ss][tb])
            
            if len(dp_pro_ctof_isr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_pro_ctof_isr_theta_phi[ss][tb][pb].RebinX(4)
            dp_pro_ctof_isr_theta_phi[ss][tb][pb].SetTitle('Proton #Delta P/P_{rec} CTOF ISR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_pro_ctof_isr_theta_phi[ss][tb][pb].Draw()
            dp_pro_ctof_isr_theta_phi[ss][tb][pb].Fit('dp_ctof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
    
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        
        temp_dp_phi_theta.append(temp_dp_phi)
        #print(temp_dp_phi_theta)
    g_dp_pro_theta_phi_ctof_isr[ss] = temp_dp_phi_theta
    
    can.Print(mon_out_file_name)


# delta p for proton with proton in CTOF and FSR
# fit the dp for proton per theta bin in each sector
g_dp_pro_theta_phi_ctof_fsr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_ctof_fsr_theta_phi[ss])

    can.Clear()
    can.Divide(8,6)
    ii=0

    temp_dp_phi_theta=[]

    for tb in range(0,max_tb):        
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')


        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range
        for pb in range(0, len(dp_pro_ctof_fsr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits=getFitLimits(dp_pro_ctof_fsr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'landau',temp_fit_limits[0], temp_fit_limits[1])  #gaus(0)*pol2(3)',temp_fit_limits[0], temp_fit_limits[1])
                      
            if len(dp_pro_ctof_fsr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_pro_ctof_fsr_theta_phi[ss][tb][pb].SetTitle('Proton #Delta P/P_{rec} CTOF FSR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_pro_ctof_fsr_theta_phi[ss][tb][pb].RebinX(4)
            dp_pro_ctof_fsr_theta_phi[ss][tb][pb].Draw()
            dp_pro_ctof_fsr_theta_phi[ss][tb][pb].Fit('dp_ctof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
        
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        temp_dp_phi_theta.append(temp_dp_phi)
    
    g_dp_pro_theta_phi_ctof_fsr[ss] = temp_dp_phi_theta

    can.Print(mon_out_file_name)

# delta p for proton with proton in FTOF and ISR
# fit the dp for proton per theta bin in each sector
g_dp_pro_theta_phi_ftof_isr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_ftof_isr_theta_phi[ss])
    
    can.Clear()
    can.Divide(8,6)
    ii=0
    
    temp_dp_phi_theta=[]
    
    for tb in range(0,max_tb):
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')

        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range
        for pb in range(0, len(dp_pro_ftof_isr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits = getFitLimits(dp_pro_ftof_isr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'gaus',temp_fit_limits[0], temp_fit_limits[1])
            dp_fit.SetParameter(1,temp_fit_limits[2])
            dp_fit.SetParameter(2,temp_fit_limits[3])
            
            if len(dp_pro_ftof_isr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_pro_ftof_isr_theta_phi[ss][tb][pb].SetTitle('Proton #Delta P/P_{rec} FTOF ISR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_pro_ftof_isr_theta_phi[ss][tb][pb].RebinX(4)
            dp_pro_ftof_isr_theta_phi[ss][tb][pb].Draw()
            dp_pro_ftof_isr_theta_phi[ss][tb][pb].Fit('dp_ftof_isr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
        
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        temp_dp_phi_theta.append(temp_dp_phi)
        
    g_dp_pro_theta_phi_ftof_isr[ss] = temp_dp_phi_theta
    can.Print(mon_out_file_name)

# delta p for proton with proton in FTOF and FSR
# fit the dp for proton per theta bin in each sector
g_dp_pro_theta_phi_ftof_fsr = {}
for ss in range(1,7):
    delta_theta = 7
    theta_min = 0
    theta_max = delta_theta
    max_tb=len(dp_pro_ftof_fsr_theta_phi[ss])

    can.Clear()
    can.Divide(8,6)
    ii=0

    temp_dp_phi_theta=[]
    
    for tb in range(0,max_tb):
        dp_phi_x, dp_phi_errx = array('d'), array('d')
        dp_fit_y, dp_fit_erry = array('d'), array('d')

        phi_max = 35.0
        phi_range = 5
        phi_low = phi_max - phi_range

        for pb in range(0, len(dp_pro_ftof_fsr_theta_phi[ss][tb])):
            can.cd(ii+1)

            temp_fit_limits=getFitLimits(dp_pro_ftof_fsr_theta_phi[ss][tb][pb],2.5)
            dp_fit = TF1('dp_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'landau',temp_fit_limits[0], temp_fit_limits[1])
            #dp_fit.SetParameter(1,temp_fit_limits[2])
            #dp_fit.SetParameter(2,temp_fit_limits[3])

            if len(dp_pro_ftof_fsr_theta_phi[ss][tb][pb]) == 0 : continue
            dp_pro_ftof_fsr_theta_phi[ss][tb][pb].SetTitle('Proton #Delta P/P_{rec} FTOF FSR Sector '+str(ss)+' Theta Bin ' + str(tb) + ' Phi Bin ' + str(pb) +' ; #Delta P/P_{rec}; counts')
            dp_pro_ftof_fsr_theta_phi[ss][tb][pb].RebinX(4)
            dp_pro_ftof_fsr_theta_phi[ss][tb][pb].Draw()
            dp_pro_ftof_fsr_theta_phi[ss][tb][pb].Fit('dp_ftof_fsr_fit_s'+str(ss)+'_tb'+str(tb)+'_pb'+str(pb),'R')        
            dp_phi_x.append( phi_low + (phi_max - phi_low)/2.0)
            dp_phi_errx.append(0.0)
            dp_fit_y.append(dp_fit.GetParameter(1))
            dp_fit_erry.append(dp_fit.GetParameter(2))

            ii+=1
            phi_max=phi_low
            phi_low-=phi_range
        
        temp_dp_phi = []
        temp_dp_phi.append(dp_phi_x)
        temp_dp_phi.append(dp_phi_errx)
        temp_dp_phi.append(dp_fit_y)
        temp_dp_phi.append(dp_fit_erry)
        temp_dp_phi_theta.append(temp_dp_phi)

    g_dp_pro_theta_phi_ftof_fsr[ss] = temp_dp_phi_theta
    
    can.Print(mon_out_file_name)



##########################################################################################
##########################################################################################
# plot the graph delta p/p for PROTON for each case per sector
# i.e. use the lists in g_dp_pro_theta_phi_*
l_mg_dp_res_phi_pro_ctof_isr = []
l_temp0=[]
for ss in range(1,len(g_dp_pro_theta_phi_ctof_isr)+1):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,9)
    ii=1

    mg_dp_res_phi_pro_ctof_isr = TMultiGraph()

    for tb in range(1, len(g_dp_pro_theta_phi_ctof_isr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_pro_theta_phi_ctof_isr[ss][tb]
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_pro_ctof_isr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_pro_ctof_isr_dep.SetTitle("Proton #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " CTOF ISR ; #phi (deg); #Delta p/p ")
        dp_res_phi_pro_ctof_isr_dep.SetMarkerStyle(21)
        dp_res_phi_pro_ctof_isr_dep.SetMarkerColor(kSpring+tb)
        store_graphs.append(dp_res_phi_pro_ctof_isr_dep)
        dp_res_phi_pro_ctof_isr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_pro_ctof_isr_dep.GetHistogram().SetMinimum(-0.15)
        mg_dp_res_phi_pro_ctof_isr.Add(dp_res_phi_pro_ctof_isr_dep)        
        l0=TLine(dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        l_temp0.append(l0)
        store_line.append(l0)

        dp_res_phi_pro_ctof_isr_dep.Draw("AP")
        ii+=1
        
    store_graphs.append(mg_dp_res_phi_pro_ctof_isr)
    can.cd(5)
    mg_dp_res_phi_pro_ctof_isr.SetTitle("Proton S" + str(ss) + " #Delta P/P per Theta Bin vs #phi CTOF ISR; #phi (deg); #Delta p/p" )
    l_mg_dp_res_phi_pro_ctof_isr.append(mg_dp_res_phi_pro_ctof_isr)

    can.Print(mon_out_file_name)

# delta p / p for proton in ftof with isr 
l_mg_dp_res_phi_pro_ftof_isr = []
l_temp1=[]
for ss in range(1,len(g_dp_pro_theta_phi_ftof_isr)+1 ):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,7)
    ii=1
    mg_dp_res_phi_pro_ftof_isr = TMultiGraph()

    for tb in range(1, len(g_dp_pro_theta_phi_ftof_isr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_pro_theta_phi_ftof_isr[ss][tb]
        print(dp_phi_dep[0])
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_pro_ftof_isr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_pro_ftof_isr_dep.SetTitle("Proton #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " FTOF ISR ; #phi (deg); #Delta p/p ")
        dp_res_phi_pro_ftof_isr_dep.SetMarkerStyle(21)
        dp_res_phi_pro_ftof_isr_dep.SetMarkerColor(kSpring+tb)
        
        store_graphs.append(dp_res_phi_pro_ftof_isr_dep)
        mg_dp_res_phi_pro_ftof_isr.Add(dp_res_phi_pro_ftof_isr_dep)
        dp_res_phi_pro_ftof_isr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_pro_ftof_isr_dep.GetHistogram().SetMinimum(-0.15)
        dp_res_phi_pro_ftof_isr_dep.Draw("AP")
        l0=TLine(dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        l_temp1.append(l0)
        store_line.append(l0)


        ii+=1
    store_graphs.append(mg_dp_res_phi_pro_ftof_isr)
    can.cd(5)
    mg_dp_res_phi_pro_ftof_isr.SetTitle("Proton S" + str(ss) + " #Delta P/P per Theta Bin vs #phi FTOF ISR; #phi (deg); #Delta p/p" )
    l_mg_dp_res_phi_pro_ftof_isr.append(mg_dp_res_phi_pro_ftof_isr)
    
    can.Print(mon_out_file_name)

# delta p / p for proton in ctof with fsr 
l_mg_dp_res_phi_pro_ctof_fsr = []
l_temp2=[]
for ss in range(1,len(g_dp_pro_theta_phi_ctof_fsr)+1 ):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,4)
    ii=1
    mg_dp_res_phi_pro_ctof_fsr = TMultiGraph()
    for tb in range(1, len(g_dp_pro_theta_phi_ctof_fsr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_pro_theta_phi_ctof_fsr[ss][tb]
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_pro_ctof_fsr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_pro_ctof_fsr_dep.SetTitle("Proton #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " CTOF FSR ; #phi (deg); #Delta p/p ")
        dp_res_phi_pro_ctof_fsr_dep.SetMarkerStyle(21)
        dp_res_phi_pro_ctof_fsr_dep.SetMarkerColor(kSpring+tb)

        store_graphs.append(dp_res_phi_pro_ctof_fsr_dep)
        dp_res_phi_pro_ctof_fsr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_pro_ctof_fsr_dep.GetHistogram().SetMinimum(-0.15)
        dp_res_phi_pro_ctof_fsr_dep.Draw("AP")
        l0=TLine(dp_phi_dep[0][0], 0.0, dp_phi_dep[0][-1], 0.0)
        l0.SetLineColor(kRed)
        store_line.append(l0)

        l_temp2.append(l0)

        ii+=1
    store_graphs.append(mg_dp_res_phi_pro_ctof_fsr)
    can.cd(5)
    mg_dp_res_phi_pro_ctof_fsr.SetTitle("Proton S " + str(ss) + " #Delta P/P per Theta Bin vs #phi CTOF FSR; #phi (deg); #Delta p/p" )
    l_mg_dp_res_phi_pro_ctof_fsr.append(mg_dp_res_phi_pro_ctof_fsr)

    can.Print(mon_out_file_name)


# delta p / p for proton in ftof with fsr 

l_mg_dp_res_phi_pro_ftof_fsr = []
for ss in range(1,len(g_dp_pro_theta_phi_ftof_fsr)+1 ):
    can.Clear()
    #can.SetWindowSize(600,2400)
    can.Divide(1,6)
    ii=1
    mg_dp_res_phi_pro_ftof_fsr = TMultiGraph()
    for tb in range(1, len(g_dp_pro_theta_phi_ftof_fsr[ss]) ):
        can.cd(ii)
        print(ii)
        dp_phi_dep = g_dp_pro_theta_phi_ftof_fsr[ss][tb]
        if len(dp_phi_dep[0]) == 0 : continue 
        dp_res_phi_pro_ftof_fsr_dep = TGraphErrors( len(dp_phi_dep[0]), dp_phi_dep[0], dp_phi_dep[2], dp_phi_dep[1], dp_phi_dep[3])
        dp_res_phi_pro_ftof_fsr_dep.SetTitle("Proton #Delta p/p vs #Phi Sector " + str(ss) + " Theta Bin " + str(tb) + " FTOF FSR ; #phi (deg); #Delta p/p ")
        dp_res_phi_pro_ftof_fsr_dep.SetMarkerStyle(21)
        dp_res_phi_pro_ftof_fsr_dep.SetMarkerColor(kSpring+tb)
        store_graphs.append(dp_res_phi_pro_ftof_fsr_dep)
        mg_dp_res_phi_pro_ftof_fsr.Add(dp_res_phi_pro_ftof_fsr_dep)

        dp_res_phi_pro_ftof_isr_dep.GetHistogram().SetMaximum(0.15)
        dp_res_phi_pro_ftof_isr_dep.GetHistogram().SetMinimum(-0.15)
        dp_res_phi_pro_ftof_fsr_dep.Draw("AP")
        ii+=1
    store_graphs.append(mg_dp_res_phi_pro_ftof_fsr)
    can.cd(5)
    mg_dp_res_phi_pro_ftof_fsr.SetTitle("Proton S " + str(ss) + " #Delta P/P per Theta Bin vs #phi FTOF FSR; #phi (deg) ; #Delta p/p" )
    l_mg_dp_res_phi_pro_ftof_fsr.append(mg_dp_res_phi_pro_ftof_fsr)
    can.Print(mon_out_file_name)



########################################################################
########################################################################
## overlay the deltap/p per theta bin vs phi in each sector

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_ele_ctof_isr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_ele_ctof_isr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_ele_ctof_isr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_ele_ctof_isr[ss].Draw("AP")
    le0[ss].Draw('same')
    ii+=1
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_ele_ftof_isr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_ele_ftof_isr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_ele_ftof_isr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_ele_ftof_isr[ss].Draw("AP")
    le1[ss].Draw('same')
    ii+=1
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_ele_ctof_fsr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_ele_ctof_fsr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_ele_ctof_fsr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_ele_ctof_fsr[ss].Draw("AP")
    le2[ss].Draw('same')
    ii+=1
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_ele_ftof_fsr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_ele_ftof_fsr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_ele_ftof_fsr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_ele_ftof_fsr[ss].Draw("AP")
    le3[ss].Draw('same')
    ii+=1
can.Print(mon_out_file_name)

#################################################
## check the same thing but for protons
can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_pro_ctof_isr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_pro_ctof_isr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_pro_ctof_isr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_pro_ctof_isr[ss].Draw("AP")
    ii+=1
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_pro_ftof_isr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_pro_ftof_isr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_pro_ftof_isr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_pro_ftof_isr[ss].Draw("AP")
    ii+=1
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_pro_ctof_fsr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_pro_ctof_fsr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_pro_ctof_fsr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_pro_ctof_fsr[ss].Draw("AP")
    ii+=1
can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(0, len(l_mg_dp_res_phi_pro_ftof_fsr) ):
    can.cd(ii+1)
    l_mg_dp_res_phi_pro_ftof_fsr[ss].GetHistogram().SetMaximum(0.60)
    l_mg_dp_res_phi_pro_ftof_fsr[ss].GetHistogram().SetMinimum(-0.60)
    l_mg_dp_res_phi_pro_ftof_fsr[ss].Draw("AP")
    ii+=1
can.Print(mon_out_file_name)



can.Print('{}]'.format(mon_out_file_name))


