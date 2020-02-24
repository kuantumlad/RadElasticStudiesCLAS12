from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLine, TLatex
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kSpring
import numpy as np

from array import array
import math
import sys,os

# script to plot, view, correct electron and proton momentum
# in each theta and phi bin

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)


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

def fit_slices(histo, x_range, x_bin_step, y_fit_range):

    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)

        fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus')
        fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
        fit.SetName(histo.GetTitle() + '_fit{}'.format(i))

        if y_fit_range:
            fit.SetParameter(1, 0.5 * (y_fit_range[0] + y_fit_range[1]))
 
        projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])

        slices.append(projec)
        fits.append(fit)

        x_low = histo.GetXaxis().GetBinCenter(x_bin)
        x_high = histo.GetXaxis().GetBinCenter(x_bin + x_bin_step)
        x_values.append(0.5 * (x_high + x_low))

    means = array('d')
    means_err = array('d')
    stds = array('d')
    stds_err = array('d')
    zeros = array('d')

    for f in fits:
        means.append(f.GetParameter(1))
        means_err.append(f.GetParError(1))
        stds.append(f.GetParameter(2))
        stds_err.append(f.GetParError(1)) ### change if want the sigma of dist.
        zeros.append(0.0)

    graph = TGraphErrors(len(x_values), x_values, means, zeros, stds_err)
    graph.SetName('g_' + histo.GetName())
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(kRed)

    return graph, slices, fits


def plot_fits(canvas, histos, x_range, x_bin_step, title_formatter,
              save_name, label, y_fit_range, y_range=None,
              title=None, xtitle=None, ytitle=None, hline=None):

    canvas.Clear()
    canvas.Divide(2,3)

    root_is_dumb = []
    for i in range(1,7):#len(histos)):
        canvas.cd(i)

        title = title_formatter.format(i)
        graph, slices, fits = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range=y_fit_range)
        graph.SetMarkerStyle(8)
        graph.SetMarkerSize(1)
        
        if y_range:
            graph.GetHistogram().SetMinimum(y_range[0])
            graph.GetHistogram().SetMaximum(y_range[1])
            graph.Draw('AP')
            root_is_dumb.append(graph)
        else:
            graph.Draw('AP')
            root_is_dumb.append(graph)
            
        #if hline:
        line = TLine(x_range[0], hline, x_range[1], hline)
        #line.SetLineStyle(8)
        line.SetLineWidth(2)
        line.SetLineColor(kRed)
        line.Draw()
        root_is_dumb.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

    # For the slices 
    slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
    slice_pdfname = title_formatter.split('_{}')[0] + '_slices.pdf'
    slice_can.Print(slice_pdfname + '[')
    for i in range(1,7):
        title = title_formatter.format(i)
        graph, slices, fits = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range)

        # Size of slices page
        nrows = 5
        ncols = int(np.ceil(len(slices) / nrows) + 1)
        slice_can.Clear() 
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices, fits)):
            slice_can.cd(j+1)
            s.Draw()
            lab.DrawLatex(0.15, 0.88, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(f.GetParameter(1), f.GetParameter(2)))
            
        slice_can.Print(slice_pdfname)
    slice_can.Print(slice_pdfname + ']')



# load data and file
datatype = sys.argv[2]
hhs={}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj


base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_kin_corr2d_'+datatype+'.pdf'
can = TCanvas('can','can')
can.SetCanvasSize(1200,1200)
can.Print('{}['.format(mon_out_file_name))


lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)

gStyle.SetOptStat(00000)


######################
## check 2D Histograms
## delta p vs p 

p_ele_dp_ele_from_ang_ctof_isr = {}
p_ele_dp_ele_from_ang_ftof_isr = {}
p_ele_dp_ele_from_ang_ctof_fsr = {}
p_ele_dp_ele_from_ang_ftof_fsr = {}

p_pro_dp_pro_from_ang_ctof_isr = {}
p_pro_dp_pro_from_ang_ftof_isr = {}
p_pro_dp_pro_from_ang_ctof_fsr = {}
p_pro_dp_pro_from_ang_ftof_fsr = {}


for ss in range(0,7):
    ################################################333
    ## get 2D histograms
    if 'p_ele_dp_ele_from_angles_CTOF_'+str(ss)+'_isr' in hhs:
        p_ele_dp_ele_from_ang_ctof_isr[ss] = hhs['p_ele_dp_ele_from_angles_CTOF_'+str(ss)+'_isr']
    if 'p_ele_dp_ele_from_angles_FTOF_'+str(ss)+'_isr' in hhs:
        p_ele_dp_ele_from_ang_ftof_isr[ss] = hhs['p_ele_dp_ele_from_angles_FTOF_'+str(ss)+'_isr']
    if 'p_ele_dp_ele_from_angles_CTOF_'+str(ss)+'_fsr' in hhs:
        p_ele_dp_ele_from_ang_ctof_fsr[ss] = hhs['p_ele_dp_ele_from_angles_CTOF_'+str(ss)+'_fsr']
    if 'p_ele_dp_ele_from_angles_FTOF_'+str(ss)+'_fsr' in hhs:
        p_ele_dp_ele_from_ang_ftof_fsr[ss] = hhs['p_ele_dp_ele_from_angles_FTOF_'+str(ss)+'_fsr']

    if 'p_pro_dp_pro_from_angles_CTOF_'+str(ss)+'_isr' in hhs:
        p_pro_dp_pro_from_ang_ctof_isr[ss] = hhs['p_pro_dp_pro_from_angles_CTOF_'+str(ss)+'_isr']
    if 'p_pro_dp_pro_from_angles_FTOF_'+str(ss)+'_isr' in hhs:
        print (' here ' )
        p_pro_dp_pro_from_ang_ftof_isr[ss] = hhs['p_pro_dp_pro_from_angles_FTOF_'+str(ss)+'_isr']
    if 'p_pro_dp_pro_from_angles_CTOF_'+str(ss)+'_fsr' in hhs:
        print( ' pro ctof here ' )
        p_pro_dp_pro_from_ang_ctof_fsr[ss] = hhs['p_pro_dp_pro_from_angles_CTOF_'+str(ss)+'_fsr']
    if 'p_pro_dp_pro_from_angles_FTOF_'+str(ss)+'_fsr' in hhs:
        p_pro_dp_pro_from_ang_ftof_fsr[ss] = hhs['p_pro_dp_pro_from_angles_FTOF_'+str(ss)+'_fsr']
        

################################
## info from electron angles
can.Clear()
can.Divide(2,3)


fit_range = [ [-0.4, 0.4],
              [-0.4, -0.05]
            ]

store_fits =[]
store_graphs = []
l_temp0 = []
g_temp0=[]


can.Clear()
can.Divide(2,3)
ii=0
for ss in range(1,7):
    can.cd(ii+1)
    p_ele_dp_ele_from_ang_ftof_isr[ss].SetTitle('El. P_{e} from #Theta_{e}#Theta_{p}, #Delta P_{e} vs P_{e} Pr (FTOF) S'+str(ss)+ ', #gamma (ISR) ; P_{e} (GeV); #Delta P_{e}')    
    p_ele_dp_ele_from_ang_ftof_isr[ss].RebinX(2)
    p_ele_dp_ele_from_ang_ftof_isr[ss].RebinY(2)
    p_ele_dp_ele_from_ang_ftof_isr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)

g_temp1=[]
l_temp1=[]
for ss in range(1,7):
    can.Clear()
    can.Divide(10,10)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_ele_dp_ele_from_ang_ftof_isr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_ele_dp_ele_from_ang_ftof_isr[ss].ProjectionY('profY_p_ele_dp_ele_from_ang_ftof_isr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_ele_dp_ele_from_ang_ftof_isr_s'+str(ss)+'_tb'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_ele_dp_ele_from_ang_ftof_isr_s'+str(ss)+'_tb'+str(xx),'R')
        h_temp.Draw()                
        temp_x.append(p_ele_dp_ele_from_ang_ftof_isr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)

    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('El. P_{e} from #Theta_{e}#Theta_{p} #Delta P_{e} vs P_{e} Pr (FTOF) S'+str(ss)+ ', #gamma (ISR) ; P_{e} (GeV); #Delta P_{e}')    
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    g_temp1.append(g_temp)
    l0 = TLine(temp_x[0], 0.0, temp_x[-1], 0.0)
    l0.SetLineColor(kRed)
    l_temp1.append(l0)

    store_graphs.append(g_temp)


can.Clear()
can.Divide(2,3)
for ss in range(0, len(g_temp1)):
    can.cd(ss+1)
    g_temp1[ss].Draw('AP')
    l_temp1[ss].Draw('same')
can.Print(mon_out_file_name)

ii=0
can.Clear()
can.Divide(2,3)
for ss in range(1,7):
    can.cd(ii+1)
    p_ele_dp_ele_from_ang_ctof_isr[ss].Draw("colz")    
    p_ele_dp_ele_from_ang_ctof_isr[ss].SetTitle("El. from #Theta_{e}#Theta_{p} #Delta P_{e} vs P_{e} Pr (CTOF) " + str(ss) +", #gamma (ISR) ; P_{e} (GeV); #Delta P_{e}")    
    ii+=1
can.Print(mon_out_file_name)

for ss in range(1,7):
    can.Clear()
    can.Divide(10,10)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_ele_dp_ele_from_ang_ctof_isr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_ele_dp_ele_from_ang_ctof_isr[ss].ProjectionY('profY_p_ele_dp_ele_from_ang_ctof_isr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = getFitFractionalHeight(h_temp, 'p_ele_dp_ele_from_ang_ctof_isr_'+str(xx), 0.4, store_fits)
        #TF1('p_ele_dp_ele_from_ang_ctof_isr_'+str(xx),'gaus',-0.4,0.4)                  
        #h_temp.Fit('p_ele_dp_ele_from_ang_ctof_isr_'+str(xx),'R')
        #h_temp.Draw()                
        temp_x.append(p_ele_dp_ele_from_ang_ctof_isr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)
    
    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('El. from #Theta_{e}#Theta_{p} #Delta P_{e} vs P_{e} Pr (CTOF) S'+str(ss)+ ', #gamma (ISR) ; P_{e} (GeV); #Delta P_{e}')    
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    l0 = TLine(temp_x[0], 0.0, temp_x[-1], 0.0)
    l0.SetLineColor(kRed)
    l_temp0.append(l0)
    store_graphs.append(g_temp)
    g_temp0.append(g_temp)

can.Clear()
can.Divide(2,3)
for ss in range(0, len(g_temp0)):
    can.cd(ss+1)
    g_temp0[ss].Draw('AP')
    l_temp0[ss].Draw('same')
can.Print(mon_out_file_name)


##############################################################
## next in principle look at delta p vs theta


###############################################################################################################################3
## Look at Delta P vs P for different theta bin
#p_ele_dp_ele_from_angles_FTOF_3_thetabin4_isr
#p_ele_dp_ele_from_angles_FTOF_tb_isr_sector = []
#p_ele_dp_ele_from_angles_FTOF_tb_isr = []
#for ss in range(1,7):
#    temp_tb =[]
#    for tb in range(0, 5):
#        if 'p_ele_dp_ele_from_angles_FTOF_'+str(ss)+'_thetabin'+str(tb)+'_isr' in hhs:
#            temp_tb.append(hhs['p_ele_dp_ele_from_angles_FTOF_'+str(ss)+'_thetabin'+str(tb)+'_isr'])
#    p_ele_dp_ele_from_angles_FTOF_tb_isr_sector.append(temp_tb)

#can.Clear()
#can.Divide(3,3)
#for temph in range(0, len(p_ele_dp_ele_from_angles_FTOF_tb_isr_sector) ):
#    for hh in temph:
#        can.cd(ss+1)
#        hh.Draw("colz")
        
plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='W vs. #theta_{e}', xtitle='#theta_{e}', ytitle='W', hline=0.0, y_fit_range=[-0.15, 0.15])

    

'''
can.Clear()
can.Divide(2,3)
ii=0
for ss in range(1,7):
    can.cd(ii+1)
    p_ele_dp_ele_from_ang_ctof_fsr[ss].SetTitle("Electron From Ang. #Delta p vs p CTOF FSR; p (GeV); #Delta p")
    p_ele_dp_ele_from_ang_ctof_fsr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)

g_temp2=[]
l_temp2=[]
for ss in range(1,7):
    can.Clear()
    can.Divide(10,10)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_ele_dp_ele_from_ang_ctof_fsr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_ele_dp_ele_from_ang_ctof_fsr[ss].ProjectionY('profY_p_ele_dp_ele_from_ang_ctof_fsr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_ele_dp_ele_from_ang_ctof_fsr_'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_ele_dp_ele_from_ang_ctof_fsr_'+str(xx),'R')
        h_temp.Draw()                
        temp_x.append(p_ele_dp_ele_from_ang_ctof_fsr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)
    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('Electron #Delta P vs P CTOF FSR S'+str(ss))
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    g_temp2.append(g_temp)
    store_graphs.append(g_temp)
    l0 = TLine(temp_x[0], 0.0, temp_x[-1], 0.0)
    l0.SetLineColor(kRed)
    l_temp2.append(l0)


can.Clear()
can.Divide(2,3)
for ss in range(0, len(g_temp1)):
    can.cd(ss+1)
    g_temp2[ss].Draw('AP')
    l_temp2[ss].Draw('same')
can.Print(mon_out_file_name)


ii=0
can.Clear()
can.Divide(2,3)
for ss in range(1,7):
    can.cd(ii+1)
    p_ele_dp_ele_from_ang_ftof_fsr[ss].SetTitle("Electron From Ang. #Delta p vs p FTOF FSR; p (GeV); #Delta p")
    p_ele_dp_ele_from_ang_ftof_fsr[ss].RebinX(5)
    p_ele_dp_ele_from_ang_ftof_fsr[ss].RebinY(5)
    p_ele_dp_ele_from_ang_ftof_fsr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)

g_temp3 =[]
l_temp3= []
for ss in range(1,7):
    can.Clear()
    can.Divide(5,5)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_ele_dp_ele_from_ang_ftof_fsr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_ele_dp_ele_from_ang_ftof_fsr[ss].ProjectionY('profY_p_ele_dp_ele_from_ang_ftof_fsr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_ele_dp_ele_from_ang_ftof_fsr_'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_ele_dp_ele_from_ang_ftof_fsr_'+str(xx),'R')
        #h_temp.Draw()                
        temp_x.append(p_ele_dp_ele_from_ang_ftof_fsr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)
    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('Electron #Delta P vs P FTOF ISR S'+str(ss))
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    store_graphs.append(g_temp)
    g_temp3.append(g_temp)
    l0 = TLine(temp_x[0], 0.0, temp_x[-1], 0.0)
    l0.SetLineColor(kRed)
    l_temp3.append(l0)


can.Clear()
can.Divide(2,3)
for ss in range(0, len(g_temp1)):
    can.cd(ss+1)
    g_temp3[ss].Draw('AP')
    l_temp3[ss].Draw('same')
can.Print(mon_out_file_name)
'''
#######################################################################################
## proton information 2D

'''
can.Clear()
can.Divide(1,1)
ii=0
for ss in range(0,0):
    can.cd(ii+1)
    p_pro_dp_pro_from_ang_ctof_isr[ss].SetTitle("Proton From Ang. #Delta p vs p CTOF ISR ; p (GeV); #Delta p")
    p_pro_dp_pro_from_ang_ctof_isr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)


for ss in range(0,0):
    can.Clear()
    can.Divide(10,10)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_pro_dp_pro_from_ang_ctof_isr.GetNbinsY()) :
        can.cd(xx)
        h_temp = p_pro_dp_pro_from_ang_ctof_isr[ss].ProjectionY('profY_p_pro_dp_pro_from_ang_ctof_isr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_pro_dp_pro_from_ang_ctof_isr_'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_pro_dp_pro_from_ang_ctof_isr_'+str(xx),'R')
        #h_temp.Draw()                
        temp_x.append(p_pro_dp_pro_from_ang_ctof_isr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)
    can.Clear()
    can.Divide(1,1)
    can.cd(1)
    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('Proton #Delta P vs P CTOF ISR S'+str(ss))
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    g_temp.Draw('AP')
    store_graphs.append(g_temp)
    can.Print(mon_out_file_name)

can.Clear()
can.Divide(2,3)
ii=0
for ss in range(1,7):
    can.cd(ii+1)
    p_pro_dp_pro_from_ang_ftof_isr[ss].SetTitle("Proton From Ang. #Delta p vs p FTOF ISR; p (GeV); #Delta p")
    p_pro_dp_pro_from_ang_ftof_isr[ss].RebinX(4)
    p_pro_dp_pro_from_ang_ftof_isr[ss].RebinY(4)
    p_pro_dp_pro_from_ang_ftof_isr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)

g_temp4=[]
l_temp4=[]
for ss in range(1,7):
    can.Clear()
    can.Divide(3,4)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_pro_dp_pro_from_ang_ftof_isr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_pro_dp_pro_from_ang_ftof_isr[ss].ProjectionY('profY_p_pro_dp_pro_from_ang_ftof_isr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_pro_dp_pro_from_ang_ftof_isr_'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_pro_dp_pro_from_ang_ftof_isr_'+str(xx),'R')
        #h_temp.Draw()                
        temp_x.append(p_pro_dp_pro_from_ang_ftof_isr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)

    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('Proton #Delta P vs P FTOF ISR S'+str(ss))
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    g_temp4.append(g_temp)
    store_graphs.append(g_temp)
    l0 = TLine(temp_x[0], 0.0, temp_x[-1], 0.0)
    l0.SetLineColor(kRed)
    l_temp4.append(l0)


can.Clear()
can.Divide(2,3)
for ss in range(0, len(g_temp1)):
    can.cd(ss+1)
    g_temp4[ss].Draw('AP')
    l_temp4[ss].Draw('same')
can.Print(mon_out_file_name)


can.Clear()
can.Divide(1,1)
ii=0
for ss in range(0,1):
    can.cd(ii+1)
    p_pro_dp_pro_from_ang_ctof_fsr[ss].SetTitle("Proton From Ang. #Delta p vs p CTOF FSR; p (GeV); #Delta p")
    p_pro_dp_pro_from_ang_ctof_fsr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)

for ss in range(0,1):
    can.Clear()
    can.Divide(10,10)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_pro_dp_pro_from_ang_ctof_fsr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_pro_dp_pro_from_ang_ctof_fsr[ss].ProjectionY('profY_p_pro_dp_pro_from_ang_ctof_fsr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_pro_dp_pro_from_ang_ctof_fsr_'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_pro_dp_pro_from_ang_ctof_fsr_'+str(xx),'R')
        #h_temp.Draw()                
        temp_x.append(p_pro_dp_pro_from_ang_ctof_fsr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)
    can.Clear()
    can.Divide(1,1)
    can.cd(1)
    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('Proton #Delta P vs P CTOF FSR S'+str(ss))
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    g_temp.Draw('AP')
    store_graphs.append(g_temp)
    can.Print(mon_out_file_name)


can.Clear()
can.Divide(2,3)
ii=0
for ss in range(1,7):
    can.cd(ii+1)
    p_pro_dp_pro_from_ang_ftof_fsr[ss].SetTitle("Proton From Ang. #Delta p vs p FTOF FSR; p (GeV); #Delta p")
    p_pro_dp_pro_from_ang_ftof_fsr[ss].RebinX(4)
    p_pro_dp_pro_from_ang_ftof_fsr[ss].RebinY(4)
    p_pro_dp_pro_from_ang_ftof_fsr[ss].Draw("colz")    
    ii+=1
can.Print(mon_out_file_name)

g_temp5=[]
l_temp5=[]
for ss in range(1,7):
    can.Clear()
    can.Divide(4,4)

    temp_x, temp_xerr, temp_y, temp_yerr = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )    
    for xx in range(1,p_pro_dp_pro_from_ang_ftof_fsr[ss].GetNbinsY()) :
        can.cd(xx)
        h_temp = p_pro_dp_pro_from_ang_ftof_fsr[ss].ProjectionY('profY_p_pro_dp_pro_from_ang_ftof_fsr_'+str(xx)+'_sector'+str(ss),xx,xx+1)    
        f_temp = TF1('p_pro_dp_pro_from_ang_ftof_fsr_'+str(xx),'gaus',-0.2,0.2)
        h_temp.Fit('p_pro_dp_pro_from_ang_ftof_fsr_'+str(xx),'R')
        #h_temp.Draw()                
        temp_x.append(p_pro_dp_pro_from_ang_ftof_fsr[ss].GetXaxis().GetBinCenter(xx))
        temp_xerr.append(0)
        temp_y.append(f_temp.GetParameter(1))
        temp_yerr.append(f_temp.GetParError(1))        
    can.Print(mon_out_file_name)
    if len(temp_x) == 0: continue 
    g_temp = TGraphErrors(len(temp_x),temp_x, temp_y, temp_xerr, temp_yerr)
    g_temp.SetTitle('Proton #Delta P vs P FTOF FSR S'+str(ss))
    g_temp.SetMarkerStyle(21)
    g_temp.SetMarkerColor(kRed)
    g_temp.GetHistogram().SetMaximum(0.1)
    g_temp.GetHistogram().SetMinimum(-0.1)
    store_graphs.append(g_temp)
    g_temp5.append(g_temp)
    l0 = TLine(temp_x[0], 0.0, temp_x[-1], 0.0)
    l0.SetLineColor(kRed)
    l_temp5.append(l0)


can.Clear()
can.Divide(2,3)
for ss in range(0, len(g_temp1)):
    can.cd(ss+1)
    g_temp5[ss].Draw('AP')
    l_temp5[ss].Draw('same')
can.Print(mon_out_file_name)
'''
can.Print('{}]'.format(mon_out_file_name))

