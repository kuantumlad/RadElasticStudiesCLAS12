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

def rebinHist(title, rebinX=1, rebinY=1, start_range = 1, end_range = 7):
    for ii in range(start_range,end_range):        
        if title.format(ii) in hhs:
            hhs[title.format(ii)].RebinX(rebinX)
            hhs[title.format(ii)].RebinY(rebinY)

def fit_graph(title,graph):
    print(' Fitting {}'.format(graph.GetTitle()))
    #graph_fit = 
    graph.Fit("poly2")

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

def fit_slices(histo, x_range, x_bin_step, y_fit_range, y_fit_special):
    print(' Status of special fit is ')
    print(y_fit_special)
    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)
        fit = None
        #fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus')
        #fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
        #fit.SetName(histo.GetTitle() + '_fit{}'.format(i))
        
        if y_fit_range and not y_fit_special:
            print(y_fit_range)
            print(y_fit_special)
            fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus')
            fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
            fit.SetName(histo.GetTitle() + '_fit{}'.format(i))        
            fit.SetParameter(1, 0.5 * (y_fit_range[0] + y_fit_range[1]))
            projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])

        if y_fit_special:
            print(" OMG USING SPECIAL FIT FUNCTION " )
            fit = getFitFractionalHeight(projec,histo.GetTitle() + 'special_fit{}'.format(i), 0.6, fits)

        #projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])

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
        means.append( round(f.GetParameter(1),3) )
        means_err.append(f.GetParError(1))
        stds.append(f.GetParameter(2))
        stds_err.append(f.GetParError(1)) ### change if want the sigma of dist.
        zeros.append(0.0)

    graph = TGraphErrors(len(x_values), x_values, means, zeros, stds_err)
    graph.SetName('g_' + histo.GetName())
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(kRed)

    htitle = histo.GetName()
    # write out to txt files here
    f_txt_out = open(htitle+'_bcsim.txt','w')
    for ii in range(0, len(x_values) ):
        f_txt_out.write(str(x_values[ii]) + ' ' + str(means[ii]) + ' ' + str(stds[ii]) + '\n' )
    f_txt_out.close()
        

    ## find the x-min and x-max values to use for fitting range limits 
    ## of the delta p or theta vs p or theta
    min_i = next((i for i, x in enumerate(means) if x), None) # x!= 0 for strict match
    max_i = [i for i, e in enumerate(means) if e != 0]
    x_fit_range_min=0
    x_fit_range_max=0    
    if False:# min_i != None and max_i != 0 and len(x_values)>0 :        
        x_fit_range_min=x_values[min_i]
        x_fit_range_max=x_values[max_i[-1]]
            
    return graph, slices, fits, x_fit_range_min, x_fit_range_max


def plot_fits(canvas, histos, x_range, x_bin_step, title_formatter,
              save_name, label, y_fit_range, y_range=None,
              title=None, xtitle=None, ytitle=None, hline=None, 
              fit_special=None, output_txt=None, start_range=1, end_range=7, manual_fit_par=None):


    print( ' plot fits function status of fit special is ' )
    print(fit_special)
    canvas.Clear()
    canvas.Divide(3,2)

    root_is_dumb = []

    #add can shift if starting at non traditional value of 1.
    can_shift=0
    if start_range < 0: can_shift=-start_range+1
    if start_range == 0: can_shift=1
    print(' need to shift canvas by {}'.format(can_shift))

    for i in range(start_range,end_range):#len(histos)):
        canvas.cd(i+can_shift)

        title_temp = title_formatter.format(i)
        graph, slices, fits, x_min, x_max = fit_slices(histos.get(title_temp, default_histo2d), x_range, x_bin_step, y_fit_range, fit_special)
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
            graph.SetTitle("")
            label.DrawLatex(0.1, 0.925, title + str(i))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

        if output_txt:
            print('>>> trying to fit graph')
            if manual_fit_par: 
                x_min, x_max = manual_fit_par[i][0], manual_fit_par[i][1]
                print(' using manual fit parameters {} {}'.format(x_min, x_max))
            poly_fit = TF1(title + str(i), "[2]*x*x + [1]*x + [0]",x_min, x_max)            
            graph.Fit(poly_fit,'R')
            poly_fit.SetLineColor(kRed)
            poly_fit.Draw('same')
            f_out = open('/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/python/fit_parameters/'+title_temp+'.txt','w')
            f_out.write('a {}'.format(poly_fit.GetParameter(0)) + ' \n')
            f_out.write('b {}'.format(poly_fit.GetParameter(1)) + ' \n')
            f_out.write('c {}'.format(poly_fit.GetParameter(2)) + ' \n')
            f_out.close()
            
            
            

    canvas.Print(save_name)

    # For the slices 
    slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
    slice_pdfname = title_formatter.split('_{}')[0] + '_slices.pdf'
    slice_can.Print(slice_pdfname + '[')
    for i in range(start_range,end_range): #usually 1,7
        title = title_formatter.format(i)
        graph, slices, fits, x_min, x_max = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range, fit_special)

        # Size of slices page
        nrows = 5
        ncols = int(np.ceil(len(slices) / nrows) + 1)
        slice_can.Clear() 
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices, fits)):
            slice_can.cd(j+1)
            s.Draw()
            lab.DrawLatex(0.15, 0.84, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(f.GetParameter(1), f.GetParameter(2)))
            
        slice_can.Print(slice_pdfname)
    slice_can.Print(slice_pdfname + ']')
    

def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False,
                     y_fit_range=None, landscape=False, x_range=None, y_range=None, vline=None,
                     hline=None, start_range=1, end_range=7):

    root_garbage_can = []
    
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    #add can shift if starting at non traditional value of 1.
    can_shift=0
    if start_range < 0: can_shift=-start_range+1
    if start_range == 0: can_shift=1
    print(' need to shift canvas by {}'.format(can_shift))

    for i in range(start_range,end_range):
        canvas.cd(i+can_shift)
        
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
                label.DrawLatex(0.15, 0.84, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
        #elif isinstance(histos[title_formatter.format(i)], TH2F):
        elif isinstance(histos.get(title_formatter.format(i), default_histo), TH2F):
            # histos[title_formatter.format(i)].Draw('colz')
            histos.get(title_formatter.format(i), default_histo).Draw('colz')
            if log:
                gPad.SetLogz() 

            if y_range:
                histos.get(title_formatter.format(i), default_histo).GetYaxis().SetRangeUser(y_range[0], y_range[1])

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
            label.DrawLatex(0.1, 0.925, title + ' ' + str(i))
            label.DrawLatex(0.75, 0.85, str( int(histos.get(title_formatter.format(i), default_histo).GetEntries())) )

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

def plot_page(canvas, histos, histo_title, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False, hline=None):
    
    canvas.Clear() 
    root_garbage_can=[]
    if isinstance(histos[histo_title], TH1F):
        histos[histo_title].Draw()
    elif isinstance(histos[histo_title], TH2F):
        histos[histo_title].Draw('colz')
        if log:
            gPad.SetLogz() 
    else:
        #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
        pass
    
    if title:
        histos.get(histo_title, default_histo).SetTitle('')
        label.DrawLatex(0.1, 0.925, title)
        
    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)
        
    if hline:
        xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
        xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 
        line = TLine(xmin, hline, xmax, hline)
        line.SetLineColor(1)
        line.SetLineStyle(1)
        line.Draw('same')
        root_garbage_can.append(line)

    canvas.Print(save_name)



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

################################
## elastic scattering
### look at DeltaP/P vs local Phi vs different theta bins
## elastic with proton in ftof

rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
#phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin1_pass_all_elastic
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S1 TB', log=True, hline=0.00001, 
                 start_range=1, end_range=6)
plot_fits(can, hhs, x_range=[-14.0, 14.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S1 TB', hline=0.00001, start_range=1, end_range=6)

rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_2_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_2_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S2', log=True, hline=0.00001, 
                 start_range=1, end_range=6)
plot_fits(can, hhs, x_range=[-14.0, 14.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_2_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S2 TB', hline=0.00001, start_range=1, end_range=6)

rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_3_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_3_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S3', log=True, hline=0.00001, 
                 start_range=1, end_range=6)
plot_fits(can, hhs, x_range=[-14.0, 14.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_3_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S3 TB', hline=0.00001, start_range=1, end_range=6)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_4_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_4_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S4', log=True, hline=0.00001, 
                 start_range=1, end_range=6)
plot_fits(can, hhs, x_range=[-14.0, 14.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_4_thetabin{}_pass_all_elastic', 
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S4 TB', hline=0.00001, start_range=1, end_range=6)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_5_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_5_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S5', log=True, hline=0.00001
                 , start_range=2, end_range=6)
plot_fits(can, hhs, x_range=[-14.0, 14.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_5_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S5 TB', hline=0.00001, start_range=1, end_range=6)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_6_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_6_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S6', log=True, hline=0.00001, 
                 start_range=1, end_range=6)
plot_fits(can, hhs, x_range=[-14.0, 14.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_6_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), S6 TB', hline=0.00001, start_range=1, end_range=6)


### elastic with proton in ctof

rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_1_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_1_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S1 TB', log=True, hline=0.00001, start_range=1, end_range=4)
plot_fits(can, hhs, x_range=[-14.0, 20.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_1_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S1 TB', hline=0.00001, start_range=1, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_2_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_2_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S2', log=True, hline=0.00001, start_range=1, end_range=4)
plot_fits(can, hhs, x_range=[-14.0, 20.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_2_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S2 TB', hline=0.00001, start_range=1, end_range=4)

rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_3_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_3_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S3', log=True, hline=0.00001, start_range=1, end_range=4)
plot_fits(can, hhs, x_range=[-14.0, 20.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_3_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S3 TB', hline=0.00001, start_range=1, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_4_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_4_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S4', log=True, hline=0.00001, start_range=1, end_range=4)
plot_fits(can, hhs, x_range=[-14.0, 20.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_4_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S4 TB', hline=0.00001, start_range=1, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_5_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_5_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S5', log=True, hline=0.00001, start_range=1, end_range=4)
plot_fits(can, hhs, x_range=[-14.0, 20.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_5_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S5 TB', hline=0.00001, start_range=1, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_6_thetabin{}_pass_all_elastic', rebinX=5, rebinY=1)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_6_thetabin{}_pass_all_elastic',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S6', log=True, hline=0.00001, start_range=1, end_range=4)
plot_fits(can, hhs, x_range=[-14.0, 20.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_6_thetabin{}_pass_all_elastic',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.05, 0.05],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='Elastic, El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), S6 TB', hline=0.00001, start_range=1, end_range=4)



###########

## delta P/P for ISR in FTOF and CTOF in theta bins as a function of phi of electron
## 
rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin{}_isr', rebinX=10, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S1 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 28.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_1_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S1 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_2_thetabin{}_isr', rebinX=10, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_2_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S2 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 28.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_2_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S2 TB', hline=0.00001, start_range=0, end_range=4)

rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_3_thetabin{}_isr', rebinX=10, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_3_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S3 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 28.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_3_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S3 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_4_thetabin{}_isr', rebinX=10, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_4_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S4 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 28.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_4_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S4 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_5_thetabin{}_isr', rebinX=10, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_5_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S5 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 28.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_5_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S5 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_FTOF_6_thetabin{}_isr', rebinX=10, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_FTOF_6_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S6 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 28.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_FTOF_6_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (FD), #gamma (ISR), S6 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_1_thetabin{}_isr', rebinX=2, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_1_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S1 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 30.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_1_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S1 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_2_thetabin{}_isr', rebinX=2, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_2_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S2 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 30.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_2_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S2 TB', hline=0.00001, start_range=0, end_range=4)

rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_3_thetabin{}_isr', rebinX=2, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_3_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S3 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 30.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_3_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S3 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_4_thetabin{}_isr', rebinX=2, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_4_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S4 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 30.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_4_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S4 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_5_thetabin{}_isr', rebinX=2, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_5_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S5 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 30.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_5_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S5 TB', hline=0.00001, start_range=0, end_range=4)


rebinHist('phi_local_ele_dpp_ele_from_angles_CTOF_6_thetabin{}_isr', rebinX=2, rebinY=2, start_range=0, end_range=4)
plot_sector_page(can,hhs,'phi_local_ele_dpp_ele_from_angles_CTOF_6_thetabin{}_isr',lab, save_name=mon_out_file_name,
                 xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S6 TB', log=True, hline=0.00001, start_range=0, end_range=4)
plot_fits(can, hhs, x_range=[-20.0, 30.0], x_bin_step=1, title_formatter='phi_local_ele_dpp_ele_from_angles_CTOF_6_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_fit_range=[-0.2, 0.2], y_range=[-0.025, 0.025],
          xtitle='#Phi_{e} (deg)', ytitle='#Delta P_{e}/P_{e}', title='El. #Delta P_{e}/P_{e} vs local #Phi_{e}, Pr. (CD), #gamma (ISR), S6 TB', hline=0.00001, start_range=0, end_range=4)


'''


######################################################################################################################################################################################
######################################################################################################################################################################################

rebinHist('p_ele_dp_ele_from_angles_FTOF_{}_isr',rebinX=6,rebinY=1)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. #Delta P_{e} vs. P_{e}, #Theta_{e}#Theta_{p}, Pr. (FTOF), #gamma (ISR) S', log=True,  hline=0.0)
        
#output txt files here
plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.15, 0.15],
          title='El. #Delta P_{e} vs. P_{e}, #Theta_{e}#Theta_{p}, Pr. (FTOF), #gamma (ISR) ', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.2, 0.2])
 
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. #Delta P_{e} vs. P_{e}, #Theta_{e}#Theta_{p}, Pr. (CTOF), #gamma (ISR) S', log=True,  hline=0.0)

# output text files here
rebinHist('p_ele_dp_ele_from_angles_CTOF_{}_isr',rebinX=6,rebinY=1)
plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.5, 0.5],
          title='El. #Delta P_{e} vs. P_{e}, #Theta_{e}#Theta_{p}, Pr. (CTOF), #gamma (ISR) ', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, fit_special=True, y_fit_range=[-0.15, 0.15])

## look at fraction of delta p/p vs p for electron
#rebi
rebinHist('p_ele_dpp_ele_from_angles_FTOF_{}_isr',rebinX=6, rebinY=1)
plot_sector_page(can,hhs, 'p_ele_dpp_ele_from_angles_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e}/P_{e}', title ='El. (FD), Pr (FD), #Delta P_{e}/P_{e} vs P_{e} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', log=True,  hline=0.0001)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dpp_ele_from_angles_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.15, 0.15],
          title='El. (FD), Pr (FD), #Delta P_{e}/P_{e} vs P_{e} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e}/P_{e} (GeV)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])

plot_sector_page(can,hhs, 'p_ele_dpp_ele_from_angles_CTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e}/P_{e}', title ='El. (FD), Pr (CD), #Delta P_{e}/P_{e} vs P_{e} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', log=True,  hline=0.0001)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dpp_ele_from_angles_CTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.15, 0.15],
          title='El. (FD), Pr (CD), #Delta P_{e}/P_{e} vs P_{e} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e}/P_{e} (GeV)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])

#### now look at #DeltaP/P vs P for different theta bins
rebinHist('theta_ele_dp_ele_from_angles_FTOF_{}_isr',rebinX=10, rebinY=4)
plot_sector_page(can, hhs, 'theta_ele_dp_ele_from_angles_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='#Theta_{e} (deg)', ytitle='#Delta P_{e} (GeV)', title='El #Delta P_{e} vs #Theta_{e}, Pr. (FTOF), #gamma (ISR) S', log=True, hline=0.0001)
plot_fits(can,hhs,x_range=[10.0, 35.0], x_bin_step=1, title_formatter='theta_ele_dp_ele_from_angles_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-1, 1],
          title='El. (FD), Pr (FD), #Delta P_{e} vs #Theta_{e} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', xtitle='#Theta_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])


rebinHist('theta_ele_dp_ele_from_angles_CTOF_{}_isr',rebinX=6, rebinY=1)
plot_sector_page(can, hhs, 'theta_ele_dp_ele_from_angles_CTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='#Theta_{e} (deg)', ytitle='#Delta P_{e} (GeV)', title='El #Delta P_{e} vs #Theta_{e}, Pr. (CTOF), #gamma (ISR) S', log=True, hline=0.0001)
plot_fits(can,hhs,x_range=[10.0, 35.0], x_bin_step=1, title_formatter='theta_ele_dp_ele_from_angles_CTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-1, 1],
          title='El. (FD), Pr (CD), #Delta P_{e} vs #Theta_{e} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', xtitle='#Theta_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])

# FSR here is useless 
plot_sector_page(can, hhs, 'theta_ele_dp_ele_from_angles_FTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 xtitle='#Theta_{e} (deg)', ytitle='#Delta P_{e} (GeV)', title='El #Delta P_{e} vs #Theta_{e}, Pr. (FTOF), #gamma (FSR) S', log=True, hline=0.0001)
plot_sector_page(can, hhs, 'theta_ele_dp_ele_from_angles_CTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 xtitle='#Theta_{e} (deg)', ytitle='#Delta P_{e} (GeV)', title='El #Delta P_{e} vs #Theta_{e}, Pr. (CTOF), #gamma (FSR) S', log=True, hline=0.0001)

#look at delta theta vs p
rebinHist('p_ele_dtheta_ele_from_angles_FTOF_{}_isr', rebinX=10, rebinY=8)
plot_sector_page(can, hhs, 'p_ele_dtheta_ele_from_angles_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta #Theta (deg)', title='El. #Delta #Theta_{e} vs P_{e}, Pr. (FTOF), #gamma (ISR) S', log=True, hline=0.0000001)
plot_fits(can,hhs,x_range=[2.0, 8.0], x_bin_step=1, title_formatter='p_ele_dtheta_ele_from_angles_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-1, 1],
          title='El. #Delta #Theta_{e} vs P_{e}, Pr. (FTOF), #gamma (ISR) S', xtitle='P_{e} (GeV)', ytitle='#Delta #Theta_{e} (deg)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])


rebinHist('p_ele_dtheta_ele_from_angles_CTOF_{}_isr', rebinX=4, rebinY=4)
plot_sector_page(can, hhs, 'p_ele_dtheta_ele_from_angles_CTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta #Theta (deg)', title='El. #Delta #Theta_{e} vs P_{e}, Pr. (CTOF), #gamma (ISR) S', log=True, hline=0.0000001)
plot_fits(can,hhs,x_range=[2.0, 8.0], x_bin_step=1, title_formatter='p_ele_dtheta_ele_from_angles_CTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-5, 5],
          title='El. #Delta #Theta_{e} vs P_{e}, Pr. (CTOF), #gamma (ISR) S', xtitle='P_{e} (GeV)', ytitle='#Delta #Theta_{e} (deg)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])

rebinHist('p_ele_dtheta_ele_from_angles_FTOF_{}_fsr',rebinX=10, rebinY=4)
plot_sector_page(can, hhs, 'p_ele_dtheta_ele_from_angles_FTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta #Theta (deg)', title='El. #Delta #Theta_{e} vs P_{e}, Pr. (FTOF), #gamma (FSR) S', log=True, hline=0.0000001)
plot_fits(can,hhs,x_range=[5.0, 8.0], x_bin_step=1, title_formatter='p_ele_dtheta_ele_from_angles_FTOF_{}_fsr',
          save_name=mon_out_file_name, label=lab, y_range=[-1, 1],
          title='El. #Delta #Theta_{e} vs P_{e}, Pr. (FTOF), #gamma (FSR) S', xtitle='P_{e} (GeV)', ytitle='#Delta #Theta_{e} (deg)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])

rebinHist('p_ele_dtheta_ele_from_angles_CTOF_{}_fsr',rebinX=8, rebinY=4)
plot_sector_page(can, hhs, 'p_ele_dtheta_ele_from_angles_CTOF_{}_fsr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta #Theta (deg)', title='El. #Delta #Theta_{e} vs P_{e}, Pr. (CTOF), #gamma (FSR) S', log=True, hline=0.0000001)
plot_fits(can,hhs,x_range=[7.0, 9.0], x_bin_step=1, title_formatter='p_ele_dtheta_ele_from_angles_CTOF_{}_fsr',
          save_name=mon_out_file_name, label=lab, y_range=[-1, 1],
          title='El. #Delta #Theta_{e} vs P_{e}, Pr. (CTOF), #gamma (FSR) S', xtitle='P_{e} (GeV)', ytitle='#Delta #Theta_{e} (deg)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])

# look at the dependence of delta p vs p in different theta bins for each sector
# not really looking at sector but it works since there are so few bins in theta

fit_min_max={ -1: [0, 0],
              0:  [6.5, 8.5],
              1:  [4.5, 8.25],
              2:  [3.0, 7.0],
              3:  [2.5, 5.5],
              4:  [1.5, 3.75],
              5:  [0, 0],
              6:  [0, 0]
            } 

fit_min_max_pr_ctof={ -1: [0, 0],
                       0: [2.2, 8.0],
                       1: [2.0, 6.0],
                       2: [2.0, 4.0],
                       3: [2.0, 3.0],
                       4: [2.0, 2.3],
                       5: [0,0]
            } 

'''

'''
rebinHist('p_ele_dp_ele_from_angles_FTOF_1_thetabin{}_isr',rebinX=6,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_1_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.1  #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_1_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.35, 0.35],
          title='El. Sect. 1 #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_1_thetabin{}_isr_corr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.35, 0.35],
          title='El. Sect. 1 #Delta P_{e}/P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) Corr, #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e}/P_{e}', hline=0.0, y_fit_range=[-0.15, 0.15],
          start_range = 0, end_range=6)


rebinHist('p_ele_dp_ele_from_angles_FTOF_2_thetabin{}_isr',rebinX=6,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_2_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.2  #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_2_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 2 #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6,manual_fit_par=fit_min_max)


rebinHist('p_ele_dp_ele_from_angles_FTOF_3_thetabin{}_isr',rebinX=6,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_3_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.3  #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0,start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_3_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 3 #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max)



rebinHist('p_ele_dp_ele_from_angles_FTOF_4_thetabin{}_isr',rebinX=6,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_4_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.4  #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_4_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 4 #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max)

rebinHist('p_ele_dp_ele_from_angles_FTOF_5_thetabin{}_isr',rebinX=6,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_5_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.5  #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_5_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 5 #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max)

rebinHist('p_ele_dp_ele_from_angles_FTOF_6_thetabin{}_isr',rebinX=6,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_FTOF_6_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.6  #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_FTOF_6_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 6 #Delta P_{e} vs. P_{e}, Pr. (FTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max)


## now for when proton is in the CTOF
rebinHist('p_ele_dp_ele_from_angles_CTOF_1_thetabin{}_isr',rebinX=1,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_1_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.1  #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_1_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 1 (FD) #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max_pr_ctof)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_1_thetabin{}_isr_corr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 1 (FD) #Delta P_{e}/P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR), P Corrected, #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e}/P_{e}', hline=0.0, 
          y_fit_range=[-0.15, 0.15], start_range = 0, end_range=6)

rebinHist('p_ele_dp_ele_from_angles_CTOF_2_thetabin{}_isr',rebinX=1,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_2_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.2  #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_2_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 2 (FD) #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max_pr_ctof)

rebinHist('p_ele_dp_ele_from_angles_CTOF_3_thetabin{}_isr',rebinX=1,rebinY=1, start_range=-1, end_range=7)
plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_3_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.3  #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_3_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 3 (FD) #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max_pr_ctof)

plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_4_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.4  #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_4_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 4 (FD) #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max_pr_ctof)

plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_5_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='El. Sect.5  #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_5_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 5 (FD) #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max_pr_ctof)

plot_sector_page(can,hhs,'p_ele_dp_ele_from_angles_CTOF_6_thetabin{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta Pe_{e} (GeV)', title='El. Sect.6  #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', log=True,  hline=0.0, start_range=0, end_range=6)

plot_fits(can, hhs, x_range=[2.0,9.0], x_bin_step=1, title_formatter='p_ele_dp_ele_from_angles_CTOF_6_thetabin{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.1, 0.1],
          title='El. Sect. 6 (FD) #Delta P_{e} vs. P_{e}, Pr. (CTOF), #gamma (ISR) #Theta bin', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15],
          output_txt=True, start_range = 0, end_range=6, manual_fit_par=fit_min_max_pr_ctof)

rebinHist('theta_ele_dtheta_ele_from_beam_pangle_FTOF_{}_fsr',rebinX=1,rebinY=2)
plot_sector_page( can, hhs, 'theta_ele_dtheta_ele_from_beam_pangle_FTOF_{}_fsr', lab, save_name=mon_out_file_name,
                  xtitle='#Theta_{e} (deg)', ytitle='#Delta #Theta_{e} (deg)', title='El. #Delta #Theta_{e} vs. #theta_{e}, Pr. (FTOF), #gamma (FSR) S', log=True,  hline=0.0001)

plot_fits(can, hhs, x_range=[10.0,18.0], x_bin_step=4, title_formatter='theta_ele_dtheta_ele_from_beam_pangle_FTOF_{}_fsr',
          save_name=mon_out_file_name, label=lab, y_range=[-2, 2],
          title='El. #Delta #Theta_{e} vs. #Theta_{e}, Pr. (FTOF), #gamma (FSR) S', xtitle='#Theta_{e} (deg)', ytitle='#Delta #Theta_{e} (deg)', hline=0.0, y_fit_range=[-1.5, 1.5])#fit_special=True)

rebinHist('theta_ele_dtheta_ele_from_beam_pangle_CTOF_{}_fsr',rebinX=2,rebinY=4)

plot_sector_page( can, hhs, 'theta_ele_dtheta_ele_from_beam_pangle_CTOF_{}_fsr', lab, save_name=mon_out_file_name,
                  xtitle='#Theta_{e} (deg)', ytitle='#Delta #Theta_{e} (deg)', title='El. #Delta #Theta_{e} vs. #theta_{e}, Pr. (CTOF), #gamma (FSR) S', log=True,  hline=0.0001)

#plot_fits(can, hhs, x_range=[6.5,12.0], x_bin_step=2, title_formatter='theta_ele_dtheta_ele_from_beam_pangle_CTOF_{}_fsr',
#          save_name=mon_out_file_name, label=lab, y_range=[-5, 5],
#          title='El. #Delta #Theta_{e} vs. #Theta_{e}, Pr. (CTOF), #gamma (FSR) S', xtitle='#Theta_{e} (deg)', ytitle='#Delta #Theta_{e} (deg)', hline=0.0, y_fit_range=[-1.5, 1.5], fit_special=True,)
'''

#####################################################################################################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################
## 
## 
## now we look at the protons 
'''
rebinHist('p_pro_dp_pro_FTOF_{}_isr',rebinX=10,rebinY=4)
plot_sector_page(can,hhs,'p_pro_dp_pro_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{p} (GeV)', ytitle='#Delta P_{p} (GeV)', title='Pr. #Delta P_{p} vs. P_{p}, Pr. (FTOF), #gamma (ISR) S', log=True,  hline=0.00001)#, y_range=[-0.4,0.4])

plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_pro_dp_pro_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.6, 0.6],
          title='Pr. #Delta P_{p} vs. P_{p}, Pr. (FTOF), #gamma (ISR) ', xtitle='P_{p} (GeV)', ytitle='#Delta P_{p} (GeV)', hline=0.0, y_fit_range=[-0.15, 0.15])

#plot_page(can,hhs,'p_pro_dp_pro_CTOF_0_isr', lab, save_name=mon_out_file_name,
#                 xtitle='P_{p} (GeV)', ytitle='#Delta P_{p} (GeV)', title='Pr. #Delta P_{p} vs. P_{p}, Pr. (CTOF), #gamma (ISR) S', log=True,  hline=0.00001)#, y_range=[-0.4, 0.4])

## proton dpp
rebinHist('p_pro_dpp_pro_FTOF_{}_isr',rebinX=6,rebinY=2)
plot_sector_page(can,hhs, 'p_pro_dpp_pro_FTOF_{}_isr', lab, save_name=mon_out_file_name,
                 xtitle='P_{p} (GeV)', ytitle='#Delta P_{p}/P_{p}', title ='El. (FD), Pr (FD), #Delta P_{p}/P_{p} vs P_{p} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', log=True,  hline=0.0001, y_range=[-0.3, 0.3])
plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_pro_dpp_pro_FTOF_{}_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.15, 0.15],
          title='El. (FD), Pr (FD), #Delta P_{p}/P_{p} vs P_{p} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', xtitle='P_{p} (GeV)', ytitle='#Delta P_{p}/P_{p} (GeV)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])


#plot_page(can,hhs, 'p_pro_dpp_pro_CTOF_0_isr', lab, save_name=mon_out_file_name,
 #                xtitle='P_{p} (GeV)', ytitle='#Delta P_{p}/P_{p}', title ='El. (FD), Pr (CD), #Delta P_{p}/P_{p} vs P_{p} using Theta_{e}#Theta_{p}, #gamma (ISR) S', log=True,  hline=0.0001)

'''
'''
plot_fits(can, hhs, x_range=[2.0,8.0], x_bin_step=1, title_formatter='p_pro_dpp_pro_FTOF_0_isr',
          save_name=mon_out_file_name, label=lab, y_range=[-0.15, 0.15],
          title='El. (FD), Pr (CD), #Delta P_{p}/P_{p} vs P_{p} using #Theta_{e}#Theta_{p}, #gamma (ISR) S', xtitle='P_{p} (GeV)', ytitle='#Delta P_{p}/P_{p} (GeV)', hline=0.0001, fit_special=True, y_fit_range=[-0.15, 0.15])
'''


can.Print('{}]'.format(mon_out_file_name))

