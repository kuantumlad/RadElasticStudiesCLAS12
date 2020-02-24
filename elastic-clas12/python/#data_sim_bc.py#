from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLine
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad, TLatex
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kSpring, kBlack
import numpy as np

from array import array
import math
import sys,os

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)

def remove_bad_points(x, mu, sig, max_errorbar):
    condition = np.logical_and(mu != 0, sig < max_errorbar)
    idx = np.where(condition)[0]
    return x[idx], mu[idx], sig[idx]


def fit_slices(histo, x_range, x_bin_step, y_fit_range, x_shift, name_add):

    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(name_add+'_'+histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)
        
        fit = TF1(name_add+'_'+histo.GetTitle() + '_fit{}'.format(i), 'gaus')
        fit.SetTitle(name_add+'_'+histo.GetTitle() + '_fit{}'.format(i))
        fit.SetName(name_add+'_'+histo.GetTitle() + '_fit{}'.format(i))        
        
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
        stds_err.append(f.GetParError(2)) ### change if want the sigma of dist.
        #zeros.append(0.0)

    #transform back to np array
    #shift values
    np_x_values=np.array(x_values)
    print(' before shift ' )
    print(np_x_values)
    if x_shift:
        print(' shifting x values ')
        np_x_values += 0.5 * (np_x_values[1]-np_x_values[0])
        print(np_x_values)

    x_good, mu_good, sig_good = remove_bad_points(np_x_values, np.array(means), np.array(stds), 0.8)
    print(x_good)
    final_x = array('d')
    final_mean = array('d')
    final_zeros = array('d')
    final_sig = array('d')
    for ii in range(0,len(x_good)):
        final_x.append(x_good[ii])
        final_mean.append(mu_good[ii])
        final_zeros.append(0.0)
        final_sig.append(sig_good[ii])
        
    #graph = TGraphErrors(len(x_values), x_values, means, zeros, stds)
    
    graph = TGraphErrors(len(final_x), final_x, final_mean, final_zeros, final_sig)
    graph.SetName(name_add+'_'+'g_' + histo.GetName())
 
    return graph, slices, fits
 

def plot_fits_data_sim(canvas, histos_sim, histos_data, x_range, x_bin_step, title_formatter,
              save_name, label, y_fit_range, y_range=None,
              title=None, xtitle=None, ytitle=None, hline=None, 
                       output_txt=None, start_range=1, end_range=7, manual_fit_par=None, xshift=False):


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
        graph_sim, slices_sim, fits_sim = fit_slices(histos_sim.get(title_temp, default_histo2d), x_range, x_bin_step, y_fit_range, True, 'sim')
        graph_sim.SetMarkerStyle(8)
        graph_sim.SetMarkerSize(1)
        graph_sim.SetMarkerColor(kRed)

        graph_data, slices_data, fits_data = fit_slices(histos_data.get(title_temp, default_histo2d), x_range, x_bin_step, y_fit_range, False, 'data')
        graph_data.SetMarkerStyle(9)
        graph_data.SetMarkerSize(1)
        graph_data.SetMarkerColor(kBlack)


        multi_graph_data_sim = TMultiGraph()
        multi_graph_data_sim.Add(graph_data)
        multi_graph_data_sim.Add(graph_sim)
        multi_graph_data_sim.Draw('AP')
        multi_graph_data_sim.GetHistogram().GetXaxis().SetRangeUser(x_range[0], x_range[1]);
        can.Update()
        
        if y_range:
            multi_graph_data_sim.GetHistogram().SetMinimum(y_range[0])
            multi_graph_data_sim.GetHistogram().SetMaximum(y_range[1])
            multi_graph_data_sim.Draw('AP')
            root_is_dumb.append(multi_graph_data_sim)
        else:
            multi_graph_data_sim.Draw('AP')
            root_is_dumb.append(multi_graph_data_sim)
            
        if hline:
            line = TLine(x_range[0], hline, x_range[1], hline)
            line.SetLineStyle(8)
            line.SetLineWidth(2)
            line.SetLineColor(kRed)
            line.Draw()
            root_is_dumb.append(line)
            
        if title:
            multi_graph_data_sim.SetTitle("")
            label.DrawLatex(0.1, 0.925, title + str(i))

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
    for i in range(start_range,end_range): #usually 1,7
        title = title_formatter.format(i)
        graph_sim, slices_sim, fits_sim, = fit_slices(histos_sim.get(title, default_histo2d), x_range, x_bin_step, y_fit_range, xshift, 'sim')

        # Size of slices page
        nrows = 5
        ncols = int(np.ceil(len(slices_sim) / nrows) + 1)
        slice_can.Clear() 
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices_sim, fits_sim)):
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
            label.DrawLatex(0.75, 0.85, str( histos.get(title_formatter.format(i), default_histo).GetEntries()))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)


datatype = sys.argv[3]
hhs_sim={}
hhs_data={}
ff_sim = TFile(sys.argv[1])
ff_data = TFile(sys.argv[2])
# load sim
for kk in ff_sim.GetListOfKeys():
    obj = kk.ReadObj()
    hhs_sim[obj.GetName()] = obj

#load data
for kk in ff_data.GetListOfKeys():
    obj = kk.ReadObj()
    hhs_data[obj.GetName()] = obj

base_sim=os.path.basename(ff_sim.GetName())
base_data=os.path.basename(ff_data.GetName())
print('Analysing file %s %s' % (base_sim, base_data) )

mon_out_file_name='monitor_data_sim_'+datatype+'.pdf'
can = TCanvas('can','can')
can.SetCanvasSize(1200,1200)
can.Print('{}['.format(mon_out_file_name))


lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)

gStyle.SetOptStat(00000)


plot_sector_page(can,hhs_data,'p_electron_delta_p_electron_{}', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='Data: El. #Delta P_{e} vs. P_{e}, Pr. (CTOF), Elastic S', log=True,  hline=0.00001, landscape=True)

plot_sector_page(can,hhs_sim,'p_electron_delta_p_electron_{}', lab, save_name=mon_out_file_name,
                 xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', title='Sim: El. #Delta P_{e} vs. P_{e}, Pr. (CTOF), Elastic S', log=True,  hline=0.00001, landscape=True)

#add the plots fere to compare the theta gamma, theta egamma, and mm2
        
#plot_fits_data_sim(can, hhs_sim, hhs_data, x_range=[8.7,9.8], x_bin_step=6, title_formatter='p_electron_delta_p_electron_{}',
#          save_name=mon_out_file_name, label=lab, y_range=[-0.35, 0.35],
#                   title='El. #Delta P_{e} vs. P_{e}, From #Theta_{e}#Theta_{p}, Pr. (CTOF), Elastic', xtitle='P_{e} (GeV)', ytitle='#Delta P_{e} (GeV)', hline=0.00001, y_fit_range=[-0.2,# 0.2], xshift=True)


can.Print('{}]'.format(mon_out_file_name))



