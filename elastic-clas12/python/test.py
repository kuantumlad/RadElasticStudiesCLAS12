from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLine
from ROOT import TGraph, TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad, TLatex
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kSpring
import numpy as np

from array import array
import math
import sys,os


x, y = array('d'), array('d')
xt = [1,2,3,4,5,6,7,8,9,10,11,12]
yt = [0,0,35,57,25,22,55,47,89,44,0,0]
for ii in range(0, len(xt)):
  x.append(xt[ii])  
  y.append(yt[ii])  

print x
print y
#g = TGraph( len(x), x, y); 
#f = TF1("f", "[2] * x * x + [1] * x + [0]"); 
#g.Fit(f); 
#g.Draw("AL");


t=[2,3,51,6,7,32,64,0,1]
print(max(t))
print(min(t))


imin = t.index(min(t))
imax = t.index(max(t))
print(imin)
print(imax)

m = min(i for i in t if i > 0)
print(m)


print(' testing next() ' )
myi = next((i for i, x in enumerate(yt) if x), None) # x!= 0 for strict match
mylasti = [i for i, e in enumerate(yt) if e != 0][-1]

print(myi)
print(mylasti)

