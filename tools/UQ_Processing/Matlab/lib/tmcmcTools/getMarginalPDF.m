%===============================================
% Compute Marginal PDF
%------------------------
% INPUT:
% mydata = curgen_db matrix
% parID  = index of parameter for which to get marginal PDF
% range  = upper and lower bound for values of parameter if parID
% KDEbins= number of bins over which to discretized kde the space
% Nbins  = number of bins over which to compute marginal PDF
%
% OUTPUT:
%  mPDF = marginal PDF for parameter parID
%  bins = array of points in which mPDF is evalueated
%===============================================



function [mPDF, bins] = getMarginalPDF(mydata,parId,range,KDEbins,Nbins)

bw = (range(2) - range(1)) / KDEbins;
w=1.*ones(1,length(mydata));  %1.*mydata(:,end)'

p = kde(mydata(:,parId)',bw,w);
bins=linspace(range(1),range(2),Nbins);
mPDF = evaluate(p,bins);

