%===============================================
% Compute Gaussian KDE of 2D PDF
%------------------------
% INPUT:
% mydata  = curgen_db matrix
% param   = indecies of 2 parameters for which to get kde PDF
% bounds  = upper and lower bound for values of parameter if parID
% KDEbins = number of bins over which to discretized kde the space
% Nbins   = number of bins over which to compute marginal PDF
%
% OUTPUT:
%  h     = Gaussian KDE of PDF, with automatic kernel size computation  ("Rule of Thumb"; Silverman '86 / Scott '92)
%  p     = Gaussian KDE of PDF of given two parameters, kernel size given by bw1, bw2 
%===============================================



function p = getKde_2D_PDF(mydata,param,bounds,KDEbins,Nbins)

data = [mydata(:,param(1)), mydata(:,param(2))];
data = data';
bw1 = (bounds(param(1),2) - bounds(param(1),1)) / KDEbins(1);
bw2 = (bounds(param(2),2) - bounds(param(2),1)) / KDEbins(2);

w=1.*ones(1,length(mydata));  %1.*mydata(:,end)'

size(data);
p = kde(data,[bw1;bw2],w);

m = ksize(p, 'rot'); % "Rule of Thumb"; Silverman '86 / Scott '92xh = hist(m, Nbins,[2,1],[ bounds(param(1),:); bounds(param(2),:) ] );
p = hist(p, Nbins,[2,1],[ bounds(param(1),:); bounds(param(2),:) ] );



