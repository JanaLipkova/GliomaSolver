%=============================================================
%
%  Compute and Plot for all parameters in param - vecetor
%  1) marginals (diagonal)
%  2) scatter plot for samlpes (above diagonal)
%  3) Gaussian KDE of samples (belov diagonal)
%-------------------------------------------------------------
% INPUT:
%  mydata      = full matrix from curgen_db file
%  param       = vector of parameters of interest
%  bounds      = bounds for parameters (from prior)
%  names       = names of parameters
%  groundTruth = vector of ground truth parameters
%  dump        = 0-1 to dump output
%  KDEbins     = number of bins over which to discretized kde the space
%  Nbins       = number of bins over which to compute marginal PDF
%  bSynthetic  = 1 for synthetic data (i.e. plot groundTruth) 0 for patient case
%  path2output = where to save the plot
%=============================================================



function plotPosteriorTMCMC(mydata,param,bounds,names,groundTruth,dump,KDEbins,Nbins,bSynthetic,path2output)

% 1) colors set up - dark for marginal line, light for filling below curve
%-------------------------------------------------------------------------

a1 = [0.0, 0.2, 0.0]; % dark green
b1 = [0.0, 0.8, 0.0]; % light green

a2 = [0.0, 0.2, 1.0];  % dark blue
b2 = [0.2, 0.8, 1.0];  % light blue

a3 = [0.0, 0.8, 1.0];  % dark cyan
b3 = [0.4, 1.0, 1.0];  % light cyan

a4 = [0.6, 0.0, 0.8]; % dark purple
b4 = [0.6, 0.6, 1.0]; % light purple

a5 = [1.0, 0.0, 0.8];  % dark pink
b5 = [1.0, 0.6, 1.0];  % light pink

a6 = [1.0, 0.4, 0.2]; % orange
b6 = [1.0, 1.0, 0.4]; % yellow

a = [a1;a2;a3;a4;a5;a6];
b = [b1;b2;b3;b4;b5;b6;];

% if more than 6 parameters, repeat colors
a = [a;a];
b = [b;b];

% further plotting set upsL
pointsize=10;
N = length(param);


%2) report statistic
%----------------------------------------------------------------
bestC = find( max(mydata(:,end)) == mydata(:,end));
best = mydata(bestC,:);
meanData = mean(mydata);
varData = var(mydata);
stdData = sqrt(varData);


% 3) Processed and plot results
%----------------------------------------------------------------
figure,
hold on
fs=16;
set(gca,'Fontsize',fs);

for i = 1:N
    for j = 1:N
        
        plotId = i + N * (j-1);
        
        if( i==j)  % diagonal -> marginal
            
            parId = param(i);
            [mPDF, bins] = getMarginalPDF(mydata,parId,bounds(parId,:),KDEbins(parId),Nbins);
            
            subplot(N,N,plotId)
            set(gca,'Fontsize',fs);
            set(gca,'Linewidth',2)
            
            hold on
            z=zeros(size(mPDF));
            [fillhandle,msg]=jbfill(bins,mPDF,z,b(parId,:),a(parId,:),0,0.5);
            grid on;box on;
            xlim(bounds(parId,:));
             title([names(parId,:)]);
  
            set(gca,'ytick',[])
            ax = gca;
            ax.TickLength = [0.02 0.035];

            
        else if( j < i ) % below diagonal -> scatter samples
                
                parId1 = param(i);
                parId2 = param(j);
                
                subplot(N,N,plotId)
                set(gca,'Fontsize',fs);
                set(gca,'Linewidth',2)
                
                hold on
                scatter(mydata(:,parId1), mydata(:,parId2),pointsize,mydata(:,end),'*','Linewidth',3);
                if (bSynthetic)
                    plot(groundTruth(parId1),groundTruth(parId2),'*k','Linewidth',2);
                end;

                axis( [bounds(parId1,1) bounds(parId1,2) bounds(parId2,1) bounds(parId2,2) ])
                box on;
                ax = gca;
                ax.XAxisLocation = 'top';
                ax.YAxisLocation = 'right';
                ax.TickLength = [0.02 0.035];
                
                if(i<N)
                    set(gca,'ytick',[])
                end;
                
                if(j>1)
                    set(gca,'xtick',[])
                end;
                                
                %compute correlation coefficient
                cc = corrcoef(mydata(:,parId1),mydata(:,parId2));
                str1 = num2str(cc(1,2));
                str2 = str1;
                v = axis;
                text( v(1) + 0.33*(v(2) - v(1)), v(3) + 0.5*(v(4) - v(3)), str2, 'FontSize',16);

                
            else  % above diagonal -> KDE
                
                parId1 = param(i);
                parId2 = param(j);
                
                X=linspace( bounds(parId1,1),bounds(parId1,2),Nbins);
                Y=linspace( bounds(parId2,1),bounds(parId2,2),Nbins);
                
                
                subplot(N,N,plotId)
                set(gca,'Fontsize',fs);
                set(gca,'Linewidth',2)
                
                hold on
                Gaus_jointPDF = getKde_2D_PDF(mydata,[parId1, parId2], bounds,KDEbins,Nbins);
                
                pcolor(X,Y,Gaus_jointPDF);
                if (bSynthetic)
                    plot(groundTruth(parId1),groundTruth(parId2),'*k','Linewidth',2)
                end;
                axis( [bounds(parId1,1) bounds(parId1,2) bounds(parId2,1) bounds(parId2,2) ])
                
                if(j<N)
                    set(gca,'xtick',[])
                else
                    xlabel( names(parId1,:) );
                end;
                
                if(i>1)
                    set(gca,'ytick',[])
                else
                    ylabel( names(parId2,:) );
                end;
                
                ax=gca;
                ax.TickLength = [0.02 0.035];
                
                shading flat;% colorbar
            end;
            
            
        end;
        
    end;
end;

set(gcf, 'Position', [100, 100, 1500, 950]);
tightfig;

if(dump)
    set(gcf,'papersize',[50,40]);
    set(gcf,'PaperPositionMode','auto')
%     print([path2output,'PosteriorPDF'],'-dpdf')
    print([path2output,'PosteriorPDF'],'-djpeg')
end
