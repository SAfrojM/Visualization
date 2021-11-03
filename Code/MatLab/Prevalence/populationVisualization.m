% Visualization: prevalence of Variant-1 and variant-2
function [s]= populationVisualization()
tau=[ 0.02, 0.03 ];
newVar=1.4;

emergenceDay=[90 120 150];
seedLoc=[51,51059,51121];
seed=["VA", "Fx", "Mo"];

chopOff=zeros(18*20,4);
ind=1;
rep=60;
for i=1:length(tau)
    for j=1:length(seedLoc)
        for k=1: length(emergenceDay)
            s=poplationEpiCurve(tau(i), newVar, seedLoc(j),seed(j), emergenceDay(k));
            chopOff((ind-1)*rep+1:(ind-1)*rep+rep, 1)=tau(i);
            chopOff((ind-1)*rep+1:(ind-1)*rep+rep, 2)=seedLoc(j);
            chopOff((ind-1)*rep+1:(ind-1)*rep+rep, 3)=emergenceDay(k);
            chopOff((ind-1)*rep+1:(ind-1)*rep+rep, 4)=s;
            ind=ind+1;
            
            
        end
    end
end

s=1
end
function [s]=poplationEpiCurve(tau, newVar, seedLoc,seed, emergenceDay)
%populationCountExp50.03_var2_1.4_seed2n_51_seed2t_90
%populationCountExp5_rep_600.03_var2_1.4_seed2n_51121_seed2t_150
s=strcat('populationCountExp5_rep_60', num2str(tau),'_var2_', num2str(newVar), '_seed2n_',num2str(seedLoc),'_seed2t_',num2str(emergenceDay),'.csv');
U3all=readtable(s);
U3=table2array(U3all);
% zeroCount=length(U3(U3(:,1)==0,1));
% U3=U3(1:size(U3,1)-zeroCount, :);
day=max(U3(:,2))+1;
iter=max(U3(:,1));
Ivar1v=U3(:,5)+U3(:,6);
Ivar1m=reshape(Ivar1v,day,iter);
Ivar2v=U3(:,8)+U3(:,9);
Ivar2m=reshape(Ivar2v,day,iter);
qIvar1m=quantile(Ivar1m, [0.25,0.5,0.75],2);
qIvar2m=quantile(Ivar2m, [0.25,0.5,0.75],2);
meanIvar1m=mean(Ivar1m,2);
meanIvar2m=mean(Ivar2m,2);


totalInfectedMedian=qIvar1m+qIvar2m;
percentageSecondVariant=qIvar2m./totalInfectedMedian;

meanPercentageSecondVariant=meanIvar2m./(meanIvar1m+meanIvar2m);

% tau=0.03;
% newVar=1.6;
% emergenceDay=150;
plotday=400;
f=figure
% yyaxis left
p1=plot(1:plotday, qIvar1m(1:plotday,2),'Color', [0.9290 0.6940 0.1250],  'LineWidth',4)
hold on
plot(1:plotday, qIvar1m(1:plotday,1),'Color', [0.9290 0.6940 0.1250],  'LineWidth',2.5, 'LineStyle',':' )
hold on
plot(1:plotday, qIvar1m(1:plotday,3),'Color', [0.9290 0.6940 0.1250],  'LineWidth',2.5,'LineStyle',':')

hold on
p2=plot(emergenceDay:plotday, qIvar2m(emergenceDay:plotday,2),'Color', [0.6350 0.0780 0.1840],  'LineWidth',4)
hold on
plot(emergenceDay:plotday, qIvar2m(emergenceDay:plotday,1),'Color', [0.6350 0.0780 0.1840],  'LineWidth',2.5, 'LineStyle',':' )
hold on
plot(emergenceDay:plotday, qIvar2m(emergenceDay:plotday,3),'Color', [0.6350 0.0780 0.1840],  'LineWidth',2.5,'LineStyle',':')
hold on
%ylim([0 5000])
title(strcat('\beta_1=',num2str(tau),',',' seedV2=', seed, ', EmergenceDay=',num2str(emergenceDay)))
xlabel('day')
ylabel('Infected population')
legend([p1 p2],{'variant-1','variant-2'},'Location','northwest')
%ylim([0 15000])

yyaxis right
p3=plot(emergenceDay:plotday, percentageSecondVariant(emergenceDay:plotday,2)*100,'Color', [0 0.4470 0.7410],  'LineWidth',2)
ylabel('Percentage of variant-2', 'Color', [0 0.4470 0.7410])
hold on
plot(emergenceDay:plotday, percentageSecondVariant(emergenceDay:plotday,1)*100,'Color', [0 0.4470 0.7410],  'LineWidth',2,'LineStyle',':')
hold on
plot(emergenceDay:plotday, percentageSecondVariant(emergenceDay:plotday,3)*100,'Color', [0 0.4470 0.7410],  'LineWidth',2,'LineStyle',':')
legend([p1 p2 p3],{'variant-1','variant-2','percentage of variant-2'})
ylim([0 70])


ax = gca; % current axes
ax.FontSize = 20;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
set(f,'Position',[10 10 730 400])
set(gca,'YColor',[0 0.4470 0.7410])

hold off
saveas(gcf,strcat('Popuday_Tau',num2str(tau),'seedLoc', num2str(seedLoc), 'emergencDay',num2str(emergenceDay),'.png'))


% figure
% 
% cIvar2m=cumsum(Ivar2m);
% cIvar1m=cumsum(Ivar1m);
% hb=bar([cIvar1m(end,:); cIvar2m(end,:)]')
% set( hb, {'DisplayName'}, {'variant-1'; 'variant-2'} );
% legend
% title(strcat('\beta_1=',num2str(tau),',',' seedV2=', seed, ', EmergenceDay=',num2str(emergenceDay)))
% ylabel('Cumulative infected at day 400')
% xlabel('EpiHiper replicate')
% ax = gca; % current axes
% ax.FontSize = 16;
% 
% saveas(gcf,strcat('BarPopuday_Tau',num2str(tau),'seedLoc', num2str(seedLoc), 'emergencDay',num2str(emergenceDay),'.png'))


s=Ivar2m(end,:);
end


