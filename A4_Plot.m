clc
clear all
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
vec = [2 5 10 15]; %[1 5 10 15]; 
Mar = ['o','+','*','x','v']; %'o','*','x','v'
city='NYC357';
mult = [0.0078 0.0156 0.0312 0.0625 0.125 0.25 0.5 1 2]; %
colors=lines;
CM=[colors(1,:);colors(180,:);colors(220,:);colors(30,:);colors(256,:);]; %

%% group together
for Delay =vec%[2 5]% 
    
    for WaitingTime= vec %[1 2 5 10 15]%vec
        Improv = [];
        objs = [];
        TrackDems = [];
        TotG = [];
        Del = [];
        for multiplicator = mult %   0.01 1.6 0.005 0.01 0.03 0.05 0.1 0.2
            A=load(strcat(city,'/Results/Ppl_4Delay',num2str(Delay),'WTime',num2str(WaitingTime),'Dem',num2str(multiplicator),'.mat'),'TrackDems_temp','Sol','Cumul_delay','TotGamma');
            Improv = [Improv  1-(A.Sol(2).obj+A.Sol(3).obj)/A.Sol(1).obj];
            objs = [objs; A.Sol.obj];
            TrackDems = [TrackDems; A.TrackDems_temp];
            TotG = [TotG; A.TotGamma(1) A.TotGamma(2) A.TotGamma(3)];
            Del = [Del; A.Cumul_delay(1) A.Cumul_delay(2) A.Cumul_delay(3)];
        end
        save(strcat(city,'/Results/New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'),'TrackDems','objs','Improv','TotG','Del')
    end
end

%%

% subplot(2,1,2)
% hold on; grid on; box on;
% 
% count2 =0;
% for Delay = vec
%     count=0;
%     count2 = count2+1;
%     
%     
%     for WaitingTime= vec
%         count = count+1;
%         
%         load(strcat(city,'/Results/New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'),'Improv')
%         set(gca,'ticklabelinterpreter','Latex','fontsize',16)
%         plot(TrackDems(:,1), Improv*100,'Color',CM(count,:),'Marker',Mar(count2),'MarkerSize',8,'LineWidth',1);
%         xlabel('Requests per Hour')
%         set(gca,'ticklabelinterpreter','Latex','fontsize',16)
%         ylabel('Improvement [\%]')
%         %ylim([0 50])
%         %xlim([0 2.2e4])
%         
%     end
% end
%%
subplot(2,1,1)
%figure()
hold on; grid on; box on;
%h(1)=plot(NaN, NaN);
h(1)=plot(NaN, NaN,'Color',CM(1,:),'LineWidth',2);
h(2)=plot(NaN, NaN,'Color',CM(2,:),'LineWidth',2);
h(3)=plot(NaN, NaN,'Color',CM(3,:),'LineWidth',2);
h(4)=plot(NaN, NaN,'Color',CM(4,:),'LineWidth',2);
%h(5)=plot(NaN, NaN,'Color',CM(5,:),'LineWidth',2);
h(6)=plot(NaN, NaN,'Marker',Mar(1),'MarkerSize',8,'Color','k','Linestyle', 'none');
h(7)=plot(NaN, NaN,'Marker',Mar(2),'MarkerSize',8,'Color','k','Linestyle', 'none');
h(8)=plot(NaN, NaN,'Marker',Mar(3),'MarkerSize',8,'Color','k','Linestyle', 'none');
h(9)=plot(NaN, NaN,'Marker',Mar(4),'MarkerSize',8,'Color','k','Linestyle', 'none');
%h(11)=plot(NaN, NaN,'Marker',Mar(5),'MarkerSize',8,'Color','k','Linestyle', 'none');
h(5)=plot(NaN, NaN,'Color','k','Linestyle', '--','LineWidth',2);

for Delay = vec
for WaitingTime= vec
    
load(strcat(city,'\Results\New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'),'TrackDems','objs','Improv')
set(gca,'ticklabelinterpreter','Latex','fontsize',13,'XTickLabel',[], 'XScale', 'log') %'XScale', 'log',
plot(TrackDems(:,1), 100- 100*TrackDems(:,2)./TrackDems(:,1) ,'Marker',Mar(find(vec==Delay)),'MarkerSize',8,'Color',CM(find(vec==WaitingTime),:),'LineWidth',1);
%xlabel('')
set(gca,'ticklabelinterpreter','Latex','fontsize',13)
ylabel('Pooled Rides [\%]')
ylim([40 100])
xlim([800 107000])
end
end

% plot ref from Santi 
kk = 1.1*10^(-4);
nn=0.92;
dems = [0:0.001:3]*10^6;
fun = 100*(kk*dems.^nn)./(1 + (kk*dems.^nn));
plot(dems/24,fun,'k','LineWidth',2,'Linestyle','--')
%
set(gca,'ticklabelinterpreter','Latex','fontsize',14)
legend(h,'$\bar{t}=2$ min','$\bar{t}=5$ min','$\bar{t}=10$ min','$\bar{t}=15$ min','Santi et al. [15]','$\bar{\delta}=2$ min','$\bar{\delta}=5$ min','$\bar{\delta}=10$ min','$\bar{\delta}=15$ min','NumColumns',2)
%legend(h,'$\bar{t}=2$ min','$\bar{t}=5$ min','$\bar{t}=10$ min','Microscopic approach [25]','$\bar{\delta}=2$ min','$\bar{\delta}=5$ min','$\bar{\delta}=10$ min','NumColumns',2)
%legend(h,'$\bar{t}=5$ min','$\bar{t}=10$ min','$\bar{t}=15$ min','$\bar{\delta}=5$ min','$\bar{\delta}=10$ min','$\bar{\delta}=15$ min','NumColumns',2)

%%
subplot(2,1,2)
hold on; grid on; box on;
for Delay = vec
for WaitingTime= vec

load(strcat(city,'\Results\New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'))
set(gca,'ticklabelinterpreter','Latex','fontsize',13)%,'XTickLabel',[]);
plot(TrackDems(:,1), Del./(sum(TotG,2)) ,'Marker',Mar(find(vec==Delay)),'MarkerSize',8,'Color',CM(find(vec==WaitingTime),:),'LineWidth',1);
%yline(Delay,'LineWidth',1.5)
xlabel('Demands per Hour','fontsize',13)
set(gca,'ticklabelinterpreter','Latex','fontsize',13, 'XScale', 'log','YScale', 'lin')
set(gca, 'ytick', 0:0.5:2);
ylabel('Average Delay [min]')
%ylim([0 2.1])
xlim([800 107000])%xlim([230 30050])
end
end
%legend(h,'$\bar{t}=1$ min','$\bar{t}=2$ min','$\bar{t}=5$ min','$\bar{t}=10$ min','$\bar{t}=15$ min','$\bar{\delta}=1$ min','$\bar{\delta}=2$ min','$\bar{\delta}=5$ min','$\bar{\delta}=10$ min','$\bar{\delta}=15$ min')
%set(gca,'ticklabelinterpreter','Latex','fontsize',16)
%%
city='NYC357'
fig=figure()
for iii=[1 4 6 8] % 7 8
Matt=[];
for Delay = vec
    temp = [];
for WaitingTime= vec
  
load(strcat(city,'\Results\New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'),'TrackDems','objs','Improv')
num=iii;
temp = [temp 100*(-objs(num,2)-objs(num,3) + objs(num,1))/objs(num,1)];
end
Matt = [Matt; temp];
end
Delay = vec;
WaitingTime= vec;
[XX,YY]=meshgrid(Delay,WaitingTime);
nexttile
Dems = roundn(TrackDems(iii,1),2);

contourf(XX,YY,Matt)
caxis([0,50]); 
title(strcat ( num2str(Dems) ,' Demands/h'),'FontSize', 22,'interpreter','latex' )
caxis([0,50]);
set(gca,'ticklabelinterpreter','Latex','fontsize',22)%,'interpreter','latex')
box on;
grid on;

end
h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on'; 
cc=colorbar(h,'ticklabelinterpreter','Latex','fontsize',18)
cc.Label.String="Improvement [\%]";
cc.Position= [0.925 0.11 0.025 0.82];
cc.Label.Interpreter = 'Latex';
caxis(h,[10,50]); 
%axis equal
ylabel(h,'Maximum Delay $\bar{\delta}$ [min]','FontSize', 32,'interpreter','latex')
xlabel(h,'Maximum Waiting Time $\bar{t}$ [min]','FontSize', 32,'interpreter','latex')
set(gca,'FontSize',28)
cc.Label.FontSize=30;
%%
figure(2)
%%
subplot(1,2,2)
bar(TotG./(sum(TotG,2)),'stacked')
grid on; box on;
ylim([0 1 ])
%ylabel('Composition','FontSize', 16,'interpreter','latex')
xlabel('Demands per Hour','FontSize', 16,'interpreter','latex')
set(gca,'FontSize',16,'XTickLabel',num2str(round(TrackDems(:,1),-1)) )
%legend('Two people','Three people','Four people','FontSize',16)
title('$\bar{\delta}= 10$min, $ \bar{t} = 15$  min')
set(gca,'FontSize',16)
%%
% h = axes(fig,'visible','off'); 
% h.Title.Visible = 'on';
% h.XLabel.Visible = 'on';
% h.YLabel.Visible = 'on'; 
