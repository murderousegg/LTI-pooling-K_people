clc
clear all
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
vec = [2 5 10 15]; %[1 5 10 15]; 
Mar = ['o','+','*','x','v']; %'o','*','x','v'
city='SF';
mult = [0.0156 0.0312 0.0625 0.125 0.25 0.5 1 2]; %
colors=lines;
CM=[colors(1,:);colors(180,:);colors(220,:);colors(30,:);colors(256,:);]; %

%% group together
for Delay = vec
    
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
h(4)=plot(NaN, NaN,'Marker',Mar(1),'MarkerSize',8,'Color','k','Linestyle', 'none');
h(5)=plot(NaN, NaN,'Marker',Mar(2),'MarkerSize',8,'Color','k','Linestyle', 'none');
h(6)=plot(NaN, NaN,'Marker',Mar(3),'MarkerSize',8,'Color','k','Linestyle', 'none');
vec  = [5 10 15];
for Delay = vec
for WaitingTime= vec
    
load(strcat(city,'\Results\New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'),'TrackDems','objs','Improv')
set(gca,'ticklabelinterpreter','Latex','fontsize',13,'XTickLabel',[], 'XScale', 'log') %'XScale', 'log',
plot(TrackDems(:,1), 100- 100*TrackDems(:,2)./TrackDems(:,1) ,'Marker',Mar(find(vec==Delay)),'MarkerSize',8,'Color',CM(find(vec==WaitingTime),:),'LineWidth',1);
%xlabel('')
set(gca,'ticklabelinterpreter','Latex','fontsize',13)
ylabel('Pooled Rides [\%]')
ylim([40 100])
xlim([234 30050])
end
end

set(gca,'ticklabelinterpreter','Latex','fontsize',14)
legend('$\bar{t}=5$ min','$\bar{t}=10$ min','$\bar{t}=15$ min','$\bar{\delta}=5$ min','$\bar{\delta}=10$ min','$\bar{\delta}=15$ min','NumColumns',2)

%%
subplot(2,1,2)
hold on; grid on; box on;
vec = [5 10 15];
for Delay = vec
for WaitingTime= vec

load(strcat(city,'\Results\New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'))
set(gca,'ticklabelinterpreter','Latex','fontsize',13)%,'XTickLabel',[]);
plot(TrackDems(:,1), sum(Del./sum(TotG,2),2) ,'Marker',Mar(find(vec==Delay)),'MarkerSize',8,'Color',CM(find(vec==WaitingTime),:),'LineWidth',1);
%yline(Delay,'LineWidth',1.5)
xlabel('Demands per Hour','fontsize',13)
set(gca,'ticklabelinterpreter','Latex','fontsize',13, 'XScale', 'log','YScale', 'lin')
set(gca, 'ytick', 0:0.5:2);
ylabel('Average Delay [min]')
%ylim([0 2.1])
xlim([234 30050])%xlim([230 30050])
end
end
%legend(h,'$\bar{t}=1$ min','$\bar{t}=2$ min','$\bar{t}=5$ min','$\bar{t}=10$ min','$\bar{t}=15$ min','$\bar{\delta}=1$ min','$\bar{\delta}=2$ min','$\bar{\delta}=5$ min','$\bar{\delta}=10$ min','$\bar{\delta}=15$ min')
%set(gca,'ticklabelinterpreter','Latex','fontsize',16)

%%
figure(2)
load('SF\Results\New_Delay10WTime2.mat')
subplot(1,2,1)
PercNRP = TrackDems(:,2)./(TrackDems(:,3) + TrackDems(:,2));
bar([PercNRP TotG./(sum(TotG,2)).*(1-PercNRP)],'stacked','FaceColor','flat')
grid on; box on;
ylim([0 1 ])
ylabel('Composition','FontSize', 16,'interpreter','latex')
%xlabel('Demands per Hour','FontSize', 16,'interpreter','latex')
set(gca,'FontSize',16,'XTickLabel',num2str(round(TrackDems(:,1),-1)) )
%legend('Single person','Two people','Three people','Four people','FontSize',16)
title('$\bar{\delta}= 10$min, $ \bar{t} = 2$  min')
set(gca,'FontSize',16)

load('SF\Results\New_Delay10WTime5.mat')
subplot(1,2,2)
PercNRP = TrackDems(:,2)./(TrackDems(:,3) + TrackDems(:,2));
bar([PercNRP TotG./(sum(TotG,2)).*(1-PercNRP)],'stacked')
grid on; box on;
ylim([0 1 ])
%ylabel('Composition','FontSize', 16,'interpreter','latex')
%xlabel('Demands per Hour','FontSize', 16,'interpreter','latex')
set(gca,'FontSize',16,'XTickLabel',num2str(round(TrackDems(:,1),-1)) )
legend('Single person','Two people','Three people','Four people','FontSize',16)
title('$\bar{\delta}= 10$min, $ \bar{t} = 5$  min')
set(gca,'FontSize',16)

figure(3)
load('SF\Results\New_Delay10WTime10.mat')
subplot(1,2,1)
PercNRP = TrackDems(:,2)./(TrackDems(:,3) + TrackDems(:,2));
bar([PercNRP TotG./(sum(TotG,2)).*(1-PercNRP)],'stacked','FaceColor','flat')
grid on; box on;
ylim([0 1 ])
ylabel('Composition','FontSize', 16,'interpreter','latex')
xlabel('Demands per Hour','FontSize', 16,'interpreter','latex')
set(gca,'FontSize',16,'XTickLabel',num2str(round(TrackDems(:,1),-1)) )
%legend('Single person','Two people','Three people','Four people','FontSize',16)
title('$\bar{\delta}= 10$min, $ \bar{t} = 10$  min')
set(gca,'FontSize',16)

load('SF\Results\New_Delay10WTime15.mat')
subplot(1,2,2)
PercNRP = TrackDems(:,2)./(TrackDems(:,3) + TrackDems(:,2));
bar([PercNRP TotG./(sum(TotG,2)).*(1-PercNRP)],'stacked')
grid on; box on;
ylim([0 1 ])
%ylabel('Composition','FontSize', 16,'interpreter','latex')
xlabel('Demands per Hour','FontSize', 16,'interpreter','latex')
set(gca,'FontSize',16,'XTickLabel',num2str(round(TrackDems(:,1),-1)) )
%legend('Single person','Two people','Three people','Four people','FontSize',16)
title('$\bar{\delta}= 10$min, $ \bar{t} = 15$  min')
set(gca,'FontSize',16)