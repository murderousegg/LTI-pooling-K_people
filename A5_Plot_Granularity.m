clc
clear all
set(gca,'ticklabelinterpreter','Latex','fontsize',16)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
Mar = ['o','+','*','v','x'];

%% group together
vector = [20 40 80 120 160 250 357]; % 120 200 250

for Delay = [2 5 10] %%!!!!!!
for WaitingTime = [2 5 10]
    multiplicator = 1;
Soll = [];
TrackDems = [];
TotG = [];
Del = [];
improv = [];
numpool = [];   
for city = vector
   
A =load(strcat('NYC',num2str(city),'/Results\Ppl_2Delay',num2str(Delay),'WTime',num2str(WaitingTime),'Dem',num2str(multiplicator),'.mat'),'TrackDems_temp','Sol','Cumul_delay','TotGamma');
Soll = [Soll; A.Sol];
TrackDems = [TrackDems; A.TrackDems_temp];
TotG = [TotG A.TotGamma];
Del = [Del A.Cumul_delay];
improv = [improv, 1-(A.Sol(2).obj+A.Sol(3).obj)/A.Sol(1).obj];
numpool = [numpool, 100-100*A.TrackDems_temp(2)./A.TrackDems_temp(1)];
end
save(strcat('Utilities\New_Delay',num2str(Delay),'WTime',num2str(WaitingTime),'.mat'),'TrackDems','Soll','TotG','Del','improv','numpool')
end
end
%%

colors=lines;
CM=[colors(1,:);colors(30,:);colors(180,:);colors(220,:);colors(256,:)];
t = tiledlayout(1,1)

for Dela = [5]
    nexttile
    hold on; grid on; box on;
    set(gca,'ticklabelinterpreter','Latex','fontsize',16)
    improv = [];
    numpool = [];
    for Wait  = [2 5 10]
        A = load(strcat('Utilities\New_Delay',num2str(Dela),'WTime',num2str(Wait),'.mat'));    
        improv = [improv; A.improv];
        numpool = [numpool; A.numpool];
    end
    plot(vector,improv(1,:)*100,'Marker',Mar(1),'MarkerSize',8,'LineWidth',1);
    hold on;
    plot(vector,improv(2,:)*100,'Marker',Mar(1),'MarkerSize',8,'LineWidth',1);
    plot(vector,improv(3,:)*100,'Marker',Mar(1),'MarkerSize',8,'LineWidth',1);
    hold off;
    if Dela ~= 5
    xticklabels({''});
    end
    ylim([0 55])
    %xlabel('Nodes') %% !!!!
    %ylabel('Improvement [\%]')
    %plot(vector, 100- 100*TrackDems(:,2)./TrackDems(:,1) ,'Marker',Mar(2),'MarkerSize',8,'LineWidth',1);
    % legend('$\bar{t} = 2$ min', '$\bar{t} = 5$ min', '$\bar{t} = 10$ min')
    set(gca,'ticklabelinterpreter','Latex','fontsize',14)
    title(strcat('$\bar{\delta} = $',num2str(Dela),'min'))
    set(gca,'ticklabelinterpreter','Latex','fontsize',14)
end
set(gca,'ticklabelinterpreter','Latex','fontsize',14)
legend('$\bar{t} = 2$ min', '$\bar{t} = 5$ min', '$\bar{t} = 10$ min')
xlabel(t,'Nodes','FontSize', 14,'interpreter','latex')
ylabel(t,'Improvement [\%]','FontSize', 14,'interpreter','latex')
set(gca,'ticklabelinterpreter','Latex','fontsize',14)