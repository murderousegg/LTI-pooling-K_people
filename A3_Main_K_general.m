clc
clear all
city='NYC357';
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);
if strcmp('SF', city)
DemandS = DemandS/24; %%daily to hourly only for SF!!!!!!!
end
OriginalDemand= full(DemandS);
DemandS = full(DemandS);
TotDems = sum(DemandS,'all');
%mkdir(strcat(city,'/Results')) uncomment if folder does not exist yet
%% Layer 2
ppl = 2;

if ppl ==2
    load(strcat(city,'/MatL2.mat'))
    size2 = size(sol2_LC,1);
    Sol2 = [sol2_LC(:,1:3),NaN(size2,2),sol2_LC(:,4:7),NaN(size2,4),sol2_LC(:,8:11),NaN(size2,4)];
    FullList = [Sol2];
elseif ppl ==3
    load(strcat(city,'/MatL2.mat'))
    size2 = size(sol2_LC,1);
    Sol2 = [sol2_LC(:,1:3),NaN(size2,2),sol2_LC(:,4:7),NaN(size2,4),sol2_LC(:,8:11),NaN(size2,4)];
    
    load(strcat(city,'/MatL3.mat'))
    size3 = size(sol3_LC,1);
    Sol3 = [sol3_LC(:,1:4),NaN(size3,1),sol3_LC(:,5:10),NaN(size3,2),sol3_LC(:,11:16),NaN(size3,2)];
    
    FullList = [Sol2;Sol3];
elseif ppl ==4
    load(strcat(city,'/MatL2.mat'))
    size2 = size(sol2_LC,1);
    Sol2 = [sol2_LC(:,1:3),NaN(size2,2),sol2_LC(:,4:7),NaN(size2,4),sol2_LC(:,8:11),NaN(size2,4)];
    
    load(strcat(city,'/MatL3.mat'))
    size3 = size(sol3_LC,1);
    Sol3 = [sol3_LC(:,1:4),NaN(size3,1),sol3_LC(:,5:10),NaN(size3,2),sol3_LC(:,11:16),NaN(size3,2)];

    load(strcat(city,'/MatL4.mat'))
    FullList = [Sol2;Sol3;sol4_LC];
end

FullList = sortrows(FullList,1);

for iiii = 1:size(FullList,1)
    vect =  FullList(iiii,6:13)';
    vectR = vect(~isnan( vect));
    vectR = reshape(vectR,2,[])';
    num = size(vectR,1);
    for iii=1:num
    if DemandS(vectR(iii,2),vectR(iii,1)) == 0
    FullList(iiii,1) = 0;
    end
    end
end
FullList(FullList(:,1) >= -0.01,:) = [];
%%
CountGamma = [];
for WaitingTime = [2] %in min 2 5 10 15
    for Delay = [10] % in min
        for mult= [0.125] %0.0078 0.0156 0.0312 0.0625 0.125 0.25 0.5 1 2
            gamma = 0;
            TotGamma2 = 0;
            TotGamma3 = 0;
            TotGamma4 = 0;
            Cumul_delay2 = 0;
            Cumul_delay3 = 0;
            Cumul_delay4 = 0;
            DemandS =  mult* OriginalDemand;
            Demands_rp = zeros(N_nodes,N_nodes);
            
            for iii = 1:size(FullList,1)
                if isnan(gamma)
                iii
                break
                end
                if isnan(FullList(iii,4)) && isnan(FullList(iii,5)) && FullList(iii,2) < Delay && FullList(iii,3) < Delay && DemandS(FullList(iii,7),FullList(iii,6)) >= 10e-5 && DemandS(FullList(iii,9),FullList(iii,8)) >= 10e-5
                    
                    jj1 = FullList(iii,6);ii1 = FullList(iii,7);jj2 = FullList(iii,8); ii2 = FullList(iii,9);
                    gamma = min([DemandS(ii1,jj1),DemandS(ii2,jj2)])*probcombN([DemandS(ii1,jj1),DemandS(ii2,jj2)],WaitingTime)/2;
                    
                    Gamma0 = zeros(N_nodes,N_nodes);
                    if FullList(iii,14:17) == [1 2 1 2]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(ii1,jj2) = 1;
                        Gamma0(ii2,ii1) = 1;
                    elseif FullList(iii,14:17) == [1 2 2 1]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(ii2,jj2) = 1;
                        Gamma0(ii1,ii2) = 1;
                    end
                    
                    matrow = [jj1, ii1; jj2, ii2];
                    multip = size(unique(matrow, 'rows', 'first'),1);
                    Demands_rp = Demands_rp +  multip* gamma* Gamma0;
                    DemandS(ii1,jj1) = DemandS(ii1,jj1) - multip*gamma;
                    DemandS(ii2,jj2) = DemandS(ii2,jj2) - multip*gamma;
                    Cumul_delay2 = Cumul_delay2 + multip*gamma* (FullList(iii,2) + FullList(iii,3)) ;
                    TotGamma2 = TotGamma2 + multip*gamma;
                    %CountGamma = [CountGamma gamma];
                    
                elseif ~isnan(FullList(iii,4)) && isnan(FullList(iii,5)) && FullList(iii,2) < Delay && FullList(iii,3) < Delay && FullList(iii,4) < Delay && DemandS(FullList(iii,7),FullList(iii,6))>= 10e-5 && DemandS(FullList(iii,9),FullList(iii,8))>= 10e-5 && DemandS(FullList(iii,11),FullList(iii,10))>= 10e-5
                    
                    jj1 = FullList(iii,6);ii1 = FullList(iii,7);jj2 = FullList(iii,8); ii2 = FullList(iii,9);jj3 = FullList(iii,10); ii3 = FullList(iii,11);
                    gamma = min([DemandS(ii1,jj1),DemandS(ii2,jj2),DemandS(ii3,jj3)])*probcombN([DemandS(ii1,jj1),DemandS(ii2,jj2),DemandS(ii3,jj3)],WaitingTime)/3;
                    Gamma0 = zeros(N_nodes,N_nodes);
                    
                    if FullList(iii,14:19) == [1 2 3 1 2 3]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(jj3,jj2) = 1;
                        Gamma0(ii1,jj3) = 1;
                        Gamma0(ii2,ii1) = 1;
                        Gamma0(ii3,ii2) = 1;
                    elseif FullList(iii,14:19) == [1 2 3 1 3 2]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(jj3,jj2) = 1;
                        Gamma0(ii1,jj3) = 1;
                        Gamma0(ii3,ii1) = 1;
                        Gamma0(ii2,ii3) = 1;
                    elseif FullList(iii,14:19) == [1 2 3 2 1 3]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(jj3,jj2) = 1;
                        Gamma0(ii2,jj3) = 1;
                        Gamma0(ii1,ii2) = 1;
                        Gamma0(ii3,ii1) = 1;
                    elseif FullList(iii,14:19) == [1 2 3 2 3 1]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(jj3,jj2) = 1;
                        Gamma0(ii2,jj3) = 1;
                        Gamma0(ii3,ii2) = 1;
                        Gamma0(ii1,ii3) = 1;
                    elseif FullList(iii,14:19) == [1 2 3 3 1 2]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(jj3,jj2) = 1;
                        Gamma0(ii3,jj3) = 1;
                        Gamma0(ii1,ii3) = 1;
                        Gamma0(ii2,ii1) = 1;
                    elseif FullList(iii,14:19) == [1 2 3 3 2 1]
                        Gamma0(jj2,jj1) = 1;
                        Gamma0(jj3,jj2) = 1;
                        Gamma0(ii3,jj3) = 1;
                        Gamma0(ii2,ii3) = 1;
                        Gamma0(ii1,ii2) = 1;
                    end
                    %%% da qui
                    matrow = [jj1, ii1; jj2, ii2; jj3,ii3];
                    multip = size(unique(matrow, 'rows', 'first'),1);
                    Demands_rp = Demands_rp + multip*gamma*Gamma0;
                    DemandS(ii1,jj1) = DemandS(ii1,jj1) - multip*gamma;
                    DemandS(ii2,jj2) = DemandS(ii2,jj2) - multip*gamma;
                    DemandS(ii3,jj3) = DemandS(ii3,jj3) - multip*gamma;
                    Cumul_delay3 = Cumul_delay3 + multip*gamma*(FullList(iii,2) + FullList(iii,3) + FullList(iii,4) );
                    TotGamma3 = TotGamma3 + gamma;
                    
                elseif FullList(iii,2) < Delay && FullList(iii,3) < Delay && FullList(iii,4) < Delay && FullList(iii,5) < Delay && DemandS(FullList(iii,7),FullList(iii,6))>= 10e-5 && DemandS(FullList(iii,9),FullList(iii,8))>= 10e-5 && DemandS(FullList(iii,11),FullList(iii,10)) >= 10e-5 && DemandS(FullList(iii,13),FullList(iii,12))>= 10e-5
                    jj1 = FullList(iii,6);ii1 = FullList(iii,7);jj2 = FullList(iii,8); ii2 = FullList(iii,9);jj3 = FullList(iii,10); ii3 = FullList(iii,11);jj4 = FullList(iii,12); ii4 = FullList(iii,13);
                    
                        gamma = min([DemandS(ii1,jj1),DemandS(ii2,jj2),DemandS(ii3,jj3),DemandS(ii4,jj4)])*probcombN([DemandS(ii1,jj1),DemandS(ii2,jj2),DemandS(ii3,jj3),DemandS(ii4,jj4)],WaitingTime)/4;
                        Gamma0 = zeros(N_nodes,N_nodes);
                        
                        if FullList(iii,14:21) == [1 2 3 4 1 2 3 4]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii1,jj4) = 1;
                            Gamma0(ii2,ii1) = 1;
                            Gamma0(ii3,ii2) = 1;
                            Gamma0(ii4,ii3) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 1 3 2 4]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii1,jj4) = 1;
                            Gamma0(ii3,ii1) = 1;
                            Gamma0(ii2,ii3) = 1;
                            Gamma0(ii4,ii2) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 1 3 4 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii1,jj4) = 1;
                            Gamma0(ii3,ii1) = 1;
                            Gamma0(ii4,ii3) = 1;
                            Gamma0(ii2,ii4) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 1 2 4 3]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii1,jj4) = 1;
                            Gamma0(ii2,ii1) = 1;
                            Gamma0(ii4,ii2) = 1;
                            Gamma0(ii3,ii4) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 1 4 2 3]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii1,jj4) = 1;
                            Gamma0(ii4,ii1) = 1;
                            Gamma0(ii2,ii4) = 1;
                            Gamma0(ii3,ii2) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 1 4 3 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii1,jj4) = 1;
                            Gamma0(ii4,ii1) = 1;
                            Gamma0(ii3,ii4) = 1;
                            Gamma0(ii2,ii3) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 2 1 3 4]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii2,jj4) = 1;
                            Gamma0(ii1,ii2) = 1;
                            Gamma0(ii3,ii1) = 1;
                            Gamma0(ii4,ii3) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 2 3 1 4]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii2,jj4) = 1;
                            Gamma0(ii3,ii2) = 1;
                            Gamma0(ii1,ii3) = 1;
                            Gamma0(ii4,ii1) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 2 1 4 3]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii2,jj4) = 1;
                            Gamma0(ii1,ii2) = 1;
                            Gamma0(ii4,ii1) = 1;
                            Gamma0(ii3,ii4) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 2 3 4 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii2,jj4) = 1;
                            Gamma0(ii3,ii2) = 1;
                            Gamma0(ii4,ii3) = 1;
                            Gamma0(ii1,ii4) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 2 4 3 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii2,jj4) = 1;
                            Gamma0(ii4,ii2) = 1;
                            Gamma0(ii3,ii4) = 1;
                            Gamma0(ii1,ii3) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 2 4 1 3]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii2,jj4) = 1;
                            Gamma0(ii4,ii2) = 1;
                            Gamma0(ii1,ii4) = 1;
                            Gamma0(ii3,ii1) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 3 1 2 4]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii3,jj4) = 1;
                            Gamma0(ii1,ii3) = 1;
                            Gamma0(ii2,ii1) = 1;
                            Gamma0(ii4,ii2) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 3 1 4 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii3,jj4) = 1;
                            Gamma0(ii1,ii3) = 1;
                            Gamma0(ii4,ii1) = 1;
                            Gamma0(ii2,ii4) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 3 2 1 4]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii3,jj4) = 1;
                            Gamma0(ii2,ii3) = 1;
                            Gamma0(ii1,ii2) = 1;
                            Gamma0(ii4,ii1) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 3 2 4 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii3,jj4) = 1;
                            Gamma0(ii2,ii3) = 1;
                            Gamma0(ii4,ii2) = 1;
                            Gamma0(ii1,ii4) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 3 4 1 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii3,jj4) = 1;
                            Gamma0(ii4,ii3) = 1;
                            Gamma0(ii1,ii4) = 1;
                            Gamma0(ii2,ii1) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 3 4 2 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii3,jj4) = 1;
                            Gamma0(ii4,ii3) = 1;
                            Gamma0(ii2,ii4) = 1;
                            Gamma0(ii1,ii2) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 4 3 2 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii4,jj4) = 1;
                            Gamma0(ii3,ii4) = 1;
                            Gamma0(ii2,ii3) = 1;
                            Gamma0(ii1,ii2) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 4 3 1 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii4,jj4) = 1;
                            Gamma0(ii3,ii4) = 1;
                            Gamma0(ii1,ii3) = 1;
                            Gamma0(ii2,ii1) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 4 2 1 3]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii4,jj4) = 1;
                            Gamma0(ii2,ii4) = 1;
                            Gamma0(ii1,ii2) = 1;
                            Gamma0(ii3,ii1) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 4 2 3 1]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii4,jj4) = 1;
                            Gamma0(ii2,ii4) = 1;
                            Gamma0(ii3,ii2) = 1;
                            Gamma0(ii1,ii3) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 4 1 3 2]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii4,jj4) = 1;
                            Gamma0(ii1,ii4) = 1;
                            Gamma0(ii3,ii1) = 1;
                            Gamma0(ii2,ii3) = 1;
                        elseif FullList(iii,14:21) == [1 2 3 4 4 1 2 3]
                            Gamma0(jj2,jj1) = 1;
                            Gamma0(jj3,jj2) = 1;
                            Gamma0(jj4,jj3) = 1;
                            Gamma0(ii4,jj4) = 1;
                            Gamma0(ii1,ii4) = 1;
                            Gamma0(ii2,ii1) = 1;
                            Gamma0(ii3,ii2) = 1;
                        end
                        %
                        matrow = [jj1, ii1; jj2, ii2; jj3,ii3;jj4,ii4];
                        multip = size(unique(matrow, 'rows', 'first'),1);
                        Demands_rp = Demands_rp + multip*gamma* Gamma0;
                        DemandS(ii1,jj1) = DemandS(ii1,jj1) - multip*gamma;
                        DemandS(ii2,jj2) = DemandS(ii2,jj2) - multip*gamma;
                        DemandS(ii3,jj3) = DemandS(ii3,jj3) - multip*gamma;
                        DemandS(ii4,jj4) = DemandS(ii4,jj4) - multip*gamma;
                        Cumul_delay4 = Cumul_delay4 + multip*gamma*(FullList(iii,2)+ FullList(iii,3)+ FullList(iii,4)+ FullList(iii,5));
                        TotGamma4 = TotGamma4 + gamma;
                    end
                               
            end
            Demands_rp = Demands_rp - diag(diag(Demands_rp));
            solBase =LTIFM_reb(mult*OriginalDemand,city);
            solNP = LTIFM_reb(DemandS,city);
            solRP = LTIFM_reb(Demands_rp,city);
            TrackDems_temp = [sum(mult*OriginalDemand,'all'), sum(DemandS,'all'),sum(Demands_rp,'all')];
            Sol = [solBase, solNP, solRP];
            TotGamma = [TotGamma2, TotGamma3,TotGamma4];
            Cumul_delay = [Cumul_delay2, Cumul_delay3, Cumul_delay4];
            save(strcat(city,'/Results/','Ppl_',num2str(ppl),'Delay',num2str(Delay),'WTime',num2str(WaitingTime),'Dem',num2str(mult),'.mat'),'TrackDems_temp','Sol','Cumul_delay','TotGamma')
            
            
        end
    end
end
