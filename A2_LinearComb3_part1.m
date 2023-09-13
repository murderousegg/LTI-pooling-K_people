clc
clear all
city='SF';
load(strcat(city,'/Graphs.mat'));
load(strcat(city,'/solPart_',city,'.mat'));
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

%% Layer 3
mkdir(strcat(city,'/L3'))
for jj1=1:N_nodes %18
    parfor ii1=1:N_nodes%15:N_nodes %1
        ii1
        sol3_LC_temp= zeros(20^5,16);
        counter=1;
        for jj2= jj1 : N_nodes
            for ii2 = ii1 : N_nodes
                for jj3 = jj2 : N_nodes
                    for ii3 = ii2 : N_nodes
                        if ~any([ii1==jj1,ii1==jj2,ii1==jj3,ii2==jj1,ii2==jj2,ii2==jj3,ii3==jj1,ii3==jj2,ii3==jj3,DemandS(ii1,jj1)==0,DemandS(ii2,jj2)==0,DemandS(ii3,jj3)==0])
                            opti = sortrows([LTIFM3_SP(jj1,ii1,jj2,ii2,jj3,ii3,solPart);
                                             LTIFM3_SP(jj2,ii2,jj1,ii1,jj3,ii3,solPart);
                                             LTIFM3_SP(jj3,ii3,jj2,ii2,jj1,ii1,solPart)]);
                            sol3_LC_temp(counter,:) = opti(1,:);
                            counter = counter +1;
                        end
                    end
                end
            end
        end
    sol3_LC_temp( sol3_LC_temp(:,1)==0,: ) = []; %filter out zero obj
    sol3_LC_temp( sol3_LC_temp(:,2) > 20,: ) = []; %filter out above 20 delay
    sol3_LC_temp( sol3_LC_temp(:,3) > 20,: ) = []; %filter out above 20 delay
    sol3_LC_temp( sol3_LC_temp(:,4) > 20,: ) = []; %filter out above 20 delay
    sol3_LC_temp = sortrows(sol3_LC_temp);
    parsave3(strcat(city,'/L3/MatL3_',num2str(jj1),'_',num2str(ii1),'.mat'),sol3_LC_temp)%, '-v7.3')
    end
end