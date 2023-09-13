clc
clear all
city='SF';
load(strcat(city,'/Graphs.mat'));
load(strcat(city,'/solPart_',city,'.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

%% Layer 4
%mkdir(strcat(city,'/L4')) %uncomment only the first time
for jj1=1%:N_nodes %18
    for ii1= 7:N_nodes %1
        parfor jj2= jj1 : N_nodes
            jj2
            sol4_LC_temp= zeros(100000,21);
            counter=1;
            for ii2 = ii1 : N_nodes
                ii2
                for jj3 = jj2 : N_nodes
                    for ii3 = ii2 : N_nodes
                        for jj4 = jj3 : N_nodes
                            for ii4 = ii3 : N_nodes
                                
                                if ~any([ii1==jj1,ii1==jj2,ii1==jj3,ii1==jj4,ii2==jj1,ii2==jj2,ii2==jj3,ii2==jj4,ii3==jj1,ii3==jj2,ii3==jj3,ii3==jj4,ii4==jj1,ii4==jj2,ii4==jj3,ii4==jj4,DemandS(ii1,jj1)==0,DemandS(ii2,jj2)==0,DemandS(ii3,jj3)==0,DemandS(ii4,jj4)==0])
                                    opti = sortrows([LTIFM4_SP(jj1,ii1,jj2,ii2,jj3,ii3,jj4,ii4,solPart);
                                        LTIFM4_SP(jj2,ii2,jj1,ii1,jj3,ii3,jj4,ii4,solPart);
                                        LTIFM4_SP(jj3,ii3,jj2,ii2,jj1,ii1,jj4,ii4,solPart);
                                        LTIFM4_SP(jj4,ii4,jj2,ii2,jj3,ii3,jj1,ii1,solPart)]);
                                    sol4_LC_temp(counter,:) = opti(1,:);
                                    counter = counter +1;
                                    
                                end
                            end
                        end
                    end
                end
            end
            sol4_LC_temp( sol4_LC_temp(:,1)==0,: ) = []; %filter out zero obj
            sol4_LC_temp( sol4_LC_temp(:,2) > 20,: ) = []; %filter out above 20 delay
            sol4_LC_temp( sol4_LC_temp(:,3) > 20,: ) = []; %filter out above 20 delay
            sol4_LC_temp( sol4_LC_temp(:,4) > 20,: ) = []; %filter out above 20 delay
            sol4_LC_temp( sol4_LC_temp(:,5) > 20,: ) = []; %filter out above 20 delay
            parsave4(strcat(city,'/L4/MatL4_',num2str(jj1),'_',num2str(ii1),'_',num2str(jj2),'.mat'),sol4_LC_temp)%, '-v7.3')
            
        end
    end
end
