clear all
city='SF';
load(strcat(city,'/Graphs.mat'));
Adj = adjacency(G_road);
Binc = incidence(G_road); 
[N_nodes,N_edges]=size(Binc);

%% Layer 4 join together last 2 digits.
sol4_LC = [];
for jj1= 1 : N_nodes
    jj1
    sol4_LC = [];
    for ii1 = 1: N_nodes
        ii1
        for jj2 = jj1 : N_nodes
            load(strcat(city,'/L4/MatL4_',num2str(jj1),'_',num2str(ii1),'_',num2str(jj2),'.mat'),'sol4_LC_temp')%, '-v7.3')
            sol4_LC = [sol4_LC; sol4_LC_temp];
        end
    end
    save(strcat(city,'/L4/MatL4_',num2str(jj1),'.mat'),'sol4_LC')%, '-v7.3')
end
%sol4_LC = sortrows(sol4_LC);

sol4_LC = [];
for jj1= 1 : N_nodes
            A = load(strcat(city,'/L4/MatL4_',num2str(jj1),'.mat'));%, '-v7.3')
            sol4_LC = [sol4_LC; A.sol4_LC];
end

sol4_LC = sortrows(sol4_LC);
save(strcat(city,'/MatL4.mat'),'sol4_LC')%, '-v7.3')
%
