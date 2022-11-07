clear all
%how to make results variable
j = matfile('twininfo_997subj.mat');;
allfiles = dir('*yeo.mat');
%load testoutput files P_test_%_yeo.mat
allfiles = natsortfiles(allfiles);
results = []
%each testoutput file has a variable called
%sumofswaps:the frequency of swaps when matching to 997
%subjects, where the first two columns are coordinates
%of the swaps and the thrid is frequency counts.
% Here I am adding all these variables together for each of the 997 subject
%to find the mean sumofswaps for each pair of 
% swap coordinates (392 regions, 392 regions)



% results takes about 3-4 minutes to generate
for ii=1:length(allfiles)
x = ii;
Chosenfiles{ii} = matfile(allfiles(x).name);
allsums = Chosenfiles{ii}.sum_swaps;



results = [results;allsums];

end


%load('networkmeansumofswap.mat', 'results')

%makescatterplotshowingswapfrequencyfornetwork
RegionSwapTable = results;
RegionSwapTable(RegionSwapTable(:, 1) == RegionSwapTable(:, 2), :) = [];
RegionSwapTable(RegionSwapTable(:, 3) == 0, :) = [];

[i,~,j] = unique(RegionSwapTable(:, [1]), 'rows');
[i, accumarray(j, RegionSwapTable(:, 3), [], @mean)];
oneRegionSwapTable = ans;
oneRegionSwapTable = ans;
oneRegionSwapTable(:,3) = (oneRegionSwapTable(:,2))/997;
yeol ={'Visual','Somatomotor','Dorsal Attn','Ventral Attn','Limbic','Fronto-parietal','Default Mode','Cerebellar','SubCortical'};
opts = detectImportOptions('CC400_Yeo7_Map.csv');
opts.SelectedVariableNames = [2];
M = readmatrix('CC400_Yeo7_Map.csv',opts)
x = (1:392)'
NetworkSwappedfrom = table(M,x);
NetworkSwappedfrom = sortrows(NetworkSwappedfrom,{'M','x'},'ascend');
cc = NetworkSwappedfrom(:,1);
cc = table2array(cc);
lix = oneRegionSwapTable;
lix = array2table(lix);
lix.Properties.VariableNames([1]) = {'keytomatch'};
NetworkSwappedfrom.Properties.VariableNames([2]) = {'keytomatch'};
NetworkSwappedfromtable = outerjoin(lix, NetworkSwappedfrom);
NetworkSwappedfromtable.Properties.VariableNames([1]) = {'regionswitchedfrom'};
NetworkSwappedfromtable = sortrows(NetworkSwappedfromtable,'regionswitchedfrom','ascend');
figure
origin = NetworkSwappedfromtable.M;
origin = string(origin);
vs = violinplot(NetworkSwappedfromtable.lix3,origin);
vs(1, 8).ViolinColor = [0,0.7235,0.6505];
vs(1, 9).ViolinColor = [0,0.85,0.905];
colormap
xlabel('Networks')
ylabel('Average Swap Frequency');
title('Network Swap Frequency')
title('Average Network Swap Frequency')
ylabel('Swap Frequency');
%makeANOVAforscatterplot
RegionsubsetforANOVAtable = RegionSwapTable(:,[1 3]);
RegionsubsetforANOVAtable = array2table(RegionsubsetforANOVAtable);
RegionsubsetforANOVAtable.Properties.VariableNames([1]) = {'keytomatch'};
%NetworkSwappedfrom.Properties.VariableNames([2]) = {'keytomatch'};
AllregionSwappedfromtable = join(RegionsubsetforANOVAtable, NetworkSwappedfrom);
AllregionSwappedfromtable.RegionsubsetforANOVAtable2 = string(AllregionSwappedfromtable.RegionsubsetforANOVAtable2);
arst = table2array(AllregionSwappedfromtable);
arst = double(arst);
[~,~,stats] = anova2(arst(:,2),arst(:,3));
[c,m,h] =  multcompare(stats,'Dimension',[1],'CType','bonferroni');

%makeheatmapforregions

mc2 = accumarray( RegionSwapTable(:,[1 2]), RegionSwapTable(:,3), [], @mean ) ;
jex = table2array(NetworkSwappedfrom(:,1));
b = NetworkSwappedfrom(:,2);
b = table2array(b);
[C, ia, ic] = unique(jex);
v = [ia;392];
vstr = movmean(v,2);
vstr = filter([1 1],2,v);
vstr = vstr(2:end);
jeen = mc2(b,b)/997;

imagesc(jeen);
hold
xticks(vstr);
yticks(vstr);
yeol ={'Visual','Somatomotor','Dorsal Attn','Ventral Attn','Limbic','Fronto-parietal','Default Mode','Cerebellar','SubCortical'};
xticklabels(yeol)
yticklabels(yeol);

hold on
for i=1:length(v)
y = v(i)
yline(y,'w','LineWidth',1)
end

for i=1:length(v)
y = v(i)
xline(y,'w','LineWidth',1)
end
colormap(turbo)
title('Average Region to Region Swap Frequency')
exportgraphics(gcf,('RegiontoRegionheatmap_latest.jpg'))
saveas(gcf,('RegiontoRegionFig_ltest'));
%%%%%%%%%%
AllnetworkSwaptable = results;
AllnetworkSwaptable(AllnetworkSwaptable(:,1) == AllnetworkSwaptable(:,2),:)= [];
AllnetworkSwaptable(AllnetworkSwaptable(:,3) == 0,:)= [];
AllnetworkSwaptable(:,4) =AllnetworkSwaptable(:,3)/997;
AllnetworkSwaptable = array2table(AllnetworkSwaptable);
AllnetworkSwaptable.Properties.VariableNames([1 2 3 4]) = {'keytomatch';'keytomatch_Region2';
'SwapNo';'SwapFreq'};
AllnetworkSwaptable = join(AllnetworkSwaptable,NetworkSwappedfrom);
AllnetworkSwaptable.Properties.VariableNames([1 2 5]) = {'Region1';'keytomatch';'Network1'};
AllnetworkSwaptable = join(AllnetworkSwaptable,NetworkSwappedfrom);
AllnetworkSwaptable.Properties.VariableNames([1 2 6]) = {'Region1';'Region2';'Network2'};
AllnetworkSwaptable.SwapFreq = string(AllnetworkSwaptable.SwapFreq);
AllnetworkSwapdouble = table2array(AllnetworkSwaptable);
AllnetworkSwapdouble = double(AllnetworkSwapdouble);
NetworktoNetwork = accumarray( AllnetworkSwapdouble(:,[5 6]), AllnetworkSwapdouble(:,4), [], @mean ) ;
figure
heatmap(NetworktoNetwork)
colormap('')
colormap;
colormap(flipud(ans));
fig = gcf
dataObjs = findobj(fig,'-property','YData')
dataObjs.XLabel = yeol
dataObjs.XLabel = {'Networks'}
dataObjs.YLabel = {'Networks'}
dataObjs.XDisplayLabels = yeol
dataObjs.YDisplayLabels = yeol

exportgraphics(gcf,('NetworktoNetworkheatmap_latest.jpg'))
saveas(gcf,('NetworktoNetworkFig_ltest'));
