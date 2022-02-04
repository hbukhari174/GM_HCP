% crystal = T.CogCrystalComp_AgeAdj.values
% fluid = T.CogFluidComp_AgeAdj.values
% total = T.CogTotalComp_AgeAdj.values

clear all
j = matfile('twininfo_997subj.mat');
Tableofallsubjects = table(j.subjects, j.gender, j.age);
Tableofallsubjects.Properties.VariableNames([1 2 3]) = {'Subject ID' 'Sex' 'Age'};

allfiles = dir('*o.mat');  
allfiles = natsortfiles(allfiles);

for ii=1:length(allfiles)
   fileData{ii} = matfile(allfiles(ii).name);
   offdata{ii} = fileData{ii}.offdiag_swap_counts;
  
end




Tableofallsubjects.OffDiagCounts = transpose(offdata);

Tableofallsubjects.Sex = char(Tableofallsubjects{:,2});
Tableofallsubjects.index = (1:height(Tableofallsubjects)).';

dist = zeros(997,997);
for i=1:997
dist(i,:) = Tableofallsubjects.OffDiagCounts{i};
end
HCPbehaviordata = readtable('hcp_behaviour.csv');
Tableofallsubjects = Tableofallsubjects(strcmp(Tableofallsubjects.Sex, "F"), :)

Tableofallsubjects.Properties.VariableNames(1) = {'Subject'};
CFCdataindex = find(string(HCPbehaviordata.Properties.VariableNames) == "CogFluidComp_AgeAdj")
CFCdatatable = HCPbehaviordata(:,[1 CFCdataindex]);
Commons = innerjoin(CFCdatatable, Tableofallsubjects, 'Keys','Subject');

Commons.CogFluidComp_AgeAdj = ceil(Commons.CogFluidComp_AgeAdj/5)*5;

Z = Commons;
Zcheck = Commons;
Commons( any(ismissing(Commons),2), :) = [];

Z( any(ismissing(Z),2), :) = [];

Commons = sortrows(Commons,'CogFluidComp_AgeAdj','ascend');
Commons.OffDiagCounts = [];
Z.OffDiagCounts = [];
Commons = repmat(Commons,532,1);
Z = repmat(Z,532,1);
Z.Properties.VariableNames([1 2 3 4 5]) = {'Subject ID 2' 'ctc' 'dex' 'age' 'indx'};
Commons = sortrows(Commons,1,'ascend');

Matchedpairs = [Commons Z];
for i=1:length(Matchedpairs.indx)
fp = Matchedpairs.index(i);
sp = Matchedpairs.indx(i);

Swapsbetweenmmpairs{i} = mat2cell((triu(dist(fp,sp))),1,1);
end

Matchedpairs.Swap = Swapsbetweenmmpairs';
Matchedpairs.Swap2 = string(Matchedpairs.Swap);

Matchedpairssubset = Matchedpairs(:,([2 7 12]));
Matchedpairssubset.Swap2 = str2double(Matchedpairssubset.Swap2);
Matchedpairssubset = table2array(Matchedpairssubset);
Matchedpairssubset(Matchedpairssubset(:, 3) == 0, :) = [];

% 
[inta,~,jinta] = unique(Matchedpairssubset(:, [1,2]), 'rows');
[inta, accumarray(jinta, Matchedpairssubset(:, 3), [], @mean)];
Uniquescorepairstable = ans;

Uniquescorepairstable(:,4) = accumarray(jinta,1);
Uniquescorepairstable(Uniquescorepairstable(:, 3) == 0, :) = [];
Uniquescorepairstable(Uniquescorepairstable(:, 1) < 70, :) = [];
Uniquescorepairstable(Uniquescorepairstable(:, 2) < 70, :) = [];
Uniquescorepairstable(Uniquescorepairstable(:, 1) > 145, :) = [];
Uniquescorepairstable(Uniquescorepairstable(:, 2) > 145, :) = [];
t = Uniquescorepairstable;
figure
p =t(:,3)*100/392
s = heatmap(reshape(p,[],16));
{'65-70';'70-75';'75-80';'80-85';'85-90';'90-95';'95-100';'100-105';'105-110';'110-115';'115-120';'120-125';'125-130';'130-135';'135-140';'140-145'}
s.XDisplayLabels = ans
s.YDisplayLabels = ans
xlabel('CFC score')
ylabel('CFC score')
title("Average Swaps by Pairs of CFC score -Male- binned by  5")
colormap(parula)
% % % xlim([1.5 16.5])
% % % ylim([1.5 16.5])
% % {'70-80','80-90','90-100','100-110','110-120','120-130','130-140','140-150'}
% % % yticklabels({'70-75','80-85','90-95','100-105','110-115','120-125','130-135','140-145'})
figure
sf = heatmap(reshape(t(:,4),[],16));
{'65-70';'70-75';'75-80';'80-85';'85-90';'90-95';'95-100';'100-105';'105-110';'110-115';'115-120';'120-125';'125-130';'130-135';'135-140';'140-145'}
sf.XDisplayLabels = ans
sf.YDisplayLabels = ans
xlabel('CFC score')
ylabel('CFC score')
title("Number of Pairs per Bin - Male")
colormap(parula)