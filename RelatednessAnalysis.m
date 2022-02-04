
%run this portion in retest results folder to obtain Self swaps
j = matfile('twininfo_997subj.mat');
Retest = 1:41
T = table(j.subjects(Retest,1), j.gender(Retest,1), j.age(Retest,1));
T.Properties.VariableNames([1 2 3]) = {'Subject ID' 'Sex' 'Age'};
allfiles = dir('*o.mat');
allfiles = natsortfiles(allfiles);
for ii=1:length(allfiles)
fileData{ii} = matfile(allfiles(ii).name);
offdata{ii} = fileData{ii}.offdiag_swap_counts(1:997,1);
end
T.OffDiagCounts = transpose(offdata);
T.Sex = char(T{:,2});
dist = zeros(997,997);

for i=1:41
dist(i,:) = T.OffDiagCounts{i};
end


row = 1:41;
col = [21,29,51,67,70,96,111,133,158,168,182,192,196,204,221,223,231,265,312,322,339,351,378,388,408,414,419,442,445,459,529,552,651,721,758,759,790,867,913,914,925];
v = []
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(c,xx);
end

v = cell2mat(v)
v = v*100/(392)
violinplot(v)
close all
violinplot(v)

selfswappercent = table(v');
selfswappercent.Relation(1:41,1) = {'Self'};
clearvars -except selfswappercent

% run this code for swap matrix for retest_results

j = matfile('twininfo_997subj.mat');
Retest = 1:41
T = table(j.subjects(Retest,1), j.gender(Retest,1), j.age(Retest,1));
T.Properties.VariableNames([1 2 3]) = {'Subject ID' 'Sex' 'Age'};
allfiles = dir('*o.mat');
allfiles = natsortfiles(allfiles);
for ii=1:length(allfiles)
fileData{ii} = matfile(allfiles(ii).name);
offdata{ii} = fileData{ii}.offdiag_swap_counts(1:997,1);
end
T.OffDiagCounts = transpose(offdata);
T.Sex = char(T{:,2});
dist = zeros(997,997);
for i=1:41
dist(i,:) = T.OffDiagCounts{i};
end
row = 1:41;
col = [21,29,51,67,70,96,111,133,158,168,182,192,196,204,221,223,231,265,312,322,339,351,378,388,408,414,419,442,445,459,529,552,651,721,758,759,790,867,913,914,925];
selfswapmapercentmatt=dist(row,col)*100/392;
imagesc(selfswapmapercentmatt)
colormap(hot);
colormap
colormap(flipud(ans);
title("Distance Matrix - Swap Numbers between Test-Retest Pairs")
xlabel('Subject')
ylabel('Subject')
exportgraphics(gcf,('DistMatrixfortest-retest.jpg'))
saveas(gcf,('DISTMATRIXfortest-retestFIG'));

clearvars -except selfswappercent selfswapmapercentmatt

% go run in test_folder, make sure not to clear variables

j = matfile('twininfo_997subj.mat');
NRmat = j.matchhalfsib+j.matchDZ+j.matchfullsib+j.matchhalfsib+j.matchMZ;
NRmat = NRmat == 0;
NRmat = NRmat.*~eye(size(NRmat));
[row,col] = find(triu(NRmat)==1);

T = table(j.subjects, j.gender, j.age);
T.Properties.VariableNames([1 2 3]) = {'Subject ID' 'Sex' 'Age'};
allfiles = dir('*o.mat');
allfiles = natsortfiles(allfiles);
for ii=1:length(allfiles)
fileData{ii} = matfile(allfiles(ii).name);
offdata{ii} = fileData{ii}.offdiag_swap_counts;
end
T.OffDiagCounts = transpose(offdata);
T.Sex = char(T{:,2});
T.index = (1:height(T)).';
dist = zeros(997,997);
for i=1:997
dist(i,:) = T.OffDiagCounts{i};
end
v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
vi = v*100/(392);
NRswaps = vi';
NRswaps = table(NRswaps);
NRswaps.Relation(1:496100,1)= {'NR'};
clearvars -except selfswappercent selfswapmapercentmatt  NRswaps dist
j = matfile('twininfo_997subj.mat');
[row,col] = find(triu(j.matchfullsib)==1)
v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
v = v*100/(392)
FSswaps = v'
FSswaps = table(FSswaps);
FSswaps.Relation(1:200,1) = {'FS'}
clearvars -except selfswappercent selfswapmapercentmatt  NRswaps dist FSswaps

j = matfile('twininfo_997subj.mat');
[row,col] = find(triu(j.matchDZ)==1)

v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
v = v*100/392
DZswaps = v'
DZswaps = table(DZswaps);
DZswaps.Relation(1:95,1) = {'DZ'}
clearvars -except selfswappercent selfswapmapercentmatt  NRswaps dist FSswaps DZswaps
j = matfile('twininfo_997subj.mat');
[row,col] = find(triu(j.matchMZ)==1);
v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
v = v*100/(392)
MZswaps = v'
MZswaps = table(MZswaps);
MZswaps.Relation(1:130,1) = {'MZ'};

[~,~,unage] = unique(j.age);
agematch=unage==unage';
agematch=agematch & eye(size(agematch))==0;
agematch = triu(agematch);
NRmat = j.matchhalfsib+j.matchDZ+j.matchfullsib+j.matchhalfsib+j.matchMZ;
NRmat = NRmat == 0;
NRmat = NRmat.*~eye(size(NRmat));

NRagemat=NRmat+agematch;
NRagemat=NRagemat==2;
NRagemat=triu(NRagemat);
NRagematswaps=dist(NRagemat)*100/392;
NRAMswaps = [array2table(NRagematswaps)];
NRAMswaps.Relation(:,1) = {'NRAM'};

NRswaps.Properties.VariableNames(1) = {'Swap'};
NRAMswaps.Properties.VariableNames(1) = {'Swap'};
MZswaps.Properties.VariableNames(1) = {'Swap'};
DZswaps.Properties.VariableNames(1) = {'Swap'};
FSswaps.Properties.VariableNames(1) = {'Swap'};
selfswappercent.Properties.VariableNames(1) = {'Swap'};



Swapsia = [NRswaps;NRAMswaps;MZswaps;DZswaps;FSswaps;selfswappercent];
origin = Swapsia.Relation;
origin = cellstr(origin);
vs = violinplot(Swapsia.Swap,origin);

[p,t,stats] = anova1(Swapsia.Swap,origin)
[c,m,h] =  multcompare(stats,'Dimension',[1],'CType','bonferroni');
exportgraphics(gcf,('multcompare.jpg'))
saveas(gcf,('multcompareforrelatedness'));


