
%run this portion in retest results folder to obtain Self swaps
j = matfile('twininfo_997subj.mat');
Retest = [21,29,51,67,70,96,111,133,158,168,182,192,196,204,221,223,231,265,312,322,339,351,378,388,408,414,419,442,445,459,529,552,651,721,758,759,790,867,913,914,925];
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
Retest = [21,29,51,67,70,96,111,133,158,168,182,192,196,204,221,223,231,265,312,322,339,351,378,388,408,414,419,442,445,459,529,552,651,721,758,759,790,867,913,914,925];
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
% generate dist matrix for test dataset and then run these!
j = matfile('twininfo_997subj.mat');
NRmat =j.matchDZ+j.matchfullsib+j.matchhalfsib+j.matchMZ;
NRmat = NRmat == 0;
NRmat = NRmat.*~eye(size(NRmat));
NRmat=NRmat;
[NRx,NRy] =find(NRmat==1);
[~,~,f] = unique(j.age);
agematch=f==f';
agematch=agematch & eye(size(agematch))==0;
NRagemat=NRmat+agematch;
NRagemat=NRagemat==2;
NRna= NRmat-NRagemat;
[row,col]=find(triu(NRna==1));


for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
vi = v*100/(392);
NRNAswaps = vi';
NRNAswaps = table(NRNAswaps);
NRNAswaps.Relation(1:length(row),1)= {'NRNA'};
clearvars -except NRNAswaps dist selfswaps
j = matfile('twininfo_997subj.mat');
[row,col] = find(triu(j.matchfullsib)==1)
v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(c, xx);
end
v = cell2mat(v);
vi = v*100/(392)
FSswaps = vi'
FSswaps = table(FSswaps);
FSswaps.Relation(1:200,1) = {'FS'}
clearvars -except NRNAswaps dist FSswaps selfswaps

j = matfile('twininfo_997subj.mat');

overlap = (j.matchDZ+j.matchMZ == 2);
[ox,oy] = find(overlap ==1);
DZmat = j.matchDZ;
DZmat(ox,oy) = 0;
[row,col] = find(triu(DZmat)==1)

v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
vi = v*(100/392)
DZswaps = vi'
DZswaps = table(DZswaps);
DZswaps.Relation(1:length(row),1) = {'DZ'}
clearvars -except  NRNAswaps dist FSswaps DZswaps selfswaps
j = matfile('twininfo_997subj.mat');
[row,col] = find(triu(j.matchMZ)==1);
v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
vi = v*100/392
MZswaps = vi'
MZswaps = table(MZswaps);
MZswaps.Relation(1:length(row),1) = {'MZ'};

clearvars -except  NRNAswaps dist FSswaps DZswaps selfswaps MZswaps j
NRmat =j.matchDZ+j.matchfullsib+j.matchhalfsib+j.matchMZ;
NRmat = NRmat == 0;
NRmat = NRmat.*~eye(size(NRmat));
NRmat=NRmat;
[NRx,NRy] =find(NRmat==1);
[~,~,f] = unique(j.age);
agematch=f==f';
agematch=agematch & eye(size(agematch))==0;
NRagemat=NRmat+agematch;
NRagemat=NRagemat==2;
[row,col]=find(triu(NRagemat==1));
v = [];
for i=1:length(row)
c = row(i);
xx = col(i);
v{i} = dist(xx, c);
end
v = cell2mat(v);
vi = v*(100/392)
NRAMswaps = vi'
NRAMswaps = table(NRAMswaps);
NRAMswaps.Relation(1:length(row),1) = {'NRAM'};
NRAMswaps.Properties.VariableNames(1) = {'Swap'};
NRNAswaps.Properties.VariableNames(1) = {'Swap'};
MZswaps.Properties.VariableNames(1) = {'Swap'};
DZswaps.Properties.VariableNames(1) = {'Swap'};
FSswaps.Properties.VariableNames(1) = {'Swap'};
selfswaps.Properties.VariableNames(1) = {'Swap'};


Swapsia = [selfswaps;MZswaps;DZswaps;FSswaps;NRNAswaps;NRAMswaps];
origin = Swapsia.Relation;
origin = cellstr(origin);
figure
vs = violinplot(Swapsia.Swap,origin);
figure
[p,t,stats] = anova1(Swapsia.Swap,origin)
figure
[c,m,h] =  multcompare(stats,'Dimension',[1],'CType','bonferroni');
exportgraphics(gcf,('multcomparefornoselfRel.jpg'))
% saveas(gcf,('multcompareforrelatednessNoselfRel'));
