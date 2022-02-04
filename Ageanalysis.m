
%run in test_results folder for age based analysis
clear all
j = matfile('twininfo_997subj.mat');

%generate dist
infofile = matfile('twininfo_997subj.mat');
Tableofallsubjects = table(infofile.subjects, infofile.gender, infofile.age);
Tableofallsubjects.Properties.VariableNames([1 2 3]) = {'Subject ID' 'Sex' 'Age'};

allfiles = dir('*o.mat');  
allfiles = natsortfiles(allfiles);

for ii=1:length(allfiles)
   fileData{ii} = matfile(allfiles(ii).name);
   offdata{ii} = fileData{ii}.offdiag_swap_counts;
  
end

Tableofallsubjects.OffDiagCounts = transpose(offdata);
Tableofallsubjects.Sex = char(Tableofallsubjects{:,2});

dist = zeros(997,997);
for i=1:997
dist(i,:) = Tableofallsubjects.OffDiagCounts{i};
end

% upper tri so each swap only considered once
jist = triu(dist);
%generate matrix for male-male
M = ones(997,997)
ismale = strcmp(j.gender, "M");

M(~(ismale | ismale),:)=false;
M(:,~(ismale | ismale))=false;
M=logical(M);

%generate non-related matrix
NRmat = j.matchhalfsib+j.matchDZ+j.matchfullsib+j.matchhalfsib+j.matchMZ;
NRmat = NRmat == 0;
NRmat = NRmat.*~eye(size(NRmat));

%generate non-related m-m swaps

NRmascmat=NRmat+M;
NRmascmat=NRmascmat==2;
NRmascmat=triu(NRmascmat);
NRmascmat=logical(NRmascmat);
NRmascswap=jist(NRmascmat);

%generate non-related age-matched m-m swaps
[~,~,f] = unique(j.age);
agematch=f==f';
agematch=agematch & eye(size(agematch))==0;
agematch = triu(agematch);
NRamMasmat=agematch+NRmat+M
NRamMasmat=NRamMasmat==3;
NRamMasmat=triu(NRamMasmat);
NRamMasmat =logical(NRamMasmat);
NRagematchmascswap=jist(NRamMasmat);

%generate non-related f-f swaps

F = ones(997,997);
isfemale = ~ismale
F(~(isfemale | isfemale),:)=false;
F(:,~(isfemale | isfemale))=false;
F = triu(F);
F=logical(F)
NRfemalemat = F+NRmat
NRfemalemat=NRfemalemat==2;
NRfemalemat=triu(NRfemalemat);
NRfemalemat=logical(NRfemalemat);
NRfemswaps=jist(NRfemalemat);

%generate non-related agematched  f-f swaps
NRamFemmat=agematch+NRfemalemat
%2 here because NRfemalemats have values of 1 (line 69)

NRamFemmat=NRamFemmat==2;
NRamFemmat=triu(NRamFemmat)
NRagematchfemswaps=jist(NRamFemmat);

%generate all non-related swaps (496000 swaps)

NRswaps=jist(logical(triu(NRmat)));

%generate all non-related agematched swaps
NRagemat=NRmat+agematch;
NRagemat=NRagemat==2;
NRagemat=triu(NRagemat);
TotalNRagematchswaps=jist(NRagemat);

%generate non-related m-f swaps

NRmasfemalemat=NRmat-(NRfemalemat+NRmascmat);
NRmasfemalemat=logical(NRmasfemalemat)
NRmasfemSwaps=jist(triu(NRmasfemalemat));

%generate non-related age-matched m-f swaps
NRamMasFemmat = NRagemat-(NRamMasmat+NRamFemmat);
NRamMasFemmat=logical(NRamMasFemmat);
NRamMasFemSwaps=jist(triu(NRamMasFemmat));



FemalesNRAM = [array2table(NRagematchfemswaps)];
FemalesNRAM.Properties.VariableNames(1) = {'Swap'};
FemalesNRAM.vs(:) = {'FNRam'};

MalesNRAM = [array2table(NRagematchmascswap)];
MalesNRAM.Properties.VariableNames(1) = {'Swap'};
MalesNRAM.vs(:) = {'MNRam'};

MalesFemalesNRAM = [array2table(NRamMasFemSwaps)];
MalesFemalesNRAM.Properties.VariableNames(1) = {'Swap'};
MalesFemalesNRAM.vs(:) = {'FandMNRam'};

TotalAMswaps = [array2table(TotalNRagematchswaps)];
TotalAMswaps.Properties.VariableNames(1) = {'Swap'};
TotalAMswaps.vs(:) = {'AMNR'};

FemalesNR = [array2table(NRfemswaps)];
FemalesNR.Properties.VariableNames(1) = {'Swap'};
FemalesNR.vs(:) = {'FNR'};

MalesFemalesNR = [array2table(NRmasfemSwaps)];
MalesFemalesNR.Properties.VariableNames(1) = {'Swap'};
MalesFemalesNR.vs(:) = {'FandMNR'};

MalesNR = [array2table(NRmascswap)];
MalesNR.Properties.VariableNames(1) = {'Swap'};
MalesNR.vs(:) = {'MNR'};

TotalNR = [array2table(NRswaps)];
TotalNR.Properties.VariableNames(1) = {'Swap'};
TotalNR.vs(:) = {'NR'};

ANOvasa = [TotalNR;TotalAMswaps;FemalesNR;FemalesNRAM;MalesNR;MalesNRAM;MalesFemalesNR;MalesFemalesNRAM];
[p,t,stats] = anova1(ANOvasa.Swap,ANOvasa.vs);
%saveas(gcf,('anovadiag'));
%saveas(gcf,('anovatable'));
[c,m,h] =  multcompare(stats,'Dimension',[1],'CType','bonferroni');
%saveas(gcf,('mutlcompareage'));
