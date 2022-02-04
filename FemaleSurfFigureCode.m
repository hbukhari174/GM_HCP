clear all
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

Tableofallsubjects.index = (1:height(Tableofallsubjects)).';

u = Tableofallsubjects(strcmp(Tableofallsubjects.Sex, "F"), :);
v = Tableofallsubjects(strcmp(Tableofallsubjects.Sex, "F"), :);
%Z = T;
u = sortrows(u,'Age','ascend');
u = repmat(u,532,1);
v = repmat(v,532,1);
u = sortrows(u,1,'ascend');
v.Properties.VariableNames([1 2 3 4 5]) = {'Subject ID 2' 'dex' 'Age 2' 'ofd' 'indx'};
Matchedpairs = [u v];
for i=1:length(Matchedpairs.indx)
fp = Matchedpairs.index(i);
sp = Matchedpairs.indx(i);
Swapsbetweenpairs{i} = mat2cell(dist(fp,sp),1,1);
end

Matchedpairs.Swap = Swapsbetweenpairs';
Matchedpairs.Swap2 = string(Matchedpairs.Swap);

Matchedpairssubset = Matchedpairs(:,([3 8 12]));
Matchedpairssubset.Swap2 = str2double(Matchedpairssubset.Swap2);
Matchedpairssubset = table2array(Matchedpairssubset);
Matchedpairssubset(Matchedpairssubset(:,3) == 0, :) = [];

[inta,~,jinta] = unique(Matchedpairssubset(:, [1,2]), 'rows');
[inta, accumarray(jinta, Matchedpairssubset(:, 3), [], @mean)];
Uniqueagepairstable = ans;

t = Uniqueagepairstable;
t(t(:, 3) ==0, :) = [];
[x, y] = meshgrid(min(t(:,1)):0.1:35, min(t(:,2)):0.1:35);
f = scatteredInterpolant(t(:,1), t(:,2), t(:,3));
z = f(x, y);
% womansagesurf = 
surf(x, y, z/392);
shading interp
colormap(turbo)

title('Pairwise Swaps of Females Over Age')
xlabel('Age')
ylabel('Age')
zlabel('Swaps')
exportgraphics(gcf,('surfforFemalesOverage_latest.JPEG'))
% saveas(gcf,('surfforFemalesOverage-ex35.fig'));
% 
% 
% 


