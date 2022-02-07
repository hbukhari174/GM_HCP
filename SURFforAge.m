clear all
infofile = matfile('twininfo_997subj.mat');
Tableofallsubjects = table(infofile.subjects, infofile.gender, infofile.age);
Tableofallsubjects.Properties.VariableNames([1 2 3]) = {'Subject ID' 'Sex' 'Age'};

allfiles = dir('*o.mat');  
allfiles = natsortfiles(allfiles);

for ii=1:length(allfiles)
   fileData{ii} = matfile(allfiles(ii).name);
   offdiagdata{ii} = fileData{ii}.offdiag_swap_counts;
  
end



Tableofallsubjects.OffDiagCounts = transpose(offdiagdata);

Tableofallsubjects.Sex = char(Tableofallsubjects{:,2});

dist = zeros(997,997);
for i=1:997
dist(i,:) = Tableofallsubjects.OffDiagCounts{i};
end

Tableofallsubjects.index = (1:height(Tableofallsubjects)).';

Z = Tableofallsubjects;
Tableofallsubjects = sortrows(Tableofallsubjects,'Age','ascend');
Tableofallsubjects = repmat(Tableofallsubjects,997,1);
Z = repmat(Z,997,1);
Tableofallsubjects = sortrows(Tableofallsubjects,1,'ascend');
Z.Properties.VariableNames([1 2 3 4 5]) = {'Subject ID 2' 'dex' 'Age 2' 'ofd' 'indx'};
matchedtables = [Tableofallsubjects Z];
for i=1:length(matchedtables.indx)
fp = matchedtables.index(i);
sp = matchedtables.indx(i);
Swapsbetweenpairs{i} = mat2cell(dist(fp,sp),1,1);
end
matchedtables.Swap = Swapsbetweenpairs';
matchedtables.Swap2 = string(matchedtables.Swap);

Matchedpairssubset = matchedtables(:,([3 8 12]));
Matchedpairssubset.Swap2 = str2double(Matchedpairssubset.Swap2);
Matchedpairssubset = table2array(Matchedpairssubset);
Matchedpairssubset(Matchedpairssubset(:,3) == 0, :) = [];

[inta,~,jinta] = unique(Matchedpairssubset(:, [1,2]), 'rows');
[inta, accumarray(jinta, Matchedpairssubset(:, 3), [], @mean)];
Uniqueagepairstable = ans;
% xval = xxsurf(1:256,1);
% xval = reshape(xval,16,16);
% yval = xxsurf(1:256,2);
% yval = reshape(yval,16,16);
% zval = xxsurf(1:256,3);
% zval = reshape(zval,16,16);
% s = surf(interp2(xval,yval,zval));
t = Uniqueagepairstable;
t(t(:, 3) ==0, :) = [];

[x, y] = meshgrid(min(t(:,1)):0.1:35, min(t(:,2)):0.1:35);
f = scatteredInterpolant(t(:,1), t(:,2), t(:,3));
z = f(x, y);
surfage = surf(x, y, z/392);
shading interp
colormap(summermap)

title('Swaps Between Individuals Increase With Age.')
xlabel('Age')
ylabel('Age')
zlabel('Swaps')

%       saveas(gcf,('surfforage.fig'))
exportgraphics(gcf,('surfforage.JPEG'))




