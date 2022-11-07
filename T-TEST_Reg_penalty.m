j = matfile('twininfo_997subj.mat');
%generate NR age matchec
Rlog =j.twinDZ_GT+j.twinMZ_GT+j.hasfullsib+j.hashalfsib;
Unrelated = Rlog == 0; 
[UNrx,~]=find(Unrelated==1);

load('Allsffreshreg_penalty.mat','-mat','dist')
%dist is swap distance matrix
%mean of all unrelated agematsubset


%create vector contaiing sex for UNrx
gender = j.gender
gender = strcmp(gender(UNrx),'F');

age = j.age;
age = age(UNrx);


[findx,~]=find(gender==1);
p=randperm(100);
findx100p = findx(p)

fage100p = age(findx100p);


for i=1:100
    l1=fage100p(i)
    ind = find(age==l1 & gender == 0);
    rx = randperm(length(ind));
    mind(i) = ind(rx(1));
end

femdist=sum(dist(findx100p,findx100p))/99;

masdist=sum(dist(mind,mind))/99;

Femdist_all = sum(dist(findx100p,[findx100p;mind']),2)/199;
Masdist_all = sum(dist(mind,[findx100p;mind']),2)/199;

Femdist_male = sum(dist(findx100p,mind))/99;
masdist_female = sum(dist(mind,findx100p))/99;

%use t-test2 function
