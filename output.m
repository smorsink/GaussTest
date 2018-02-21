disp('Output results to file!')

%name='/home/kitung/Spot/Data/Data-2018-01-18-B'

%!mv '/home/kitung/Spot/FerretData' $name

%cd '/home/kitung/Spot/Data/Data-2017-11-24-X'
cd '/home/kitung/GaussTest/FerretData'

%cd name

load('OptimalSolutions.mat');
%load('MergedSolutions.mat');


para = transpose(OptimalSolutions.X);
dlmwrite('optimals.txt',para,'delimiter', '\t');

fitness = transpose(OptimalSolutions.F);
dlmwrite('fitness.txt',fitness,'delimiter', '\t');

rank = transpose(OptimalSolutions.rank);
dlmwrite('rank.txt',rank,'precision',6);

ndim = size(para)

dlmwrite('dimensions.txt',ndim,'delimiter','\t','precision',6);


