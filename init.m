function extPar = init  %defining an external parameter sent to the executable; equivalent of command line arguments to main


spotDir=fileparts(which('init.m'));

% ------------
% Compile the functions...
%
% If you see bizarre errors when mex evaluates..
% There is a file called mexopts.sh that needs to be modified.
% It is located in /home/<your user name>/.matlab/R2010a/mexopts.sh (may depend on MATLAB version)
% look in there and you will see about four places where "-ansi" is specified as a flag.
% Simply remove the "-ansi" and it should compile.
%

disp('----  Hello ------');
cd(spotDir);
mex LogLMex.cpp -L/home/kitung/GaussTest -cxx bayes.cpp -cxx nrutil.c 


disp('------ Init ------');

% ------------
% extPar is external parameters
%
%
extPar.fixed.mean=0.0;
extPar.fixed.sigma=1.0;
% 10^4 data points
extPar.fixed.bins=100000;
extPar.fixed.numdatabins=10000;
extPar.fixed.data = load('dataGausstest/data1d1e4.dat');
% 10^2
extPar.fixed.bins=10000;
extPar.fixed.numdatabins=100;
extPar.fixed.data = load('dataGausstest/data1d1e2.txt');
%extPar.fixed.data = datafile(:);

disp(extPar.fixed.data(1));



