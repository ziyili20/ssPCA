
% Example code to use Fused/Grouped sPCA 

% set working directory
work_directory = '/your/working/directory/to/ssPCA_matlab_code/';
cd(work_directory);
addpath('./functions')

% loading data set
load './exampledata.mat';

% add path to your cvx folder
cvx_directory = '/your/directory/to/your/cvx/folder';
cd(cvx_directory);
cvx_setup;

cd(work_directory);

% call Fused sparse PCA to analyze GBMD_affy_processed_data

ngrid=20; % number of searching grid points = 20
mygamma=2;  % when Grouped is used, mygamma is set to 8 in our implementations
eta=0.5;
r=2;  % generate the first two PCs
cumuvar='F';    % generate the cumulative variation explained results
method='Fused';    % method is Fused sPCA, you can also use 'Grouped'
    
% BIC to select the best tuning parameter with training data
% this may take a while, especially when you use Grouped sPCA
% decrease ngrid can speed up the program
[BICout]=ssPCA_BIC(Xtrain,r,ngrid,edgesX,weightsX,mygamma,eta,method);


% instead of using BIC, you can also specify tuning parameter
% yourself.  Just use
% cvout.optTau = [*,*]  * can be any tuning parameter you want.

% generate the first two PC loadings (saved in rPCresults.rPCload)
cumuvar='T';
[rPCresults]=generateRPC(r,Xtrain,Xtest,BICout.optTau,edgesX,weightsX,mygamma,eta,method,cumuvar);

%  if you want PC loadings, rPCresults.rPCload is the variable you are
%  looking for.

save './exmpleoutput.mat';  % output results
