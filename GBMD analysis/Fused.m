
% mastersheet for using Fused sPCA to analyze GBMD affy data

% set working directory
% Input your directory to ssPCA_matlab_code below:
work_directory = '/your/directory/to/ssPCA_matlab_code/';
cd(work_directory);
addpath('./functions')

% loading data set
load './GBMD analysis/GBMD_affy_processed_data.mat';

% Input your directory to cvx folder
cvx_directory = '/your/directory/to/cvx';
cd(cvx_directory);
cvx_setup;

cd(work_directory);

% call Fused sparse PCA to analyze GBMD_affy_processed_data

rng(12345); % set random seed so that the results are reproducible

Fusedresults=[];
ngrid=50; % number of searching grid points = 50
mygamma=2;
eta=0.5;
r=2;  % generate the first two PCs
cumuvar='F';    % generate the cumulative variation explained results
method='Fused';    % method is Fused sPCA
    
% Use BIC to select the best tuning parameter with training data    
[out]=ssPCA_BIC(Xtrain,r,ngrid,edgesX,weightsX,mygamma,eta,method);

% generate the first two PC loadings (saved in rPCresults.rPCload)
cumuvar='T';
[rPCresults]=generateRPC(r,Xtrain,Xtest,out.optTau,edgesX,weightsX,mygamma,eta,method,cumuvar);

% check and save the sparseness of the first tow PC loadings
Fusedresults.firstloading=sum(rPCresults.rPCload(:,1)~=0);
Fusedresults.secondloading=sum(rPCresults.rPCload(:,2)~=0);

% check and save the cumulative proportions of variation explained in
% testing data set
% Fusedresults.cumuvarprop is the cumulative variation explained by testing
% datasets.
Fusedresults.cumuvarprop=rPCresults.cumuvarprop;

% generate and save the cumulative proportions of variation explained in
% training data set
% Fusedresults.traincumuvarprop is the cumulative variation explained by
% training datasets.
PCcombine=rPCresults.rPCload;
[~,corrX,~,~]=mynonspcaall(Xtrain);
[varexplained,varprop,cumuvarprop]=deal(zeros(1,r));
        k=1;
        varexplained(k)=PCcombine(:,k)'*corrX*PCcombine(:,k);
        varprop(k)=varexplained(k)/trace(corrX);
        cumuvarprop(k)=varprop(k);
        for k=2:r
            varexplained(k)=PCcombine(:,k)'*corrX*PCcombine(:,k);
            varprop(k)=varexplained(k)/trace(corrX);
            cumuvarprop(k)=cumuvarprop(k-1)+varprop(k);
        end

Fusedresults.traincumuvarprop=cumuvarprop;

save './GBMD analysis/results/GBMD_affy_Fusedresults.mat';


%%%%%%%%%  This is extra step to present results!  %%%%%%%%%%

% plot the tuning parameter versus BIC values

plot(1:50,out.check(1:50,2))
saveas(gcf,'./GBMD analysis/results/Fused_BIC_plot.jpeg')

% save PC loadings to a separate file
csvwrite('./GBMD analysis/results/Fused_PCloadings.txt',rPCresults.rPCload);

gscatter(rPCresults.rPCload(:,1),rPCresults.rPCload(:,2))

% plot all subjects using the first two PCs grouped by their subtype
twoPC=Xuse*rPCresults.rPCload;
twoPC(:,3)=SubtypeRec(:,2);

csvwrite('./GBMD analysis/results/Fused_TwoPCs.txt',twoPC); % save the first two PCs to a separate file

