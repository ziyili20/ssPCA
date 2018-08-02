%  Function Name: ssPCA_BIC.m

%  Purpose: This code use BIC criteria to obtain optimal tuning parameters 
%   for the first r principle component vector using the proposed methods

%  Inputs: 
%  X is a n x p input dataset (for example, n can be the number of
%   patients, p is the number of variables);
%  r is the number of desired principal components;
%  edgesX is a M x 2 matrix of indices of edges so that [1,2] means node 1
%   is connected to node 2, [5,3] means node 5 is connected to node 3 etc;
%  weightsX is p x 1 vector of weights for p variables;
%  mygamma is a scalar parameter for Grouped Method;
%  eta is the proportion parameter for weighting structural information and
%   l1 panelty in Fused/Grouped Methods;
%  method is the selection of computing methods, either Grouped or Fused

%  Output:
%  out is the output variable including several components:
%   (1) out.optTau is the optinal selection of tuning parameter
%   (2) out.check output a matrix, first column is the observation number
%   of tuning parameter, last column is the BIC values for each tuning
%   parameter.

%  Author: Ziyi Li (ziyi.li@emory.edu)

%  Date: 4/27/2016


function out=ssPCA_BIC(X,r,ngrid,edgesX,weightsX,mygamma,eta,method)

[n,p]=size(X);

Taurange= mytunerange(X,method,r);

Tauvec=[];
for l=1:size(Taurange,1)
 Tauvec(l,:)=10.^(linspace(log(Taurange(l,1))./log(10), log(Taurange(l,2))./log(10), ngrid+1));
end

cumuvar='F';

strcat('You choose number of grid points as', ngrid)
summaryBIC=[];
    for jj=1:(length(Tauvec)-1);
         strcat('Current jj =',num2str(jj))
         Tau=(Tauvec(:,jj));
            Z= bsxfun(@minus, X,mean(X,1));
            [rPCresults]=generateRPC(r,Z,Z,Tau,edgesX,weightsX,mygamma,eta,method,cumuvar);
            PCcombine=rPCresults.rPCload;
            %apply estimated canonical vectors on testing folds
             myUtest=X*(PCcombine*PCcombine');

             myfrobTest=log(norm(X-myUtest,'fro')/n/p/2)+log(n*p)/n/p*sum(sum(PCcombine~=0));
             %myfrobTest=log(norm(X-myUtest,'fro')/n/p/2)+log(n*p)/n/p*sum(sum(PCcombine~=0));
             %myfrobTest=log(norm(X-myUtest,'fro')/(n*p-sum(sum(PCcombine~=0))-1))+log(n*p)/n/p*sum(sum(PCcombine~=0));
             %myfrobTest=log(norm(X-myUtest,'fro')/(n*p-sum(sum(PCcombine~=0))-1))+8*log(n*p)/n/p*sum(sum(PCcombine~=0));
             myfrobTest(isnan(myfrobTest))=999;
             
         summaryBIC=[summaryBIC;[jj myfrobTest]];
    end

    %choose lambda yielding minimum mean difference
   if(isempty(summaryBIC))
        optTau=Tauvec(:,1);
    else
    row= find(summaryBIC(:,2)==nanmin(summaryBIC(:,2)),1,'last');
    optTau=Tauvec(:,row); 
    end

out.optTau=optTau;
out.method=method;
out.check=summaryBIC;
end