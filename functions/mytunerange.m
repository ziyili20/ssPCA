%  Function Name: mytunerange.m

%  Purpose: This function generates the range of tuning parameters for
%  Fused and Grouped sPCA

%  Inputs: 
%  X is a n x p input dataset (for example, n can be the number of
%   patients, p is the number of variables);
%  method is the selection of computing methods, either Grouped or Fused;
%  Number is the number of desired principal components.  

%  Output:
%  tunerange is range of tuning parameters for selected method.

%  Author: Ziyi Li (ziyi.li@emory.edu)

%  Date: 4/27/2016

function [Taurange]= mytunerange(X,method,Number)

%obtain first nonsparse input
[~,corrX,myalpha,~]=mynonspcaall(X);
myalpha=myalpha(:,1:Number);

[n,p]=size(X);
if(strcmp(method, 'Grouped'))
    mymaxlambdaalpha(1,:)=max(abs(corrX*myalpha));
    Taurange=[sqrt(log(p)/n)*ones(1,Number)'  (mymaxlambdaalpha)']; % a matrix with noofcanvectors
elseif(strcmp(method,'Fused'))
    Taurange=[sqrt(log(p)/n)*ones(1,Number)' 100*ones(1,Number)'];
end