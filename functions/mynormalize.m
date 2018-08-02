%  Function Name: mynormalize.m

%  Purpose: This function normalize the input data so that each variable
%  has mean 0 and variance 1.

%  Inputs: 
%  mydata is any arbitrary data matrix that you want to normalize.

%  Output:
%  X is the normalized data matrix

%  Author: Ziyi Li (ziyi.li@emory.edu)

%  Date: 4/27/2016


function X=mynormalize(mydata);
[N,p]=size(mydata);
mycenter=mydata-repmat(mean(mydata,1),N,1);
mystddata=mycenter./repmat(std(mydata,0),N,1);
X=mystddata;