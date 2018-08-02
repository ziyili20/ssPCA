%  Function Name: mynonspcaall.m

%  Purpose: This function obtain non-sparse solution to PCA problem.

%  Inputs: 
%  X is any arbitrary data matrix that you want to apply PCA on.

%  Output:
%  tildealpha is the eigenvectors of correlation matrix;
%  tildelambda is the eigenvalues of correlation matrix;
%  Atildealpha is the result of multiplying correlation matrix by
%    tildealpha devided by the (n-1), this quantity is calculated for
%    later computation convenience;
%  corrX is the correlation matrix for X, this quantity is calculated for 
%    later computation convenience as well.

%  Author: Ziyi Li (ziyi.li@emory.edu)

%  Date: 4/27/2016

function[Atildealpha,corrX, tildealpha, tildelambda]=mynonspcaall(X)

[n,p]=size(X);

X=X';
[U,D,V]=svd(X,'econ');
R=D*V';

cr = bsxfun(@minus, R,mean(R,2));

[Vr,Dr,~]=eigs((cr*cr')/(n-1), size(cr*cr',1)-1);
tildealpha=U*Vr;
tildelambda=diag(Dr);

Atildealpha=U*(cr*cr')*U'*tildealpha/(n-1);

corrX=U*(cr*cr')*U'/(n-1);
end