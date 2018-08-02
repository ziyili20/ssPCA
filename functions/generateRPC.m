%  Function Name: generateRPC.m

%  Purpose: This function generates the first r sparse PC using the proposed methods

%  Inputs: 
%  r is the number of desired principal components;
%  X is a n x p input (training) dataset (for example, n can be the number of
%   patients, p is the number of variables);
%  Xtest is a m x p input testing dataset;  If you want to use the same
%   testing dataset as training dataset, just specify the same input matrix
%   for Xtest as X.  
%  Tau is the tuning parameter used in the program.  You can either specify
%  one or choose it using "cvssPCA_BIC" function;
%  edgesX is a M x 2 matrix of indices of edges so that [1,2] means node 1
%   is connected to node 2, [5,3] means node 5 is connected to node 3 etc;
%  weightsX is p x 1 vector of weights for p variables;
%  mygamma is a scalar parameter for Grouped Method;
%  eta is the proportion parameter for weighting structural information and
%   l1 panelty in Fused/Grouped Methods;
%  method is the selection of computing methods, either Grouped or Fused;
%  cumuvar is the option to generate cumulative variation explained or not.
%   Please specify 'T' or 'F'.

%  Output:
%  rPCresults is the output variable including several components:
%   (1) rPCresults.rPCload is the r estimated PC loadings by selected
%   methods;
%   (2) rPCresults.varprop is the vector of proportions of variation 
%   explained by each PC loading;
%   (3) rPCresults.cumuvarprop is the vector of cumulative proportions of
%   variation explained by the each PC loading. (rPCresults.cumuvarprop indeed 
%   include the same information as rPCresults.varprop.  It might be a bit redundancy, 
%   but depends on what variable you need, you can choose to use varprop or
%   cumuvarprop).

%  Author: Ziyi Li (ziyi.li@emory.edu)

%  Date: 4/27/2016


function rPCresults=generateRPC(r,X,Xtest,Tau,edgesX,weightsX,mygamma,eta,method,cumuvar)
    
 
   p=size(X,2);
    [Atildealpha,~,~,mylambdaoldall]=mynonspcaall(X);
    [~,corrX,~,~]=mynonspcaall(Xtest);
    PCsummary=[];
    [~,idx]=sort(edgesX(:,1));
    edges2=edgesX(idx,:);
    
    for i =1:r
        cvx_begin quiet
         variable alphai(p)
        
         alphai_w = alphai./ weightsX.^ (1/mygamma) ;
         sumedges = 0;
         
         if(strcmp(method,'Grouped'))
                     for k = 1 : size(edgesX,1),
                         sumedges = sumedges + norm(alphai_w(edgesX(k,:)),mygamma);
                     end
         elseif(strcmp(method,'Fused'))
                    %this avoids a for loop
                    sumedges=sum(abs(alphai(edges2(:,1))./weightsX(edges2(:,1))-alphai(edges2(:,2))./weightsX(edges2(:,2))));
         end
         
         minimize((1-eta)*sumedges + eta*norm(alphai(setdiff(1:p,unique(edgesX))),1))
         
         subject to
        
         norm( Atildealpha(:,i)-mylambdaoldall(i)*alphai,Inf)<=Tau(i);   

          for j=2:i;
               norm(PCsummary(:,j-1)'*alphai,1)<=eps;
          end;
        cvx_end
        
        alphai(abs(alphai)<=0.001)=0;
        
        if(sum(abs(alphai))~=0)
            myalphai=alphai./norm(alphai,2);
        else
             myalphai=alphai;
        end
        
       PCsummary = [PCsummary myalphai];
    end

    rPCresults.rPCload = PCsummary;

    if (strcmp(cumuvar,'T'))
        [varexplained,varprop,cumuvarprop]=deal(zeros(1,r));
        k=1;
        varexplained(k)=PCsummary(:,k)'*corrX*PCsummary(:,k);
        varprop(k)=varexplained(k)/trace(corrX);
        cumuvarprop(k)=varprop(k);
        for k=2:r
            varexplained(k)=PCsummary(:,k)'*corrX*PCsummary(:,k);
            varprop(k)=varexplained(k)/trace(corrX);
            cumuvarprop(k)=cumuvarprop(k-1)+varprop(k);
        end
        rPCresults.varprop=varprop;
        rPCresults.cumuvarprop=cumuvarprop;
    end

end
