function gene=crossover(mom,dad)

% gene = crossover(mom,dad)
%
% make children genes from parent genes
%
% NRR 2008

prob=((rand(size(mom,1),size(mom,2)).*2)>1); 
% list of rand numbers size of 
% * 2 , where is rand > 1. now we have a logical with 1,0 to get half T/F 
I=find(prob); % find indicies where prob == 1 
gene=mom;
gene(I)=dad(I); % now you half half of mom's genes and half dad 

