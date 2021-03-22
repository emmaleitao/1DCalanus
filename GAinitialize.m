function P=GAinitialize(params,npop,pinc,ppow,ndigit)
% P = GAinitialize(params,npop,p)
%
% Initialize a population for use in the genetic algorithm
%
%  params is a vector of initial parameter values to be perturbed.
%  NaN values will be filled with random params (not functional yet)
%
%  npop is the population size
%
% NRR 2008

% Default probability
if(~exist('pinc'))
	pinc=0.01;
end
if(~exist('ppow'))
	ppow=0.01;
end
if(~exist('ndigit'))
	ndigit=3;
end

gene=params;
P(1,:)=gene;
for i=2:npop
	P(i,:)=mutate(gene,pinc,ppow,ndigit);
end

