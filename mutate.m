function gene=mutate(gene,pinc,ppow,ndigit)
% gene = mutate(gene,p)
%
% function for mutating gene
%
% gene is an array of numbers
% pinc is the probability of incrementally mutating for each allele
% ppow is the probability of incrementally mutating the power of 10 for each allele
% ndigit is the number of digits to consider
%
% NRR 2008

if(max(pinc)>1 | min(pinc)<0)
	error('pinc must be in [0 1]');
end
if(max(ppow)>1 | min(ppow)<0)
	error('ppow must be in [0 1]');
end

pinc=(rand(size(gene))<=pinc); % just 0 and 1s 
m1=floor(log10(abs(gene))); % integer of log of the abs values of gene - power of your gene 
m2=floor(rand(size(gene))*ndigit); % [0 ndigit] random integer
m3=floor(rand(size(gene))*11)-5; % [-5 5] rand int
pinc=pinc.*m3.*10.^(m1-m2); % unlogging the log having removed the digits
gene=gene+pinc; 

ppow=(rand(size(gene))<=ppow); % same as line 20
m1=10.^((rand(size(gene))<.5)*2-1); % rand 1 or -1 - trying to change big or smaller
ppow=ppow.*m1; 
I=find(ppow==0);ppow(I)=1; 
gene=abs(gene.*ppow);
