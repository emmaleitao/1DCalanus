function GAoutput=GArun(params,GAparams,fhandle_fitness,station,condition)
% PP=GArun(params,ngen,npop,nelite,nrep,p,func)
%
% Run the genetic algorithm
%
%  params: initial gene of parameters, 1-by-(number of parameters)
%    ...may also be npop-by-(number of parameters)
%  GAparams: a struct with the following fields
%    .ngen: number of generations
%    .npop: population size
%    .nelite: number of elite (to keep each generation)
%    .nrep: number of reproducing chromosomes (to fill after elite)
%    .pinc: the probability of mutating by increment
%    .ppow: the probability of mutating by power of 10
%    .ndigit: the number of sig. figs. to act on
%    .wt: a selection weight
%       wt = 1 selects parents randomly for reproduction
%       wt > 1 favors those with higher fitness
%    (defaults are ngen=10, npop=20, ..., pinc=0.5, ppow=0, ndigit=2, wt=1)
%    .fitness_goal: if included, the algorithm will stop once this fitness
%       level is reached
%    .naverage: if included and ==1, the last member of the population will
%       be the averaged values of the entire population
%    .nadapt: if included, the number of generations of repetition allowed
%       before the algorithm adapts, increasing pinc and ppow temporarily
%    .paramsmax
%    .paramsmin
%  F_fitness: a function that evaluates the fitness of a gene (in char form)
%
% The output PP has three dimensions:
%   npop by length(params) by ngen
%
% NRR 2008
%
% Get genetic algorithm parameters, and set defaults if necessary
%
starttime = datetime(now,'ConvertFrom','datenum');

if isfield(GAparams,'ngen')
    ngen=GAparams.ngen;
else
    ngen=10;
end
if isfield(GAparams,'npop')
    npop=GAparams.npop;
else
    npop=20;
end
if isfield(GAparams,'nelite')
    nelite=GAparams.nelite;
else
    nelite=round(npop/4);
end
if isfield(GAparams,'nrep')
    nrep=GAparams.nrep;
else
    nrep=npop-nelite;
end
if isfield(GAparams,'pinc')
    pinc=GAparams.pinc;
else
    pinc=0.5;
end
if isfield(GAparams,'ppow')
    ppow=GAparams.ppow;
else
    ppow=0.0;
end
if isfield(GAparams,'ndigit')
    ndigit=GAparams.ndigit;
else
    ndigit=2;
end
if isfield(GAparams,'wt')
    wt=GAparams.wt;
else
    wt=1;
end
if isfield(GAparams,'fitness_goal')
    fitness_goal=GAparams.fitness_goal;
end
if isfield(GAparams,'naverage')
    naverage=GAparams.naverage;
end
if isfield(GAparams,'nadapt')
    nadapt=max(1,round(GAparams.nadapt));
end
if isfield(GAparams,'paramsmax')
    paramsmax=GAparams.paramsmax;
end
if isfield(GAparams,'paramsmin')
    paramsmin=GAparams.paramsmin;
end

disp('--- Genetic Algorithm parameters ---')
disp([' ngen = ',num2str(ngen)])
disp([' npop = ',num2str(npop)])
disp([' nelite = ',num2str(nelite)])
disp([' nrep = ',num2str(nrep)])
disp([' pinc = ',num2str(pinc)])
disp([' ppow = ',num2str(ppow)])
disp([' ndigit = ',num2str(ndigit)])
disp([' wt = ',num2str(wt)])
if exist('fitness_goal')
    disp([' fitness_goal = ',num2str(fitness_goal)])
end
if exist('naverage')
    disp(['One gene is average of population'])
end
if exist('nadapt')
    disp([' nadapt = ',num2str(nadapt)])
end
if exist('paramsmax')
    disp([' Using parameter upper bounds'])
    disp(['  ',num2str(paramsmax)])
end
if exist('paramsmin')
    disp([' Using parameter lower bounds'])
    disp(['  ',num2str(paramsmin)])
end
disp('------------------------------------')

% Initialize your set of paramosomes with GAinitialize.m
% creates and array with the size npop x length(params). With defaults for
% the calanus model, that would be 20 x 10. 

if(size(params,1)==1)
    P=GAinitialize(params,npop,pinc,ppow,ndigit);
else
    P=params;
end

inew=1;
for i=1:ngen % possibly change to while loop? ngen max - condition to quit at a certain fitness? 

	% Evaluate the fitness for each gene in the population
	
    for j=inew:npop
		param_temp=P(j,:);
        f=feval(fhandle_fitness,param_temp,station); 
        if isempty(f) == 1 % with the current fitness assessment, some f values return empty. Temporary fix...
            f = 999;
        end
		F(j,1)=f;
    end
    
    inew=nelite+1; % set eval after first gen 

	% Sort the population based on fitness
	%
	F(:,1)=F(:,1).*isfinite(F(:,1)); % turn NaNs/Infs to zero % if F == 0 set to 999 - this value was infinite and should be discarded 
    
    [F,I]=sort(F); % sort the fitness (ascending) - I is the index of the old F, before sorting
    Fx(i)=F(1); % grab the best fitness ( ie: the lowest value) This array is the best fitnesses each round - 2D array? 
	P=P(I,:); % Arrange paramosomes according to the new sorted F array. First row corresponds to the highest fitness
    
	% Save this generation
	%
	PP(:,:,i)=P; % save this generation. To see final best, P(1,:,end)
    FF(:,i) = squeeze(F); % check this  F(:,1)

	% Next generation
	%
	Ptemp(1:nelite,:)=P(1:nelite,:); % top 5 parameters saved in Ptemp
	Imom=(floor(rand(1,npop-nelite).^wt*nrep)+nelite+1)'; %parents. If wt == 1, random sexual selection, if not, favors higher fitness inds
	Idad=(floor(rand(1,npop-nelite).^wt*nrep)+nelite+1)'; % don't select the elite: add nelite +1 to only select random numbers between nelite +1 : npop
    
    if exist('nadapt') & nadapt<i & Fx(i)==Fx(i-nadapt)%min(PP(:,:,i)==PP(:,:,i-nadapt))==1
        pinc=1-(1-pinc)*.999;
        %ppow=1-(1-ppow)*.99;
        disp(['NADAPT! Using pinc = ',num2str(pinc),...
            ', ppow = ',num2str(ppow)]);
    else
        pinc=GAparams.pinc;
        ppow=GAparams.ppow;
    end
    %Ptemp(nelite+1:npop,:)=...
       % mutate(crossover(P(Imom,:),P(Idad,:)),pinc,ppow,ndigit); % cross over mom and dad, mutate the new gene 
	Ptemp(nelite+1:npop,:) = P(nelite+1:npop,:);
    if exist('naverage') & naverage==1 % save last gene to be an 'avg' gene 
        Ptemp(end,:)=mean(Ptemp,1);
    end
    P=Ptemp;
	if exist('paramsmax') %cap the parameter with a maximum
        for j=1:npop
            P(j,:) = min(P(j,:),paramsmax);
        end
    end
    if exist('paramsmin') % same but min
        for j=1:npop
            P(j,:) = max(P(j,:),paramsmin); % define an array with a max and min for each parameter 
        end
    end
    
    disp(['Fitness = ',num2str(min(F))])
	if(floor(10*i/ngen)==10*i/ngen & ngen>1)
		disp(['................................',num2str(100*i/ngen),'%'])
    end
    if(exist('fitness_goal') & min(F)<=fitness_goal)
        break
    end
    
    if(length(i)>=100 && i <= 1000)
         check = all(~diff(PP(1,i:i-99,:)));
         if(sum(check) == length(params))
             break
         end 
    end 
    
    % Save important things 
    
    statname{1} = 'AG';
    statname{2} = 'HL2';
    conditionname{1} = 'Warm';
    conditionname{2} = 'Cool';
    
    GAoutput.testname = strcat(conditionname(condition), " ", statname(station));
    GAoutput.starttime = starttime;
    GAoutput.endtime = datetime(now,'ConvertFrom','datenum');
    GAoutput.topfits = Fx;
    GAoutput.allparams = PP;
    GAoutput.bestfit = Fx(end);
    GAoutput.bestparam = PP(1,:,end);
    GAoutput.allfit = FF;
    
    save('GAoutput')
end







