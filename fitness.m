function fitval = fitness(station,nc13,nc45,nfem,cmean,lmean)

%fitval = fitness(station,nc13,nc45,nc5,nc5d,nfem,carb,lip)
%
% Evaluates the goodness of fit of the model results to the data.
% Better fitness means lower value
%
% Fitness score currently based on :
%
% 1- Must have dormant individuals
% 2- Timing of maximum relative proportion of young copepodites & females
%    and of minimum of C5
% 3- Abundance of (a) young copepodites
%                 (b) advanced copepodites
%                 (c) females
% 4- Difference of advanced copepodites abundance between the beginning
%    and the end of the year
% 5- Levels of body C & lipids in C5
%
% F.Maps 2009


%%% Compute the timing & amplitude discrepancies in the min & max of 
%%% relative abundance of copepodite and carbon & lipid contents of C5


% Getting observations

% Abundance

if station < 2

    load obs_AG;

elseif station >= 2
    % GoM deep climatology (Wilkinson Bassin).
    % COOA station WB7 2004-2009 for environmental forcing; 
    % 2005-2007 for finmarks.
    % Updated 10/2010
    load obs_WB7;

end


% Computing the timing score

% 2- Timing of relative abundance of fem, C13 & C45

tc13 = abs( find(nc13==nanmax(nc13),1) - find(tsc13obs==nanmax(tsc13obs),1) );
tc13 = min(tc13,365-tc13) ./ 182;

tc45 = abs( find(nc45==nanmin(nc45),1) - find(tsc45obs==nanmin(tsc45obs),1) );
tc45 = min(tc45,365-tc45) ./ 182;

tfem = abs( find(nfem==nanmax(nfem),1) - find(tsfemobs==nanmax(tsfemobs),1) );
tfem = min(tfem,365-tfem) ./ 182;


% Computing the amplitude score

% 3- Abundance

ac13 = abs( nanmax(nc13)-nanmax(tsc13obs) ) ./ nanmax(tsc13obs);

ac45 = abs( nanmax(nc45)-nanmax(tsc45obs) ) ./ nanmax(tsc45obs);

afem = abs( nanmax(nfem)-nanmax(tsfemobs) ) ./ nanmax(tsfemobs);

% Compare the beginning and the end of the year simulated abundances in C4-C5 
% to check for closure of life-cycle
yc45 = abs( nc45(1)-nc45(365) ) ./ nc45(1);
  
% 4- Carbon & lipid

acarb =  abs( nanmax(cmean)-nanmax(tscobs(:,1)) ) ./ nanmax(tscobs(:,1))...
       + abs( nanmin(cmean)-nanmin(tscobs(:,1)) ) ./ nanmin(tscobs(:,1));

alip  =  abs( nanmax(lmean)-nanmax(tslobs(:,1)) ) ./ nanmax(tslobs(:,1))...
       + abs( nanmin(lmean)-nanmin(tslobs(:,1)) ) ./ nanmin(tslobs(:,1));


% "fitness" is the weighted mean of the scores so that the three blocks
% "timings", "amplitudes" & "carbon+lipids" weight the same.

fitval = 4/3*(tc13 + tc45 + tfem) + ac13 + ac45 + afem + yc45 + acarb + alip;

if ~isfinite(fitval)
    fitval = 999;
end

end