function fitval = calanusmain_GA(paramosome,station)

% Program calmainmodel modified for use with the GA procedure
% FMaps 2009

statname{1} = 'AG';
statname{2} = 'HL2';
conditionname{1} = 'Warm';
conditionname{2} = 'Cool';
outfile = char(strcat(statname(station),'_Optimization')); % save file name

tic

%%% Parameters from the paramosome tuned during the GA procedure

% Coefficient for exponential mortality function
am        = paramosome(1);

% Eponent for exponential mortality function
bm        = paramosome(2);

% Reference temperature for mortality scaling to temperature
reftemp   = paramosome(3);

% Half saturation coefficient for physiological processes
% (growth, development & egg production)
% Different for nauplii and copepodites|females
ks        = [repmat(paramosome(4),1,7) repmat(paramosome(5),1,6)];

% Coefficient to attenuate the sensisitivity to food of development
% compared to growth.
devks     = paramosome(6);

% Coefficient of lipid accumulation. Fixed for C3, variable for C4 & C5
lipup     = [zeros(1,9) 0.2 paramosome(7:8) 0];

% Proportion of lipid in C4 triggering dormancy 
sleepcrit = paramosome(9);

% Proportion of lipid in C5 triggering exit from dormancy
wakecrit  = paramosome(10);

%%%

% Load fixed parameters values

load cfin_parameters


% Load physical setup
% Physical and numerical parameters

deltat = 30; %minutes

time = ndays*24*60;

%%% These climatological series include "tzchl" which has
%%% chl as f(depth) as well, matched to edges

if station == 1

    load mld_temp_forcing_AG

    load food_forcing_AG

    totdepth = 300;

    diapdepth = 150;

    diapdepx = 200;

    lat = 49.72;

    lon = -66.25;

    tzone = 4;

    % Monthly sunrise/sunset in min
    %sunrise = [ 7*60+15 ...
    %            6*60+30 ...
    %            5*60+45 ...
    %            4*60+30 ...
    %            3*60+45 ...
    %            3*60+15 ...
    %            3*60+30 ...
    %            4*60+15 ...
    %            5*60+00 ...
    %            5*60+45 ...
    %            6*60+30 ...
    %            7*60+15 ];

    %sunset  = [ 15*60+45 ...
    %            16*60+45 ...
    %            17*60+30 ...
    %            18*60+15 ...
    %            19*60+00 ...
    %            19*60+30 ...
    %            19*60+30 ...
    %            18*60+45 ...
    %            17*60+45 ...
    %            16*60+30 ...
    %            15*60+45 ...
    %            15*60+30 ];

elseif station == 2
    % GoM deep climatology (Wilkinson Bassin).
    % COOA station WB7 2004-2009
    % Updated 11/09

    load mld_temp_forcing_WB7

    load food_forcing_WB7

    totdepth = 250;

    diapdepth = 150;

    diapdepx = 200;

    lat = 42.86;

    lon = -69.86;

    tzone = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% TEST WARM SURFACE SCENARIO
%%% temp*1.2 above the mixed layer.

    %tztemp(:,1:3) = tztemp(:,1:3).*1.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% TEST COLD SCENARIO (Loder et al. 2001)

    %tztemp(:,1:6) = tztemp(:,1:6).*0.9;

    %tztemp(:,7:end) = tztemp(:,7:end).*0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % Monthly sunrise/sunset in min
    %sunrise = [ 7*60+15 ...
    %            6*60+45 ...
    %            6*60+00 ...
    %            5*60+00 ...
    %            4*60+15 ...
    %            4*60+00 ...
    %            4*60+15 ...
    %            4*60+45 ...
    %            5*60+15 ...
    %            5*60+45 ...
    %            6*60+30 ...
    %            7*60+00 ];

    %sunset  = [ 16*60+30 ...
    %            17*60+15 ...
    %            17*60+45 ...
    %            18*60+15 ...
    %            19*60+00 ...
    %            19*60+15 ...
    %            19*60+15 ...
    %            18*60+45 ...
    %            17*60+45 ...
    %            17*60+00 ...
    %            16*60+15 ...
    %            16*60+00 ];

end

wcedges = zedges(totdepth,tmld);

% Adjusting migrating depth to topography
if totdepth < 200
    migdepth = mean([min(tmld) totdepth]);
else
    migdepth = 100;
end

% Ensures the pseudo-random numbers generator always produce the same 
% sequence between subsequent runs

rng default


%**************************************************************************

poddeltat = deltat/(24*60); % timestep for converting rates to proper units

%%% Description of the main state variable.
%%% Re-use some during development to save space.
%
% "pod" is a matrix (metrics,individuals) holding all the pod related items
% pod(1,:)  =  Stage (MFC) | for females = C when molting
% pod(2,:)  =  Age
% pod(3,:)  =  State : 0|1; for C5 0 = diapause, 1 = active; for "females" 0 = male, and >1 = female (= their C/N ratio)
% pod(4,:)  =  pod C, structure + lipid stores
% pod(5,:)  =  pod C which is lipid stores
% pod(6,:)  =  pod N, structure only (no N in lipid) | for females = MFC for 1 generation time -> maximum life span.
% pod(7,:)  =  Target depth (m)
% pod(8,:)  =  vertical position (m)
% pod(9,:)  =  Food concentration at the individual location 
% pod(10,:) =  Temperature at the individual location
% pod(11,:) =  Demographic weight of the individual when resampling population
% pod(12,:) =  vector ID
%
%%%

%%% Initializing a winter population (with diapause)
%

pod = zeros(12,initnum);

if station < 2

    %%% Initial conditions for AG

    c5d = round(initnum.*0.95);

    pod = zeros(12,initnum);

    % In case we use reduced development in C5d
    pod(1,1:c5d) = 11.5;

    pod(1,c5d+1:end) = 12.5;

    pod(2,c5d+1:end) = 30;

    pod(3,c5d+1:end) = 1;

    pod(6,:) = 30;

    pod(4,1:c5d) = pod(6,1:c5d).*9;

    pod(4,c5d+1:end) = pod(6,c5d+1:end).*c2n(13);

    pod(5,:) = pod(4,:).*0.5;

    pod(8,:) = diapdepth + (rand(1,length(c5d))-0.5).*25;

    pod(9,:) = interp1(wcedges(1,:),tzchl(1,:),pod(8,:));

    pod(10,:) = interp1(wcedges(1,:),tztemp(1,:),pod(8,:));

    pod(11,:) = 30;

    pod(12,:) = 1:initnum;

    idcount = initnum;

elseif station >= 2

    %%% Initial conditions for WB7

    c5d = round(initnum.*0.95);

    pod = zeros(12,initnum);

    % In case we use reduced development in C5d
    pod(1,1:c5d) = 11.8;

    pod(1,c5d+1:end) = 12.5;

    pod(2,c5d+1:end) = 50;

    pod(3,c5d+1:end) = 1;

    pod(6,:) = 20;

    pod(4,1:c5d) = pod(6,1:c5d).*8;

    pod(4,c5d+1:end) = pod(6,c5d+1:end).*c2n(13);

    pod(5,:) = pod(4,:).*0.5;

    pod(8,:) = diapdepth + (rand(1,length(c5d))-0.5).*25;

    pod(9,:) = interp1(wcedges(1,:),tzchl(1,:),pod(8,:));

    pod(10,:) = interp1(wcedges(1,:),tztemp(1,:),pod(8,:));

    pod(11,:) = 30;

    pod(12,:) = 1:initnum;

    idcount = initnum;

end
%%%

% Interval for saving the output data (in ndays)
saveinterval = 1;

% Matrices for saving details about the copepods
podsave = nan(12,maxnum,round(ndays/saveinterval)+1);


% Mortality function & mortality related values

% First compute rates decreasing exponentially for each of 13 stages
mortini = 0.01 + am.*exp(-bm.*[1:13]);

% Critical development lengths for each stage ~ die of old age
dtcrit = ones(1,13);

dtcrit(1) = aD(1)*(mintemp+cD)^bD;

for i = 2:12
    dtcrit(i) = dtcrit(i-1) + aD(i)*(mintemp+cD)^bD;
end

% Critical carbon mass for each stage ~ die if too skinny 
% < body C of preceeding stage at 12C (Campbell et al. 2001)
wgtcrit(1:9) = 0;

wgtcrit(10) = 9.39-0.276.*12; % first C3

wgtcrit(11) = 28.9-1.09.*12;  % C4

wgtcrit(12) = 94-4.71.*12;    % C5

wgtcrit(13) = wgtcrit(12);    % Females 

% Temperature scaling of mortality
tmort = (reftemp+cD)^bD;


% Setting up variables for later use
dayselapsed = 0; % counter for doing the interval to save pod data

saves = 1; % counter for incrementing the index for saving pod data


% Migration related variables
night = 1;

migok = 0;


%**************************************************************************
% This is the start of the first main loop, 
% units are in minutes for the running of the model
%

for i = 0:deltat:time

    nowtime = i/(24*60);

    day = mod( nowtime, 365 );

    iday = floor(day+1);

    hr = mod( i./60, 24 );

    % Sunrise/sunset calculation (NOAA "simple equations")

    if hr == 0

        % Four eq below in radian
        gamma = 2.*pi./365 .* ( iday -1 + (hr-12)./24 );

        eqtime = 229.18 .* (  0.000075 + 0.001868.*cos(gamma) ...
                            - 0.032077.*sin(gamma) ...
                            - 0.014615.*cos(2.*gamma) ...
                            - 0.040849.*sin(2.*gamma) );

        decl =  0.006918 - 0.399912.*cos(gamma) + 0.070257.*sin(gamma) ...
              - 0.006758.*cos(2.*gamma) + 0.000907.*sin(2.*gamma) ...
              - 0.002697.*cos(3.*gamma) + 0.00148.*sin(3.*gamma);

        % !!! ATTENTION: mix of degree & radian !!!
        ha = acos( cosd(90.833)./(cosd(lat).*cos(decl)) - tand(lat).*tan(decl) );

        % Convert from radian to degree
        ha = ha.*180./pi;

        sunrise = mod( round( 720+4.*(lon-ha)-eqtime )+tzone*60, 1440 );

        sunset  = mod( round( 720+4.*(lon+ha)-eqtime )+tzone*60, 1440 );

    end

    if hr>sunrise/60-1 && hr<sunset/60-1 && night==1
        migok = 1;
        night = 0;
    elseif hr>sunset/60-1 && night==0
        migok = -1;
        night = 1;
    end

    % Simple 12/12 photoperiod
    %if hr>6 && hr<18 && night==1
    %    migok = 1;
    %    night = 0;
    %elseif hr>18 && night==0
    %    migok = -1;
    %    night = 1;
    %end

    % F.Maps 2009. Use vertically structured food field

    zchl = zeros(1,13);

    if structfood == 0
        zchl(2) = max(tzchl(iday,:));
        zchl(1) = 0.7.*zchl(2);
        zchl(3) = 0.9.*zchl(2);
    else
        zchl(1:4) = tzchl(iday,1:4);
    end

%**************************************************************************
% Copepod section


    %%% Retrieve the stage distribution at this point

    stg = min(13,ceil(pod(1,:)));

    mig = stg>10 & pod(3,:)>0;

    egg = stg==1;

    c5d = stg==12 & pod(3,:)==0;

    fem = stg==13;


    %%% Everyone gets older

    pod(2,:) = pod(2,:) + poddeltat; % units of poddeltat need to be in ndays

    % Females do not age this way...
    pod(2,fem) = pod(2,fem) - poddeltat;

    % Max life span of females is 2 generations at the temperature they develop to
    pod(2,fem) = pod(2,fem) + poddeltat./( 2.*sum(aD(1:end-1)).*(pod(10,fem)+cD).^bD ); 


    %%% Update food and temperature at the cop location before computing
    %%% physiological rates

    pod(9,~mig) = interp1(wcedges(iday,:),zchl,pod(8,~mig));

    % Update food only at night when DVM -> remember [food] during day at depth
    if night==1
        pod(9,mig) = interp1(wcedges(iday,:),zchl,pod(8,mig));
    end

    food = max( 0., pod(9,:)-0.3 );

    pod(10,:) = interp1(wcedges(iday,:),tztemp(iday,:),pod(8,:));


    %%% Development

    % Belehradek relationship with temperature
    % * Campbell based development reducer at low food
    % 0.3 mg Chla m-3 x C/Chla=50 --> ~ 15 mg C m-3

    dev =   1./( aD(stg).*(pod(10,:)+cD).^bD )...
         .* food.^2./(food.^2+(devks.*ks(stg)).^2);

    pod(1,:) = pod(1,:) + dev.*poddeltat;


    %%% Growth
    % Limitation of growth by food is more pronounced than for development

    % C growth
    Cmax =   max( 0., aC(stg).*pod(10,:).^2+bC(stg).*pod(10,:)+cC(stg) )...
          .* food.^2./(food.^2+ks(stg).^2)...
          .* ceil(pod(3,:)); % C5d do not grow

    pod(4,:) = pod(4,:) .* (1.+Cmax.*poddeltat);

    % Lipids
    pod(5,:) = pod(5,:) + lipup(stg).*pod(4,:).*Cmax.*poddeltat.*ceil(pod(3,:)); % C5d use lipids;

    % Function for lipid metabolism in C5d �gC.�gN-1.day-1
    meta = max( 0., 24.*0.74.*12.011.*10.^(0.0442.*pod(10,c5d)+2.117).*1e-6 );

    pod(4,c5d) = pod(4,c5d) - meta.*pod(6,c5d).*poddeltat;

    pod(5,c5d) = pod(5,c5d) - meta.*pod(6,c5d).*poddeltat;

    % N growth
    Nmax =   max( 0., aN(stg).*pod(10,:).^2+bN(stg).*pod(10,:)+cN(stg) )...
          .* food.^2./(food.^2+ks(stg).^2)...
          .* ceil(pod(3,:)); % C5d do not grow

    pod(6,:) = pod(6,:) .* (1.+Nmax.*poddeltat); 


    %%%%%%%%%%%%%%%%%%%    LIPID & DORMANCY CONTROL    %%%%%%%%%%%%%%%%%%%%

    % Send C4 on the dormant path if they have stored enough lipids

    sleep = pod(5,:)>sleepcrit.*pod(4,:) & stg==11;

    pod(3,sleep) = 0.5; 

    % Enter dormancy when the C5 development is completed

    sleep = pod(1,:)>11.99 & pod(3,:)==0.5;

    pod(3,sleep) = 0;

    % Wake up when reaching lower lipid threshold or minimum C/N

    wake =  ( stg==12 & pod(3,:)==0 & pod(5,:)<wakecrit.*pod(4,:) )...
          | ( stg==12 & pod(3,:)==0 & pod(4,:)./pod(6,:)<4 );

    pod(1,wake) = 12;

    pod(3,wake) = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Egg production

    if sum(fem)>0

        % No egg made from lipid storage

        % Hirche 1997
        femmax = 0.029 + 0.00727.*pod(10,fem);

        % Have to deal with variable C/N ratios and the possibility for females
        % to grow heavier, especially after exit from dormancy.
        % I consider the distance between the females' actual c2n and a
        % maximum c2n = 6 (Campbell et al. 2001).
        % If their c2n is below, part of the EPR is converted in growth for
        % the female. This portion decreases asymptotically as c2n approaches 6.

        % Compute the N content of the matter ingested and transfered into eggs
        cn = c2n(1)+c2n(13)-pod(4,fem)./pod(6,fem);

        if sum(cn<1.05*c2n(1))>0
            cn(cn<1.05*c2n(1)) = c2n(1);
        end

        Cmax = pod(4,fem).*femmax.*food(fem).^2./(food(fem).^2+ks(stg(fem)).^2).*poddeltat;

        pod(4,fem) = pod(4,fem) + Cmax;

        Nmax = Cmax./cn;

        pod(6,fem) = pod(6,fem) + Nmax;

        % Store the N accumulated since last clutch into stage holder
        pod(1,fem) = pod(1,fem) + Nmax;

        % Subtracts the weight of the newly made eggs right before going down
        negg = pod(1,fem).*0;

        if migok == 1

            % Check if they have enough N to make at least 1 egg
            negg = floor((pod(1,fem)-12.5)./eggwgt.*c2n(1));

            pod(1,fem) = pod(1,fem).*(1-real(negg==0))+12.5.*real(negg==0);

            pod(4,fem) = pod(4,fem) - negg.*eggwgt;

            pod(6,fem) = pod(6,fem) - negg.*eggwgt./c2n(1);

            % Inherit the depth, food, temp and demographic weight of the female
            inherit = pod(8:11,pod(1,:)>12.5);

            pod(1,pod(1,:)>12.5) = 12.5;

        end

    end


    %%% Swimming

    %Leising's function
    %step = bl(stg).*(1-pod(9,:)./(ks(stg)+pod(9,:)));

    % The current time step is too long to take full advantage of the
    % statisticall nature of his function.
    % I used mine which is more determinitic and better adapted to our
    % constraints (in abscence of more detailed knowledge of swimming
    % behavior...)
    % Daily change in optimal depth.

    if migok == -1

        % Trade-off between growth and development
        % Target depth = depth which would produce the largest individual

        stgi = unique(stg);

        for j = 1:length(stgi)

            if stgi(j)<5
                % For N1-3 optimize DT (no growth)
                best =   ( zchl.^2./(zchl.^2+(devks.*ks(stgi(j))).^2) )...
                      ./ ( aD(stgi(j)).*(tztemp(iday,:)+cD).^bD );

                pod(7,stg==stgi(j)) = wcedges(iday, min(4,find(best==max(best),1)) );

            elseif stgi(j)<13

                % Others optimize growth
                best =   aD(stgi(j)).*(tztemp(iday,:)+cD).^bD...
                      .* max( 0., aC(stgi(j)).*tztemp(iday,:).^2+bC(stgi(j)).*tztemp(iday,:)+cC(stgi(j)) )...
                      .* zchl.^2./(zchl.^2+ks(stgi(j)).^2);

                pod(7,stg==stgi(j)) = wcedges(iday, find(best==max(best),1) );

            else

               % Females optimize EPR
                best =   (0.029+0.00727.*tztemp(iday,:))...
                      .* zchl.^2./( zchl.^2+ks(stgi(j)).^2 );

                pod(7,stg==stgi(j)) = wcedges(iday, find(best==max(best),1) );

            end

        end

    end

    step = deltat.*60.*migspeed.*bl(stg);

    pod(7,c5d) = diapdepth;

    % General case: normal distribution around target depth
    pod(8,:) = max( 1, pod(7,:)+randn(1,size(pod,2)).*step );

    % Case of falling eggs (~10 m.d-1)
    if sum(egg)>0
        pod(8,egg) = min( totdepth, pod(8,egg)+abs(randn(1,sum(egg)).*step(egg)) );
    end


    %%% Mortality (make room now for eggs)

    % Base mortality + temperature effect + age effect (except for c5d) + size effect
    mort = poddeltat.*(   mortini(stg).*tmort./(pod(10,:)+cD).^bD...
                       +  floor(pod(2,:)./dtcrit(stg)).*ceil(pod(3,:))...
                       +  floor(wgtcrit(stg)./pod(4,:)) );

    mort(c5d) = 0.1.*mort(c5d);

    kill = zeros(1,size(pod,2));

    % For meta-individuals with only 1 individual
    if sum(pod(11,:)==1)>0
        kill(pod(11,:)==1) = min(1, floor( rand(1,sum(pod(11,:)==1))+mort(pod(11,:)==1) ) );
    end

    % For meta-individuals with multiple individuals
    kill(pod(11,:)>1) = pod(11,pod(11,:)>1).*mort(pod(11,:)>1);

    pod(11,:) = pod(11,:)-kill;

    pod(:,pod(11,:)<1) = [];

    % If all the meta-individuals disapear, return.
    if sum(size(pod))==0 && negg==0

        fitval = 999;

        save(outfile,'fitval','paramosome','station','totdepth','ndays',...
             'wcedges','podsave');

        return

    end


    %%% Take care of the sex-ratio = get rid of the males in this step
    % (Procedure similar to mortality)

    mature = pod(1,:)>=12 & pod(1,:)<12.5;

    % Set their "stage" as 12.5 at maturity
    pod(1,mature) = 12.5;

    % Reset their age to 0 at maturity
    pod(2,mature) = 0;

    male = zeros(1,size(pod,2));

    male(mature) = sexratio;

    % For meta-individuals with only 1 individual
    if sum(pod(11,:)==1)>0
        male(pod(11,:)==1) = min(1, floor( rand(1,sum(pod(11,:)==1))+male(pod(11,:)==1)-1e-6 ) );
    end

    % For meta-individuals with multiple individuals
    male(pod(11,:)>1) = pod(11,pod(11,:)>1).*male(pod(11,:)>1);

    pod(11,:) = pod(11,:)-male;

    pod(:,pod(11,:)<1) = [];

    % If all the meta-individuals disapear, return.
    if sum(size(pod))==0

        fitval = 999;

        save(outfile,'fitval','paramosome','station','totdepth','ndays',...
             'wcedges','podsave');

        return

    end


    %%%%%%%%%%%%%%%%%%%%%%%   Newly produced eggs   %%%%%%%%%%%%%%%%%%%%%%%

    if migok==1 && ~isempty(negg)

        % Total number of meta-individuals    
        alive = size(pod,2);

        % Total number of eggs
        ne = sum(negg);

        if ne > 0

            % Subsampling if > maxnum and adding eggs

            % Define the lineage of the eggs
            lineage = zeros(4,ne);

            dummy = negg(negg>0);

            for j = 1:length(dummy)

                jj = sum(dummy(1:j-1));

                lineage(:,jj+1:jj+dummy(j)) = repmat(inherit(:,j),1,dummy(j));

            end

            % Check if subsampling is needed
            if alive+ne <= maxnum

                % Add newly produced eggs
                % Stage
                pod(1,alive+1:alive+ne) = .0001;

                % Age
                pod(2,alive+1:alive+ne) = poddeltat;

                % State
                pod(3,alive+1:alive+ne) = 1.;

                % Carbon
                pod(4,alive+1:alive+ne) = eggwgt;

                % Lipids
                pod(5,alive+1:alive+ne) = 0.;

                % Nitrogen
                pod(6,alive+1:alive+ne) = eggwgt./c2n(1);

                % Depth, food, temp & demographic weight
                % inherited from the female
                pod(8:11,alive+1:alive+ne) = lineage;

                % Vector ID
                pod(12,alive+1:alive+ne) = idcount+[1:ne];
                idcount = idcount+ne;

            else

                % Record pre-sampling stage composition in order to
                % ensure that no stage would disappear because of subsampling
                stg_in = min(13,ceil(pod(1,:)));

                % Check if some stage should not be subsampled
                sample_stg = 1:alive;

                num_in = zeros(1,13);

                for j = 1:13

                    if (sum(stg_in==j)<maxnum*1e-2)
                        sample_stg(stg_in==j) = 0;
                    end

                    num_in(j) = sum(pod(11,stg_in==j));

                end

                sample_stg(sample_stg==0) = [];

                % Define subsample size according to maxnum
                sample = floor( length(sample_stg) * (maxnum-alive+length(sample_stg))/(length(sample_stg)+ne) );

                % Randomly select individuals
                permut = randperm(length(sample_stg));

                % Keep selected individuals
                pod(:,sample_stg(permut(sample+1:end))) = [];

                % Update number of meta-individuals
                alive = size(pod,2);

                % Randomly add newly produced eggs
                permut = randperm(ne);

                pod(1,alive+1:maxnum) = .0001;

                pod(2,alive+1:maxnum) = poddeltat;

                pod(3,alive+1:maxnum) = 1.;

                pod(4,alive+1:maxnum) = eggwgt;

                pod(5,alive+1:maxnum) = 0.;

                pod(6,alive+1:maxnum) = eggwgt./c2n(1);

                % Depth, food, temp & demographic weight
                % inherited from the female
                pod(8:11,alive+1:maxnum) = lineage(:,permut(1:maxnum-alive));

                % Vector ID
                pod(12,alive+1:maxnum) = idcount+[1:maxnum-alive];
                idcount = idcount+maxnum-alive;

                % Modify the demographic weight of the remaining individuals
                stg_out = min(13,ceil(pod(1,:)));

                pod(11,stg_out==1) = pod(11,stg_out==1).*(num_in(1)+sum(lineage(4,:)))/sum(pod(11,stg_out==1));

                for j = 2:13
                    pod(11,stg_out==j) = pod(11,stg_out==j).*num_in(j)/sum(pod(11,stg_out==j));
                end

            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alive = size(pod,2);

    % Saving pods just before they migrate down
    if migok == 1;

        migok = 0;

        if dayselapsed == saveinterval

            dayselapsed = 0;

            saves = saves+1;

            podsave(:,1:alive,saves) = pod;

        end

        % Now update the depth of cops performing DVM
        stg = min(13,ceil(pod(1,:)));

        mig = stg>10 & pod(3,:)>0;

        pod(7,mig) = migdepth;

        pod(8,mig) = pod(7,mig);

    elseif migok == -1

        migok = 0;

    end

    if hr == 0
        dayselapsed = dayselapsed+1;        
    end

% End copepod section
%**************************************************************************

if floor(10*i/time)==(10*i/time)
    disp(['.................................... ' num2str(100*i/time) ' %'])
end

end %time loop


% Evaluate fitness

podfit = podsave(:,:,end-364:end);

stg = min(13,ceil(podfit(1,:,:)));

c13 = stg>=8 & stg<11;

c45 = stg>=11 & stg<13;

c5  = stg>=12 & stg<13;

c5d = stg>=12 & stg<13 & podfit(3,:,:)==0;

fem = stg==13;


% Abundances of relevent stages

nc5d = zeros(size(squeeze(podfit(1,:,:))));

nc5d(c5d) = podfit(11,c5d);

nc5d = sum(nc5d);

% First fitness condition: are there dormant copepodites?

if sum(nc5d)==0

    fitval = 999;

    save(outfile,'fitval','paramosome','station','totdepth','ndays',...
         'wcedges','podsave');

    return

end

nc13 = zeros(size(squeeze(podfit(1,:,:))));

nc13(c13) = podfit(11,c13);

nc13 = sum(nc13);

nc45 = zeros(size(squeeze(podfit(1,:,:))));

nc45(c45) = podfit(11,c45);

nc45 = sum(nc45);

nfem = zeros(size(squeeze(podfit(1,:,:))));

nfem(fem) = podfit(11,fem);

nfem = sum(nfem);


% Carbon content

nc5 = zeros(size(squeeze(podfit(1,:,:))));

nc5(c5) = podfit(11,c5);

carb = zeros(size(squeeze(podfit(1,:,:))));

carb(c5) = podfit(4,c5);

cmean = sum(carb.*nc5)./sum(nc5);


% Lipid content

lip = zeros(size(squeeze(podfit(1,:,:))));

lip(c5) = podfit(5,c5)./podfit(4,c5);

lmean = sum(lip.*nc5)./sum(nc5);


% Call fitness function

fitval = fitness(station,nc13,nc45,nfem,cmean,lmean);


save(outfile,'fitval','paramosome','station','totdepth','ndays',...
     'wcedges','podsave');

toc

end
