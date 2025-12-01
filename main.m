% Example based on Aiyagari (1994) with calibration from Fehr and Kindermann 
% See https://github.com/fabiankindermann/ce-fortran/tree/main/code-book/prog09/prog09_09
clear,clc,close all
% Toolkit path home pc and/or personal laptop
addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab'))
% Toolkit path university laptop 
%addpath(genpath('C:\Users\dinolaa\Documents\GitHub\VFIToolkit-matlab'))

%% Set folders
res_dir_discrete = fullfile('results','discrete'); % pure discretization
res_dir_interp   = fullfile('results','interp');   % interpolation

%% Set parameters
[n_d,d_grid,n_a,a_grid,n_z,z_grid,pi_z,Params,vfoptions,simoptions,heteroagentoptions] = set_params();

if vfoptions.gridinterplayer==0
    res_dir = res_dir_discrete;
elseif vfoptions.gridinterplayer==1
    res_dir = res_dir_interp;
end

%% Set return function
DiscFactorNames={'beta'};

ReturnFn=@(aprime,a,z,alpha,delta,r,crra) f_returnfn(aprime,a,z,alpha,delta,r,crra);

%% Create functions to be evaluated
FnsToEval.A = @(aprime,a,z) a;
FnsToEval.taxrev = @(aprime,a,z,alpha,delta,r,tau_a) f_taxrev(aprime,a,z,alpha,delta,r,tau_a);

%% General equil conditions
GEPriceNames = {'r'};
Params.r     = 0.04;
GeEqmEqns.CapitalMkt = @(r,A,alpha,delta,L_agg) r-(alpha*(A^(alpha-1))*(L_agg^(1-alpha))-delta);

%% Test value function, distribution and aggregate stats

fprintf('Grid sizes are: %d points for assets, and %d points for exogenous shock \n', n_a,n_z)

vfoptions.preGI = 1;

disp('Start VFI...')
tic
[V1,Policy1]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscFactorNames,[],vfoptions);
time_vfi=toc;

fprintf('vfoptions.gridinterplayer = %d \n',vfoptions.gridinterplayer)
fprintf('vfoptions.preGI           = %d \n',vfoptions.preGI)
fprintf('Run time VFI              = %f \n',time_vfi)

vfoptions.preGI = 0;

disp('Start VFI...')
tic
[V0,Policy0]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscFactorNames,[],vfoptions);
time_vfi=toc;

errV = max(abs(V0-V1),[],"all");
errP = max(abs(Policy0-Policy1),[],"all"); 

fprintf('vfoptions.gridinterplayer = %d \n',vfoptions.gridinterplayer)
fprintf('vfoptions.preGI           = %d \n',vfoptions.preGI)
fprintf('Run time VFI              = %f \n',time_vfi)

fprintf('Difference pre-GI post-GI, value function = %f \n',errV)
fprintf('Difference pre-GI post-GI, policy function = %f \n',errP)

disp('Start distribution...')
StatDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);

disp('Start EvalFnOnAgentDist_AggVars_Case1...')
AggVars=EvalFnOnAgentDist_AggVars_Case1(StatDist,Policy,FnsToEval,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

AggVars.A.Mean
AggVars.taxrev.Mean

%% General equilibrium

fprintf('Calculating price vector corresponding to the stationary general eqm \n')
[p_eqm,~,GeEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEval,GeEqmEqns,Params,DiscFactorNames,[],[],[],GEPriceNames,heteroagentoptions,simoptions,vfoptions);

disp('GE found, now recalculate V,Policy and Distribution at GE price')
Params.r=p_eqm.r; % Put the equilibrium interest rate into Params

[V,Policy] = ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscFactorNames,[],vfoptions);
StatDist   = StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);

%% Compute policy functions and values for other variables
% This will give you the policy in terms of values rather than index
PolicyValues = PolicyInd2Val_Case1(Policy,n_d,n_a,n_z,d_grid,a_grid,vfoptions); 
PolicyValues = reshape(PolicyValues,[n_a,n_z]);

% Additional functions to evaluate
FnsToEval.consumption = @(aprime,a,z,alpha,delta,r,tau_a) f_consumption(aprime,a,z,alpha,delta,r,tau_a);
FnsToEval.labearnings = @(aprime,a,z) z;

% Aggregate variables and values on grid
AllStats     = EvalFnOnAgentDist_AllStats_Case1(StatDist,Policy,FnsToEval,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
ValuesOnGrid = EvalFnOnAgentDist_ValuesOnGrid_Case1(Policy,FnsToEval,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

%% Display results

L_agg = Params.L_agg;
K_agg = AllStats.A.Mean;
Y_agg = K_agg^Params.alpha*L_agg^(1-Params.alpha);
C_agg = AllStats.consumption.Mean;
I_agg = Params.delta*K_agg;

coefvar.labearnings = AllStats.labearnings.StdDeviation/AllStats.labearnings.Mean;
coefvar.consumption = AllStats.consumption.StdDeviation/AllStats.consumption.Mean;
coefvar.wealth      = AllStats.A.StdDeviation/AllStats.A.Mean;

disp('TIMING')
fprintf('gridinterplayer: %d   \n',vfoptions.gridinterplayer);
fprintf('Run time for GE: %.3f \n',time_ge);

disp('MODEL RESULTS')
fprintf('GE Capital Market:            %.3f \n',GeEqmCondn.CapitalMkt);
fprintf('Capital-output ratio     K/Y: %.3f \n',K_agg/Y_agg)
fprintf('Interest rate (in %%)      r:  %.3f \n',100*p_eqm.r)
fprintf('Consumption-output ratio C/Y: %.3f \n',C_agg/Y_agg)
fprintf('Investment-output ratio  I/Y: %.3f \n',I_agg/Y_agg)

fprintf('CV labor earnings:            %.3f \n',coefvar.labearnings)
fprintf('CV of consumption:            %.3f \n',coefvar.consumption)
fprintf('CV of wealth:                 %.3f \n',coefvar.wealth)

% Open a text file for writing (this will overwrite if it already exists)
fid = fopen(fullfile(res_dir,'results.txt'),'w');

fprintf(fid, 'TIMING');
fprintf(fid,'gridinterplayer: %d   \n',vfoptions.gridinterplayer);
fprintf(fid,'Run time for GE: %.3f \n',time_ge);

fprintf(fid, 'MODEL RESULTS');
fprintf(fid, 'GE Capital Market:            %.3f \n', GeEqmCondn.CapitalMkt);
fprintf(fid, 'Capital-output ratio     K/Y: %.3f \n', K_agg/Y_agg);
fprintf(fid, 'Interest rate (in %%)      r:  %.3f \n', 100*p_eqm.r);
fprintf(fid, 'Consumption-output ratio C/Y: %.3f \n', C_agg/Y_agg);
fprintf(fid, 'Investment-output ratio  I/Y: %.3f \n', I_agg/Y_agg);

fprintf(fid, 'CV labor earnings:            %.3f \n', coefvar.labearnings);
fprintf(fid, 'CV of consumption:            %.3f \n', coefvar.consumption);
fprintf(fid, 'CV of wealth:                 %.3f \n', coefvar.wealth);

% Always close the file
fclose(fid);


%% Plots

StatDist_a = sum(StatDist,2);

figure
plot(a_grid,a_grid,'--')
hold on
plot(a_grid,PolicyValues(:,1))
hold on
plot(a_grid,PolicyValues(:,n_z))
legend('45 line','z_1','z_{nz}')
title('Policy for assets')
xlabel('Assets today')
ylabel('Assets tomorrow')

figure
plot(a_grid,ValuesOnGrid.consumption(:,1))
hold on
plot(a_grid,ValuesOnGrid.consumption(:,n_z))
legend('z_1','z_{nz}')
title('Policy for consumption')
xlabel('Assets today')
ylabel('Consumption')

figure
plot(a_grid,StatDist_a)
title('Stationary distribution')
xlabel('Assets today')

%% Plots as in Kindermann chapter
a_lim  = find(a_grid>=30, 1 );
a_lim2 = find(a_grid>=0.5, 1 );

figure
plot(a_grid(1:a_lim),StatDist_a(1:a_lim),'linewidth',2)
title('Stationary distribution, zoomed in')
xlabel('Wealth a')
ylabel('Share of households')
xlim([0,30])
print(fullfile(res_dir,'stat_dist_a'),'-dpng')

figure
plot(a_grid(1:a_lim2),PolicyValues(1:a_lim2,1),'linewidth',2)
hold on
plot(a_grid(1:a_lim2),PolicyValues(1:a_lim2,2),'linewidth',2)
hold on
plot(a_grid(1:a_lim2),PolicyValues(1:a_lim2,3),'linewidth',2)
hold on
plot(a_grid(1:a_lim2),PolicyValues(1:a_lim2,4),'linewidth',2)
hold on
title('Policy function for a'', zoomed in')
legend('z_1','z_2','z_3','z_4','Location','northwest')
xlabel('Wealth a')
ylabel('Savings for tomorrow a'' ')
xlim([0,0.5])
print(fullfile(res_dir,'pol_aprime'),'-dpng')

