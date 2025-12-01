function [n_d,d_grid,n_a,a_grid,n_z,z_grid,pi_z,Params,vfoptions,simoptions,heteroagentoptions] = set_params()
% Note: the standard calibration for Aiyagari is in Kindermann book

%% Size of the grids
n_a = 1000;
n_z = 7;

%% Set toolkit options

vfoptions.verbose           = 0;    % on-screen display in VFI (default is 0)
vfoptions.tolerance         = 1e-8; % Tolerance VFI (default is 1e-9)
vfoptions.maxiter           = 1000; % Default is inf
vfoptions.howards           = 80;   % Default: 150
vfoptions.maxhowards        = 500;  % Default: 500
vfoptions.howardsgreedy     = 0;    % Default: 0 (leave it at 0!!)
vfoptions.howardssparse     = 0;    % Use sparse in Howards iterations. Default: 0
vfoptions.lowmemory         = 0;    % Default: 0, lowmemory in Return Matrix
vfoptions.separableReturnFn = 0;    % Default: 0
vfoptions.divideandconquer  = 0;    % Default: 0
vfoptions.gridinterplayer   = 1;
vfoptions.ngridinterp       = 15;
vfoptions.multigridswitch   = 10000;

simoptions.verbose          = 0;       % On-screen display in Distribution, stats
simoptions.tolerance        = 1e-8;    % Tolerance distribution
simoptions.gridinterplayer  = vfoptions.gridinterplayer;
simoptions.ngridinterp      = vfoptions.ngridinterp;
simoptions.whichstats       = zeros(7,1);
simoptions.whichstats(1)  = 1; % mean
simoptions.whichstats(3)  = 1; % standard deviation

heteroagentoptions.verbose=1;
heteroagentoptions.maxiter=70; % Max no. of iterations in GE loop (def 1000)
heteroagentoptions.toleranceGEprices=1e-12; % Accuracy of general eqm prices
heteroagentoptions.toleranceGEcondns=1e-6; % Accuracy of general eqm eqns
heteroagentoptions.fminalgo = 0; % 0=fzero; 1=fminsearch (default);
% 4=CMA-ES algorithm; 5=shooting;
% 7=fsolve; 8=lsqnonlin

%% Set economic parameters
Params.beta  = 0.96; % Discount factor
Params.crra  = 2.0;  % CRRA utility of consumption
Params.delta = 0.08; % Depreciation rate
Params.alpha = 0.36; % Capital share in Cobb-Douglas production
% --- AR(1) worker ability eps
uncmean_eps = 0;
rho_eps     = 0.6;
var_eps     = 0.04*(1-rho_eps^2);
sig_eps     = sqrt(var_eps);

%% Grids for a and z

% --- Productivity shock, z
[z_grid,pi_z]=discretizeAR1_Rouwenhorst(uncmean_eps,rho_eps,sig_eps,n_z);
z_grid = exp(z_grid);

% Normalize so that z has mean one
aux = pi_z^1000;
prob_z = aux(1,:);
prob_z = prob_z/sum(prob_z);
% mean_z = dot(z_grid,prob_z);
% z_grid = z_grid/mean_z;

% Aggregate labor is the integral of z and is exogenously set to 1
Params.L_agg = dot(z_grid,prob_z);

% --- Assets, a
a_min   = 0;
a_max   = 150;%z_grid(n_z)*100;
a_curv  = 3;
a_grid  = a_min+(a_max-a_min)*(linspace(0,1,n_a).^a_curv)';

%% Grid for d variable
n_d=0;
d_grid = [];

end %end function