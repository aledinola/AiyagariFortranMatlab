function [n_d,d_grid,n_a,a_grid,n_z,z_grid,pi_z,Params,vfoptions,simoptions,heteroagentoptions] = set_params()

% Note: the standard calibration for Aiyagari is in Kindermann book


%% Size of the grids
n_a = 1000;
n_z = 7;

do_superstar = 0;

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
vfoptions.gridinterplayer   = 0;
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
Params.beta  = 0.96;%0.924;     % Discount factor
Params.crra  = 2.0; %1.5      % CRRA utility of consumption
Params.delta = 0.08;%0.06;    % Depreciation rate
Params.alpha = 0.36;%0.33;    % Capital share in non-entre sector
% --- AR(1) worker ability eps
uncmean_eps = 0;
rho_eps     = 0.6;%0.94;
var_eps     = 0.04*(1-rho_eps^2);%0.02;
sig_eps     = sqrt(var_eps);
eps_super       = 11.2; % Value (on grid) of super star shock
prob_super      = 0.0085; % Prob of becoming superstar
prob_super_back = 0.4; % Prob of falling back

% --- Taxes
Params.tau_a = 0.0;

%% Grids for a and z

% --- Productivity shock, z

if do_superstar==1
    eps_grid = zeros(n_z,1);
    pi_eps   = zeros(n_z,n_z);
    [eps_grid_help,pi_eps_help]=discretizeAR1_Rouwenhorst(uncmean_eps,rho_eps,sig_eps,n_z-1);
    
    aux = pi_eps_help^1000;
    p_unc_eps_help = aux(1,:)';
    p_unc_eps_help = p_unc_eps_help/sum(p_unc_eps_help);

    % Add superstar epsilon to eps grid
    eps_grid(1:n_z-1) = exp(eps_grid_help);
    eps_grid(n_z)     = eps_super;

    % Construct the transition prob P_eps, including transitions to superstar shock
    pi_eps(1:n_z-1,1:n_z-1) = (1-prob_super)*pi_eps_help;
    pi_eps(1:n_z-1,n_z) = prob_super;
    pi_eps(n_z,1:n_z-1) = prob_super_back*p_unc_eps_help';
    pi_eps(n_z,n_z)     = 1-prob_super_back;

    % Switch to toolkit notation
    z_grid = eps_grid;
    pi_z   = pi_eps;

elseif do_superstar==0
    [z_grid,pi_z]=discretizeAR1_Rouwenhorst(uncmean_eps,rho_eps,sig_eps,n_z);
    z_grid = exp(z_grid);
end

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