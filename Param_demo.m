%% Param Demo
% Runs short simulations with different scaling parameters
%   * Varying static weights G
%   * Varying input Win and Winp
%   * Varying feedback Q

%% - Add paths
f = filesep;
main_dir = 'Param_demo_thalamus_scaled';
mkdir(['Output' f main_dir])

%% - Save files
savename = 'VarG';
savefolder_G = ['.' f 'Output' f main_dir f savename];

savename = 'VarWin';
savefolder_Win = ['.' f 'Output' f main_dir f savename];

savename = 'VarWinp';
savefolder_Winp = ['.' f 'Output' f main_dir f savename];

savename = 'VarQ';
savefolder_Q = ['.' f 'Output' f main_dir f savename];

%% - Global network parameters
% Set the parameters for all three simulations
N = 2000;        % number of neurons
N_th = 200;     % number of thalamus neurons
N_train = 2;    % number of training trials
N_test = 10;     % number of validation trials
N_total = 1;    % number of epochs (only functions properly wiht one)
alpha = 0.05;    % learning rate

% logicals
FORCE = false;   % apply FORCE learning during trials
makespikes = false; % make the trial spiking structures 

% Dale's law
Pexc = 0;

%% -- Varying G
%{
% varying parameters
Win = [0];       % scales the input weights
% scales the static weights
G = [0 1 10 20 30 40 50 60 70 80 90 100];
%G = [1:15];
Q = [0];         % scales the feedback weights
Winp = [0];      % network sparsity

%% --- Run varying G
run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, savefolder_G);
%}
%% -- Varying Win
% varying parameters
G = [10];        % scales the static weights
Q = [0];         % scales the feedback weights
Win = [0:0.1:1];       % scales the input weights
Winp = [1];      % network sparsity

%% --- Run varying Win
run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, savefolder_Win);

%% -- Varying Winp

% varying parameters
G = [10];        % scales the static weights
Q = [0];         % scales the feedback weights
Win = [0.5];       % scales the input weights
Winp = [0:0.1:1];      % network sparsity

%% --- Run varying Winp
run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, savefolder_Winp);

%% -- Varying Q
%{
% varying parameters
G = [10];        % scales the static weights
Q = [1:10];         % scales the feedback weights
Win = [0.05];       % scales the input weights
Winp = [1];      % network sparsity
FORCE = true;
%% --- Run varying Q

run_sim_PSTH(N, N_th, N_train, N_test, N_total, Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, savefolder_Q);
%}
%% -- Plots

% Plot varying G
%{
figure
subplot(2,2,1)
param_name = 'G';
scale_param_plot(savefolder_G, param_name, N_total, N_test)

subplot(2,2,2)
cv_G_plot(savefolder_G, N_total, N_test)

subplot(2,2,3)
cv_hist_plot(savefolder_G, N_total, N_test)
%}

% Plot varying Win
figure 
subplot(1,2,1)
param_name = 'Win';
scale_param_plot(savefolder_Win, param_name, N_total, N_test)

% Plot varying Winp
subplot(1,2,2)
param_name = 'Winp';
scale_param_plot(savefolder_Winp, param_name, N_total, N_test)

% Plot varying Q
%{
figure
cdf_Q_plot(savefolder_Q, N_total, N_test)
%}
%% Helper functions

function colormat = make_colors(Nvar, colorscheme)

cmap = colormap(colorscheme);

colorindex = floor(linspace(1, size(cmap,1), Nvar));

%colormat = zeros(Nvar, size(cmap, 2));
colormat = cmap(colorindex, :);
end

function [cs,index] = sort_nat(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Set default value for mode if necessary.
if nargin < 2
	mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
	error('sort_nat:sortDirection',...
		'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
	num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
	num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[unused,index] = sortrows(comp);
if is_descend
	index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end

%% Plot functions

function scale_param_plot(input_folder, param_name, epoch, trials)
%SCALE_PARAM_PLOT Plots the firing rate and Coefficient of variation
% over different values of a scaling parameter
%   Detailed explanation goes here

f = filesep;
%% Get a list of the structs in the folder
structs = dir(fullfile(input_folder, '*.mat'));
struct_names = sort_nat({structs.name}, 'ascend');


%% Load the structs and save the relevant information
scale_param = zeros(1, length(structs));
Cv = [];
A_t = [];
scale_param_var = [];


% loop through the files in the folder and save the relevant information
for i = 1 : length(struct_names)
    
    filenname = [input_folder f struct_names{i}];
    input = load(filenname);
    
    scale_param(i) = getfield(input.scale_param, param_name);
    
    output_Cv = 0;
    output_A_t = 0;
    
    for j = 1 : trials
        output_Cv = output_Cv + input.training_output(epoch).test_output.stats{1, j}.Cv';
        output_A_t = output_A_t + input.training_output(epoch).test_output.stats{1, j}.A_t;
    end
    
    output_Cv = output_Cv/trials;
    output_A_t = output_A_t/trials;
    
    Cv = [Cv output_Cv];
    A_t = [A_t output_A_t];
    scale_param_var = [scale_param_var (scale_param(i) * ones(1, length(output_Cv)))];
    
end

Cv(isnan(Cv)) = 0;

%% Create the plot
if strcmp(param_name, 'G')
    scatter3(Cv, A_t, scale_param_var, 15, scale_param_var, 'filled');
    xlim([0 6])
    xlabel('CV')
    ylim([0 40])
    ylabel('firing rate in Hz')
else
    scatter3(A_t, Cv, scale_param_var, 15, scale_param_var)
    xlim([0 25])
    xlabel('firing rate in Hz')
    ylim([0 3])
    ylabel('CV') 
end
colormap(cool)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = param_name;
set(colorTitleHandle ,'String',titleString);
title('CV vs firing rate')
end

function cv_G_plot(input_folder, epoch, trials)
f = filesep;
%% Get a list of the structs in the folder
structs = dir(fullfile(input_folder, '*.mat'));
struct_names = sort_nat({structs.name}, 'ascend');

%% Plot CV vs G
mean_Cv = zeros(1, length(struct_names));
Gvar = zeros(1, length(struct_names));

% loop through the files in the folder and save the relevant information
for i = 1 : length(struct_names)
    
    filenname = [input_folder f struct_names{i}];
    input = load(filenname);
    
    Gvar(i) = input.scale_param.G;
    
    output_Cv = 0;
    new_Cv = 0;
    
    for j = 1 : trials
        new_Cv = input.training_output(epoch).test_output.stats{j}.Cv';
        new_Cv(isnan(new_Cv)) = 0;
        output_Cv = output_Cv + new_Cv;
    end
    
    output_Cv = output_Cv/trials;
    
    
    mean_Cv(i) = mean(output_Cv);
end


plot(log10(Gvar), mean_Cv,'k','linewidth',1.1)
xlabel('log G')
ylabel('CV')
end

function cv_hist_plot(savefolder, epoch, trials)
f = filesep;
structs = dir(fullfile(savefolder, '*.mat'));

Cv = [];
G_var = [];
Gs = zeros(1, length(structs));


for i = 1 : length(structs)
    filenname = [savefolder f structs(i).name];
    input = load(filenname);
    
    % get the value for G
    Gs(i) = input.scale_param.G;
    
    % get the mean Cv over all the trials
    output_Cv = 0;
    
    colormat = make_colors(length(structs), cool);
    for j = 1 : trials
        output_Cv = output_Cv +...
            input.training_output(epoch).test_output.stats{1, j}.Cv';
    end
    output_Cv = output_Cv/trials;
    
    % append values to the arrays
    Cv = [Cv output_Cv];
    G_var = [G_var (Gs(i) * ones(1, length(output_Cv)))];
    
    hold on
    edges = 0:0.1:6;
    histogram(output_Cv, edges, 'displaystyle', 'stairs',...
        'EdgeColor', colormat(i, :), 'LineWidth' , 1.5);
    colormap(cool)
    colorbar
    xlim([0 6])
end
caxis([min(Gs) max(Gs)])
hold off

% make the matrix to plot
%{
histogram2(Cv, G_var,'DisplayStyle', 'tile', 'EdgeAlpha', 0, 'FaceAlpha', 1);
title('Tile histogram CV neurons')
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '# neurons';
set(colorTitleHandle ,'String',titleString);
xlabel('Cv')
ylabel('G')
%}
end

function cdf_Q_plot(savefolder, epoch, trials)

f = filesep;
structs = dir(fullfile(savefolder, '*.mat'));

struct_names = sort_nat({structs.name}, 'ascend');
colormat = make_colors(length(structs), cool);
Qs = zeros(1, length(structs));

for i = 1 : length(structs)
    
    filenname = [savefolder f struct_names{i}];
    input = load(filenname);
    
    % get the value for G
    Qs(i) = input.scale_param.Q;
    
    % get the mean Cv over all the trials
    avg_fire_rate = 0;
    
    for j = 1 : trials
        avg_fire_rate  = avg_fire_rate  +...
            input.training_output(epoch).test_output.stats{1, j}.A_t;
    end
    avg_fire_rate  = avg_fire_rate /trials;

    hold on
    nBins = 150;
    histogram(avg_fire_rate, nBins, 'displaystyle', 'stairs',...
        'EdgeColor', colormat(i, :),'Normalization', 'cdf'...
        ,'LineWidth' , 1.5);
    ylim([0 1])
    colormap(cool)
    colorbar
    
end
caxis([min(Qs) max(Qs)])
hold off
    
end






