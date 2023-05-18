%% Header
% Trial change for github check
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["timeInS", "i_aInA", "i_bInA", "i_cInA", "u_aInV", "u_bInV", "u_cInV", "epsilon_elInRad"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
ReducedTorqueData = readtable("./ReducedTorqueData.csv", opts);


%% Clear temporary variables
clear opts
%% Given Parameters
n_me = 16.67;  % Hz aka 1000 min^-1
U_DC = 300  ;  % Volts DC
Rs = 18e-3  ;  % Ohms
Ld = 370e-6 ;  % Henry
Lq = 1200e-6;  % Henry
psi_p = 66e-3; % Vs
i_dqmax = 240; % Amps

%%
% 
% 
% $$T = \frac{3}{2}p(\psi_di_q - \psi_qi_d) $$
%
% $$\psi_d = L_di_d + \psi_p$$
%
% $$\psi_q - L_qi_q$$
%

p = 3; % Pole Pair Number
dt = 1e-6; % time step in seconds
index = 500000;
time = [0:dt:.5-dt];
for i = 1:index
    i_d(i) = 2/3*[cos(ReducedTorqueData.epsilon_elInRad(i)) -cos(ReducedTorqueData.epsilon_elInRad(i) + pi/3) -cos(ReducedTorqueData.epsilon_elInRad(i) - pi/3)]*[ReducedTorqueData.i_aInA(i); ReducedTorqueData.i_bInA(i);ReducedTorqueData.i_cInA(i)];
    i_q(i) = 2/3*[-sin(ReducedTorqueData.epsilon_elInRad(i)) sin(ReducedTorqueData.epsilon_elInRad(i) + pi/3) sin(ReducedTorqueData.epsilon_elInRad(i) - pi/3)]*[ReducedTorqueData.i_aInA(i); ReducedTorqueData.i_bInA(i);ReducedTorqueData.i_cInA(i)];
    torque(i) = 3/2*p *((Ld*i_d(i) + psi_p)*i_q(i) - (Lq*i_q(i)*i_d(i)));
end
plot(time, torque)


% m = [Timeseries120rpm.timeInS(1:index) Timeseries120rpm.i_aInA(1:index) Timeseries120rpm.i_bInA(1:index) Timeseries120rpm.i_cInA(1:index) Timeseries120rpm.u_aInV(1:index) Timeseries120rpm.u_bInV(1:index) Timeseries120rpm.u_cInV(1:index) Timeseries120rpm.epsilon_elInRad(1:index)];
% writematrix(m,'ReducedTorqueData.csv') 
% 
% for i = 1:index
%     i_d(i) = 2/3*[cos(Timeseries120rpm.epsilon_elInRad(i)) -cos(Timeseries120rpm.epsilon_elInRad(i) + pi/3) -cos(Timeseries120rpm.epsilon_elInRad(i) - pi/3)]*[Timeseries120rpm.i_aInA(i); Timeseries120rpm.i_bInA(i);Timeseries120rpm.i_cInA(i)];
%     i_q(i) = 2/3*[-sin(Timeseries120rpm.epsilon_elInRad(i)) sin(Timeseries120rpm.epsilon_elInRad(i) + pi/3) sin(Timeseries120rpm.epsilon_elInRad(i) - pi/3)]*[Timeseries120rpm.i_aInA(i); Timeseries120rpm.i_bInA(i);Timeseries120rpm.i_cInA(i)];
%     torque(i) = 3/2*p *((Ld*i_d(i) + psi_p)*i_q(i) - (Lq*i_q(i)*i_d(i)));
% end
% 

