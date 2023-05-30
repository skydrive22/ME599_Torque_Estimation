%% Torque Estimation
% ME 599 - Data-Driven Modeling
clc;clear; 

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
ReducedTorqueData = readtable("../ReducedTorqueData.csv", opts);


%% Clear temporary variables
clear opts
%%
%% Given Parameters
n_me = 16.67;  % Hz aka 1000 min^-1
U_DC = 294  ;  % Volts DC
Rs = 18e-3  ;  % Ohms
Ld = 370e-6 ;  % Henry
Lq = 1200e-6;  % Henry
psi_p = 66e-3; % Vs
i_dqmax = 240; % Amps

%%
% To start we want to first import our data and analyze the torque graph
% with time to obatin a clear idea of the quality of data. This also allows
% us to compare the data with the averages given the torque dataset
% provided by the publishers to verify the values are the same.
%
% $$T = \frac{3}{2}p(\psi_di_q - \psi_qi_d) $$
%
% where T = torque, p = pole pair number, $\psi_d$ and $\psi_q$ = , and
% $i_d$ and $i_q$ = current in after a current transform from the $i_{abc}$
% current readings. Next, we can calculate both $\psi_d$ and $\psi_q$ with
% the following equations:
%
% $$\psi_d = L_di_d + \psi_p$$
%
% $$\psi_q = L_qi_q$$
%
% We also want to find $u_{dq}$ with the following equation:
%
% $$u_{dq} = Q(\epsilon_{el})\frac{U_{DC}}{3}\left[\begin{array}{ccc}
% 1 & -\frac{1}{2} & -\frac{1}{2}\\
% 0 & \frac{\sqrt{3}}{2} & - \frac{\sqrt{3}}{2}
% \end{array}\right]v_n$ Where $$Q(\epsilon_{el}) =
% \left[\begin{array}{ccc}
% cos(\epsilon_{el}) & sin(\epsilon_{el})\\
% -sin(\epsilon_{el}) & cos(\epsilon_{el})
% \end{array}\right]$

p     = 3                        ; % Pole Pair Number
dt    = 1e-6                     ; % Time step in seconds
index = size(ReducedTorqueData,1); %
time  = [0:dt:.5-dt]             ;
vn    = [1;1;1]                  ; % vn = switch state of sa, sb, sc, but
                                   % the research group only used switch
                                   % state 1
%%
angleCos = [cos(ReducedTorqueData.epsilon_elInRad),...
        -cos(ReducedTorqueData.epsilon_elInRad + pi/3),...
        -cos(ReducedTorqueData.epsilon_elInRad - pi/3)];
angleSin = [-sin(ReducedTorqueData.epsilon_elInRad),...
        sin(ReducedTorqueData.epsilon_elInRad + pi/3),...
        sin(ReducedTorqueData.epsilon_elInRad - pi/3)];
iabc = [ReducedTorqueData.i_aInA ReducedTorqueData.i_bInA ReducedTorqueData.i_cInA]';
vabc = [ReducedTorqueData.u_aInV ReducedTorqueData.u_bInV ReducedTorqueData.u_cInV]';
Q = [cos(ReducedTorqueData.epsilon_elInRad) sin(ReducedTorqueData.epsilon_elInRad)
        -sin(ReducedTorqueData.epsilon_elInRad) cos(ReducedTorqueData.epsilon_elInRad)];

i_d = zeros(size(ReducedTorqueData,1),1);
i_q = i_d; u_d = i_d; u_q = i_d;
%%
for i = 1:length(time)
    i_d(i) = 2/3*angleCos(i,:)*iabc(:,i);
    i_q(i) = 2/3*angleSin(i,:)*iabc(:,i);
    u_d(i) = 2/3*angleCos(i,:)*vabc(:,i);
    u_q(i) = 2/3*angleSin(i,:)*vabc(:,i);
end
%%
i_dPass = lowpass(i_d,.00005);
i_qPass = lowpass(i_q,.00005);
i_dAve = smoothdata(i_d);
i_qAve = smoothdata(i_q);
torque = 3/2*p *((Ld*i_d + psi_p).*i_q - (Lq*i_q.*i_d));
torquePass = 3/2*p *((Ld*i_dPass + psi_p).*i_qPass - (Lq*i_qPass.*i_dPass));
torqueAve = 3/2*p *((Ld*i_dAve + psi_p).*i_qAve - (Lq*i_qAve.*i_dAve));
averagedTorque = smoothdata(torque);

%%
subplot(2,1,1)
plot(time, torque,time,torquePass)
legend('Raw Torque Data','Low Pass Filter - \pi/20000 rad/sample')
subplot(2,1,2)
plot(time,averagedTorque,time,torqueAve)
xlabel('Time')
ylabel('Torque (Nm)')
legend('Averaged Torque Data','Averaged Current')
title('Torque Vs. Time')
%%
% Set up library functions
sinEps = sin(ReducedTorqueData.epsilon_elInRad);
cosEps = cos(ReducedTorqueData.epsilon_elInRad);
x = [i_d i_q u_d u_q sinEps cosEps];
%x = [i_dPass i_qPass u_d u_q sinEps cosEps];

M = size(x,2);
for i=3:length(x)-3
    for k = 1:M
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));
    end
end
%%
x = [x(3:end-3,:)];
polyorder = 1;
usesine   = 0;
Theta = poolData(x,M,polyorder,usesine);
lambda = .5; % 

Xi = sparsifyDynamics(Theta,dx,lambda,6);
A = poolDataLIST({'id','iq','ud','uq','sin(e)','cos(e)'},Xi,M,polyorder,usesine);
%%
r = 5;
options = odeset('RelTol',1e-3,'AbsTol',1e-6*ones(1,r+1));
x0 = x(1,:);%
[tD,xD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),time,x0,options);

torqueSINDy = 3/2*p*((Ld*xD(:,1) + psi_p).*xD(:,2) - (Lq*xD(:,2).*xD(:,1)));
figure(2)
n = 500;
plot(tD,torqueSINDy,time,averagedTorque,'r',time(1:n:end),torquePass(1:n:end))
xlabel('Time')
ylabel('Torque (Nm)')
legend('SINDy Estimated Torque','Averaged Torque Data','Low Pass Torque Data')
title('Torque Vs. Time')