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
m     = 1                       ; % Number of reduced datapoints
dt    = 1e-6 * m                 ; % Time step in seconds
index = size(ReducedTorqueData,1)/m; %
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

angleCos = angleCos(1:m:end,:);
angleSin = angleSin(1:m:end,:);
sinEps = sin(ReducedTorqueData.epsilon_elInRad(1:m:end));
cosEps = cos(ReducedTorqueData.epsilon_elInRad(1:m:end));
iabc     = iabc(:,1:m:end);
vabc     = vabc(:,1:m:end);
Q        = Q(1:m:end);
i_d      = zeros(size(ReducedTorqueData,1)/m,1);
i_q = i_d; u_d = i_d; u_q = i_d;
%%
for i = 1:length(time)
    i_d(i) = 2/3*angleCos(i,:)*iabc(:,i);
    i_q(i) = 2/3*angleSin(i,:)*iabc(:,i);
    u_d(i) = 2/3*angleCos(i,:)*vabc(:,i);
    u_q(i) = 2/3*angleSin(i,:)*vabc(:,i);
end
%%
% i_dPass = lowpass(i_d,1/1e5);
% i_qPass = lowpass(i_q,1/1e5);
% u_dPass = lowpass(u_d,1/1e5);
% u_qPass = lowpass(u_q,1/1e5);
% sinEpsPass = lowpass(sinEps,1/1e5);
% cosEpsPass = lowpass(cosEps,1/1e5);
method = 'movmean';
i_dAve = smoothdata(i_d,method);
i_qAve = smoothdata(i_q,method);
u_dAve = smoothdata(u_d,method);
u_qAve = smoothdata(u_q,method);
sinEpsAve = smoothdata(sinEps,method);
cosEpsAve = smoothdata(cosEps,method);
torque = 3/2*p *((Ld*i_d + psi_p).*i_q - (Lq*i_q.*i_d));
% torquePass = 3/2*p *((Ld*i_dPass + psi_p).*i_qPass - (Lq*i_qPass.*i_dPass));
torqueAve = 3/2*p *((Ld*i_dAve + psi_p).*i_qAve - (Lq*i_qAve.*i_dAve));
averagedTorque = smoothdata(torque,method);

%%
figure(1)
plot(time, torque,time,torquePass,time,averagedTorque)
legend('Raw Torque Data','Low Pass Filter - \pi/20000 rad/sample','Averaged Torque Data')
xlabel('Time')
ylabel('Torque (Nm)')
title('Torque Vs. Time')
%%
% Set up library functions

%x = [torque i_d i_q u_d u_q sinEps cosEps];
%x = [torquePass i_dPass i_qPass u_dPass u_qPass sinEpsPass cosEpsPass];
x = [torqueAve i_dAve i_qAve u_dAve u_qAve sinEpsAve cosEpsAve];

M = size(x,2);
for i=3:length(x)-3
    for k = 1:M
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));
    end
end

x = [x(3:end-3,:)];
polyorder = 3;
usesine   = 0;
Theta = poolData(x,M,polyorder,usesine);
lambda = .05; % 

Xi = sparsifyDynamics(Theta,dx,lambda,6);
poolDataLIST({'torque', 'id','iq','ud','uq','sin(e)','cos(e)'},Xi,M,polyorder,usesine);
%%
r = 6;
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,r+1));
x0 = x(1,:);%
[tD,xD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),time,x0,options);

torqueSINDy = 3/2*p*((Ld*xD(:,2) + psi_p).*xD(:,3) - (Lq*xD(:,3).*xD(:,2)));
%torqueSINDy = xD(:,1);
figure(2)
plot(time,torqueAve,'r-',tD,xD(:,1),'b')
xlabel('Time')
ylabel('Torque (Nm)')
legend('Average Torque Data','SINDy Estimated Average Torque')
title('SINDy Torque Vs. Time')