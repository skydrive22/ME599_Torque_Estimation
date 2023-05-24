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
ReducedTorqueData = readtable("./ReducedTorqueData.csv", opts);


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

for i = 1:index
    i_d(i) = 2/3*[cos(ReducedTorqueData.epsilon_elInRad(i)) -cos(ReducedTorqueData.epsilon_elInRad(i) + pi/3) -cos(ReducedTorqueData.epsilon_elInRad(i) - pi/3)]*[ReducedTorqueData.i_aInA(i); ReducedTorqueData.i_bInA(i);ReducedTorqueData.i_cInA(i)];
    i_q(i) = 2/3*[-sin(ReducedTorqueData.epsilon_elInRad(i)) sin(ReducedTorqueData.epsilon_elInRad(i) + pi/3) sin(ReducedTorqueData.epsilon_elInRad(i) - pi/3)]*[ReducedTorqueData.i_aInA(i); ReducedTorqueData.i_bInA(i);ReducedTorqueData.i_cInA(i)];
    torque(i) = 3/2*p *((Ld*i_d(i) + psi_p)*i_q(i) - (Lq*i_q(i)*i_d(i)));

    Q = [cos(ReducedTorqueData.epsilon_elInRad(i)) sin(ReducedTorqueData.epsilon_elInRad(i))
        -sin(ReducedTorqueData.epsilon_elInRad(i)) cos(ReducedTorqueData.epsilon_elInRad(i))];
    u_d(i) = 2/3*[cos(ReducedTorqueData.epsilon_elInRad(i)) -cos(ReducedTorqueData.epsilon_elInRad(i) + pi/3) -cos(ReducedTorqueData.epsilon_elInRad(i) - pi/3)]*[ReducedTorqueData.u_aInV(i); ReducedTorqueData.u_bInV(i);ReducedTorqueData.u_cInV(i)];
    u_q(i) = 2/3*[-sin(ReducedTorqueData.epsilon_elInRad(i)) sin(ReducedTorqueData.epsilon_elInRad(i) + pi/3) sin(ReducedTorqueData.epsilon_elInRad(i) - pi/3)]*[ReducedTorqueData.u_aInV(i); ReducedTorqueData.u_bInV(i);ReducedTorqueData.u_cInV(i)];
end
%%
% Could also pull Q and u_dq outside the loop and have (:) instead of (i)
% in ...elInRad(:), then multiply u_d(:);u_q(:) with Q(:)
averagedTorque = smoothdata(torque);
plot(time, torque,time,averagedTorque)
xlabel('Time')
ylabel('Torque (Nm)')
legend('Torque Data', 'Averaged Torque Data')
%%
% Set up library functions
x = [i_d' i_q' u_d' u_q' sin(ReducedTorqueData.epsilon_elInRad(:))...
    cos(ReducedTorqueData.epsilon_elInRad(:))];
M = size(x,2);

for i=3:length(x)-3
    for k = 1:M
        dx(i-2,k) = (1/(12*dt))*(-x(i+2,k)+8*x(i+1,k)-8*x(i-1,k)+x(i-2,k));
    end
end

x = [x(3:end-3,:)];
polyorder = 1;
usesine   = 0;
Theta = poolData(x,M,polyorder,usesine);
lambda = dt; % 

Xi = sparsifyDynamics(Theta,dx,lambda,3);
poolDataLIST({'id','iq','ud','uq','sin(e)','cos(e)'},Xi,M,polyorder,usesine);

%%
r = 5;
options = odeset('RelTol',1e-8,'AbsTol',1e-8*ones(1,r+1));
x0 = x(1,:);%
[tD,xD]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),time,x0,options);

% i_dApprox = xD(:,1);
% i_qApprox = xD(:,2);
% torqueApprox = 3/2*p.*((Ld.*i_dApprox(:) + psi_p).*i_qApprox(:) -...
%     (Lq.*i_qApprox(:).*i_dApprox(:)));
figure(2)
plot(tD,xD(:,1))
%%
% Next we want to try and figure out the variables we need. We know we want
% to find $\frac{d}{dt}i_{dq} = f(i_{dq})$