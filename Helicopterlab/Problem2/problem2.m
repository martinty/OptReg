%% Problem 10.2

%% Initialization and model definition
init04;
h = 0.25;                             

% Discrete time system model. x = [lambda r p p_dot]'
A1 = [1 h           0             0;
      0 1      -h*K_2             0;
      0 0           1             h;
      0 0 -h*K_1*K_pp -h*K_1*K_pd+1];
B1 = [0 0           0    h*K_1*K_pp]';

% Number of states and inputs
mx = size(A1,2);                        
mu = size(B1,2);                        

% Initial values
x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x0 = [x1_0 x2_0 x3_0 x4_0]';            % Initial values

% Time horizon and initialization
N  = 100;                               
M  = N;                                
z  = zeros(N*mx+M*mu,1);                
z0 = z;                                

% Bounds
ul 	    = -30*pi/180;        
uu 	    = 30*pi/180;             

xl      = -Inf*ones(mx,1);       
xu      = Inf*ones(mx,1);           
xl(3)   = ul;                           
xu(3)   = uu;                           

% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N,M,xl,xu,ul,uu); 
vlb(N*mx+M*mu)  = 0;                    
vub(N*mx+M*mu)  = 0;                   

% Generate the matrix Q and the vector c 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                           
Q1(2,2) = 0;                            
Q1(3,3) = 0;                            
Q1(4,4) = 0;                            
P1 = 1;                                 
Q = 2*genq2(Q1,P1,N,M,mu);              
c = zeros(N*mx+M*mu,1);                 

%% Generate system matrixes for linear model
Aeq = gena2(A1,B1,N,mx,mu);             
beq = zeros(1,N*mx)';                   
beq(1:mx) = A1*x0;                     

%% Solve QP problem with linear model
tic
[z,lambda] = quadprog(Q, [], [], [], Aeq, beq, vlb, vub);
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; 

x1 = [x0(1);z(1:mx:N*mx)];              
x2 = [x0(2);z(2:mx:N*mx)];              
x3 = [x0(3);z(3:mx:N*mx)];            
x4 = [x0(4);z(4:mx:N*mx)];              

num_variables = 5/h;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u   = [zero_padding; u; zero_padding];
x1P1 = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3P1  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

%% Plotting
t = 0:h:h*(length(u)-1);
time = outputP1.time;

pitchP01 = outputP01.signals.values(:,2);
pitchP1 = outputP1.signals.values(:,2);
pitchP10 = outputP10.signals.values(:,2);
travelP01 = outputP01.signals.values(:,1);
travelP1 = outputP1.signals.values(:,1);
travelP10 = outputP10.signals.values(:,1);

figure('Name','Problem 10.2')
s11 = subplot(321);
plot(t,x1P01, 'r--');
hold on;
grid on;
plot(time,travelP01, 'b');
legend('\lambda*', '\lambda');
ylabel('Angle [rad]'), xlabel('Time [s]');
axis([0 20 -2 3.5]);
s12 = subplot(322);
plot(t,x3P01, 'r--');
hold on;
grid on;
plot(time,pitchP01, 'b');
legend('p*', 'p');
ylabel('Angle [rad]'), xlabel('Time [s]');
axis([0 20 -0.6 0.6]);
s21 = subplot(323);
plot(t,x1P1, 'r--');
hold on;
grid on;
plot(time,travelP1, 'b');
legend('\lambda*', '\lambda');
ylabel('Angle [rad]'), xlabel('Time [s]');
axis([0 20 -2 3.5]);
s22 = subplot(324);
plot(t,x3P1, 'r--');
hold on;
grid on;
plot(time,pitchP1, 'b');
legend('p*', 'p');
ylabel('Angle [rad]'), xlabel('Time [s]');
axis([0 20 -0.6 0.6]);
s31 = subplot(325);
plot(t,x1P10, 'r--');
hold on;
grid on;
plot(time,travelP10, 'b');
legend('\lambda*', '\lambda');
ylabel('Angle [rad]'), xlabel('Time [s]');
axis([0 20 -2 3.5]);
s32 = subplot(326);
plot(t,x3P10, 'r--');
hold on;
grid on;
plot(time,pitchP10, 'b');
legend('p*', 'p');
ylabel('Angle [rad]'), xlabel('Time [s]');
axis([0 20 -0.6 0.6]);

title(s11,'Travel response, q = 0.1');
title(s12,'Pitch response, q = 0.1');
title(s21,'Travel response, q = 1');
title(s22,'Pitch response, q = 1');
title(s31,'Travel response, q = 10');
title(s32,'Pitch response, q = 10');

% figure('Name','States')
% subplot(511)
% stairs(t,u),grid
% ylabel('u')
% subplot(512)
% plot(t,x1,'m',t,x1,'mo'),grid
% ylabel('\lambda')
% subplot(513)
% plot(t,x2,'m',t,x2','mo'),grid
% ylabel('r')
% subplot(514)
% plot(t,x3,'m',t,x3,'mo'),grid
% ylabel('p')
% subplot(515)
% plot(t,x4,'m',t,x4','mo'),grid
% xlabel('tid (s)'),ylabel('pdot')