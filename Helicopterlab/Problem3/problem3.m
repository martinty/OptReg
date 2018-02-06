%% Problem 10.3

%% Initialization and model definition
init; 
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
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x = [x1 x2 x3 x4];

%% LQ
Q_LQ = diag([4,2,0,0]);
R_LQ = 0.1;
K_LQ = dlqr(A1, B1, Q_LQ, R_LQ);

%% Plotting
t = 0:h:h*(length(u)-1);
time = output.time;

figure('Name', 'Problem 10.3');
s1 = subplot(211);
y5 = output.signals.values(:,1);
plot(time,y5, 'b');
hold on;
grid on;
plot(t,x1,'r--');
xlabel('Time [s]'),ylabel('Angle [rad]');
legend('\lambda', '\lambda*');
axis([0 20 -0.5 3.5])
s2 = subplot(212);
y6 = output.signals.values(:,2);
plot(time,y6, 'b');
hold on;
grid on;
plot(t,x3,'r--');
xlabel('Time [s]'),ylabel('Angle [rad]');
legend('p', 'p*');
axis([0 20 -0.6 0.6])

title(s1, 'Travel respone');
title(s2, 'Pitch respone');

% figure('Name','Weight = 1')
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
% 
% figure('Name', 'Travel');
% y5 = output.signals.values(:,1);
% plot(time,y5, 'b');
% hold on;
% plot(t,x1,'r--');
% xlabel('Time [s]'),ylabel('Angle [rad]');
% legend('\lambda', '\lambda*');
% axis([0 25 -0.5 3.5])
% 
% figure('Name', 'Pitch');
% y6 = output.signals.values(:,2);
% plot(time,y6, 'b');
% hold on;
% plot(t,x3,'r--');
% xlabel('Time [s]'),ylabel('Angle [rad]');
% legend('p', 'p*');
% axis([0 25 -0.6 0.6])