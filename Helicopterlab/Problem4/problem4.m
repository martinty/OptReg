%% Problem 10.4

%% Initialization and model definition
init;
h = 0.25;

global N mx beta alpha lambda_t

% Discrete time system model. x = [lambda r p p_dot e e_dot]'
A1 = [1 h           0            0           0            0;
      0 1      -h*K_2            0           0            0;
      0 0           1            h           0            0;
      0 0 -h*K_1*K_pp 1-h*K_1*K_pd           0            0;
      0 0           0            0           1            h;
      0 0           0            0 -h*K_3*K_ep 1-h*K_3*K_ed];
  
B1 = [0 0 0 h*K_1*K_pp 0          0;
      0 0 0          0 0 h*K_3*K_ep]';

% Number of states and inputs
mx = size(A1,2); 
mu = size(B1,2); 

% Initial values
x1_0 = pi;                              % lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';  % Initial values

% Time horizon and initialization
N  = 40;                                
M  = N;                                 
z  = zeros(N*mx+M*mu,1);                
z0 = z;                                 

% Bounds
ul 	    = -30*pi/180*[1; inf];          
uu 	    = 30*pi/180*[1; inf];          

xl      = -Inf*ones(mx,1);              
xu      = Inf*ones(mx,1);               
xl(3)   = ul(1);                        
xu(3)   = uu(1);                        
beta = 20;
alpha = 0.2;
lambda_t = 2*pi/3;


% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N,M,xl,xu,ul,uu); 
vlb(N*mx+M*mu)  = 0;                    
vub(N*mx+M*mu)  = 0;                    

% Generate the matrix Q and the vector c 
r1 = 1;
r2 = 0.5;
Q = zeros(mx,mx);
Q(1,1) = 1;                            
Q(2,2) = 0;                            
Q(3,3) = 0;                           
Q(4,4) = 0;                            
Q(5,5) = 0;                           
Q(6,6) = 0;                           
R = diag([r1 r2]);                     
G = 2*genq2(Q,R,N,M,mu);              
c = zeros(N*mx+M*mu,1);                

%% Generate system matrixes for nonlinear model
Aeq = gena2(A1,B1,N,mx,mu);            
beq = zeros(1,N*mx)';                 
beq(1:mx) = A1*x0;                    

%% Objective
f = @(Z) 0.5*Z'*G*Z;

%% Solve QP problem with nonlinear model
tic
opt = optimoptions(@fmincon, 'MaxFunEvals', 100000);
[Z, fval] = fmincon(f, z0, [], [], Aeq, beq, vlb, vub, @constraint, opt);
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+G(i,i)*Z(i)*Z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u_1 = [Z(N*mx+1:2:N*mx+M*mu); Z(N*mx+M*mu-1)];
u_2 = [Z(N*mx+2:2:N*mx+M*mu); Z(N*mx+M*mu)];
  
  
travel = [x0(1);Z(1:mx:N*mx)];         
travel_rate = [x0(2);Z(2:mx:N*mx)];    
pitch = [x0(3);Z(3:mx:N*mx)];         
pitch_rate = [x0(4);Z(4:mx:N*mx)];     
elevation = [x0(5);Z(5:mx:N*mx)];      
elevation_rate = [x0(6);Z(6:mx:N*mx)]; 

num_variables = 5/h;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u_1   = [zero_padding; u_1; zero_padding];
u_2   = [zero_padding; u_2; zero_padding];
u = [u_1 u_2];

travel  = [pi*unit_padding; travel; zero_padding];
travel_rate  = [zero_padding; travel_rate; zero_padding];
pitch  = [zero_padding; pitch; zero_padding];
pitch_rate  = [zero_padding; pitch_rate; zero_padding];
elevation  = [zero_padding; elevation; zero_padding];
elevation_rate  = [zero_padding; elevation_rate; zero_padding];
x = [travel travel_rate pitch pitch_rate elevation elevation_rate];

%% LQ
Q_LQ = diag([4,2,0,0,6,2]);
R_LQ = diag([1,0.4]);
K_LQ = dlqr(A1, B1, Q_LQ, R_LQ);

%% Plotting
t = 0:h:h*(length(u_1)-1);
time = output.time;

figure('Name','Problem 10.4.4')
s1 = subplot(311);
x1 = output.signals.values(:,1);
plot(time,x1, 'b');
hold on;
grid on;
plot(t,travel,'r--');
xlabel('Time [s]'),ylabel('Angle [rad]');
legend('\lambda', '\lambda*');
s2 = subplot(312);
x3 = output.signals.values(:,2);
plot(time,x3, 'b');
hold on;
grid on;
plot(t,pitch,'r--');
xlabel('Time [s]'),ylabel('Angle [rad]');
legend('p', 'p*');
s3 = subplot(313);
x5 = output.signals.values(:,3);
plot(time,x5, 'b');
hold on;
grid on;
plot(t,elevation,'r--');
xlabel('Time [s]'),ylabel('Angle [rad]');
legend('e', 'e*');

title(s1, 'Travel respone');
title(s2, 'Pitch respone');
title(s3, 'Elevation respone');

% figure('Name','Weight = 1')
% subplot(711)
% stairs(t,u_1),grid
% ylabel('u')
% subplot(712)
% plot(t,travel,'m',t,travel,'mo'),grid
% ylabel('\lambda')
% subplot(713)
% plot(t,travel_rate,'m',t,travel_rate','mo'),grid
% ylabel('r')
% subplot(714)
% plot(t,pitch,'m',t,pitch,'mo'),grid
% ylabel('p')
% subplot(715)
% plot(t,pitch_rate,'m',t,pitch_rate','mo'),grid
% ylabel('$\dot{\textrm{p}}$','interpreter','latex', 'fontsize', 15)
% subplot(716)
% plot(t,elevation,'m',t,elevation','mo'),grid
% ylabel('e')
% subplot(717)
% plot(t,elevation_rate,'m',t,elevation_rate','mo'),grid
% xlabel('tid (s)'),ylabel('$\dot{\textrm{e}}$','interpreter','latex', 'fontsize', 15)
% 
%
% figure('Name', 'Travel');
% x1 = output.signals.values(:,1);
% plot(time,x1, 'b');
% hold on;
% plot(t,travel,'r--');
% xlabel('Time [s]'),ylabel('Angle [rad]');
% legend('\lambda', '\lambda*');
% 
% figure('Name', 'Pitch');
% x3 = output.signals.values(:,2);
% plot(time,x3, 'b');
% hold on;
% plot(t,pitch,'r--');
% xlabel('Time [s]'),ylabel('Angle [rad]');
% legend('p', 'p*');
% 
% figure('Name', 'Elevation');
% x5 = output.signals.values(:,3);
% plot(time,x5, 'b');
% hold on;
% plot(t,elevation,'r--');
% xlabel('Time [s]'),ylabel('Angle [rad]');
% legend('e', 'e*');

