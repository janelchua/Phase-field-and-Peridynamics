
clc
clear all
close all

%%
%Why it stress_store not storing at all?
%%
%Defining Constants
global E1 E2 E3 c2 L1 LL alpha gamma_M gamma_m
rho = 1;
gamma_M = 1;
gamma_m = 2;
E1= 10; 
E2 = 6 ;  c2 = 16;
E3 = 2;

L1 = 0.1; % Goes into H function in stress response and in plotting
LL = 0.3; % Goes into initial condition for Phi
alpha = 1E1; %% INTERFACE SHOULD BE 10x delta_x
mil = sqrt(E3/rho);
nu = 1;
M = 10; %whacking speed
c = sqrt(E3/rho);

delta_x = 0.1;
delta_t = 0.00001;  
x0 = 0;
L = 200; % Domain size
t0 = 0;
t_tot = 20;
nx = (L/delta_x)+1
nt = round((t_tot/delta_t)+1)

%%Defining variables
store = 6
U = nan(nx,store);
V = nan(nx,store);
Phi = nan(nx,store);
stress_store = nan(nx,store);
n=1; % Counter for storing the values

U_old = nan(nx,1);
V_old = nan(nx,1);
Phi_old = nan(nx,1);

U_new = nan(nx,1);
V_new = nan(nx,1);
Phi_new = nan(nx,1);

x_store = nan(nx,1);
%t_store = nan(1,nt);

%% Discretization
for i = 1:nx
    x = x0+(i-1)*delta_x;
    x_store(i) = x;
end

%% Extra functions needed
H = @(X) 0.5*(1+tanh(X/LL));

%% Initial Conditions for u(x) and phi(x)
position = (5*L/10)
u0 = zeros(1,nx);

left = -10;
right = (((10*4)^0.5)/2) ; %3.16

CCC = left*position;
CC =  left*position - right*position;
for i = 1:nx
    if i <= round(5*nx/10)
        u0(i) = left*x_store(i); %- CCC;
    else
        u0(i) = right*x_store(i) + CC; %- CCC;
    end
end

phi0 = zeros(1,nx);
for i = 1:nx
    phi0(i) = H(x_store(i)-position);
end

% Initial condition for v
v0 = 0;

%% Setting initial conditions
U_old(:) = u0(:);
V_old(:) = v0;
Phi_old(:) = phi0(:);



%% Rest of the timesteps
for j = 2:nt
    j;
    t = t0+(j-1)*delta_t;
    %% Strain, Phi_x and Phi_xx
    Ux = nan(nx,1);
    Phi_x = nan(nx,1);
    Phi_xx = nan(nx,1);
    for i = 2:nx-1
        Ux(i) = (U_old(i+1) - U_old(i-1))/(2*delta_x);
        Phi_x(i) = (Phi_old(i+1) - Phi_old(i-1))/(2*delta_x);
        Phi_xx(i) = (Phi_old(i+1) - 2*Phi_old(i) + Phi_old(i-1))/(delta_x^2);
    end
    Ux(1) = (U_old(2)-U_old(1))/delta_x;       
    Ux(nx) = (U_old(nx)-U_old(nx-1))/delta_x;  
    
    Phi_x(1) = (Phi_old(2)-Phi_old(1))/delta_x;
    Phi_x(nx) = (Phi_old(nx)-Phi_old(nx-1))/delta_x;
    
    Phi_xx(1) = (Phi_old(3) - 2*Phi_old(2) + Phi_old(1))/(delta_x)^2;
    Phi_xx(nx) = (Phi_old(nx)-2*Phi_old(nx-1)+Phi_old(nx-2))/(delta_x)^2;
       
    %% Stress
    stress = nan(nx,1);
    
    for i = 1:nx
        stress(i) = stress_response100(Ux(i),Phi_old(i));
     end
    
    %%Boundary conditions, i=1
    i=1;
    U_new(i) = 0;
    V_new(i) = V_old(i) + delta_t*((1/rho)*(stress(i+1)-stress(i))/(delta_x)...
                                   + nu*(V_old(i+2) - 2*V_old(i+1) + V_old(i))/(delta_x^2));
    Phi_new(i) = 0;
        
    %%Boundary conditions, i=nx
    i = nx;  
    %U_new(i) = u0(i) + M*c*t;
    U_new(i) = u0(i) ;
    V_new(i) = V_old(i) + delta_t*((1/rho)*(stress(i)-stress(i-1))/(delta_x)...
                                   + nu*(V_old(i) - 2*V_old(i-1) + V_old(i-2))/(delta_x^2));
    Phi_new(i) = 1;
    
    %%Interior nodes  %%%%%%%%%%%%% PAUSED HERE
    for i = 2:nx-1 % space    

        U_new(i) = U_old(i) + delta_t*V_old(i);
        V_new(i) = V_old(i) + delta_t*((1/rho)*(stress(i+1)-stress(i-1))/(2*delta_x)...
                                        + nu*(V_old(i+1) - 2*V_old(i) + V_old(i-1))/(delta_x^2));
        Phi_new(i) = Phi_old(i) + delta_t*(1/mil)*(f_driving100(Ux(i),Phi_old(i),Phi_xx(i))...
                                                    + G_nucleation100(Ux(i),Phi_old(i))); %%Continue editing this
        
        

    end
    %V_new(:) = 0.95*V_new(:); %% TURN THIS OFF WHEN YOU PUT IN f
    
        %% Storing the values in U, V and Phi
    if j==2 || j==round(nt*1/(store-1)) || j==round(nt*2/(store-1)) || j==round(nt*3/(store-1)) || j==round(nt*4/(store-1)) || j==nt-1 
        U(:,n) = U_old(:);
        V(:,n) = V_old(:);
        Phi(:,n) = Phi_old(:);
        stress_store(:,n) = stress(:);

        n = n+1;
    end
    
    U_old(:) = U_new(:);
    V_old(:) = V_new(:);
    Phi_old(:) = Phi_new(:);
end

%% POST PROCESSING

%% Ux %% CLEAN UP THIS PART
Ux = zeros(nx,store); 
x_Ux = zeros(nx,1);
Phix = zeros(nx,store); 
x_Phix = zeros(nx,1);
for j = 1:store
    for i = 2:nx-1
        Ux(i,j) = (U(i+1,j)-U(i-1,j))/(2*delta_x);
        x_Ux(i,1) = x_store(i,1);
        
        Phix(i,j) = (Phi(i+1,j)-Phi(i-1,j))/(2*delta_x);
        x_Phix(i,1) = x_store(i,1);
    end
end
Ux;
x_Ux;

%%
figure(1)
grid on
box on
plot(x_store,u0(1,1:nx),'m-','LineWidth',2)
ylabel('initial Condition [u]')
xlabel('position, x')

figure(2)
grid on
box on
plot(x_store,phi0(1,1:nx),'m-','LineWidth',2)
ylabel('initial Condition [phi]')
xlabel('position, x')

%%
range = 1:store;
figure(3)
hold on
plot(x_store,U(:,range),'LineWidth',1.5)
ylabel('U(x)')
xlabel('x')
hold off
% 
figure(4)
hold on
plot(x_Ux(2:nx-1),Ux(2:nx-1,range),'LineWidth',1.5)
ylabel('Ux(x)')
xlabel('x')
hold off
% 
figure(5)
hold on
plot(x_store,Phi(:,range),'LineWidth',1.5)
ylabel('Phi')
xlabel('x')
hold off
% 
% figure(6)
% hold on
% plot(x_Phix(2:nx-1),Phix(2:nx-1,1:50),'LineWidth',1.5)
% ylabel('Phix(x)')
% xlabel('x')
% hold off
% 

% Extra functions needed
H_plot = @(X) 0.5*(1+tanh(X/L1));
figure(7)
hold on
plot(x_store,H_plot(Phi(:,range)-0.5),'LineWidth',1.5) 
ylabel('H(phi-0.5)')
xlabel('x')
hold off
% 
figure(8)
hold on
plot(x_store,stress_store(:,range),'LineWidth',1.5) 
ylabel('stress')
xlabel('x')
hold off

%when phi is the argument, use 0.1 and when x is the argument use 0.3 or
%whatever:)