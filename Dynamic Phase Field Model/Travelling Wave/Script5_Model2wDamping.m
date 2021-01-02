% clear all;
% clc;
% close all;

%notes
%Solving for Ux

%%
%Defining Constants
global rho gamma_M gamma_m E1 E2 E3 c2 c L1 LL alpha Nu G0
rho = 1;
gamma_M = 1;
gamma_m = 2;
E1= 10;
E3 = 2;
E2 = -(E1*gamma_M - E3*gamma_m)/(gamma_M - gamma_m) ;  c2 = (gamma_M*gamma_m*(E3 - E1))/(gamma_M - gamma_m);
c = sqrt(E3/rho);

L1 = 0.1; % Goes into H function in stress response and in plotting
LL = 0.3; % Goes into initial condition for Phi
alpha = 20; %% INTERFACE SHOULD BE 10x delta_x
Nu = 0.2%0.1;%0.01;%0.1;
G0 = 1;

%% Defining Constants that will be changed
Kappa = 0.02%[0.3,1,15]%0.3: black, 1: green, 15: red %This kappa is from v = Kappa*f
%red:0.5, blue:0.2, pink:0.05, black:0.01, green:0.05 different IC

% M3 = [-0.9,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,0.9,1.05]%linspace(-1,1.5,51)% s_dot/c3
% M3 = [-0.75,-0.6,-0.45]
% M3 = [-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,0.9,1.05] %10,10 and ((10*4)^0.5)/(10) ;((10*4)^0.5)/(2) ;

M3 = [1.2]%0.15%0%-0.15%[-0.3]

Length = [40,1,1,1,1]; %red:30,10; blue:20,20, magenta: 10,20; green: 20,10
EE = [10,1,1,1,1];
NN = (Length.*max(EE))+1 ;%largest number of elements needed is 500 so *5 is enough
UxStore = zeros(length(M3),max(NN));
residueStore = zeros(length(M3),1);
ff = zeros(3,length(M3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(Kappa)
    K = Kappa(k);%viscosity
    for p = 1:length(M3)
        %% Defining Grid
        if M3(p)<1
            if Kappa(k)==Kappa(1)
                L = Length(1)
                E = EE(1)
            elseif Kappa(k)==Kappa(2)
                L = Length(1)
                E = EE(1)
            elseif Kappa(k)==Kappa(3)
                L = Length(1)
                E = EE(1)
            end
        elseif M3(p) >= 1
            if Kappa(k)==Kappa(1)
                L = Length(1)
                E = EE(1)
            elseif Kappa(k)==Kappa(2)
                L = Length(1)
                E = EE(1)
            elseif Kappa(k)==Kappa(3)
                L = Length(1)
                E = EE(1)
            end
        end
        x0 = 0
        e = L*E %Number of elements
        nx = e+1 %Number of gridpoints
        delta_x = L/(e);
        
        
        %% Discretizing
        xi = zeros(1,nx);
        for j = 1:nx
            xi(j) = x0 + (delta_x)*(j-1); %xi is the value of x at each node
        end
        
        %% Extra functions needed
        H = @(X) 0.5*(1+tanh(X/LL));
        %% Initial Conditions for U(x) and phi(x)
        position = (5*L/10);
        ux0 = zeros(1,nx);
        
        left = -0.3%((10*4)^0.5)/(10) ;
        right = 2%((10*4)^0.5)/(2) ;
        
        CCC = left*position;
        CC =  left*position - right*position;
        for i = 1:nx
            if i <= round(5*nx/10)
                ux0(i) = left%*xi(i); %- CCC;
            else
                ux0(i) = right%*xi(i) + CC; %- CCC;
            end
        end
        
        phi0 = zeros(1,nx);
        for i = 1:nx
%             phi0(i) = 0.5*H(xi(i)-position)+0.2;
            phi0(i) = H(xi(i)-position);
        end
        
        ux04 = [ux0,0];
        phi04 = [phi0,0];
        initialValues = [ux04;phi04];
        
        %% Solving
        [F] =  @(initialVals) travellingWaveEqn5(nx,delta_x,K,M3(p),initialVals)
        %% fmincon
        A = [];
        zcol = zeros(nx+1,1);
        
        A1 = -1*eye(nx+1,nx+1);
        A11 = [A1,zcol,zcol];
        A2 = eye(nx+1,nx+1);
        A22 = [zcol,A2,zcol];
        A = A11+A22;
        b = zeros(nx+1,1);
        
        Aeq = zeros(nx+1,nx+1);
        %Aeq(:,ceil((nx+1)/2)) = 1;
        Aeq(1,:) = 1;
        beq = zeros(nx+1,1);
        
        options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',10000000,'MaxIterations',10000);
        [Ux_Phi,FF] = fmincon(F,initialValues,[],[],[],[],[],[],[],options);
        %         %[U_Phi,FF] = fmincon(F,initialValues,[],[],Aeq,beq,[],[],[],options);
        %% fminunc
        %         options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',6000000,'MaxIterations',1000); %MaxFunctionEvaluations refers to F... It will call F that many times max
        %         [U_Phi,FF] = fminunc(F,initialValues, options);
        %% fminsearch
        %         options = optimset('MaxFunEvals',500000,'MaxIter',10000); %MaxFunctionEvaluations refers to F... It will call F that many times max
        %         [U_Phi] = fminsearch(F,initialValues,options);
        %%
        UxStore(p,1:nx+1) = Ux_Phi(1,1:nx+1);
        residueStore(p,1) = FF
        %% Post-Processing
        % Uxx
        Uxx = zeros(nx);
        %         x_Uxx = zeros(nx);
        for i = 2:nx-1
            Uxx(i) = (Ux_Phi(1,i+1)-Ux_Phi(1,i-1))./(2*delta_x);
            %             x_Uxx(i) = xi(i);
        end
        Uxx(1) = (Ux_Phi(1,2) - Ux_Phi(1,1))./(delta_x);
        Uxx(nx) = (Ux_Phi(1,nx) - Ux_Phi(1,nx-1))./(delta_x);
        
        %% Plotting Results
        color = 'b--'
        figure(1)
        hold on
        box
        plot(xi(1,1:nx),ux0(1,1:nx),color,'LineWidth',2)
        ylabel('initial Condition, Ux(x)')
        xlabel('position, x')
        hold off
        figure(2)
        hold on
        box
        plot(xi(1,1:nx),phi0(1,1:nx),color,'LineWidth',2)
        ylabel('initial Condition, Phi(x)')
        xlabel('position, x')
        hold off
        
        figure(16)
        hold on
        box
        if Kappa(k)== Kappa(1)
            plot(xi(1,1:nx),Ux_Phi(1,1:nx),color,'LineWidth',2)
        elseif Kappa(k)== Kappa(2)
            plot(xi(1,1:nx),Ux_Phi(1,1:nx),'g-','LineWidth',2)
        elseif Kappa(k)== Kappa(3)
            plot(xi(1,1:nx),Ux_Phi(1,1:nx),'r-','LineWidth',2)
        end
        ylabel('Ux(x)')
        xlabel('position, x')
        
        
        figure(18)
        hold on
        box
        if Kappa(k)== Kappa(1)
            plot(xi(1,1:nx),Ux_Phi(2,1:nx),color,'LineWidth',2)
        elseif Kappa(k)== Kappa(2)
            plot(xi(1,1:nx),Ux_Phi(2,1:nx),'g-','LineWidth',2)
        elseif Kappa(k)== Kappa(3)
            plot(xi(1,1:nx),Ux_Phi(2,1:nx),'r-','LineWidth',2)
        end
        ylabel('Phi(x)')
        xlabel('position, x')
        
        %% H_l
        % Extra functions needed
        H_plot = @(X) 0.5*(1+tanh(X/L1));
        figure(19)
        hold on
        if Kappa(k)== Kappa(1)
            plot(xi(1,1:nx),H_plot(Ux_Phi(2,1:nx)-0.5),color,'LineWidth',2)
        elseif Kappa(k)== Kappa(2)
            plot(xi(1,1:nx),H_plot(Ux_Phi(2,1:nx)-0.5),'g-','LineWidth',2)
        elseif Kappa(k)== Kappa(3)
            plot(xi(1,1:nx),H_plot(Ux_Phi(2,1:nx)-0.5),'r-','LineWidth',2)
        end
        ylabel('H(phi-0.5)')
        xlabel('x')
        hold off
        
        % Plotting (rho*s_dot^2 - E(phi))Ux
        E_Phi = @(Phi) (1-H_plot(Phi-0.5))*E1 + H_plot(Phi-0.5)*E3;
        Temp_Func = @(Ux,Phi) (rho.*(M3(p).*c).^2 - E_Phi(Phi)).*Ux;
        figure(20)
        hold on
        if Kappa(k)== Kappa(1)
            plot(xi(1,1:nx),Temp_Func(Ux_Phi(1,1:nx),Ux_Phi(2,1:nx)),color,'LineWidth',2)
        elseif Kappa(k)== Kappa(2)
            plot(xi(1,1:nx),Temp_Func(Ux_Phi(1,1:nx),Ux_Phi(2,1:nx)),'g-','LineWidth',2)
        elseif Kappa(k)== Kappa(3)
            plot(xi(1,1:nx),Temp_Func(Ux_Phi(1,1:nx),Ux_Phi(2,1:nx)),'r-','LineWidth',2)
        end
        ylabel('(rho*s_dot^2 - E(phi))Ux')
        xlabel('x')
        hold off
        
        % Plotting (rho*s_dot^2 - E(phi))Ux+Uxx(nu*rho*s_dot)
        s_dot = M3(p).*c
        E_Phi = @(Phi) (1-H_plot(Phi-0.5))*E1 + H_plot(Phi-0.5)*E3;
        Temp_Func2 = @(Ux,Uxx,Phi) (rho.*(s_dot).^2 - E_Phi(Phi)).*Ux + Uxx.*(Nu*rho*s_dot);
        figure(6)
        hold on
        if Kappa(k)== Kappa(1)
            plot(xi(1,1:nx),Temp_Func2(Ux_Phi(1,1:nx),Uxx(1:nx),Ux_Phi(2,1:nx)),color,'LineWidth',2)
        elseif Kappa(k)== Kappa(2)
            plot(xi(1,1:nx),Temp_Func2(Ux_Phi(1,1:nx),Uxx(1:nx),Ux_Phi(2,1:nx)),'g-','LineWidth',2)
        elseif Kappa(k)== Kappa(3)
            plot(xi(1,1:nx),Temp_Func2(Ux_Phi(1,1:nx),Uxx(1:nx),Ux_Phi(2,1:nx)),'r-','LineWidth',2)
        end
        ylabel('(rho*s_dot^2 - E(phi))Ux+Uxx(nu*rho*s_dot)')
        xlabel('x')
        hold off

        %% Plotting the f versus M curve
        f = (E1-E3)*(gamma_M*gamma_m - Ux_Phi(1,2).*Ux_Phi(1,nx-1))/2
%         f = (E1-E3)*(gamma_M*gamma_m - 0.2*1.8)/2
        color2 = 'g+-'
        figure(3)
        hold on
        if Kappa(k)== Kappa(1)
            plot(M3(p),f,color2,'linewidth',2)
        elseif Kappa(k)== Kappa(2)
            plot(M3(p),f,'g+','linewidth',2)
        elseif Kappa(k)== Kappa(3)
            plot(M3(p),f,'r+','linewidth',2)
        end
        ylabel('f')
        xlabel('M [s-dot/c3]')
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Storing all the f values for all three omega's together
        if Kappa(k)==Kappa(1)
            ff(1,p) = (E1-E3)*(gamma_M*gamma_m - Ux_Phi(1,2).*Ux_Phi(1,nx-1))/2;
        elseif Kappa(k)==Kappa(2)
            ff(2,p) = (E1-E3)*(gamma_M*gamma_m - Ux_Phi(1,2).*Ux_Phi(1,nx-1))/2;
        elseif Kappa(k)==Kappa(3)
            ff(3,p) = (E1-E3)*(gamma_M*gamma_m - Ux_Phi(1,2).*Ux_Phi(1,nx-1))/2;
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Admissible region
gamma_M = 1;
gamma_m = 2;
E1 = 10; c1 = 0;
E2 = -6 ; c2 = 16;
E3 = 2; c3 = 0;%-12
f_func = @(gamma_neg,gamma_pos) (1/2)*(E1-E3)*(gamma_M*gamma_m - gamma_neg*gamma_pos)
% s_dot [0,-c3]
c3 = sqrt(E3/rho);
c1 = sqrt(E1/rho);
s_dot_a1 = linspace(0,-0.999*c3,100);
M_a1 = s_dot_a1./c3;
f_a1 = zeros(1,length(s_dot_a1));
g_neg = gamma_M;
% g_pos = [(E1/E3)*gamma_M,inf]
syms g_pos
for j = 1:length(s_dot_a1)
    g_POS = solve((s_dot_a1(j))^2 == (c3^2*g_pos - c1^2*g_neg)/(g_pos - g_neg) ,g_pos,'Real',true);
    f_a1(j) = f_func(g_neg,g_POS);
end
%%
% s_dot [0,c*]
c_star = sqrt((c1^2 + gamma_m*c3^2)/(1+gamma_m));
s_dot_a2 = linspace(0,0.999*c_star,50);
M_a2 = s_dot_a2./c3;
f_a2 = zeros(1,length(s_dot_a2));
% g_neg = [(E1/E3)*gamma_m,-1]
g_positive = gamma_m
syms g_neg
for j = 1:length(s_dot_a2)
    g_NEG = solve((s_dot_a2(j))^2 == (c3^2*g_positive - c1^2*g_neg)/(g_positive - g_neg) ,g_neg,'Real',true);
    f_a2(j) = f_func(g_NEG,g_positive);
end
%% Plotting f versus M curve, everything together
JJ = 22;
figure(JJ)
box on
grid on
hold on
plot(M_a1,f_a1,'k-','linewidth',1)
plot(M_a2,f_a2,'k-','linewidth',1)
plot(linspace(-100,100,10),zeros(10),'k--','linewidth',1)
plot(zeros(10),linspace(-100,100,10),'k--','linewidth',1)

plot(M3,ff(1,:),color2,'linewidth',1,'MarkerSize',5)
% plot(M3,ff(2,:),'g*-','linewidth',1,'MarkerSize',5)
% plot(M3,ff(3,:),'ro-','linewidth',1,'MarkerSize',5)
ylabel('f')
xlabel('M')
ylim([-25,19])
xlim([-1,1.25])
%legend('','','','','omega = 0.3','omega = 1','omega = 15')
