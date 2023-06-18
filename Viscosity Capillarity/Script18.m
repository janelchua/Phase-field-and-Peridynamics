clear all;
clc;
%close all;

%notes
%Non-dimensionalize the equation Properly...
%Minimizing the equation with respect to dd (constant of integration)


%% Defining Constants
omega = [0.3,1,15]%[0.3,1,15]%[0.3,1,15]%0.3: black, 1: green, 15: red
lambda = 0.1%strain gradient coefficient, The greater the value of lambda the less steep the slope

rho = 1%density %when this is too small you get perturbations in your soluton
dd =0 %dd constant of integration should be minimized too

M3 = linspace(-1,1.5,51)%[1,0.95,0.9,0.85,0.8,0.75,0.7]%[-0.8,-0.6,-0.5,-0.4,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.8] % s_dot/c3

LL = [2,10,10,50,150]
EE = [60,10,20,10,1]
NN = (LL.*5)+1 %largest number of elements needed is 500 so *5 is enough
gammaStore = zeros(length(M3),max(NN))
ff = zeros(3,length(M3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(omega)
    v = (2*lambda^(0.5))/(omega(k))%viscosity
    for p = 1:length(M3)
        %% Defining Grid
        if M3(p)<1
            if omega(k)==0.3
                L = LL(4)
                E = EE(4)
            elseif omega(k)==1
                L = LL(2)
                E = EE(2)
            elseif omega(k)==15
                L = LL(2)
                E = EE(2)
            end
        elseif M3(p) >= 1
            if omega(k)==0.3
                L = LL(3)
                E = EE(3)
            elseif omega(k)==1
                L = LL(1)
                E = EE(1)
            elseif omega(k)==15
                L = LL(5)
                E = EE(5)
            end
        end
        x0 = 0
        e = L*E %Number of elements %100 works %500
        N = e+1 %Number of gridpoints
        delta_x = L/(e);
        
        %% Discretizing
        xi = zeros(1,N);
        for j = 1:N
            xi(j) = x0 + (delta_x)*(j-1); %xi is the value of x at each node
        end

        
        %% Initial Conditions
        u03 = zeros(1,N);
        if omega(k)==15 && M3(p) >= 1
            for j = 1:N
                if j <= round(N/2)
                    u03(j) = 0%((10*4)^0.5)/10%1-(p-1)*0.1; %xi is the value of x at each node 0.5,3.5
                else
                    u03(j) = 2.055%((10*4)^0.5)/2%3.5%3+(p-1)*0.05;
                end
            end
        else
            for j = 1:N
                if j <= round(N/2)
                    u03(j) = ((10*4)^0.5)/10%1-(p-1)*0.1; %xi is the value of x at each node 0.5,3.5
                else
                    u03(j) = ((10*4)^0.5)/2%3.5%3+(p-1)*0.05;
                end
            end
        end

        u04 = [u03,0]
        %% Solving
        [F] =  @(u) travellingWaveEqn18(N,delta_x,lambda,v,M3(p),rho,dd,u)
        A = [];
        zcol = zeros(N-1,1);
        A1 = -1*eye(N-1,N-1);
        A11 = [A1,zcol,zcol];
        A2 = eye(N-1,N-1);
        A22 = [zcol,A2,zcol];
        A = A11+A22;
        b = zeros(N-1,1);
        options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',100000000,'MaxIterations',5000);
        [gamma,FF] = fmincon(F,u04,[],[],[],[],[],[], [], options);
        gammaStore(p,1:N) = gamma(1:N);
        
        %% Plotting stress-strain curve
        Points = 200;
        figure(3)
        hold on
        A = 1; % When this A is changed, the 'slope' of the curve will change
        u2 = linspace(-1,4,Points);
        sigma = 4.*u2.*((u2.^2)-A);
        plot(u2,sigma,'m-','Linewidth',1.5)
        ylabel('sigma-hat')
        xlabel('gamma')
        grid on
        hold off
        
        figure(4)
        hold on
        gamma_M = 1;
        gamma_m =+2;
        E1 = 10; c1 = 0;
        E2 = -6 ; c2 = 16;
        E3 = 2; c3 = 0;%-12
        
        sigma2 = zeros(1,Points);
        for n=1:Points
            if u2(n) < gamma_M
                E = E1;
                sigma2(n) = E*u2(n) + c1;
            elseif u2(n) > gamma_M && u2(n) < gamma_m
                E = E2;
                sigma2(n) = E*u2(n) + c2;
            else
                E = E3;
                sigma2(n) = E*u2(n) + c3;
            end
        end
        plot(u2,sigma2,'b-','Linewidth',2)
        ylabel('sigma-hat')
        xlabel('strain')
        grid on
        hold off
        
        %% Plotting Results
        figure(1)
        hold on
        box
        plot(xi(1,1:N),u03(1,1:N),'m-','LineWidth',2)
        ylabel('initial Condition')
        xlabel('position, x')
        hold off
        
        figure(16)
        hold on
        box
        if omega(k)== 0.3
            plot(xi(1,1:N),gamma(1,1:N),'k-','LineWidth',2)
        elseif omega(k)== 1
            plot(xi(1,1:N),gamma(1,1:N),'g-','LineWidth',2)
        elseif omega(k)== 15
            plot(xi(1,1:N),gamma(1,1:N),'r-','LineWidth',2)
        end
        ylabel('gamma(x)')
        xlabel('position, x')
        
        %% Plotting the f versus M curve
        f = (E1-E3)*(gamma_M*gamma_m - gamma(1)*gamma(N))/2
        f02 = (E1-E3)*(gamma_M*gamma_m)/2
        M3(p)
        f/f02
        figure(2)
        hold on
        if omega(k)== 0.3
            plot(M3(p),f,'k+','linewidth',2)
        elseif omega(k)== 1
            plot(M3(p),f,'g+','linewidth',2)
        elseif omega(k)== 15
            plot(M3(p),f,'r+','linewidth',2)
        end
        ylabel('f')
        xlabel('M [s-dot/c3]')
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Storing all the f values for all three omega's together
    if omega(k)==0.3
        ff(1,p) = (E1-E3)*(gamma_M*gamma_m - gammaStore(p,1)*gammaStore(p,N))/2;
    elseif omega(k)==1
        ff(2,p) = (E1-E3)*(gamma_M*gamma_m - gammaStore(p,1)*gammaStore(p,N))/2;
    elseif omega(k)==15
        ff(3,p) = (E1-E3)*(gamma_M*gamma_m - gammaStore(p,1)*gammaStore(p,N-10))/2;
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
JJ = 2;
figure(JJ)
box on
grid on
hold on
plot(M_a1,f_a1,'k-','linewidth',1) %Line at the boundary of the admissible region, negative M region.
plot(M_a2,f_a2,'k-','linewidth',1) %Line at the boundary of the admissible region, positive M region.
plot(linspace(-100,100,10),zeros(10),'k--','linewidth',1)
plot(zeros(10),linspace(-100,100,10),'k--','linewidth',1)

plot(M3,ff(1,:),'k+-','linewidth',1,'MarkerSize',5)
plot(M3,ff(2,:),'g*-','linewidth',1,'MarkerSize',5)
plot(M3,ff(3,:),'ro-','linewidth',1,'MarkerSize',5)
ylabel('f')
xlabel('M')
ylim([-25,19])
xlim([-1,1.25])
%legend('','','','','omega = 0.3','omega = 1','omega = 15')
