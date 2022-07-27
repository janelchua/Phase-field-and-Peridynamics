clc
clear all
close all
%% Defining Grid
L0 = 1;

II = 200*L0; % the interval over which R is computed
L = II + 5*L0; % 0.5 of entire interval
xx = L/2; % where the strain discontinuity is
delta_x = L0/8; %small_L_0/10
e = ((L)/delta_x) % number of elements between 0 to L
N = e + 1 % number of gridpoints for 0 to L, N=2B+1+I
B = ((L-II)/(2*delta_x))
I = ((II)/(delta_x))+1 % number of gridpoints in II

%% Defining Constants
rho = 1%density %when this is too small you get perturbations in your soluton
v = 0.5 %20*L0,L0/20 doesnt work for M3 = 0.15
M3 = [0.9]%linspace(-1.05,1.50,18)%[0,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1]%linspace(-1,1.2,12)%linspace(-1,1,41)%[1,0.95,0.9,0.85,0.8,0.75,0.7]%[-0.8,-0.6,-0.5,-0.4,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.8] % s_dot/c3
length(M3)

gamma_M = 1;
gamma_m = 2;
E1 = 10; c1 = 0;
E2 = -6 ; c2 = 16;
E3 = 2; c3 = 0;%-12

%% Discretizing
x = zeros(1,N);
for j = 1:N
    x(j) = 0 + (delta_x)*(j-1); %x is the value of x at each node from -L to L
end
%%
x
displStore1 = zeros(length(M3),N)%v = 1
gammaStore1 = zeros(length(M3),N)
R1 = zeros(length(M3),1)
displStore03 = zeros(length(M3),N)%v = 0.3
gammaStore03 = zeros(length(M3),N)
R03 = zeros(length(M3),1)
displStore003 = zeros(length(M3),N)%v = 0.03
gammaStore003 = zeros(length(M3),N)
R003 = zeros(length(M3),1)

ff = zeros(3,length(M3)) %3 cause 3 different values of v

%%
for k = 1:length(v)   
    for p = 1:length(M3)
    
        %% Initial Conditions
        u05 = zeros(1,N);
        left = (((10*4)^0.5)/10) ;
        right = (((10*4)^0.5)/2) ;
        CCC = left*(L/2); 
        CC =  left*(L/2) - right*(L/2); 
        for j = 1:N
            if j <= round(N/2)
                u05(j) = left*x(j)% - CCC; 
            else
                u05(j) = right*x(j) + CC %- CCC ;   
            end
        end
        u05

        %% Solving
        [F] =  @(u) peri103_3(N,B,I,delta_x,M3(p),rho,L0,u,x,v(k),E1,E2,E3,c1,c2,c3,gamma_M,gamma_m)
        A = [];
        zcol = zeros(N,1);
        A1 = -1*eye(N,N);
        A11 = [A1,zcol,zcol];
        A2 = eye(N,N);
        A22 = [zcol,A2,zcol];
        A = A11+A22;
        b = zeros(N,1);
        
        Aeq = zeros(N,N);
        beq = zeros(N,1);
        
        options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',6000000,'MaxIterations',5000); %MaxFunctionEvaluations refers to F... It will call F that many times max
        [disp,FF] = fmincon(F,u05,[],[],[],[],[],[],[], options);

        if v(k)==v(1)
            displStore1(p,1:N) = disp(1:N);
            R1(p,1) = FF;
        elseif v(k)==v(2)
            displStore03(p,1:N) = disp(1:N);
            R03(p,1) = FF;
        elseif v(k)==v(3)
            displStore003(p,1:N) = disp(1:N);
            R003(p,1) = FF;
        end

        gamma1 = zeros(1,size(disp,2)); %gamma
        gamma03 = zeros(1,size(disp,2)); %gamma
        gamma003 = zeros(1,size(disp,2)); %gamma
        x_gamma = zeros(1,size(disp,2));
        for i = 2:size(disp,2)-1
            gamma1(i) = (displStore1(p,i+1)-displStore1(p,i-1))/(2*(x(i+1)-x(i)));
            gammaStore1(p,i) = gamma1(i);
            
            gamma03(i) = (displStore03(p,i+1)-displStore03(p,i-1))/(2*(x(i+1)-x(i)));
            gammaStore03(p,i) = gamma03(i);
            
            gamma003(i) = (displStore003(p,i+1)-displStore003(p,i-1))/(2*(x(i+1)-x(i)));
            gammaStore003(p,i) = gamma003(i);
            x_gamma(i) = x(i);
        end
        gamma1
        x_gamma
        %end
        
        %% Plotting Results
        figure(1)
        hold on
        grid on
        box on
        plot(x(1,1:N),u05(1,1:N),'b-','LineWidth',2)
        ylabel('initial Condition')
        xlabel('position, x')
        hold off
        
        figure(2)
        hold on
        box
        plot(x(1,1:N),displStore1(p,1:N),'b-','LineWidth',2)
        plot(x(1,1:N),displStore03(p,1:N),'b-','LineWidth',2)
        plot(x(1,1:N),displStore003(p,1:N),'r-','LineWidth',2)
        ylabel('u(x)')
        xlabel('position, x')
        hold off
        
        figure(3)
        hold on
        box on
        plot(x_gamma(1,2:N-1),gammaStore1(p,2:N-1),'b-','LineWidth',2)
        plot(x_gamma(1,2:N-1),gammaStore03(p,2:N-1),'b-','LineWidth',2)
        plot(x_gamma(1,2:N-1),gammaStore003(p,2:N-1),'r-','LineWidth',2)
        ylabel('gamma(x)')
        xlabel('position, x')
        hold off
        
        % Plotting the f versus M curve
        points =  10
        figure(6)
        hold on
        ylabel('f')
        xlabel('M [s-dot/c3]')
        f = zeros(3,length(M3));
        if v(k)== v(1)
            top = mean(gammaStore1(p,B+round((6/8)*I):points+B+round((6/8)*I)));
            bottom = mean(gammaStore1(p,B+round((2/8)*I)-points:B+round((2/8)*I)));
            %top = gammaStore1(p,N-B);
            %bottom = gammaStore1(p,B);
            f(1,p) = (E1-E3)*(gamma_M*gamma_m - bottom*top)/2;
            plot(M3(p),f(1,p),'g+-','linewidth',2) %v = 1
            
        elseif v(k)== v(2)
            top = mean(gammaStore03(p,B+round((6/8)*I):points+B+round((6/8)*I)));
            bottom = mean(gammaStore03(p,B+round((2/8)*I)-points:B+round((2/8)*I)));
            %top = gammaStore03(p,N-B-5);
            %bottom = gammaStore03(p,B+5);
            f(2,p) = (E1-E3)*(gamma_M*gamma_m - bottom*top)/2;
            plot(M3(p),f(2,p),'k+-','linewidth',2)%v = 0.3
            
        elseif v(k)== v(3)
            top = mean(gammaStore003(p,B+round((6/8)*I):points+B+round((6/8)*I)));
            bottom = mean(gammaStore003(p,B+round((2/8)*I)-points:B+round((2/8)*I)));
            f(3,p) = (E1-E3)*(gamma_M*gamma_m - bottom*top)/2;
            plot(M3(p),f(3,p),'r+-','linewidth',2)%v = 0.03
        end
        
        
        
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checking Stuff
test = 1

figure(100)
hold on
grid on
box on
plot(x(1,1:N),displStore1(test,1:N),'k-','LineWidth',2)
ylabel('u(x)')
xlabel('position, x')

figure(101)
hold on
grid on
box on
plot(x_gamma(1,2:N-1),gammaStore1(test,2:N-1),'k-','LineWidth',2)
ylabel('gamma(x)')
xlabel('position, x')

%% Plotting the f versus M curve

% Storing all the f values for all three v's together
points = 100;
ff = zeros(3,length(M3))
for k = 1:length(v)
    for p = 1:length(M3)
        if v(k)== v(1)
            %top = mean(gammaStore1(p,B+round((6/8)*I):points+B+round((6/8)*I)))
            %bottom = mean(gammaStore1(p,B+round((2/8)*I)-points:B+round((2/8)*I)))
            %
            top = 1.747%gammaStore1(p,N-B);
            bottom = 0.08135%gammaStore1(p,B);
            ff(1,p) = (E1-E3)*(gamma_M*gamma_m - bottom*top)/2;
        elseif v(k)== v(2)
            top = mean(gammaStore03(p,B+round((6/8)*I):points+B+round((6/8)*I)))
            bottom = mean(gammaStore03(p,B+round((2/8)*I)-points:B+round((2/8)*I)))
            ff(2,p) = (E1-E3)*(gamma_M*gamma_m - bottom*top)/2;
        elseif v(k)== v(3)
            top = mean(gammaStore003(p,B+round((6/8)*I):points+B+round((6/8)*I)))
            bottom = mean(gammaStore003(p,B+round((2/8)*I)-points:B+round((2/8)*I)))
            ff(3,p) = (E1-E3)*(gamma_M*gamma_m - bottom*top)/2;
        end
    end
end
figure(1)
hold on
box on
grid on
plot(M3,ff(1,:),'r+','linewidth',2)%v = 1
%plot(M3,ff(2,:),'b+-','linewidth',2)%v = 0.3
%plot(M3,ff(3,:),'r+-','linewidth',2)%v = 0.03
ylabel('fftest')
xlabel('M [s-dot/c3]')
%ylim([-8,6])
%xlim([-1,0.2])


