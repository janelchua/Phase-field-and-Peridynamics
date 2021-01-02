function [F,gradFF] = travellingWaveEqn18(N,delta_x,lambda,v,M,rho,dd,a)
%syms u [1 M]
u=a(1:(N+1)); % u(N+1)is the constant of integration dd

%F = []; %System of linear equations M equations M unknowns
F1 = 0; F2 = 0; F3 = 0;
gradF = zeros(N+1,1);
%% Equation for Energy for that number of elements

gamma_M = 1;
gamma_m = 2;
E1 = 10; c1 = 0;
E2 = -6 ; c2 = 16;
E3 = 2; c3 = 0;%-12

for n=1:N
    if u(n) <= gamma_M
        E = E1;
        sigma(n) = E*u(n) + c1;
    elseif u(n) > gamma_M && u(n) <= gamma_m
        E = E2;
        sigma(n) = E*u(n)+ c2;
    else
        E = E3;
        sigma(n) = E*u(n) + c3;
    end
end

C3 = sqrt(E3/rho);
c  = sqrt(E3/rho);
s_dot = M*C3;% Check what this is, velocity of the moving strain discontinuity %%Put into the function

for n = 1:N % All the nodes
    if n ==1 %u(n-1) = u(n)
        F1 = ( (lambda/c^2)*(1/delta_x^2)*(u(n+1)-2*u(n)+u(n)) ...
            + ((v*s_dot)/c^2)*(1/(2*delta_x))*(u(n+1)-u(n))...
            + (M^2)*u(n) ...
            - (1/(rho*c^2))*u(N+1) ...
            - (1/(rho*c^2))*sigma(n) )^2;
            %- sigma(n) )^2;
                    
    elseif n == N %u(n+1) = u(n)
        F3 = ( (lambda/c^2)*(1/delta_x^2)*(u(n)-2*u(n)+u(n-1)) ...
            + ((v*s_dot)/c^2)*(1/(2*delta_x))*(u(n)-u(n-1))...
            + (M^2)*u(n) ...
            - (1/(rho*c^2))*u(N+1) ...
            - (1/(rho*c^2))*sigma(n) )^2;
            %- sigma(n) )^2;
                
    else % the inbetween
        F2 = F2 + ( (lambda/c^2)*(1/delta_x^2)*(u(n+1)-2*u(n)+u(n-1)) ...
            + ((v*s_dot)/c^2)*(1/(2*delta_x))*(u(n+1)-u(n-1))...
            + (M^2)*u(n) ...
            - (1/(rho*c^2))*u(N+1) ...
            - (1/(rho*c^2))*sigma(n) )^2;
            %- sigma(n) )^2;

    end
end
F = double(F1+F2+F3)
%F = F1+F2+F3;
gradFF = double(gradF);
end
