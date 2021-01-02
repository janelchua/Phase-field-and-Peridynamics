function [TotalResidual,gradFF] = travellingWaveEqn5(nx,delta_x,K,M,a)
%% Model2wDamping
%%Solving for Ux
%% Global Constants
global rho gamma_M gamma_m E1 E2 E3 c2 c L1 LL alpha Nu G0
s_dot = M*c;

Ux = a(1,1:nx+1); % U(nx+1)is the constant of integration that has to be minimized too
Phi = a(2,1:nx+1);

%TotalResidual = 0; 
Ux_Residual = 0; Phi_Residual = 0;
gradF = zeros(nx+1,1);

%% Defining Variables: Strain, Phi_x and Phi_xx
U_xx = nan(nx,1);
Phi_x = nan(nx,1);
Phi_xx = nan(nx,1);
for i = 2:nx-1
    U_xx(i) = (Ux(i+1) - Ux(i-1))./(2*delta_x);
    Phi_x(i) = (Phi(i+1) - Phi(i-1))./(2*delta_x);
    Phi_xx(i) = (Phi(i+1) - 2*Phi(i) + Phi(i-1))./(delta_x^2);
end
U_xx(1) = 0;%(Ux(2) - Ux(1))./(delta_x);
U_xx(nx) = 0;%(Ux(nx) - Ux(nx-1))./(delta_x); %you want the stuff in

Phi_x(1) = (Phi(2)-Phi(1))./delta_x;
Phi_x(nx) = (Phi(nx)-Phi(nx-1))./delta_x;
Phi_xx(1) = (Phi(3) - 2*Phi(2) + Phi(1))./(delta_x)^2;
Phi_xx(nx) = (Phi(nx) - 2*Phi(nx-1) + Phi(nx-2))./(delta_x)^2;


%%
for n = 1:nx % All the nodes

    % f
    W0_Phi = -(4.*Ux(n).^2 - 8)./(L1.*(cosh((2.*Phi(n) - 1)./L1) + 1));
    f = (-(W0_Phi) + alpha.*Phi_xx(n));
    
    % Nucleation Term
    H = @(x) 0.5.*(1+tanh(x./L1));
    transitionPoint = (gamma_M+gamma_m)/2;
    G = G0.*(-(1-H(Ux(n)-transitionPoint)).*(H(Phi(n)-0.5)) + (H(Ux(n)-transitionPoint)).*(1-H(Phi(n)-0.5)));
    
    % E_Phi
    E_Phi = (1-H(Phi(n) - 0.5)).*E1 + H(Phi(n) - 0.5).*E3;
    
    Ux_Residual = Ux_Residual + ...
        ( ((M^2)*rho - (1/(c^2)).*E_Phi).*Ux(n)...
        + ((Nu*rho)/(c)).*M.*U_xx(n)...
        - (1/c^2).*Ux(nx+1) )^2;
    
    Phi_Residual = Phi_Residual + ...
        ( -s_dot.*Phi_x(n)...
        - abs(Phi_x(n)).*K.*(f)...
        - G )^2;
end
TotalResidual = Ux_Residual + Phi_Residual;

gradFF = double(gradF);
end
