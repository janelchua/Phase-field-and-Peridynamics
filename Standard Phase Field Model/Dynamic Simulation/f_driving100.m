function [f] =f_driving100(Ux,Phi,Phi_xx)
global E1 E2 E3 c2 gamma_m gamma_M L1 alpha

Phi_stiffness = 1000;
firstTerm = c2*gamma_m - c2*gamma_M + Phi_stiffness*(Phi^2)*(2*Phi - 2) + (E1*gamma_M^2)/2 + ...
    (E2*gamma_M^2)/2 - (E2*gamma_m^2)/2 - ...
    (E3*gamma_m^2)/2 + Phi_stiffness*2*Phi*(Phi - 1)^2 + ...
    (Ux^2*(tanh((Phi - 1/2)/L1)^2 - 1)*(E1 - E3))/(4*L1);


f = (-(firstTerm) + alpha*Phi_xx);



end


