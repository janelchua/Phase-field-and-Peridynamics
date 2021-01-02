function [f] =f_driving100_1(Ux,Phi,Phi_xx)
global L1 alpha

W0_Phi = -(4*Ux^2 - 8)/(L1*(cosh((2*Phi - 1)/L1) + 1));
f = (-(W0_Phi) + alpha*Phi_xx);


end
