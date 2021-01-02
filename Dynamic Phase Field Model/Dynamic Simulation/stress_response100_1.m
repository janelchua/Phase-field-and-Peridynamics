function [stress] = stress_response100_1(Ux,Phi)
global E1 E3 L1

H_energy = @(X) 0.5*(1+tanh(X/L1));

stress = E3*Ux*(H_energy(Phi-0.5)) + E1*Ux*(1-H_energy(Phi-0.5));

%stress = (Ux*(E1 + E3 - E1*tanh((2*Phi - 1)/(2*L1)) + E3*tanh((2*Phi - 1)/(2*L1))))/2;

end
