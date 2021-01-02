function [stress] = stress_response100(Ux,Phi)
global E1 E3 L1

H_energy = @(X) 0.5*(1+tanh(X/L1));

stress = E3*Ux*(H_energy(Phi-0.5)) + E1*Ux*(1-H_energy(Phi-0.5));


end
