function [G] =G_nucleation100(Ux,Phi)
global E1 E2 E3 c2 gamma_m gamma_M L1 alpha

H = @(x) 0.5*(1+tanh(x/L1));

G0 = 20000;
transitionPoint = (gamma_M+gamma_m)/2;
G = G0*(-(1-H(Ux-transitionPoint)).*(H(Phi-0.5)) + (H(Ux-transitionPoint)).*(1-H(Phi-0.5))); %Loading

% if ux0 ==0
%     G = G0*(-(1-H(Ux-gamma_M)).*(H(Phi-0.5)) + (H(Ux-gamma_M)).*(1-H(Phi-0.5))); %Loading
% else
%     G = G0*(-(1-H(Ux-gamma_m)).*(H(Phi-0.5)) + (H(Ux-gamma_m)).*(1-H(Phi-0.5))); %Unloading
% end
end


