function [I_syn, drdt] = inhibitory(topo, V, g_gaba, r)
% Constants
E_Cl = -80;  % mV
alpha_r = 5;  % mM^-1ms^-1
beta_r = 0.18;  % ms^-1
T_max = 1.5;  % mM
K_p = 5;  % mV
V_p = 7;  %mV

V_post = ones(length(V),1) * V';
V_pre = V * ones(1,length(V));

I_syn = g_gaba .* r .* (E_Cl - V_post) .* topo;
T = T_max ./ (1 + exp(-(V_pre - V_p)/K_p)) .* topo;
drdt = (alpha_r * T * (1-r) - beta_r * r) .* topo;
end
