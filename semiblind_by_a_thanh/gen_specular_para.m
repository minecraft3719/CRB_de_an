function [fading,delay,DOA] = gen_specular_para(M,N_t)
% Number of multipaths : M
% Number of TX         : N_t


% Specular (Some components are stronger than others)
% Simulated case: M = 4, N_t = 2;
if M  == 4
    fading = [0.9 0.8 0.04 0.02 ; 0.9 0.7 0.05 0.03]';
    delay  = [0.1 0.2 0.7 0.8 ; 0.2 0.3 0.7 0.9]';
    DOA    = pi./[2 4 6 8;3 5 7 9]';
elseif M == 5
    % Simulated case: M = 5, N_t = 2;
    fading = [0.9 0.8 0.7 0.02 0.01; 0.9 0.7 0.5 0.3 0.1]';
    delay  = [0.1 0.2 0.3 0.7 0.8; 0.2 0.3 0.4 0.9 0.8]';
    DOA    = pi./[2 4 6 8 9 ;3 5 7 9 8]';
end
% non-specular
% fading = rand(M,N_t);
% delay  = rand(M,N_t);
% DOA    = pi * rand(M,N_t);
end