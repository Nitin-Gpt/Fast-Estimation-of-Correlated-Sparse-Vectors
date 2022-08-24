function [phi,W] = dictionary_generation(N,M, array_response)
% dictionary_generation: generates the dictionary matrix
% Here the input arguments are just N and M.
W=(1/sqrt(N))*exp(1i*2*pi*rand(M,N));
phi=W*array_response;
end