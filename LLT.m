% AERO3260 Assignment 2 Section 1.2
% Lifting Line Solution 
% Author: Thomas Ryan

% Start with a clear workspace
clear; clc

% Define geometry and plow properties
S = 6.92*2; % From VLM program
Cla = 2*pi;
alpha_0 = deg2rad(-3.5);
alpha = deg2rad(5); 
cr = 1.58;
b = 5.84;
AR = (2*b)^2/S;
lambda = 0.5;
N = 8; 
c = linspace(1,lambda,N)'*cr;
n = 1:2:(2*N-1);
theta = 0.5*pi*(1:N)'/N;

% Calculate sigma coefficient and a star to simplify calcs
sigma = Cla*c/(8*b);
alpha_star = alpha - alpha_0;

% Right hand side matrix 
RHS = sigma.*sin(theta).*alpha_star;
coefs = sin(theta.*n).*(sigma.*n+repmat(sin(theta),1,N));

% Calculate A coefficients 
A_vec = coefs\RHS;
gamma_vec = 4*sin(theta*n)*A_vec;

% Calculate aerodynamic coefficients
CL = pi*AR*A_vec(1);
CDi = pi*AR*sum(n'.*A_vec.^2);
delta = sum(n(2:end)'.*A_vec(2:end).^2)/A_vec(1)^2;
e =  1/(1 + delta);

% Print out results
disp(strcat('CL  = ',string(CL)))
disp(strcat('Cdi = ',string(CDi)))
disp(strcat('e   = ',string(e)))


