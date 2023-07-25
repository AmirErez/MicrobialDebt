clear;clc
%Test script for dilution codes. For every variant of the model, we created
%a separate set of scripts.

Bt = 1000;
Ct = 1000;
m = 5;
p = 2;
s = transpose(linspace(0,1,m));
s = [s , 1-s];
P = [0.7; 0.3];
K = 1000;
bo = zeros(m,1) + 1/m;
chi = [0 0; 0.5 0];
spike = 0.05;
batch_limit = 100; 
P_cycle = [0 1; 1 0];
errtype = 1;
variable_K = [1000;2000];
variable_Y = [2;1];
variable_rc = 1e-1;
max_batch = 1000;
plt = 1;

%Simulate
tic
[b,Nr] = serialdilMADAPT(Bt,Ct,s,P,m,p,plt,bo,errtype,K);
%[b,Nr] = serialdilMADAPT_spike(Bt,Ct,s,P,m,p,plt,bo,K,spike,batch_limit);
%[b,Nr] = serialdilMADAPT_troph(Bt,Ct,s,P,m,p,plt,bo,errtype,K,chi);
%[b,Nr] = serialdilMADAPT_bolus_cycling(Bt,Ct,s,P_cycle,m,p,plt,bo,K);
%[b,Nr] = serialdilMADAPT_variable_K_Y(Bt,Ct,s,P,m,p,plt,bo,errtype,variable_K,variable_Y);
%[b,Nr] = serialdil_variable_end(Bt,Ct,s,P,m,p,plt,bo,errtype,K,variable_rc,max_batch);
toc

%Verification 
growth = s*transpose(Nr); %Compute the fitnesses
mult = exp(growth); %Verifying that the growth factors are all the same
g = log((Ct+Bt)/Bt);  %Theoretical fitness prediction
