function [psi,M] = go_to_work(A,psi,M,RG_vec)

A = A.*(RG_vec.*RG_vec'); % only green are allowed to immigrate
psi = psi + (psi./M)*A - (psi./M)*A'; % immigration in and out
M = sum (psi);