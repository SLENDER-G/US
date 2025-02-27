% Author: Sobhan Goudarzi
% IMage Processing And Characterization of Tissue (IMPACT) Group
% Concordia University
% email address: sobhan.goudarzi@concordia.ca
% August 2022
clc
clear all
close all
warning off
%%
count = 0;
Phi = sparse(425984,484137);% update the matrix size based on your selections
for i=1:8
    load(['B_k_',num2str(i),'.mat'])
    Phi = Phi+B_k;
end
save('Phi','Phi')
