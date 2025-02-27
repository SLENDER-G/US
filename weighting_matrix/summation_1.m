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
for i=1:8
    B_k = sparse(425984,484137);% update the matrix size based on your selections
    for j=1:16
        count = count+1;
        load(['A_k_',num2str(count),'.mat'])
        B_k = B_k+A_k;
    end
    save(['B_k_',num2str(i)],'B_k')
end