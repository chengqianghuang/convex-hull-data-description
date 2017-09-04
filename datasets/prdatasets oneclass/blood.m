% Blood Transfusion Service Center Data Set 
%
% http://archive.ics.uci.edu/ml/datasets/Blood+Transfusion+Service+Center
%
% A = blood
%
% Attribute Information (5 attributes):
% 
% R (Recency - months since last donation), 
% F (Frequency - total number of donation), 
% M (Monetary - total blood donated in c.c.), 
% T (Time - months since first donation), and 
% a binary variable representing whether he/she donated blood in March 2007 (1 stand for donating blood; 0 stands for not donating blood). 

function x = blood

load('transfusion.mat');

x = pr_dataset(transfusion(:,1:4),transfusion(:,5));
x = setname(x,'Blood Transfusion Service Center Data Set');