%Climate Model Simulation Crashes Data Set
%
% http://archive.ics.uci.edu/ml/datasets/Climate+Model+Simulation+Crashes#
%
% A = climate
%
% Attribute Information:
% 
% The goal is to predict climate model simulation outcomes (column 21, fail or succeed) given scaled values of climate model input parameters (columns 3-20). 
% 
% Column 1: Latin hypercube study ID (study 1 to study 3)   (ignored)
% 
% Column 2: simulation ID (run 1 to run 180)                (ignored)
% 
% Columns 3-20: values of 18 climate model parameters scaled in the interval [0, 1] 
% 
% Column 21: simulation outcome (0 = failure, 1 = success)

function x = climate

load('popfailures.mat');

x = prdataset(popfailures(:,1:18),popfailures(:,19));
x = setname(x,'Climate Model Simulation Crashes Data Set');