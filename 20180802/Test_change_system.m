%% initiallize workspace
clear 
close all
%function [y1,y2,y3,y4,xhat1,xhat2,xhat3,xhat4,t_s] = code20180724_2()
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 10;
seed = 8;
%rng('shuffle');
%seed = randi(1000,1,1);
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
%Y = [0,1];
r_c = 0.1;
n_ori = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% Remove Node
n_ori.remove_node([1])
n_ori.plot()
%% Add Node
n_ori.add_node(m, d, b, Y, r_c);
n_ori.plot()
%% Add Edge
n_ori.add_edge([11,5],[0,1]);
n_ori.plot()
%% Remove Edge
n_ori.remove_edge([7,5]);
n_ori.plot()