clear all
close all
%% data load
load('original');
load('iddata_ideal');
load('iddata_plus_1_same');
load('iddata_plus_2_same');
load('iddata_plus_3_same');

%% Ideal
data0 = iddata(v0,w0,Ts);
data1 = iddata(v1,w1,Ts);
data2 = iddata(v2,w2,Ts);
data0.interSample = 'foh';
data1.interSample = 'foh';
data2.interSample = 'foh';
%% polynoiminal
%model_arx = arx(data0,[6,6,0]);
%model_iv4 = iv4(data0,[6,6,0]);
%model_armax = armax(data0,[6,6,6,0]);
%model_oe_0 = oe(data0,[6,6,0]);
%model_oe_1 = oe(data1,[6,6,0]);
%model_oe_2 = oe(data2,[6,6,0]);
%model_bj = bj(data0,[6,6,6,6,0]);
%model_pem= pem(data0,[6,6,6,6,6,0]);

my_pem_opt0 = my_pem_opt(w0,v0,[6,6,0],t);
my_pem_opt1 = my_pem_opt(w1,v1,[6,6,0],t);
my_pem_opt2 = my_pem_opt(w2,v2,[6,6,0],t);
%% state space
%model_n4sid_cva = n4sid(data0,6,'N4weight','CVA','Feedthrough',1,'DisturbanceModel','estimate','Form','canonical');
%model_n4sid_arx = n4sid(data0,6,'N4weight','SSARX','Feedthrough',1,'DisturbanceModel','estimate','Form','canonical');
%model_n4sid_esp = n4sid(data0,6,'N4weight','MOESP','Feedthrough',1,'DisturbanceModel','estimate','Form','canonical');

%% Idealbode
figure
bode(sys_env)
hold on;
%bode(model_arx);
%bode(model_iv4);
%bode(model_armax);
%bode(model_oe_0);
%bode(model_oe_1);
%bode(model_oe_2);
%bode(model_bj);
%bode(model_pem);
bode(my_pem_opt0);
bode(my_pem_opt1);
bode(my_pem_opt2);