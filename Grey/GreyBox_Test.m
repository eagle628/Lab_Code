clear
close all
load('IDDATA_N30');

sys_env = c2d(sys_env,Ts);
sys_local = c2d(sys_local,Ts);

data = iddata(zeros(length(v),1),[v,w],Ts);
data.InputName = {'v','w'};
data.OutputName = {'y'};