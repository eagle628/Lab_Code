close all
clear

%%
data.net_seed = 1;
data.node_number = 6;
data.net_local = 1;
data.Ts = 0.1;
data.belief_N = 24;
data.initial_controller_seed = 240;
data.controller_state_number = 4;
data.value_deep_net_seed = 10;
data.value_fixed_apx_funciton_enable = true;
data.value_initial_lr = 0.1;
data.value_constraint_enable = false;
data.value_trigger_enable = false;
data.value_trigger_period = 1;
data.value_lr_inf_lim = 0.1;
data.policy_pi_sigma = 1;
data.policy_initial_lr = 0.1;
data.policy_constraint_enable = true;
data.policy_pi_grad_enable = false;
data.policy_trigger_enable = false;
data.policy_lr_inf_lim = 1e-6;
data.policy_trigger_period = 1;
data.train_seed = 128;
data.train_initial_seed = 1024;
data.train_Te = 50;
data.train_max_episode = 1000;
data.train_snapshot = 100;
data.train_fixed_apx_fucntion_period = 1;
data.train_gamma = 1;
data.train_render_enable = false;
data.test_initial_seed = 64;
data.test_Te = 100;

save tmp data

%%
value_initial_lrs = [1, 0.1, 0.01];
train_seeds = [32, 64, 128, 256, 512, 1024, 2048, 4096];


for value_initial_lr = value_initial_lrs
    parfor itr1 = 1 : length(train_seeds)
        S = load('tmp', 'data');
        data = S.data;
        data.value_initial_lr = value_initial_lr;
        data.train_seed = train_seeds(itr1);
        test_rl_net_function(data);
    end
end
