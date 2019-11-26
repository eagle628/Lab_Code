close all
clear
%%
net_seed1 = 1;
Node_number = 6;
belief_Ns = [1, 3, 6, 8, 10];
deepnet_seed = 10;
train_policys = [1, 2, 3, 4, 5, 6, 7, 8];
policy_initial_lrs = 1;
value_initial_lrs = [1, 0.1, 0.01];

for belief_N = belief_Ns
    for policy_initial_lr = policy_initial_lrs
        for value_initial_lr = value_initial_lrs
            parfor train_policy = train_policys
                filename = strcat(...
                                'net_seed_',num2str(net_seed1),'_',...
                                'net_node_',num2str(Node_number),'_',...
                                'belief_N_',num2str(belief_N),'_',...
                                'deep_seed_',num2str(deepnet_seed),'_',...
                                'train_seed',num2str(train_policy),'_',...
                                'policy_initial_lr',num2str(policy_initial_lr),'_',...
                                'value_initial_lr',strrep(num2str(value_initial_lr),'.',''),'_'...
                                );
                rl_test_func(net_seed1,Node_number,belief_N,deepnet_seed,train_policy,policy_initial_lr,value_initial_lr,filename)
            end
        end
    end
end