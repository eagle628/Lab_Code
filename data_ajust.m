% % % % close all
% % % clear
% % % 
% % % % dirname = '.\data_test_global\controller_enable\single_dataset\';
% % % dirname = '.\data_test_global\controller_enable\many_dataset\';
% % % % dirname = '.\data_test_global\controller_disable\many_dataset\';
% % % % dirname = '.\data_test_global\controller_disable\a_single_dataset\';
% % % sss = dir(dirname);
% % % 
% % % dataset_ga = {};
% % % dataset_multi = {};
% % % dataset_slqr = {};
% % % dataset_elqr = {};
% % % 
% % % sys_set1 = {};
% % % sys_set2 = {};
% % % for itr1 = 1 : length(sss)
% % %     if~sss(itr1).isdir
% % %         load(fullfile(dirname,sss(itr1).name))
% % %         if exist('x_multi', 'var')
% % %             f_pso_set = f_multi_set;
% % %             x_pso = x_multi;
% % %         end
% % %         dataset_ga = [dataset_ga, f_ga_set];
% % %         dataset_multi = [dataset_multi, f_pso_set'];
% % %         dataset_slqr = [dataset_slqr, f_slqr_set'];
% % %         dataset_elqr = [dataset_elqr, f_elqr_set'];
% % %         
% % % %         ss_model.set_params(x_pso)
% % % %         [tmpa,tmpb,tmpc,tmpd] = ss_model.gen_ss.get_ss();
% % % %         sys_set1 = [sys_set1, {ss(tmpa,tmpb,tmpc,tmpd,model.Ts)}];
% % % %         
% % % %         ss_model.set_params(x_ga)
% % % %         [tmpa,tmpb,tmpc,tmpd] = ss_model.gen_ss.get_ss();
% % % %         sys_set2 = [sys_set2, {ss(tmpa,tmpb,tmpc,tmpd,model.Ts)}];
% % %         
% % %     end
% % % end
% % % 
% % % % figure
% % % % bode(sys_env,'r')
% % % % hold on
% % % % cellfun(@(x)bode(x,'b'),sys_set1)
% % % % cellfun(@(x)bode(x,'g'),sys_set2)
% % % 
% % % dataset_ga = cell2mat(dataset_ga);
% % % dataset_multi = cell2mat(dataset_multi);
% % % dataset_slqr = cell2mat(dataset_slqr);
% % % dataset_elqr = cell2mat(dataset_elqr);
% % % 
% % % figure
% % % boxplot([dataset_ga',dataset_multi',dataset_slqr',dataset_elqr'],{'GA','Multi','SLQR','ELQR'})

%%
clear

store1 = load('tmp1');
store2 = load('tmp2');
% store1 = load('tmp3');
% store2 = load('tmp4');

max_length = max([length(store1.dataset_ga),length(store1.dataset_multi),length(store2.dataset_ga),length(store2.dataset_multi)]);

dataset_ga1 = nan(1, max_length);
dataset_multi1 = nan(1, max_length);
dataset_ga2 = nan(1, max_length);
dataset_multi2 = nan(1, max_length);
dataset_slqr = nan(1, max_length);
dataset_elqr = nan(1, max_length);

dataset_ga1(1:length(store1.dataset_ga)) = store1.dataset_ga;
dataset_multi1(1:length(store1.dataset_multi)) = store1.dataset_multi;
dataset_ga2(1:length(store2.dataset_ga)) = store2.dataset_ga;
dataset_multi2(1:length(store2.dataset_multi)) = store2.dataset_multi;
dataset_slqr(1:length(store1.dataset_slqr)) = store1.dataset_slqr;
dataset_elqr(1:length(store1.dataset_elqr)) = store1.dataset_elqr;

boxplot([dataset_ga1',dataset_multi1',dataset_ga2',dataset_multi2',dataset_slqr',dataset_elqr'],{'GA1','Multi1','GA2','Multi2','SLQR','ELQR'})
