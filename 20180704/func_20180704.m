function name = func_20180704(n,c_n,n_n,id_in_p,noise_p,id_method,ITR,location)

    %% Figure Name
    ob_v  = {strcat('v_node',num2str(c_n))};
    ob_w  = {strcat('w_node',num2str(c_n))};
    ob_y  = {strcat('y_node',num2str(c_n))};
    ID_in = {strcat('d_node',num2str(c_n))};
    Noise = cell(1,numel(n_n));
    for i = 1 : numel(n_n)
        Noise(i) = {strcat('d_node',num2str(n_n(i)))};
    end
    name_c = '';
    for i = 1 : numel(c_n)
        name_c = strcat(name_c,num2str(c_n(i)),'+');
    end
    name_n = '';
    for i = 1 : numel(n_n)
        name_n = strcat(name_n,num2str(n_n(i)),'+');
    end
    str = replace(num2str(noise_p),'.','%');
    name = strcat(id_method,'_N',num2str(n.N),'_Cnode',name_c,'_Nnode',name_n,'_IDinput',num2str(id_in_p),'_Noise',str);
    %% add I/O port for identificaiton
    sys_org = n.get_sys();
    for i =  1 : n.N
        sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
    end
    min = -10;
    max =  10;
    [~, sys_env] = n.get_sys_local(c_n);
    [mag_env,phase_env,omega] = bode(sys_env,{10^min,10^max});
    mag_env=reshape(mag_env(1,1,:),[],1);
    phase_env=reshape(phase_env(1,1,:),[],1);

    %% identification simulation
    %ITR = 100;
    rng(28);
    rand_s = randi(1000,2,ITR);
    Mag = zeros(length(omega),ITR);
    Pha = zeros(length(omega),ITR);

    % Model dimension
    dim = 6;
    % OE option
    opt = oeOptions;
    opt.Display = 'on';
    opt.SearchOption.MaxIter = 50;
    opt.SearchOption.Tolerance = 1e-2;
    opt.SearchMethod = 'lm';
    
    %ID input Number
    N = 100000;
    Ts = 0.01;
    t = (0:N-1)'*Ts;
    % simulation
    lsim_type = 'foh';
    % resample number
    if string(id_method) == 'OE'
        R = 1;
    else
        R = 1;
    end
    model_set = cell(1,ITR);
    for i = 1 : ITR
        % ID input
        cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(c_n),'RandomStream','mt19937ar with seed','Seed',rand_s(1,i));
        d1 = cn1();
        d1 = d1*id_in_p;
        d1(:,1) = zeros(N,1);
        cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(n_n),'RandomStream','mt19937ar with seed','Seed',rand_s(2,i));
        d2 = cn2();
        d2 = d2*noise_p;
%         for itr = 1:numel(n_n)
%             d2(:,itr*2-1) = zeros(N,1);
%         end
        % Response of v&w 
        v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
        w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
        data_o = iddata(v,w,Ts);
        data = resample(data_o,1,R);
        if string(id_method) == 'OE'
            init_sys = oe(data,[dim,dim,0],opt);
        elseif string(id_method) == 'ARMAX'
            %init_sys = pem(data,[dim,dim+1,,0],opt);
            init_sys = armax(data,[dim,dim+1,dim,0],opt);
        elseif string(id_method) == 'ARX'
            init_sys = arx(data,[dim,dim+1,0],opt);
        end
        model_p = d2c(init_sys);
        model = ss(model_p);
        model_set(i) = {model};

        [mag,phase] = bode(model,omega);
        Mag(:,i) =reshape(mag(1,1,:),[],1);
        Pha(:,i) =reshape(phase(1,1,:),[],1);
    end

    %% Bode Plot
    fig1 = figure_config.set_figure_bode(name);
    figure_config.plot_bode(fig1.axes,omega,mag_env,phase_env,{'r-','linewidth',3});
    figure_config.plot_bode(fig1.axes,omega,Mag,Pha,{'b-'});
    %% Add Controller  (Retro Fit)
    Q = kron(eye(1*numel(c_n)),diag([1,1000]));
    R = kron(eye(1*numel(c_n)),diag([1e-3]));
    n.controllers = {};
    n.add_controller(c_n,ss([],[],[],0), Q, R);
    controlled_sys_simple = n.get_sys_controlled(sys_org);
    controlled_sys_simple_K = n.controllers{1}.sys_K;
    %% Add Controller (Extend Retro fit)
    Q = kron(eye(1*numel(c_n)),diag([1,1000]));
    R = kron(eye(1*numel(c_n)),diag([1e-3]));
    % Ideal 
    n.controllers = {};
    %n.add_controller(c_n,sys_env, Q, R);
    n.add_controller(c_n, Q, R);
    controlled_sys_I = n.get_sys_controlled(sys_org);
    controlled_sys_I_K = n.controllers{1}.sys_K;

    % Identification
    controlled_sys_set = cell(1,ITR);
    controlled_sys_K_set = cell(1,ITR); 
    for i = 1 : ITR
        n.controllers = {};
        n.add_controller(c_n,model_set{i}, Q, R);
        controlled_sys_set(i) = {n.get_sys_controlled(sys_org)};
        controlled_sys_K_set(i) = {n.controllers{1}.sys_K};
    end

    %% Response simulation
    s_t = 300;%simulation_time
    time = 0:0.01:s_t-0.01;
    time = time';
    % Original Response
    yO = impulse(sys_org(ob_y,ID_in), time);
    % Ideal Response
    yI = impulse(controlled_sys_I(ob_y,ID_in), time);
    yI_de = impulse(controlled_sys_I_K('y','x'),time);
    % Retro Responsse
    y_simple = impulse(controlled_sys_simple(ob_y,ID_in), time);
    y_simple_de = impulse(controlled_sys_simple_K('y','x'),time);
    % Identification Response
    y_id_set = cell(1,ITR);
    y_id_de_set = cell(1,ITR); 
    for i = 1 : ITR
        sys_re = controlled_sys_set{i};
        sys_de = controlled_sys_K_set{i};
        y_id_set(i) = {impulse(sys_re(ob_y,ID_in), time)};
        y_id_de_set(i) = {impulse(sys_de('y','x'),time)};
    end

    %% Response plot
    check = 1;
    state = 2;
    fig2 = figure('Name',strcat('Response_',name));
    % All response
    ax1 = subplot(3,1,1);
    hold on;
    for i = 1 : ITR
        id = plot(time,y_id_set{i}(:,state,check),'b');
    end
    ori = plot(time,yO(:,state,check),'g','linewidth',1.0);
        % Extend Ideal
    ideal = plot(time,yI(:,state,check),'r','linewidth',1.5);
    legend([ori,ideal,id],'original','Simple','Extend');
        % Retro
    %retro = plot(time,y_simple(:,state,check),'r','linewidth',1.5);
    %legend([ori,retro,id],'original','Retro','Identification');
    % Controlled Response
    ax2 = subplot(3,1,2);
    hold on
        % identification
    for i = 1 : ITR
        plot(time,y_id_de_set{i}(:,state,check),'b');
    end
        % Extend Ideal
    plot(time,yI_de(:,state,check),'r','linewidth',1.5);
        % retro
    %plot(time,y_simple_de(:,state,check),'r','linewidth',1.5);
    % Error Response
    ax3 = subplot(3,1,3);
    hold on
        % Identification
    for i = 1 : ITR
        plot(time,y_id_set{i}(:,state,check)-y_id_de_set{i}(:,state,check),'b');
    end
        % Extend Ideal
    plot(time,yI(:,state,check)-yI_de(:,state,check),'r','linewidth',1.5);
        % retro
    %plot(time,y_simple(:,state,check)-y_simple_de(:,state,check),'r','linewidth',1.5);
    
    ax1.Title.Interpreter = 'latex';
    ax1.Title.FontSize = 12;
    ax1.Title.String = '$y$';
    ax1.Box = 'on';
    ax2.Title.Interpreter = 'latex';
    ax2.Title.FontSize = 12;
    ax2.Title.String = '$\hat{y}$';
    ax2.Box = 'on';
    ax3.Title.FontSize = 12;
    ax3.Title.Interpreter = 'latex';
    ax3.Title.String = '$\tilde{y}$';
    ax3.Box = 'on';

    %% save
    if ~isempty(location)
        savefig(fig1.f,strcat(location,'\',fig1.f.Name),'compact')
        savefig(fig2,strcat(location,'\',fig2.Name),'compact')
        close all
        clear fig1 fig2
        save(strcat(location,'\',name),'-v7.3')
    end
end
