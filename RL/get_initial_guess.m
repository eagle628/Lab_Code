function get_initial_guess(P, recorder, gen_ss)

x0 = gen_ss.get_params();
opt = optimoptions('fminunc', 'Display', 'iter');
fminunc(@(x) eval_fun(x, P, recorder, gen_ss), x0, opt);

end

function val = eval_fun(x, P, recorder, gen_ss)
    Ts = P.ts;
    controller = gen_ss.get_sys(Ts, x);
%     loop = loopsens(P, controller*recorder);
    sys_loop = feedback(P, -controller*recorder);
    val = norm(sys_loop);

end

