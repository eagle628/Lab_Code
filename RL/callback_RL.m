function callback_RL(episode, t, apx_x_all, cost_history, reward_history)
    subplot(3,1,1)
    plot(t, apx_x_all)
    title(['Episode-',num2str(episode)])
    grid on
    subplot(3,1,2)
    plot(nonzeros(cost_history),'-r')
    ylabel('Cost')
    subplot(3,1,3)
    plot(nonzeros(reward_history),'-b')
    ylabel('Culumative Reward')
    drawnow
end
