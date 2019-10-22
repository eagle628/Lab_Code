"""This code is checking operations for chainer."""
from chainer.datasets import mnist
from chainer.datasets import split_dataset_random
from chainer import iterators
from chainer import optimizers
from chainer.dataset import concat_examples
from chainer.cuda import to_cpu
import chainer.links as L
import chainer.functions as F
import random
import numpy as np
import chainer
import matplotlib.pyplot as plt
import myenv
import gym


def reset_seed(seed=0):
    """Set random seed."""
    random.seed(seed)
    np.random.seed(seed)
    if chainer.cuda.available:
        chainer.cuda.cupy.random.seed(seed)


class simple_net(chainer.Chain):
    """Make simple network."""

    def __init__(self, n_mid_units=100, n_out=1):
        """
        Set network Specification.

        ('Parameters')
        n_mid_units : hidden layers number.
        n_out : output number.
        """
        super(simple_net, self).__init__()
        with self.init_scope():
            self.l1 = L.Linear(None, n_mid_units)
            self.l2 = L.Linear(n_mid_units, n_mid_units)
            self.l3 = L.Linear(n_mid_units, n_out)

    def __call__(self, x):
        """Network foward calculation."""
        h1 = F.relu(self.l1(x))
        h2 = F.relu(self.l2(h1))
        return self.l3(h2)

class TD_one_step_AC(object):
    """One-step Actor critic by TD-error"""

    def __init__(self, env, policy_opt, value_opt, policy_sigma=1, gamma=0.9, max_epoch=200, max_episode=200):
        self.env = env
        self.policy_opt = policy_opt
        self.value_opt = value_opt
        self.policy_sigma = policy_sigma
        self.gamma = gamma
        self.max_epoch = max_epoch
        self.max_episode = max_episode

    def train(self):
        average_reward_epoch = []
        for epoch in range(self.max_epoch):
            for episode in range(self.max_episode):
                episode_reward = 0
                y = [self.env.reset()]
                done = False
                while not done:
                    action = self.policy_opt.target(y[-1].reshape(1,4).astype(np.float32))
                    real_action = action + np.random.randn()
                    next_step = env.step(real_action.data[0][0])
                    y.append(next_step[0])
                    reward = chainer.Variable(np.array([next_step[1],], np.float32))
                    done = next_step[2]
                    with chainer.no_backprop_mode() :
                        V_k1 = self.value_opt.target(y[-1].reshape(1,4).astype(np.float32))
                    V_k0 = self.value_opt.target(y[-2].reshape(1,4).astype(np.float32))
                    TD_error = reward + self.gamma*V_k1 - V_k0

                    policy_loss = (real_action-action)*TD_error
                    self.policy_opt.target.cleargrads()
                    policy_loss.backward()
                    policy_opt.update()

                    value_loss = F.mean_squared_error(TD_error, chainer.Variable(np.array([[0]], np.float32)))
                    self.value_opt.target.cleargrads()
                    value_loss.backward()
                    value_opt.update()

                    episode_reward += reward

            average_reward_epoch.append(episode_reward/self.max_episode)
            print(average_reward_epoch[-1])



if __name__ == '__main__':
    reset_seed(12)

    env = gym.make('CPC-v0')
    policy_net = simple_net()
    policy_opt = chainer.optimizers.SGD(lr=0.01).setup(policy_net)
    value_net = simple_net()
    value_opt = chainer.optimizers.SGD(lr=0.01).setup(value_net)

    train = TD_one_step_AC(env, policy_opt, value_opt)
    train.train()
