from gym.envs.registration import register

register(
    id='CPC-v0',
    entry_point='myenv.cartpolecontinuous:CartPoleContinuousEnv'
)
