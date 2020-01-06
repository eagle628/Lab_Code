"""Define chainer Deep Network model."""
import chainer
import chainer.links as L
import chainer.functions as F
import numpy as np
import random


def reset_seed(seed=0):
    """Set random seed."""
    random.seed(seed)
    np.random.seed(seed)
    if chainer.cuda.available:
        chainer.cuda.cupy.random.seed(seed)


class simple_net(chainer.Chain):
    """Make Simple network."""

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

    def predict(self, x, enable_backprop):
        """Network foward calculation."""
        with chainer.using_config('enable_backprop', enable_backprop):
            h1 = F.relu(self.l1(x))
            h2 = F.relu(self.l2(h1))
            return self.l3(h2)

    def reset_state(self):
        """Reset mid unit states."""
        pass


class test_clss():
    """test class."""

    def __init__(self, a, b):
        """
        Set test Specification.

        ('Parameters')
        a : property1
        b : property2
        """
        self.a = a
        self.b = b

    def func1(self, c):
        """Test funciton1."""
        print('test')
        d = self.a+c
        e = self.b+c
        return (d, e, type(c))

    def func2(self, a):
        """Print args class name."""
        return type(a)


class simple_linear5_net(chainer.Chain):
    """Make linear 5layers network."""

    def __init__(self, n_mid_units=100, n_out=1):
        """
        Set network Specification.

        ('Parameters')
        n_mid_units : hidden layers number.
        n_out : output number.
        """
        super(simple_linear5_net, self).__init__()
        with self.init_scope():
            self.l1 = L.Linear(None, n_mid_units)
            self.l2 = L.Linear(n_mid_units, n_mid_units)
            self.l3 = L.Linear(n_mid_units, n_mid_units)
            self.l4 = L.Linear(n_mid_units, n_mid_units)
            self.l5 = L.Linear(n_mid_units, n_out)

    def __call__(self, x):
        """Network foward calculation."""
        h1 = F.relu(self.l1(x))
        # h2 = F.dropout(F.relu(self.l2(h1)), ratio=0.25)
        h2 = F.dropout(F.relu(self.l2(h1)), ratio=0)
        h3 = F.relu(self.l3(h2))
        # h4 = F.dropout(F.relu(self.l4(h3)), ratio=0.25)
        h4 = F.dropout(F.relu(self.l4(h3)), ratio=0)
        return self.l5(h4)

    def predict(self, x, enable_backprop):
        """Network foward calculation."""
        with chainer.using_config('enable_backprop', enable_backprop):
            return self(x)

    def reset_state(self):
        """Reset mid unit states."""
        pass


class Base_Net(chainer.Chain):
    """Base network class."""

    def predict(self, x, enable_backprop):
        """Network foward calculation."""
        with chainer.using_config('enable_backprop', enable_backprop):
            return self(x)

    def reset_state(self):
        """Reset mid unit states."""
        pass


class linear3_net(Base_Net):
    """Net of lieanr 3 lasyers."""

    def __init__(self, n_mid_units=100, n_out=1):
        """
        Set network Specification.

        ('Parameters')
        n_mid_units : hidden layers number.
        n_out : output number.
        """
        super(linear3_net, self).__init__()
        with self.init_scope():
            self.l1 = L.Linear(None, n_mid_units)
            self.l2 = L.Linear(n_mid_units, n_mid_units)
            self.l3 = L.Linear(n_mid_units, n_out)
            self.bn1 = L.BatchNormalization(n_mid_units)
            self.bn2 = L.BatchNormalization(n_mid_units)

    def __call__(self, x):
        """Network foward calculation."""
        h1 = F.relu(self.bn1(self.l1(x)))
        h2 = F.relu(self.bn2(self.l2(h1)))
        return self.l3(h2)


class linear5_net(Base_Net):
    """Net of linear 5 layers."""

    def __init__(self, n_mid_units=100, n_out=1, ratio=0.25):
        """
        Set network Specification.

        ('Parameters')
        n_mid_units : hidden layers number.
        n_out : output number.
        """
        super(linear5_net, self).__init__()
        with self.init_scope():
            self.l1 = L.Linear(None, n_mid_units)
            self.l2 = L.Linear(n_mid_units, n_mid_units)
            self.l3 = L.Linear(n_mid_units, n_mid_units)
            self.l4 = L.Linear(n_mid_units, n_mid_units)
            self.l5 = L.Linear(n_mid_units, n_out)
        self.dropout_ratio = ratio

    def __call__(self, x):
        """Network foward calculation."""
        h1 = F.relu(self.l1(x))
        h2 = F.dropout(F.relu(self.l2(h1)), ratio=self.dropout_ratio)
        h3 = F.relu(self.l3(h2))
        h4 = F.dropout(F.relu(self.l4(h3)), ratio=self.dropout_ratio)
        return self.l5(h4)
