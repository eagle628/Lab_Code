"""Define chainer Deep Network model."""
import chainer
import chainer.links as L
import chainer.functions as F


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

    def predict(self, x):
        """Network foward calculation."""
        h1 = F.relu(self.l1(x))
        h2 = F.relu(self.l2(h1))
        return self.l3(h2)


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
