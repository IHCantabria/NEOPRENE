# PSOMdJ.py

"""
Library to perform Particle Swarm Optimization.

Authors: 
    + Manuel del Jesus
"""

import numpy as np

class Worker(object):
    __slots__ = ["_position", "_velocity", "_best", "_best_position", "_parent", "fvalue"]
    OMEGA = 0.9
    PHIW = 0.9
    PHIS = 0.9
    
    def __init__(self, parent):
        self._position = np.random.rand(parent._D) * parent._inc + parent._inf
        self._velocity = np.random.rand(parent._D) * parent._inc
        self._best = parent._Function(self._position)
        self._best_position = self._position.copy()
        self._parent = parent
        self.fvalue = None

    def update(self, swarmBest):
        self._position += self._velocity
        self._position = np.maximum(np.minimum(self._position, self._parent._limits[:,1]), self._parent._limits[:,0])
        self._velocity = self.OMEGA * self._velocity \
            + self.PHIW * np.random.rand(self._parent._D) * (self._best_position - self._position) \
            + self.PHIS * np.random.rand(self._parent._D) * (swarmBest - self._position)
        self.fvalue = self._parent._Function(self._position)
        self.testValue()

    def testValue(self):
        raise NotImplementedError("You are calling the parent class.")

class Worker_min(Worker):
    def testValue(self):
        if self.fvalue < self._best:
            self._best = self.fvalue
            self._best_position = self._position.copy()
        self.fvalue = None

class Worker_max(Worker):
    def testValue(self):
        if self.fvalue > self._best:
            self._best = self.fvalue
            self._best_position = self._position.copy()
        self.fvalue = None

class Swarm(object):
    def __init__(self, D, limits, n, F, optimum, verbose=False):
        self._D = D
        self._limits = limits
        self._inf = limits[:, 0]
        self._inc = limits[:, 1] - limits[:, 0]
        self._Function = F
        self._workers = [Worker_min(self) if optimum == 'min' else Worker_max(self) for _ in range(n)]
        self._bestList = [worker._best for worker in self._workers]
        self._best = np.min(self._bestList) if optimum == 'min' else np.max(self._bestList)
        self._best_position = self._workers[np.argmin(self._bestList)]._best_position.copy() if optimum == 'min' else self._workers[np.argmax(self._bestList)]._best_position.copy()
        if verbose:
            print("Initial best position ", self._best, self._best_position)

    def getBest(self):
        return (self._best, self._best_position)

    def setBest(self, position):
        best = self._Function(position)
        if best < self._best:
            self._best = best
            self._best_position = position

    def step(self):
        for worker in self._workers:
            worker.update(self._best_position)
        self._bestList = [worker._best for worker in self._workers]
        tentative = np.min(self._bestList)
        if tentative < self._best:
            self._best = tentative
            self._best_position = self._workers[np.argmin(self._bestList)]._best_position.copy()

def testFunction(x):
    return np.sum((x - 0.83) ** 2)

if __name__ == "__main__":
    a, b = 0.8, 1.2
    S = Swarm(6, np.array([[a, b]] * 6), 100, testFunction, 'min', verbose=True)
    with open("results_prueba.txt", 'wt') as fid:
        for _ in range(100):
            for _ in range(100):
                S.step()
            b, bp = S.getBest()
            fid.write(f'{b:6.5e}\t{bp}\n')
