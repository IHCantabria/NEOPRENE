'''
Library to perform Particle Swarm Optimization.

	Authors: 
        + Manuel del Jesus
    
'''

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
        self._position = np.maximum(np.minimum(self._position, self._parent._limits[:,1]),self._parent._limits[:,0])
        #print "++ ", swarmBest, self._best_position, self._position
        self._velocity = self.OMEGA * self._velocity \
            + self.PHIW * np.random.rand(self._parent._D) * (self._best_position - self._position) \
            + self.PHIS * np.random.rand(self._parent._D) * (swarmBest - self._position)
        #self._velocity /= np.sqrt(np.sum(self._velocity**2))
        #print np.sqrt(np.sum(self._velocity**2))
        #print self._position, swarmBest
        self.fvalue = self._parent._Function(self._position)
        self.testValue()

    def testValue(self):
        raise NotImplementedError("You are calling the parent class.")


class Worker_min(Worker):
    def testValue(self):
        if self.fvalue < self._best:
            #In case the current value is better than the old one
            self._best = self.fvalue
            self._best_position = self._position.copy()
        self.fvalue = None

class Worker_max(Worker):
    def testValue(self):
        if self.fvalue > self._best:
            #In case the current value is better than the old one
            self._best = self.fvalue
            self._best_position = self._position.copy()
        self.fvalue = None

class Swarm(object):
    def __init__(self, D, limits, n, F, optimum, verbose=False):
        self._D = D
        self._limits = limits
        self._inf = limits[:,0]
        self._inc = limits[:,1] - limits[:,0]
        self._Function = F
        if optimum=='min':
            self._workers = [Worker_min(self) for i in range(n)]
        elif optimum=='max':
            self._workers = [Worker_max(self) for i in range(n)]
        else:
            raise SystemError("Optimum not known...\toptimum = %s"%optimum)
        self._bestList = [i._best for i in self._workers]
        if optimum=='min':
            self._best = np.min(self._bestList)
            self._best_position = self._workers[np.argmin(self._bestList)]._best_position.copy()
        elif optimum=='max':
            self._best = np.max(self._bestList)
            self._best_position = self._workers[np.argmax(self._bestList)]._best_position.copy()
        if verbose:
            print("Initial best position ", self._best, self._best_position)

    def getBest(self):
        #print("Best value : %6.5e"%self._best, end=" ")
        #print("Position : ", end=" ")
        #print("\t", self._best_position)
        return (self._best, self._best_position)

    def setBest(self, posicion):
        best = self._Function(posicion)
        try:
            if best < self._best:
                self._best = best
                self._best_position = posicion
        except:
            self._best = best
            self._best_position = posicion

    def step(self):
        [i.update(self._best_position) for i in self._workers]
        self._bestList = [i._best for i in self._workers]
        tentative = np.min(self._bestList)
        if tentative < self._best:
            self._best = tentative
            self._best_position = self._workers[np.argmin(self._bestList)]._best_position.copy()

def testFunction(x):
    return np.sum((x-0.83)**2)
    #return (1-x[0])**2 + 100.0*(x[1]-x[0]**2)**2
    #return np.sum(100.0*(x[1:] - x[:-1]**2)**2 + (x[:-1] - 1.0)**2)

if __name__ == "__main__":
    a = 0.8
    b = 1.2
    S = Swarm(6, np.asarray([[a,b],[a,b],[a,b],[a,b],[a,b],[a,b]]), 100, testFunction, 'min')

    with open("results_prueba.txt",'wt', 0) as fid:
        for i in xrange(100):
            for j in xrange(100):
                S.step()
            (b, bp) = S.getBest()
            fid.write('%6.5e\t%s\n'%(b, str(bp)))


