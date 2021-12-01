import numpy as np

class HMM:
  def __init__(self, ro, mu):
    self.ro = ro
    self.mu = mu
    self.states

def forward(P, C, X, s, pi):
  fwd = np.zeros((s, len(X)))
  fwd[:, 0] = pi * C[:, int(X[0]) - 1]

  # Very similar to Viterbi in the sense that we update fwd[j, i]
  # with the probabilities of getting to state j at time i. However,
  # instead of taking the max from the previous time, we take the sum
  # of all transitions into state j from time i - 1
  for i in range(1, len(X)):
    for j in range(s):
      fwd[j, i] = sum(fwd[:, i - 1] * P[:, j] * C[j, int(X[i]) - 1])

  return fwd

def main():
  s = 3
  n = 10
  C = np.array([[.9, .1],
                [.2, .8],
                [.5, .5]])
  P = np.array([[.2, .7, .1],
                [.4, 0.0, .6],
                [0.0, 1.0, 0.0]])
  pi = [0.3335, 0.3333, 0.3332]

  x = '2111212212'

  # Get probability for this X
  print('Forward says P(x_1:10) = ')
  print(sum(forward(P, C, x, s, pi)[:, len(x) - 1]))
