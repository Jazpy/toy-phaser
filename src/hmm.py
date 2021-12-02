import numpy as np

class HMM:
  def __init__(self, ro, mu, H_i):
    self.ro = ro
    self.mu = mu
    self.s  = len(H_i)
    self.transitions = np.zeros((len(H_i), len(H_i)))
    self.emissions = np.array([[1 - self.mu, self.mu],
                               [self.mu, 1 - self.mu]])
    self.pi = np.full(len(H_i), 1 / len(H_i))
    self.H_i = np.array(H_i)

    if len(H_i) > 1:
      ind_recomb_rate = self.ro / (len(H_i) - 1)
    else:
      ind_recomb_rate = 0

    for j in range(len(H_i)):
      for k in range(len(H_i)):
        if j == k:
          self.transitions[j][k] = 1 - ind_recomb_rate
        else:
          self.transitions[j][k] = ind_recomb_rate

  def forward(self, h0, h1):
    hs = [h0, h1]
    n  = len(h0)
    ret = []

    for h in hs:
      fwd = np.zeros((self.s, n))
      f_col_em = [self.emissions[x, h[0]] for x in self.H_i[:, 0]]
      fwd[:, 0] = self.pi * f_col_em

      # Very similar to Viterbi in the sense that we update fwd[j, i]
      # with the probabilities of getting to state j at time i. However,
      # instead of taking the max from the previous time, we take the sum
      # of all transitions into state j from time i - 1
      for i in range(1, n):
        for j in range(self.s):
          fwd[j, i] = sum(fwd[:, i - 1] * self.transitions[:, j] *
                          self.emissions[self.H_i[j][i], h[i]])

      ret.append(sum(fwd[:, n - 1]))

    return ret
