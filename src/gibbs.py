import haplotypes
import random

from hmm import HMM

class Gibbs:
  def __init__(self, genotypes):
    self.G = genotypes
    self.n = len(genotypes)
    self.S = [haplotypes.haplo_space(g) for g in self.G]
    self.H = [random.choice(s) for s in self.S]
    self.E = [[] for _ in range(self.n)]
    self.p = []

  def run(self, iters, burn=100):
    # Gibbs sampling depends on some burn-in iterations and then
    # averaging over the next iterations
    for curr_iter in range(burn + iters):
      # Some basic logging
      if curr_iter % 50 == 0:
        if curr_iter < burn:
          print(f'Burn-in iteration {curr_iter} / {burn}...')
        else:
          print(f'Main iteration {curr_iter - burn} / {iters}...')

      new_H = []

      # Sample from P(h0, h1 | H_i) for all i
      for i, h in enumerate(self.H):
        # Construct H_i, where we exclude the i-th haplo pair
        H_i = []
        for j, h_j in enumerate(self.H):
          if j == i:
            continue

          h_j_0 = [x[0] for x in h_j]
          h_j_1 = [x[1] for x in h_j]
          H_i.append(h_j_0)
          H_i.append(h_j_1)

        # Construct HMM for this H_i
        model = HMM(2e-8, 2e-8, H_i)

        # Get P(h0 | H_i) and P(h1 | H_i) for each h0, h1 possible
        weights = []
        for pair in self.S[i]:
          h0 = [x[0] for x in pair]
          h1 = [x[1] for x in pair]
          # Calculate forward algorithm with h0, h1 observations
          h_p = model.forward(h0, h1)

          # Get P(H = h0, h1 | H_i)
          d = 1 if h0 == h1 else 0
          weights.append((2 - d) * h_p[0] * h_p[1])

        # Update haplotype space based on sampling from weights
        update_idx = random.choices(list(range(len(self.S[i]))),
                                    weights=weights, k=1)[0]
        self.H[i] = self.S[i][update_idx]

        # If we're past burn-in iters, add update_idx to update list
        if curr_iter >= burn:
          self.E[i].append(update_idx)

    # Get expected value of i-th haplotype
    for i, indices in enumerate(self.E):
      expected = round(sum(indices) / len(indices))
      self.H[i] = self.S[i][expected]

  def get_haplotypes(self):
    return self.H

  def get_random_haplotypes(self):
    return [random.choice(s) for s in self.S]

  def print(self):
    print('* Current haplotypes *')
    for h in self.H:
      print(' '.join([f'{x[0]}|{x[1]}' for x in h]))
