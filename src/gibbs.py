import haplotypes
import random

class Gibbs:
  def __init__(self, genotypes):
    self.G = genotypes
    self.n = len(genotypes)
    self.S = [haplotypes.haplo_space(g) for g in self.G]
    self.H = [[random.choice(s) for s in self.S]]
    self.p = []

  def run(self):
    
