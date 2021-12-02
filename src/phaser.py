import dna_sim
import haplotypes
from gibbs import Gibbs

if __name__ == "__main__":
  # Generate some haplotypes, get genotypes and haplotypes
  n = 50
  l = 50000
  genotypes, true_haplotypes = dna_sim.gen_haplotypes(n, l)
  print(f'Phasing {n} individuals over {len(genotypes[0])} SNPs...')

  # Run our statistical phaser on the genotypes
  iters = 150
  GS = Gibbs(genotypes)
  GS.run(iters)
  phased_haplotypes = GS.get_haplotypes()

  # Compare our phased haplotypes to the true haplotypes.
  # Also get switch error for random haplotypes and true haplotypes
  # as sanity checks
  phased_sw_errs = []
  random_sw_errs = []
  true_sw_errs   = []

  for true, phased in zip(true_haplotypes, phased_haplotypes):
    phased_sw_errs.append(haplotypes.get_switch_err(true, phased))
  for true, rand in zip(true_haplotypes, GS.get_random_haplotypes()):
    random_sw_errs.append(haplotypes.get_switch_err(true, rand))
  for true, other_true in zip(true_haplotypes, true_haplotypes):
    true_sw_errs.append(haplotypes.get_switch_err(true, other_true))

  print('***********')
  print('* RESULTS *')
  print('***********\n')

  print('True haplotype SWE:')
  print(f'{sum(true_sw_errs) / len(true_sw_errs):.3f}%')
  print('Random haplotype SWE:')
  print(f'{sum(random_sw_errs) / len(random_sw_errs):.3f}%')
  print('Phased haplotype SWE:')
  print(f'{sum(phased_sw_errs) / len(phased_sw_errs):.3f}%')
