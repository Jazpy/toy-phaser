import dna_sim
import haplotypes
import multiprocessing
from gibbs import Gibbs

def run_replicate(packed):
  n, l, iters = packed
  # Generate some haplotypes, get genotypes and haplotypes
  genotypes, true_haplotypes = dna_sim.gen_haplotypes(n, l)
  print(f'Phasing {n} individuals over {len(genotypes[0])} SNPs...')

  # Run our statistical phaser on the genotypes
  GS = Gibbs(genotypes)
  GS.run(iters)
  phased_haplotypes = GS.get_haplotypes()

  # Compare our phased haplotypes to the true haplotypes.
  # Also get switch error for random haplotypes
  phased_sw_errs = []
  random_sw_errs = []

  for true, rand in zip(true_haplotypes, GS.get_random_haplotypes()):
    random_sw_errs.append(haplotypes.get_switch_err(true, rand))
  for true, phased in zip(true_haplotypes, phased_haplotypes):
    phased_sw_errs.append(haplotypes.get_switch_err(true, phased))

  r_avg = sum(random_sw_errs) / len(random_sw_errs)
  p_avg = sum(phased_sw_errs) / len(phased_sw_errs)
  return [r_avg, p_avg]

def main():
  replicates = 10
  l = 4000
  iters = 150
  results = []

  # Repeat tests for n = 5, 10, 15, ..., 50
  for n in range(5, 51, 5):
    r_avg = 0.0
    p_avg = 0.0

    # Run 10 times to get an average
    pool = multiprocessing.Pool(replicates)
    t_results = pool.map(run_replicate, [(n, l, iters)] * replicates)

    for t in t_results:
      r_avg += t[0]
      p_avg += t[1]

    r_avg /= replicates
    p_avg /= replicates
    results.append([n, r_avg, p_avg])
    print(results)

  with open('results.txt', 'w') as out_f:
    for r in results:
      out_f.write(' '.join([str(x) for x in r]) + '\n')

  print('All done!')

if __name__ == "__main__":
  main()
