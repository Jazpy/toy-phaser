import msprime
import random

def gen_haplotypes(n, length):
  dem_model = msprime.Demography()
  dem_model.add_population(name='P0', initial_size=10_000)

  sample_list = [
    msprime.SampleSet(n, population='P0', time=0,  ploidy=2),
  ]

  trees = msprime.sim_ancestry(samples=sample_list, demography=dem_model,
                               sequence_length=length, recombination_rate=2e-8)

  mutations = msprime.sim_mutations(trees, rate=2e-8)

  haps = [[] for _ in range(n * 2)]
  for var in mutations.variants():
    for i, gt in enumerate(var.genotypes):
      haps[i].append(gt)

  gts = []
  hps = []
  for i in range(0, len(haps), 2):
    h0 = haps[i]
    h1 = haps[i + 1]
    gt = []
    hp = []

    for x, y in zip(h0, h1):
      if x != y:
        gt.append([0, 1])
      else:
        gt.append([x, y])
      hp.append([x, y])

    gts.append(gt)
    hps.append(hp)

  return gts, hps
