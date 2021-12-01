def read_genotypes(filepath):
  genotypes = []

  with open(filepath, 'r') as in_f:
    for line in in_f:
      snps = line.split()
      genotype = []

      for snp in snps:
        alleles = snp.split('/')
        genotype.append((alleles[0], alleles[1]))

      genotypes.append(genotype)

  return genotypes
