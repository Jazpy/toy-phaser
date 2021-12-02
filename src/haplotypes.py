def __segsites(g):
  ret = []
  first = True

  for i, snp in enumerate(g):
    if snp[0] != snp[1]:
      if first:
        first = False
      else:
        ret.append(i)

  return ret

def haplo_space(g):
  segsites = __segsites(g)
  ret = [g.copy() for _ in range(2 ** len(segsites))]

  chunk_size = 2 ** (len(segsites) - 1)

  for i in segsites:
    curr_phase = [0, 1]
    j = 0

    for hap in ret:
      hap[i] = curr_phase
      j += 1

      if j == chunk_size:
        curr_phase = curr_phase[::-1]
        j = 0

    chunk_size /= 2

  return ret

def __switch_aux(true, phased, flip):
  het_sites = 0
  misses    = 0

  for t, p in zip(true, phased):
    # Skip homozygous sites
    if t[0] == t[1]:
      continue
    else:
      het_sites += 1

    t0 = t[0]
    p0 = p[0] if not flip else p[1]

    if p0 != t0:
      misses += 1
      flip = not flip

  if not het_sites:
    sw_err = 0.0
  else:
    sw_err = misses / het_sites * 100.0

  return sw_err

def get_switch_err(true, phased):
  return min(__switch_aux(true, phased, False),
             __switch_aux(true, phased, True))
