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
    curr_phase = ('0', '1')
    j = 0

    for hap in ret:
      hap[i] = curr_phase
      j += 1

      if j == chunk_size:
        curr_phase = curr_phase[::-1]
        j = 0

    chunk_size /= 2

  return ret
