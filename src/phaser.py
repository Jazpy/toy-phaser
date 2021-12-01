import sys
import haplo_io
from gibbs import Gibbs

if __name__ == "__main__":
  genotypes = haplo_io.read_genotypes(sys.argv[1])
  GS = Gibbs(genotypes)
  GS.run()
