import argparse
from dna.Recuit import recuit_main
from dna.Traditionnal import traditionnal_main
from dna.Genetic import algo_genetique as genetic_main
from dna.Genetic import stats

# Gestion des arguments -> voir le README.md
parser = argparse.ArgumentParser()
parser.add_argument(
    "-m", "--mode", nargs='?',
    help="Choose mode : 'traditional'[default] , 'recuit'[training] or 'genetic'[training]",
    default='traditional')
parser.add_argument(
    "-d", "--dna", nargs='?', help="input filename of DNA sequence",
    default='data/plasmid_8k.fasta')
parser.add_argument("-j", "--json", nargs='?',
                    help="input filename of JSON file", default='dna/table.json')
parser.add_argument("-i", "--max-iters", nargs='?',
                    help="max iterations for recuit mode", default=100, type=int)
parser.add_argument("-p", "--pop-size", nargs='?', default=10,
                    type=int, help="population size for genetic mode")
parser.add_argument("-s", "--stat", action='store_true',
                    help="best scores by population for genetic mode")
parser.parse_args()
args = parser.parse_args()


def main():
    """Fonction principale, qui redirige vers les fonctions de l'algorithme choisi"""
    if args.mode == "recuit":
        seqs = []
        for filename in ("data/plasmid_8k.fasta", "data/plasmid_180k.fasta"):
            lineList = [line.rstrip('\n') for line in open(filename)]
            seq = ''.join(lineList[1:])
            seqs.append(seq)
        recuit_main(seqs, args.json, args.max_iters)
    elif args.mode == "genetic":
        lineList = [line.rstrip('\n') for line in open(args.dna)]
        seq = ''.join(lineList[1:])
        if args.stat:
            stats(seq)
        else:
            genetic_main(seq, args.pop_size)
    elif args.mode == "traditional":
        lineList = [line.rstrip('\n') for line in open(args.dna)]
        seq = ''.join(lineList[1:])
        traditionnal_main(seq, args.dna, args.json)


if __name__ == "__main__":
    main()
