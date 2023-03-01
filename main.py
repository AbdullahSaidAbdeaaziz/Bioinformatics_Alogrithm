

def product(*args: tuple, repeat=1):
    """
    :param args: tuple from start to end
    :param repeat: no of repeated for
    :return: list of all possible motifs from start to end of args
    """

    pools = [tuple(pool) for pool in args] * repeat  # [[0...n-k+1], [0..n-k+1],....t_list]]
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool]
    for prod in result:
        yield prod


def cal_score(motif: list[int], DNA: list[str], k_mers: int):
    """
    :param motif: Possible Motif
    :param DNA: String of DNA
    :param k_mers: K-string to find on DNA
    :return: Calculate of Score
    """
    consensus = ''
    score = 0
    aline_mat = [DNA[j][i:i + k_mers] for j, i in zip(range(len(DNA)), motif)]

    print(aline_mat, motif)

    for i in range(len(aline_mat[0])):
        counts = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }
        for j in range(len(aline_mat)):
            counts[aline_mat[j][i]] += 1
        score += max(counts.values())
        consensus += max(counts, key=counts.get)
    return score, consensus


def BruteForceMotifSearch(DNA: list[str], t_row: int, n_col: int, k_motif: int):
    """
    :param DNA: Strings for DNA
    :param t_row: no.of row of DNA
    :param n_col: length of DNA
    :param k_motif: k-string to find
    :return: best consensus string of DNA
    """
    best_consensus = None
    best_score = 0
    best_motif = None
    prob_s = product(range(n_col - k_motif + 1), repeat=t_row)
    for i in prob_s:
        motif = i
        score, consensus = cal_score(motif, DNA, k_motif)
        if score > best_score:
            best_score = score
            best_motif = motif
            best_consensus = consensus

    return best_motif, best_score, best_consensus


def main():
    DNA = ['GTACAG', 'TGACCG', 'AAAACG', 'CTAGCG']
    t = len(DNA)
    n = len(DNA[0])
    k = 4
    best_motif, best_score, consensus = BruteForceMotifSearch(DNA, t, n, k)

    print(f"{'#' * 40}")
    print(f" Best motif : {best_motif} \n"
          f" Consensus motif: {consensus} ")
    print(f" With best score: {best_score}")
    print(f"{'#' * 40}")


if __name__ == '__main__':
    main()
