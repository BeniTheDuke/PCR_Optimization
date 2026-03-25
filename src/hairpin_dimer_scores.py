def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp[b] for b in reversed(seq))

def dimer_score(p1, p2, min_run=3, end_zone=5):
    rc2 = reverse_complement(p2)
    n1, n2 = len(p1), len(rc2)
    score = 0

    # slide p1 and rc(p2) past each other in every alignment
    for offset in range(-n2 + 1, n1):
        run = 0
        for i in range(max(0, offset), min(n1, offset + n2)):
            j = i - offset
            if p1[i] == rc2[j]:
                run += 1
                # only count if the match involves the 3' end of either primer
                if i >= n1 - end_zone or j >= n2 - end_zone:
                    score = max(score, run)
            else:
                run = 0
    return score
