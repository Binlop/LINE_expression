import pytest
import subprocess


def check_KIF7():
    command = ['samtools', 'view', '../../alter_LINE_expr/ERR4231012_sorted.bam', 'chr15:89630789-89630849']
    res = subprocess.run(command, check=True, stdout=subprocess.PIPE).stdout.decode()
    reads = res.split('\n')
    only_seqs = [read.split('\t')[9] for read in reads if read != '']
    # print(only_seqs)
    seed = 'GTCTGCTGGGTAAGCTTCTACAgatgggga'
    x = 0
    # print(len(only_seqs))
    for seq in only_seqs:
        # print(seq)
        # print(seed)
        # if blast_read_and_seeding(seq, seed) is True:
        #     x += 1
        if blast_read_and_seeding(seq, seed) is True:
            x += 1
    print(x)
    # get_local_alignment_and_percent_of_matches(seq, seed)

def max_matches(string, substring):
    max_matches = 0
    for i in range(len(string) - len(substring) + 1):
        matches = sum(1 for j in range(len(substring)) if string[i + j].lower() == substring[j].lower())
        max_matches = max(max_matches, matches)
    
    return True if max_matches/30 >= 0.7 else False

def blast_read_and_seeding(seq: str, subseq: str):
    """Сравнивает попарно рид и затравку LINE со смещением каждый раз на 1 нуклеотид к концу рида, определяет наибольший % совпадений рида с затравкой"""
    matches = []
    for i in range(len(seq) - len(subseq) + 1):
        matches.append(count_mathes(seq[i:i+(len(subseq))], subseq))

    max_match = max(matches)
    return True if max_match/30 >= 0.7 else False


def count_mathes(seq1: str, seq2: str) -> int:
    """Возвращает отношение совпадений к общей длине последовательностей"""
    m = 0
    for i, nuc in enumerate(seq1):
        if nuc.lower() == seq2[i].lower():
            m += 1

    return m

if __name__ == '__main__':
    check_KIF7()

