import random
import string
import collections
import math
import matplotlib.pyplot as plt
def s_1(N, N_sequence):
    one = ["1"] * N
    zero = ["0"] * (N_sequence - N)
    seq = one + zero
    random.shuffle(seq)
    os = ''.join(seq)
    return os

def s_2(N_sequence):
    surname = 'Великий'
    N1 = len(surname)
    ls = list(surname)
    lz = ['0'] * (N_sequence - N1)
    seq = ls + lz
    os_2 = ''.join(seq)
    return os_2

def s_3(N_sequence):
    surname = 'Великий'
    N1 = len(surname)
    ls = list(surname)
    lz = ['0'] * (N_sequence - N1)
    seq = ls + lz
    random.shuffle(seq)
    os_3 = ''.join(seq)
    return os_3

def s_4(N_sequence):
    surname = 'Великий'
    letters = list(surname + '529')
    let = len(letters)
    rep = N_sequence // let
    rem = N_sequence % let
    l = letters * rep
    l += letters[:rem]
    os_4 = ''.join(map(str, l))
    return os_4

def s_5():
    w = "Ве529"
    char = list(w)
    gs = char * 20
    random.shuffle(gs)
    os_5 = ''.join(gs)
    return os_5

def s_6(N_sequence):
    let, dig = list("Ве"), list("529")
    letters, digits = int(0.7 * N_sequence), int(0.3 * N_sequence)
    l_100 = []
    for i in range(letters):
        l_100.append(random.choice(let))
    for i in range(digits):
        l_100.append(random.choice(dig))
    random.shuffle(l_100)
    os_6 = ''.join(l_100)
    return os_6

def s_7(N_sequence):
    elements = string.ascii_lowercase + string.digits
    l_100 = [random.choice(elements) for _ in range(N_sequence)]
    os_7 = ''.join(l_100)
    return os_7

def s_8(N_sequence):
    os_8 = "1" * N_sequence
    return os_8
def save_sequence(os, seq_size, seq_alphabet, symbols, probability, uniformity, entropy, source):
        with open("results_sequence.txt", "a", encoding="utf8") as file:
            file.write(f"Послідовність: {os}\n")
            file.write(f"Розмір послідовності: {seq_size} байт\n")
            file.write(f"Розмір алфавіту: {seq_alphabet}\n")
            file.write(f"Ймовірність появи символів: {symbols}\n")
            file.write(f"Середнє арифметичне ймовірностей: {probability}\n")
            file.write(f"Ймовірність розподілу символів: {uniformity}\n")
            file.write(f"Ентропія: {entropy}\n")
            file.write(f"Надмірність джерела: {source}\n" + "-" * 100 + "\n")
def sequence_txt(seq):
        with open("sequence.txt", "a", encoding="utf8") as file:
            file.write(f"{seq}\n")

if __name__ == "__main__":
    N = 5
    N_sequence = 100
    oss = [s_1(N, N_sequence), s_2(N_sequence), s_3(N_sequence),s_4(N_sequence), s_5(), s_6(N_sequence),s_7(N_sequence), s_8(N_sequence)]
    results = []
    for seq in oss:
        counts = collections.Counter(seq)
        pro = {symbol: count / N_sequence for symbol, count in counts.items()}
        probability = sum(pro.values()) / len(pro)
        equal = all(abs(prob - probability) < 0.05 * probability for prob in pro.values())
        if equal:
            uniformity = "рівна"
        else:
            uniformity = "нерівна"
        entropy = -sum(p * math.log2(p) for p in pro.values())
        if len(set(seq)) > 1:
            source= 1 - entropy / math.log2(len(set(seq)))
        else:
            source = 1
        probability_str = ', '.join([f"{symbol}={prob:.4f}" for symbol, prob in pro.items()])
        sequence_txt(seq)
        save_sequence(seq, len(seq), len(set(seq)), probability_str, probability, uniformity, entropy, source)
        results.append([len(set(seq)), round(entropy, 2), round(source, 2), uniformity])
    fig, ax = plt.subplots(figsize=(14 / 1.54, 8 / 1.54))
    headers = ['Розмір алфавіту', 'Ентропія', 'Надмірність', 'Ймовірність']
    row = [f'Послідовність {i}' for i in range(1, 9)]
    ax.axis('off')
    table = ax.table(cellText=results, colLabels=headers, rowLabels=row, loc='center', cellLoc='center')
    table.set_fontsize(14)
    table.scale(0.8, 2)
    fig.savefig("Характеристики сформованих послідовностей.png", dpi=600)