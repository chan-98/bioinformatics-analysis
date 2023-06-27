#!/usr/bin/env python

"""
This is a code to return the list of all the most frequent patterns of a given length that appear in a given sequence
Usage: python3 freq_words.py [length of pattern] [sequence text]
"""
import sys

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for kmer in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[kmer] == m:
            words.append(kmer)
            words.sort()
    return words
  
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if (bool(freq.get(Pattern))):
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1;
    return freq

k = sys.argv[1]
sequence = sys.argv[2]

print(FrequentWords(sequence, k))
