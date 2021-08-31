#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys

from . import sp800_22_monobit_test

def read_bits_from_file(filename,bigendian):
    bitlist = list()
    if filename == None:
        f = sys.stdin
    else:
        f = open(filename, "rb")
    while True:
        bytes = f.read(16384)
        if bytes:
            for bytech in bytes:
                if sys.version_info > (3,0):
                    byte = bytech
                else:
                    byte = ord(bytech) 
                for i in range(8):
                    if bigendian:
                        bit = (byte & 0x80) >> 7
                        byte = byte << 1
                    else:
                        bit = (byte >> i) & 1
                    bitlist.append(bit)    
        else:
            break
    f.close()
    return bitlist




# X 3.1  Frequency (Monobits) Test
# X 3.2  Frequency Test within a Block
# X 3.3  Runs Test
# X 3.4  Test for the Longest Run of Ones in a Block
# X 3.5  Binary Matrix Rank Test
# X 3.6  Discrete Fourier Transform (Specral) Test
# X 3.7  Non-Overlapping Template Matching Test
# X 3.8  Overlapping Template Matching Test
# X 3.9  Maurers Universal Statistical Test
# X 3.10 Linear Complexity Test
# X 3.11 Serial Test
# X 3.12 Approximate Entropy Test
# X 3.13 Cumulative Sums Test
# X 3.14 Random Excursions Test
# X 3.15 Random Excursions Variant Test 


testlist = [
        'monobit_test',
        'frequency_within_block_test',
        'runs_test',
        'longest_run_ones_in_a_block_test',
        # 'binary_matrix_rank_test',
        'dft_test',
        'non_overlapping_template_matching_test',
        # 'overlapping_template_matching_test',
        # 'maurers_universal_test',
        # 'linear_complexity_test',
        'serial_test',
        'approximate_entropy_test',
        'cumulative_sums_test',
        'random_excursion_test',
        'random_excursion_variant_test',
        ]

print("Tests of Distinguishability from Random")

def runTestRoutines(filename, bigendian = False):
    print("[sp800_22] runTestRoutines()")
    bits = read_bits_from_file(filename, bigendian)    
    gotresult=False

    results = list()
    
    num_tests_passed = 0
    total_tests = 0
    for testname in testlist:
        print("[sp800_22] TEST: %s" % testname)
        # print(sys.path)
        import sp800_22_monobit_test

        m = __import__ ("sp800_22_"+testname)
        func = getattr(m,testname)
        
        (success,p,plist) = func(bits)

        summary_name = testname
        if success:
            print("[sp800_22]   PASS")
            summary_result = "PASS"
            num_tests_passed += 1
        else:
            print("[sp800_22]   FAIL")
            summary_result = "FAIL"
        
        if p != None:
            print("[sp800_22]   P="+str(p))
            summary_p = str(p)
            
        if plist != None:
            for pval in plist:
                print("[sp800_22] P="+str(pval))
            summary_p = str(min(plist))
        
        results.append((summary_name,summary_p, summary_result))
        total_tests +=1

    success_rate = num_tests_passed / total_tests * 100.
    print()
    print("[sp800_22] SUMMARY")
    print("[sp800_22] -------")
    
    for result in results:
        (summary_name,summary_p, summary_result) = result
        print("[sp800_22] ", summary_name.ljust(40),summary_p.ljust(18),summary_result)
    print("[sp800_22] -------")
    print("[sp800_22] PASSED {:}/{:} Tests. {:3.2f}%".format(num_tests_passed, total_tests, success_rate))
    return num_tests_passed, total_tests, success_rate, results
