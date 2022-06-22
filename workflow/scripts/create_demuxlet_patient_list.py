#!/usr/bin/env python3

import argparse, sys


def main(raw_args=None):
    parser = argparse.ArgumentParser(description="Create patients file for each sample given input csv")
    parser.add_argument('file_name', metavar='demuxlet.csv', action = "store",
        type=str, help="CSV file containing the sample name and individual associated patients")

    args = parser.parse_args(raw_args)


    f = open(args.file_name)
    samples = next(f).strip().split(',')
    metadata = { i: set() for i in samples}
    for line in f:
      tabs = line.strip().split(',')
      for i in range(len(samples)):
        if tabs[i] != '':
          metadata[samples[i]].add(tabs[i])

    f.close()

    for i in samples:
      f = open(f'demuxlet/output/{i}/patient_list', 'w')
      patients = list(metadata[i])
      patients.sort()
      f.write('\n'.join(patients))
      f.close()



if __name__ == '__main__':
    main()
