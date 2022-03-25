import glob, sys

def main(path):
    files = glob.glob('*/outs/molecule_info.h5')
    aggrFile = open('AggregatedDatasets.csv', 'w')

    aggrFile.write('sample_id,molecule_h5\n')
    files.sort()
    for i in files:
        aggrFile.write('%s,%s/%s\n' % (i.split('/')[0], path, i))
    aggrFile.close()

if __name__ == "__main__":
    main(sys.argv[1])
