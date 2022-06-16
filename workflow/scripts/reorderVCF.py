import subprocess, re, argparse

def reorderVCF(vcf, bam, output = 'output.vcf'):
    bam_out = subprocess.Popen("module load samtools; samtools idxstats %s | awk '{print $1}' | grep chr" % bam, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    chromosomes = bam_out.communicate()[0].decode("utf-8").strip().split('\n')

    if 'chrM' in chromosomes:
      chromosomes.remove('chrM')

    f = open(vcf)

    g = open(output, 'w')
    line = next(f)
    while "##contig=<ID" not in line:
      g.write(line)
      line = next(f)

    chr_dict = dict()
    while "##contig=<ID" in line:
      match = re.search(r"(?<=##contig=<ID=)chr[\dXY]+,", line)
      if match != None:
        chr_dict[match[0].rstrip(',')] = ','.join(line.split('>')[0].split(',')[:2]) + '>'
      line = next(f)

    for chr in chromosomes:
      g.write(chr_dict[chr] + '\n')

    g.write(line)

    for line in f:
      if line[0] == '#':
        if '#CHROM' in line:
          break
        else:
          g.write(line)

    g.write(line)

    chr_contents = {chr: list() for chr in chromosomes}
    for line in f:
      if line.split('\t')[0] in chromosomes:
        chr_contents[line.split('\t')[0]].append(line.strip())

    g.write('\n'.join(['\n'.join(chr_contents[chr]) for chr in chromosomes]))
    g.write('\n')
    
    f.close()
    g.close()



def main(raw_args=None):
    global args
    parser = argparse.ArgumentParser(description="""Rename, reorder, and filter chromosomes in vcf file based on provided bam file""")
    parser.add_argument("-v", "--vcf", metavar="orig.vcf", dest="vcf",
                        action="store", type=str,
                        help="Original vcf file")
    parser.add_argument("-b", "--bam", metavar="sample.bam", dest="bam",
                        action="store", type=str,
                        help="Sample bam file to order by")
    parser.add_argument("-o", "--output", metavar="output.vcf", dest="output",
                        action="store", type=str, required=False, default="output.vcf",
                        help="Name of output vcf file")

    args = parser.parse_args(raw_args)

    reorderVCF(args.vcf, args.bam, args.output)

if __name__=='__main__':
    main()
