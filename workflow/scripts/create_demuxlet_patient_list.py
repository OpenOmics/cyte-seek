f = open('demuxlet.csv')
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
