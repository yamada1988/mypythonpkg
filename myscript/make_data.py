N = 0
slvfe = []
for j in range(50, 49, -1):
    inpf = "MELT_{0:03d}/ERmods/loguv.tt".format(j)
    with open(inpf, 'rt') as f:
        aveuv = [line.split()[1] for line in f if not line.startswith('#')]
    inpf = "MELT_{0:03d}/ERmods/log.tt".format(j)
    try:
        with open(inpf, 'rt') as f:
            slvfe0 = [line[9:].split() for line in f if not line.startswith('#')]
            for line in slvfe0:
                for slvfe_ in line:
                    slvfe.append(slvfe_)
    except:
        continue 
   
    inpf = 'MELT_{0:03d}/ERmods/histogram.tt'.format(j)
    hists_k = []
    for k in range(1, 52+1):
        hists = []
        for i in range(1, 20+1):
            fname="MELT_{0:03d}/ERmod_{1:02d}/soln/engsln.{2:02d}".format(j, k, i)
            try:
                with open(fname, 'rt') as f:
                    histo = [float(line.split()[2]) for line in f]   
                    print("i:", i)
                    hist_ = sum(histo[447:495])
                    hists.append(hist_)
            except:
                continue

        for h in hists:
            hists_k.append(h)
    print(len(aveuv))
    print(len(slvfe))
    print(len(hists_k))

    pairs = []
    for i in range(len(slvfe)):
        print('i:{0:d}'.format(i))
        pair = (aveuv[i], slvfe[i], hists_k[i])
        pairs.append(pair)
        print(pair)
#print(pairs)
    outf = "MELT_{0:03d}/ERmods/pairs.tt".format(j)
    with open(outf, 'wt') as f:
        f.write('# aveuv\tslvfe\thists\n')
        for i in range(len(slvfe)):
            f.write('{0:5.3f}\t{1:5.3f}\t{2:8.6f}\n'.format(float(pairs[i][0]), float(pairs[i][1]), float(pairs[i][2])))
