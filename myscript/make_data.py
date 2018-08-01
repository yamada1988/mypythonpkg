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
    hists_0 = []
    hists_1 = []
    for k in range(1, 52+1):
        hists0 = []
        hists1 = []
        for i in range(1, 20+1):
            fname="MELT_{0:03d}/ERmod_{1:02d}/soln/engsln.{2:02d}".format(j, k, i)
            try:
                with open(fname, 'rt') as f:
                    histo = [float(line.split()[2]) for line in f]   
                    print("i:", i)
                    hist_0 = sum(histo[447:495])
                    hist_1 = sum(histo[:446])
                    hists0.append(hist_0)
                    hists1.append(hist_1)
            except:
                continue

        for h0 in hists0:
            hists_0.append(h0)
        for h1 in hists1:
            hists_1.append(h1)
    print(len(aveuv))
    print(len(slvfe))
    print(len(hists_0))
    print(len(hists_1))

    pairs = []
    for i in range(len(slvfe)):
        print('i:{0:d}'.format(i))
        pair = (aveuv[i], slvfe[i], hists_0[i], hists_1[i])
        pairs.append(pair)
        print(pair)
#print(pairs)
    outf = "MELT_{0:03d}/ERmods/pairs.tt".format(j)
    with open(outf, 'wt') as f:
        f.write('# aveuv\tslvfe\thist0\thist1\n')
        for i in range(len(slvfe)):
            f.write('{0:5.3f}\t{1:5.3f}\t{2:7.5f}\t{3:7.5f}\n'.format(float(pairs[i][0]), float(pairs[i][1]), float(pairs[i][2]), float(pairs[i][3])))
