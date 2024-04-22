f = open('all_nominal_outs.txt')
lines = f.readlines()
cleanlines = [x.split('_upout')[0] for x in lines]
samps = list(set(cleanlines))
samps = sorted(samps)

g = open('list_all_background_mc_samples.txt','a')

for line in samps:
    print(line)
    g.write(line+'\n')
    
g.close()
f.close()
