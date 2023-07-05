import subprocess
import pickle
import glob

normfiles = glob.glob('analysis_output_ZpAnomalon/alphaMethodNorms/Run2_161718*shiftednormsforuncs*dynormalization_alphatest_systnominal_kfnom_btagnom_muidnom_elidnom_elreconom_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.pkl')

normsToTest = {}
for f in normfiles:
    direc = f.split('shiftednormsforuncs_')[-1].split('_dynormalization')[0]
    direcstr = ''
    name = ''
    if int(direc) == -1:
        direcstr = "down"
    if int(direc) == 1:
        direcstr = "up"
    fo = open(f,'rb')
    nraw = pickle.load(fo)
    for key in nraw.keys():
        if 'nominal' in key:
            continue
        if "env" not in key:
            name = key+direcstr
        else:
            name = key
        normsToTest[name] = nraw[key]

print(normsToTest)

for key in normsToTest.keys():
    print("Running the alpha ratio fits with the normalization vaired ",key)
    subprocess.run(["python","doAlphaRatioFits_systematics.py","-m","75.0","-z","100.0","-j","300.0","-wp","0.8","-sn","nominal","-sd","nom","-normsyst",key,"-normsystval",str(normsToTest[key])])
