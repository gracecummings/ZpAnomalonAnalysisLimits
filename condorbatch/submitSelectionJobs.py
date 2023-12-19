import json
import sys
import os
import argparse
import subprocess
from datetime import date

parser = argparse.ArgumentParser()

toplevelsyst = [
    'jecup',
    'jecdwn',
    'jerup',
    'jerdwn',
    'unclup',
    'uncldwn'
]

if __name__=='__main__':
    parser.add_argument("-j","--sampleJson",help=".json file with the topiaries to do selections on")
    parser.add_argument("-s","--sample",type=str,help="specific sample")
    parser.add_argument("-c","--channel",type=str,help="channel: mumu,ee,emu")
    parser.add_argument("-syst","--syststr",nargs="*",help="the systematic flag you want to do. Can be specifc up/downs, or general, like 'jec' if you want both")
    args = parser.parse_args()

    #inputsraw = args.syststr.strip("[]").split(",")
    #insysts = list(map(lambda x:x.strip('""'),inputsraw))
    insysts = args.syststr

    #Check json
    if (args.sampleJson is None):
        print("No list of skims provided, please provide appropriate json")
        fsjson = {}
    else:
        if (not os.path.exists(args.sampleJson)):
            print("Invalid json")
        fstotop = open(args.sampleJson,'r')
        fsjson = json.load(fstotop)

    #Make log directories
    if not os.path.exists("condorMonitoringOutput/"+str(date.today())+"/"):
        os.makedirs("condorMonitoringOutput/"+str(date.today())+"/")

    syststr = ""
    if args.syststr:
        #syst = args.syststr.split(" ")[0]
        syststr = ''.join(insysts)
        print("Running with systematics ",str(syststr))
    else:
        print("Running with all systematics in json - could be none, could be all topiary level ones.")

    #Tar the working area
    print("Creating tarball of working area")
    tarballName = "selectionsForCondor.tar.gz"
    if not os.path.exists(tarballName):
        os.system("tar -hcf "+tarballName+" selection_jobs")
    else:
        print('FOUND A TARBALL -- USING, BE CAREFULL!!!!!')

    #Where do you want to save it
    eosForOutput = "root://cmseos.fnal.gov//store/user/lpcboostres/leptophobic_selections-"+syststr+"_"+str(date.today())
    eosOnlypath  = eosForOutput.split("root://cmseos.fnal.gov/")[-1]
    eosFinalDir  = eosOnlypath.split("/")[-1]

    #check eos for premade directories
    if os.system("eos root://cmseos.fnal.gov/ ls /store/user/lpcboostres/ | grep leptophobic_selections") == 0:#if there are eos directories
        eos_conflict_dirs = subprocess.check_output("eos root://cmseos.fnal.gov/ ls /store/user/lpcboostres/ | grep leptophobic_selections",shell=True).decode(sys.stdout.encoding).split()
        print("Desired EOS directory already exits.")
        if not any(eosFinalDir in path for path in eos_conflict_dirs): 
            print("Desired EOS output directory does not exist - creating it!")
            os.system("eos root://cmseos.fnal.gov mkdir {}".format(eosOnlypath))
    else:
        print("Desired EOS output directory does not exist - creating it!")
        os.system("eos root://cmseos.fnal.gov mkdir {}".format(eosOnlypath))


    #Get the samples that are in json
    gensets = list(set([f.split("_topiary")[0] for f in fsjson.keys()]))
    
    #Submit the jobs
    #Get list grouped by sample
    jobcnt = 0
    if len(fsjson.keys()) > 0:
        sampgrouped = {}
        for samp in gensets:
            if args.sample:
                if args.sample not in samp:
                    continue
                else:
                    print("Only submitting jobs to process ",samp)
            sampgrouped = [fsjson[f]["path"]+"/"+fsjson[f]["topiaries"][0] for f in fsjson.keys() if samp in f]
            topstoanalyze = []
            if insysts:
                for syst in insysts:
                    if any(syst in s for s in toplevelsyst):
                        #print("Adding files with {} to list to be selected".format(syst))
                        #print([f for f in sampgrouped if syst in f])
                        topstoanalyze += [f for f in sampgrouped if syst in f] 
                    else:
                        #print("systematics string not a toiary level syst, be careful about input .json")
                        topstoanalyze = sampgrouped
            else:
                topstoanalyze = sampgrouped

            
            sampleName = samp
            fullpaths = ["root://cmseos.fnal.gov/"+x for x in topstoanalyze]
            samplistasstr = str(fullpaths).replace(' ','')
            #print(samplistasstr)#how to pass to bash script

            
            print("Root files with selections from {0} are being written to {1}".format(sampleName,eosOnlypath))
            print("Running over {0} topiaries in this job.".format(len(topstoanalyze)))
            
            #Args to pass
            argu = "Arguments = {0} {1} {2}\n".format(eosForOutput,samplistasstr,args.channel)
            
            #Make the jdl for each sample
            #jdlName = "selections_"+sampleName+"_"+syststr+"_"+str(date.today())+".jdl"
            #jdl = open(jdlName,"w")
            #jdl.write("universe = vanilla\n")
            #jdl.write("Should_Transfer_Files = YES\n")
            #jdl.write("WhenToTransferOutput = ON_EXIT\n")
            #jdl.write("Transfer_Input_Files = "+tarballName+"\n")
            #jdl.write("Output = condorMonitoringOutput/{0}/{1}_out.stdout\n".format(str(date.today()),sampleName+"_"+syststr))
            #jdl.write("Error = condorMonitoringOutput/{0}/{1}_err.stder\n".format(str(date.today()),sampleName+"_"+syststr))
            #jdl.write("Log = condorMonitoringOutput/{0}/{1}_log.log\n".format(str(date.today()),sampleName+"_"+syststr))
            #jdl.write("Executable = selections.sh\n")
            #jdl.write(argu)
            #jdl.write('+DESIRED_Sites="T3_US_Baylor,T2_US_Caltech,T3_US_Colorado,T3_US_Cornell,T3_US_FIT,T3_US_FNALLPC,T3_US_Omaha,T3_US_JHU,T3_US_Kansas,T2_US_MIT,T3_US_NotreDame,T2_US_Nebraska,T3_US_NU,T3_US_OSU,T3_US_Princeton_ICSE,T2_US_Purdue,T3_US_Rice,T3_US_Rutgers,T3_US_MIT,T3_US_NERSC,T3_US_SDSC,T3_US_FIU,T3_US_FSU,T3_US_OSG,T3_US_TAMU,T3_US_TTU,T3_US_UCD,T3_US_UCSB,T2_US_UCSD,T3_US_UMD,T3_US_UMiss,T2_US_Vanderbilt,T2_US_Wisconsin"')
            #jdl.write("\n")
            #jdl.write("Queue 1\n")#Not sure about this one
            #jdl.close()
            
            #submit the jobs
            #os.system("condor_submit {0}".format(jdlName))

            jobcnt += 1
    print("Submitted {} jobs".format(jobcnt))
