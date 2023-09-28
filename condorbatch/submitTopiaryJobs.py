import json
import sys
import os
import argparse
from datetime import date

#def getSampleName(sampstr):
#    name = "badsamp"
#    if "WZ" in sampstr:
#        name = sampstr.split(")#

parser = argparse.ArgumentParser()

if __name__=='__main__':
    parser.add_argument("-j","--sampleJson",help=".json file with the skims to be made into topiary")
    parser.add_argument("-s","--sample",type=str,help="specific sample you want make a topiary of")
    parser.add_argument("-c","--channel",type=str,help="channel: mumu,ee,emu")
    parser.add_argument("-syst","--syststr",type=str,default="none",help="the systematic flag you want topiary to make")
    args = parser.parse_args()

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
    if "none" not in args.syststr:
        #syst = args.syststr.split(" ")[0]
        syststr = args.syststr
        print("Running with systematics ",syststr)
    else:
        print("No topiary-level systematics being applied.")

    #Tar the working area
    print("Creating tarball of working area")
    tarballName = "topiaryForCondor.tar.gz"
    os.system("tar -hcf "+tarballName+" topiary_jobs")

    #Where do you want to save it
    eosForOutput = "root://cmseos.fnal.gov//store/user/lpcboostres/topiaries_systematics-"+syststr+"_"+str(date.today())
    eosOnlypath  = eosForOutput.split("root://cmseos.fnal.gov/")[-1]
    os.system("eos root://cmseos.fnal.gov mkdir {}".format(eosOnlypath))

    #Submit the jobs
    if len(fsjson.keys()) > 0:
        for key in fsjson.keys():
            for samp in fsjson[key]:
                if args.sample:
                    if args.sample not in samp:
                        continue
                #print(fsjson[key][samp])
                sampleName = samp.split("_13TeV")[0]
                #print("    ",sampleName)
                fullpaths = ["root://cmseos.fnal.gov/"+key+"/"+x for x in fsjson[key][samp]]
                samplistasstr = str(fullpaths).replace(' ','')
                #print(samplistasstr)

           
                print("Topiaries from {0} are being written to {1}".format(sampleName,eosOnlypath))

                #Args to pass
                argu = "Arguments = {0} {1} {2} {3} {4}\n".format(eosForOutput,sampleName,samplistasstr,args.channel,args.syststr)

                #Make the jdl for each sample
                jdlName = "topiary_"+sampleName+"_"+str(date.today())+".jdl"
                jdl = open(jdlName,"w")
                jdl.write("universe = vanilla\n")
                jdl.write("Should_Transfer_Files = YES\n")
                jdl.write("WhenToTransferOutput = ON_EXIT\n")
                jdl.write("Transfer_Input_Files = "+tarballName+"\n")
                jdl.write("Output = condorMonitoringOutput/{0}/{1}_out.stdout\n".format(str(date.today()),sampleName+"_"+syststr))
                jdl.write("Error = condorMonitoringOutput/{0}/{1}_err.stder\n".format(str(date.today()),sampleName+"_"+syststr))
                jdl.write("Log = condorMonitoringOutput/{0}/{1}_log.log\n".format(str(date.today()),sampleName+"_"+syststr))
                jdl.write("Executable = topiary.sh\n")
                jdl.write(argu)
                jdl.write('+DESIRED_Sites="T3_US_Baylor,T2_US_Caltech,T3_US_Colorado,T3_US_Cornell,T3_US_FIT,T3_US_FNALLPC,T3_US_Omaha,T3_US_JHU,T3_US_Kansas,T2_US_MIT,T3_US_NotreDame,T2_US_Nebraska,T3_US_NU,T3_US_OSU,T3_US_Princeton_ICSE,T2_US_Purdue,T3_US_Rice,T3_US_Rutgers,T3_US_MIT,T3_US_NERSC,T3_US_SDSC,T3_US_FIU,T3_US_FSU,T3_US_OSG,T3_US_TAMU,T3_US_TTU,T3_US_UCD,T3_US_UCSB,T2_US_UCSD,T3_US_UMD,T3_US_UMiss,T2_US_Vanderbilt,T2_US_Wisconsin"')
                jdl.write("\n")
                jdl.write("Queue 1\n")#Not sure about this one
                jdl.close()
                
                #submit the jobs
                os.system("condor_submit {0}".format(jdlName))
