#!/usr/bin/env python3

#####################################################
#
#  Batch job manager for distributed calculations.
#  It's primary duty has been to manage distributed
#  higher-order RDM calculations, but it could be
#  adapted for any generic task.
#  (C) 2017 Leon Freitag, ETH Zurich
#
#  This script requires Python >= 2.7 to run
#  otherwise there will be syntax errors
#####################################################

import os
import sys
import getopt
import io
import configparser
from shutil import copyfile


# init global variables. Maybe it makes sense to create a class -> TODO...
submitScriptName=""
templateName=""
inputName=""

class JobList:
  def __init__(self,slist=[],rlist=[],njobs=0):
    # jobs already run successfully
    self.successlist=slist
    # jobs to be run
    self.runlist=rlist
    # total number of jobs
    self.num_jobs=njobs
    # dictionary mapping job names to queue jobids
    self.jobids={}

def usage():
   print(''' Usage: %s [ -r successlist faillist / -3 / -n ]
	-r successlist faillist: restart with a list of successful and failed jobs. Optional.
        -3: create inputs for transition 3-RDM instead of the 4-RDM
    -n: do not submit, only create directories. (works also without drmaa)
   ''')
def setup_file_names(t3rdm):
  "Setup script/template file names depending on the type of the RDM we want to calculate"
  global submitScriptName, templateName, inputName
  if (t3rdm):
    submitScriptName="submit-3rdm.sh"
    templateName="meas-3rdm-template.in"
    inputName="meas-3rdm.in"
  else:
    submitScriptName="submit-4rdm.sh"
    templateName="meas-4rdm-template.in"
    inputName="meas-4rdm.in"

def restart(successname,failname):
  """Restart the calculation by reading in the faillist (and the successlist)
  Specification of the successlist:
  1st line: <job count> jobs in total
  subsequent lines: successful job identifiers
  """
  try:
    with open(successname) as f:
      successlines=f.read().splitlines()
      if 'jobs' in successlines[0]:
        # read in total # of jobs
        try:
          num=int(successlines[0].split()[0])
          print("Restarting jobs")
          print(num, "jobs in total")
        except ValueError:
          print("Poorly formatted successfile!")
          sys.exit(2)
        # delete the first line, other lines should be job identifiers that are passed onto success list
        del successlines[0]
      else:
        print("Poorly formatted successfile -- header is wrong!")
        sys.exit(2)
    with open(failname) as f:
      faillines=f.read().splitlines()
    print(len(faillines), "jobs to be restarted")
    if len(faillines)==0:
      print("All jobs were successful, no need to restart!")
      sys.exit(0)
  except IOError:
    print("Problem reading file %s" % successname) ## TODO: fix
    sys.exit(1)
  return JobList(successlines,faillines,num)

def create_job_scripts_dummy():
  """Creates dummy jobs. Useful for debugging purposes
  """
  print("No restart option specified, will prepare a new calculation.")
  runlist=[]
  num=0
  L=2 # for testing. actually this should be read in from meas-rdm.template
  for i in range(L):
    for j in range(i+1):
      for k in range(j+1):
        for l in range(k+1):
          if not (((i==j) and (i==k) ) or ((i==j) and (i==l)) or ((i==k) and (i==l)) or ((j==k) and (j==l))):
            name="part-" + str(i) + "-" + str(j) + "-" + str(k) + "-" + str(l)
            print("Preparing " + name)
            # Create the path and make a directory for every path
            namepath = os.getcwd() + "/parts/" + name
            try:
              os.makedirs(namepath,mode=0o755)
            except OSError:
              print("Directory " + name + " could not be created -- maybe it already exists? Continuing.")
            scriptname=namepath+"/"+submitScriptName
            fromscriptname=os.getcwd() + "/template/" + submitScriptName
            copyfile(fromscriptname,scriptname)
            os.chmod(scriptname,0o755)
            runlist.append(name)
            num+=1

  return JobList([],runlist,num)

def create_job_scripts_4rdm():
  """Prepare the directory structure for the job submission.
  This function does so for the calculation of the 4-RDM, but it can be freely substituted
  by any other job preparation procedure.
  """
  print("No restart option specified, will prepare a new calculation.")
  runlist=[]
  num=0

  # open and the template and read the number of sites from it
  # we use ConfigParser which expects a configfile of the syntax "[Section]\n Key=Value"
  # QCMaquis inputs do not have sections, so we had to add one on the fly, which we do
  # ConfigParser needs a file-like object to read the config, which we get from StringIO
  templatepath = os.getcwd() + "/template/" + templateName
  try:
    with open(templatepath) as templatefile:
      templateString = "[dummy]\n" + templatefile.read()
      templateIO = io.StringIO(templateString)
      templateConfig = configparser.ConfigParser(allow_no_value=True)
      templateConfig.optionxform = str # disable converting the keys in the config file to lowercase
      templateConfig.readfp(templateIO)
      try:
        L=templateConfig.getint("dummy",'L')
      except ValueError:
        print("Cannot read L from the template file")
        sys.exit(2)
  except IOError:
    print("Cannot open the template file:", templatepath)
    sys.exit(2)
  for i in range(L):
    for j in range(i+1):
      for k in range(j+1):
        for l in range(k+1):
          if not (((i==j) and (i==k)) or ((i==j) and (i==l)) or ((i==k) and (i==l)) or ((j==k) and (j==l))):
            name="part-" + str(i) + "-" + str(j) + "-" + str(k) + "-" + str(l)
            print("Preparing " + name)
            # Create the path and make a directory for every path
            namepath = os.getcwd() + "/parts/" + name
            try:
              os.makedirs(namepath,mode=0o755)
            except OSError:
              print("Directory " + name + " could not be created -- maybe it already exists? Continuing.")
            scriptname=namepath+"/" + submitScriptName

            # set the indexes in the Maquis input file and write it back to the config
            templateConfig.set("dummy","MEASURE[4rdm]","'p4:p3:p1:p2@" + str(l) + "," + str(k) + "," + str(i) + "," + str(j) + "'")
            # clear the templateIO string buffer. maybe there's a better way to do so?
            templateIO.close()
            templateIO = io.StringIO()
            templateConfig.write(templateIO)
            # now strip the dummy section header
            with open(namepath+"/"+ inputName,'w') as inputfile:
              templateString=templateIO.getvalue()
              templateString=templateString[templateString.find('\n')+1:templateString.rfind('\n')]
              inputfile.write(templateString)

            fromscriptname=os.getcwd() + "/template/" + submitScriptName
            copyfile(fromscriptname,scriptname)
            os.chmod(scriptname,0o755)

            runlist.append(name)

            num+=1

            templateIO.close()

  return JobList([],runlist,num)

def create_job_scripts_3rdm():
  """Prepare the directory structure for the job submission of the transition 3-RDM
     This version creates 4-index slices for the 3-RDM (there will be (L-1)*L slices)
     QCMaquis supports also 3-index slices but it's not implemented as of now
  """
  print("No restart option specified, will prepare a new calculation.")
  runlist=[]
  num=0

  # open and the template and read the number of sites and the bra checkpoint name from it
  # we use ConfigParser which expects a configfile of the syntax "[Section]\n Key=Value"
  # QCMaquis inputs do not have sections, so we had to add one on the fly, which we do
  # ConfigParser needs a file-like object to read the config, which we get from StringIO
  # when we write the config back, the dummy section header gets removed again
  templatepath = os.getcwd() + "/template/" + templateName
  try:
    with open(templatepath) as templatefile:
      templateString = "[dummy]\n" + templatefile.read()
      templateIO = io.StringIO(templateString)
      templateConfig = configparser.ConfigParser(allow_no_value=True)
      templateConfig.optionxform = str # disable converting the keys in the config file to lowercase
      templateConfig.readfp(templateIO)
      try:
        L=templateConfig.getint("dummy",'L')
      except ValueError:
        print("Cannot read L from the template file")
        sys.exit(2)
      # we must read the bra checkpoint name to re-insert it later
      try:
        brachkname=templateConfig.get("dummy","MEASURE[trans3rdm]").split(';')[0].strip('"')
      except:
        print("Cannot read the bra checkpoint file from the template file")
        sys.exit(2)

  except IOError:
    print("Cannot open the template file:", templatepath)
    sys.exit(2)
  for i in range(L-1):
    for j in range(i,L):
         name="part-" + str(i) + "-" + str(j)
         print("Preparing " + name)
         # Create the path and make a directory for every path
         namepath = os.getcwd() + "/parts/" + name
         try:
           os.makedirs(namepath,mode=0o755)
         except OSError:
           print("Directory " + name + " could not be created -- maybe it already exists? Continuing.")
         scriptname=namepath + "/" + submitScriptName

         # set the indexes in the Maquis input file and write it back to the config
         templateConfig.set("dummy","MEASURE[trans3rdm]","'"+ brachkname + ";p1:p2@" + str(j) + "," + str(i) + "'")
         # clear the templateIO string buffer. maybe there's a better way to do so?
         templateIO.close()
         templateIO = io.StringIO()
         templateConfig.write(templateIO)
         # now strip the dummy section header
         with open(namepath+"/"+ inputName,'w') as inputfile:
           templateString=templateIO.getvalue()
           templateString=templateString[templateString.find('\n')+1:templateString.rfind('\n')]
           inputfile.write(templateString)

         fromscriptname=os.getcwd() + "/template/" + submitScriptName
         copyfile(fromscriptname,scriptname)
         os.chmod(scriptname,0o755)

         runlist.append(name)

         num+=1

         templateIO.close()

  return JobList([],runlist,num)

def collect_results(s,joblist):
  """Collect the job results and write successfile and failfile.
  """
  try:
    with open("successlist",'w') as successfile, open("faillist",'w') as failfile:
      # write total # of jobs as the header for the successfile
      successfile.write(str(joblist.num_jobs) + " jobs in total\n")
      for jobid in joblist.jobids:
         print("Collecting job " + jobid + " " + joblist.jobids[jobid])
         retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
         if (s.jobStatus(retval.jobId) == drmaa.JobState.DONE) and (retval.exitStatus == 0) :
           # success! write the job to the success list
           if jobid != retval.jobId:
             print("Warning, jobid ", jobid, "doesn't match retval.jobid ", retval.jobId)
           print("Job " + jobid + " " + joblist.jobids[retval.jobId] + " succeeded")
           successfile.write(joblist.jobids[retval.jobId]+"\n")
         else:
           if (retval.exitStatus != 0) or (s.jobStatus(retval.jobId) == drmaa.JobState.FAILED):
             if jobid != retval.jobId:
               print("Warning, jobid ", jobid, "doesn't match retval.jobid ", retval.jobId)
             # fail! write the job to the fail list
             print("Job " + retval.jobId + " " + joblist.jobids[retval.jobId] + " failed")
             failfile.write(joblist.jobids[retval.jobId]+"\n")
           else:
             # this should never happen
             raise AppError("Unexpected job state " + str(s.jobStatus(retval.jobId)) + "for job " + joblist.jobids[retval.jobId] + " with id " + retval.jobId)
  except IOError as e:
    print("Cannot create successlist or faillist: %s" % e.strerror)

def launch_and_wait(joblist):
  """Launch jobs and wait for their completion.
  """

  ## init DRMAA stuff
  s = drmaa.Session()
  s.initialize()
  ## run jobs
  for job in joblist.runlist:
    jt = s.createJobTemplate()
    # remote script name. WARNING: you can NOT use custom LSF commands in comments of the script! (It'd be nice to have but I still have to figure out how!)
    jt.workingDirectory = os.getcwd() + "/parts/" + job
    jt.remoteCommand = "./" + submitScriptName
    # jobname
    jt.jobName = job
    jt.errorPath = ":" + os.getcwd() + "/parts/" + job + "/error.err"
    print('Submitting job ', job)
    ### Wall time + CPU # specifications. Modify here as needed
    jt.hardWallclockTimeLimit = "24:00:00"
    jt.nativeSpecification = "-n 1"

    joblist.jobids[s.runJob(jt)] = job
    s.deleteJobTemplate(jt)

  ## here comes the very fun part: wait for all jobs to finish
  s.synchronize([i for i in joblist.jobids], drmaa.Session.TIMEOUT_WAIT_FOREVER, False)

  ## collect the job results
  collect_results(s,joblist)

  s.exit()

def main():

  assert sys.version_info >= (2,7)

  restartmode=False
  trans3rdmmode=False
  nosubmit=False
  successlistname=""
  faillistname=""

  successlist=[]
  faillist=[]

  if len(sys.argv)>1:
     try:
       opts,args = getopt.getopt(sys.argv[1:], "3rn")
       if ('-r','') in opts:
         restartmode=True
         successlistname,faillistname=args
       if ('-3','') in opts:
         trans3rdmmode=True
       if ('-n','') in opts:
         nosubmit=True
     except getopt.GetoptError as err:
       usage()
       sys.exit(1)
  setup_file_names(trans3rdmmode)
  jl = restart(successlistname,faillistname) if restartmode else (create_job_scripts_3rdm() if trans3rdmmode else create_job_scripts_4rdm())
  if not nosubmit:
    try:
        import drmaa
        launch_and_wait(jl)
    except ImportError:
        print("""This script requires python-drmaa to run. Please install it from https://pypi.org/project/drmaa/
                or ask your system administrator kindly to do it.""")
        sys.exit(3)



if __name__=='__main__':
    main()

