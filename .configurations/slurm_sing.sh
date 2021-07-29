# Configure Cromwell to submit jobs to Slurm, including support for Singularity
#
# Author:
#   - Michael Franklin <michael.franklin@unimelb.edu.au>
#
# History:
#   - 2020-07-23 - Initial + minor fixes to generalise the script
#
# Quickstart:
#   - Replace <location> with a location for a singularity cache. I'd recommend following this GitHub thread for some information about this: https://github.com/broadinstitute/cromwell/pull/5515
#   - [Optional] Add a queue after `String? queue` to be `String? queue = "yourqueue"
#
# About
#   
#   - Transform some job information and path to get a reasonable slurm job name including shard + cpu/mem (easier to track)
#   - We submit a 'wrap' job (currently it's only implemented for submit-docker) to catch times where SLURM kills the job
#   - The regular 'submit' just submits the variables as required
#   - For "submit-docker", we use a cache location to pull images to.
#   - 'duration' is in seconds, and can be passed from your WDL runtime (it's not currently a recognised K-V)
#       - [OpenWDL #315](https://github.com/openwdl/wdl/pull/315)
#       - Cromwell doesn't (/ didn't) support ToolTimeRequirement for CWL

akka: {
  "actor.default-dispatcher.fork-join-executor": {
    "parallelism-max": 3
  }
}

system: {
  "job-shell": "/bin/sh"
}

backend: {
  "default": "slurm-singularity",
  "providers": {
    "slurm-singularity": {
      "actor-factory": "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory",
      "config": {
        "filesystems": {
          "local": {
            "localization": [
              "hard-link",
              "cached-copy"
            ],
            "enabled": true,
            "caching": {
              "duplication-strategy": [
                "hard-link",
                "cached-copy",
                "copy",
                "soft-link"
              ],
              "hashing-strategy": "fingerprint"
            }
          }
        },
        "runtime-attributes": """
Int duration = 86400
Int? cpu = 1
Int memory_mb = 3500
String? docker
String? queue
String cacheLocation = "<location>"
""",

        "submit": """
jobname='${sub(sub(cwd, ".*call-", ""), "/", "-")}-cpu-${cpu}-mem-${memory_mb}'
sbatch \
    -J $jobname \
    -D ${cwd} \
    -o ${out} \
    -e ${err} \
    -t 0:${duration} \
    ${"-p " + queue} \
    ${"-n " + cpu} \
    --mem=${memory_mb} \
    --wrap "/usr/bin/env ${job_shell} ${script}" 
""",
        "submit-docker": """

docker_subbed=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})
image=${cacheLocation}/$docker_subbed.sif
lock_path=${cacheLocation}/$docker_subbed.lock

if [ ! -f "$image" ]; then
  singularity pull $image docker://${docker}
fi

# Submit the script to SLURM
jobname=${sub(sub(cwd, ".*call-", ""), "/", "-")}-cpu-${cpu}-mem-${memory_mb}
JOBID=$(sbatch \
    --parsable \
    -J $jobname \
    --mem=${memory_mb} \
    --cpus-per-task ${select_first([cpu, 1])} \
    ${"-p " + queue} \
    -D ${cwd} \
    -o ${cwd}/execution/stdout \
    -e ${cwd}/execution/stderr \
    -t '0:${duration}' \
    --wrap "singularity exec --bind ${cwd}:${docker_cwd} $image ${job_shell} ${docker_script}") \
     && NTOKDEP=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afternotokay:$JOBID --wrap '[ ! -f rc ] && (echo 1 >> ${cwd}/execution/rc) && (echo "A slurm error occurred" >> ${cwd}/execution/stderr)') \
    && echo Submitted batch job $JOBID""",
        "kill": "scancel ${job_id}",
        "check-alive": "scontrol show job ${job_id}",
        "job-id-regex": "Submitted batch job (\\d+).*"
      }
    }
  }
}
call-caching: {
  "enabled": true
}