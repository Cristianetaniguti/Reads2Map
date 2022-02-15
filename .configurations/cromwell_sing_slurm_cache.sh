backend {
  default = slurm

  providers {
    slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        filesystems {
           local {
             localization: [
                   "hard-link", "copy"
             ]
           }
        },
        max-concurrent-workflows = 1
        concurrent-job-limit = 5
        runtime-attributes = """
        String? time
        String? tasks
        String? mem
        String? node
        String? docker
        """

        submit = """
            sbatch \
              --wait \
              -J ${job_name} \
              -D ${cwd} \
              -o ${out} \
              -e ${err} \
              --time ${time} \
              ${tasks} \
              ${node} \
              ${mem} \
              --wrap "/bin/bash ${script}"
        """

         submit-docker = """

            export SINGULARITY_CACHEDIR=$SCRATCH/.singularity
            export SINGULARITY_TMPDIR=$SCRATCH/temp

            module load WebProxy # load the Singularity

            # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
            # based on the users home.
            if [ -z $SINGULARITY_CACHEDIR ];
                then CACHE_DIR=$SCRATCH/.singularity
                else CACHE_DIR=$SINGULARITY_CACHEDIR
            fi

            # Build the Docker image into a singularity image
            DOCKER_NAME=$(sed -e 's/[^A-Za-z0-9._-]/_/g' <<< ${docker})

            # The image will live together with all the other images to force "caching" of the .sif files themselves - note, always use docker hub tags!!!

            IMAGE=$SINGULARITY_CACHEDIR/$DOCKER_NAME.sif

            if [ ! -f $IMAGE ]; then  # If we already have the image, skip everything
                singularity pull $IMAGE docker://${docker}
            fi

            echo -e "#!/bin/sh\nsingularity exec --containall --bind ${cwd}:${docker_cwd} $IMAGE ${job_shell} ${docker_script}" > ${cwd}/script_${job_name}.sh                                                                                      
            
            # Submit the script to SLURM
            sbatch \
              --wait \
              -J ${job_name} \
              -D ${cwd} \
              -o ${cwd}/execution/stdout \
              -e ${cwd}/execution/stderr \
              --time ${time} \
              ${tasks} \
              ${node} \
              ${mem} \
              --export=NONE \
              ${cwd}/script_${job_name}.sh
        """

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"
      }
    }
  }
}

database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://$FIRST_NODE@127.0.0.1:/cromwell?rewriteBatchedStatements=true&useSSL=false"
    user = "root"
    password = "1234"
    connectionTimeout = 5000
  }
}

call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}