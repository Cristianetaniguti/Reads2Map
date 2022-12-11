version 1.0

task RunPedigreeSimulator {
  input {
    File mapfile
    File founderfile
    File chromfile
    File parfile
  }

  Int disk_size = ceil(size(mapfile, "GiB") * 2 + size(founderfile, "GiB") + size(chromfile, "GiB") * 2 + size(parfile, "GiB") * 2) 
  Int memory_size = 3000

  command <<<
    cp ~{parfile} parfile.txt
    set -e
    sed -i 's+chromosome.txt+~{chromfile}+g' parfile.txt
    sed -i 's+mapfile.txt+~{mapfile}+g' parfile.txt
    sed -i 's+founderfile.txt+~{founderfile}+g' parfile.txt
    java -jar /usr/jars/PedigreeSim.jar parfile.txt

  >>>

  runtime {
    docker: "cristaniguti/java-in-the-cloud:0.0.1"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "RunPedigreeSimulator"
    mem:"~{memory_size}M"
    time:"05:00:00"
  }

  meta {
    author: "Cristiane Taniguti"
    email: "chtaniguti@tamu.edu"
    description: "Runs [PedigreeSim](https://www.wur.nl/en/show/Software-PedigreeSim.htm)."
  }

  output {
    File genotypes_dat = "sim_genotypes.dat"
  }

}