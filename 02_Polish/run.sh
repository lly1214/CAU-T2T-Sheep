set -ex
date
hostname
/opt/jdk-11.0.6/bin/java  -Xmx20g   -Dsystem.job-shell=/bin/bash     -jar /export/pipeline/RNASeq/Software/Cromwell/cromwell-80.jar   run  main.wdl    -i run.json
date
touch done
