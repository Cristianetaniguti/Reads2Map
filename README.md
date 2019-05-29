## Linkage map simulation

This workflow provides linkage map simulation using markers from RAD-seq methodology. The main goal of the workflow is to mesuare the impact of different error probabilities coming from SNP calling methodologies in the linkage map building in OneMap.

## Quickstart

```
# Gerar as imagens
bash .scripts/build_images.sh

# Adequar o caminho dos inputs no F2.json

# Executar o workflow
java -jar -Dconfig.file=.configurations/cromwell.conf -jar cromwell.jar run -i F2.json F2.wdl
```
