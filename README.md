Alguma descrição

## Quickstart

```
# Gerar as imagens
bash .scripts/build_images.sh

# Adequar o caminho dos inputs no F2.json

# Executar o workflow
java -jar -Dconfig.file=.configurations/cromwell.conf -jar cromwell.jar run -i F2.json F2.wdl
```
