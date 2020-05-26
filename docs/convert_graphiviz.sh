
for i in $(cat workflows); do
    filename=$(basename -- "$i")
    echo $filename
    filename2=${filename%.*}
    echo $filename2
    java -jar ~/womtool-49.jar graph $i > images/$filename2.gv
    dot -Tsvg images/$filename2.gv > images/$filename2.svg
done
