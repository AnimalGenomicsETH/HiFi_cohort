
## steps

download the cpg islands

```
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/bbi/GCF_002263795.3_ARS-UCD2.0.cpgIslandExtUnmasked.bb
```

download alias file

```
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/002/263/795/GCF_002263795.3/GCF_002263795.3.chromAlias.txt
```

and rename

```
awk 'NR==FNR{a[$1]=$2;b[$1]=$3;next}$1 in a{$1=($1~/NC/?a[$1]:b[$1])}1' GCF_002263795.3.chromAlias.txt GCF_002263795.3_ARS-UCD2.0.cpgIslandExtUnmasked.bed > GCF_002263795.3_ARS-UCD2.0.cpgIslandExtUnmasked.NEW.bed
```
