# HaploFrequencies 

The latest version of the AADR table can be downloaded from the [David Reich Lab](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data). 

The first argument specifies the DNA type (mtDNA/YDNA), the second argument specifies the range of a population in years (e.g. 100/500/1000), the third argument specifies the spatial grouping method (Country/Locality) and the last argument specifies the level of basal haplogroup (1/2). HaploFrequencies can for example be run like this:

```
./HaploFrequencies.R mtDNA 500 Locality 2
```

Alternatively, a previously created output file can be visualized without recalculating haplogroup frequencies:

```
./HaploFrequencies.R frequencies.xlsx
```
