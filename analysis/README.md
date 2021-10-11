# Assembly summaries

##### hifiasm results
```
sum = 938048665, n = 1367, ave = 686209.70, largest = 46801135
N50 = 14953593, n = 20
N60 = 12536664, n = 27
N70 = 9770061, n = 36
N80 = 7795771, n = 47
N90 = 2567837, n = 65
N100 = 3559, n = 1367
N_count = 0
Gaps = 0
```

##### purge_haplotigs results

```
sum = 869214813, n = 217, ave = 4005598.22, largest = 46801135
N50 = 16234331, n = 18
N60 = 13465811, n = 24
N70 = 11076384, n = 31
N80 = 8788749, n = 40
N90 = 7028129, n = 51
N100 = 10153, n = 217
N_count = 0
Gaps = 0
```

##### Bionano Round 1 results
```
sum = 862310619, n = 16, ave = 53894413.69, largest = 122473895
N50 = 74699160, n = 5
N60 = 70180950, n = 6
N70 = 55717035, n = 8
N80 = 48331264, n = 9
N90 = 33677862, n = 11
N100 = 368418, n = 16
N_count = 6412036
Gaps = 1388
```

##### Hi-C round 1 results
```
sum = 861332878, n = 10, ave = 86133287.80, largest = 122473895
N50 = 88736738, n = 5
N60 = 82690113, n = 6
N70 = 79495078, n = 7
N80 = 74699160, n = 8
N90 = 70180950, n = 9
N100 = 63045887, n = 10
N_count = 6414013
Gaps = 1391
```

##### Bionano round 2 results
```
sum = 858656212, n = 10, ave = 85865621.20, largest = 101681876
N50 = 88736738, n = 5
N60 = 82690113, n = 6
N70 = 81806108, n = 7
N80 = 78761988, n = 8
N90 = 74699160, n = 9
N100 = 70176875, n = 10
N_count = 6591978
Gaps = 1449
```

##### Racon (6 rounds) and DENTIST gap closure
```
sum = 857417049, n = 10, ave = 85741704.90, largest = 101371656
N50 = 88645700, n = 5
N60 = 82636096, n = 6
N70 = 81616834, n = 7
N80 = 78696137, n = 8
N90 = 74632051, n = 9
N100 = 70031136, n = 10
N_count = 6229811
Gaps = 48
```


# BUSCO results for hifiasm genome
```
        C:94.0%[S:88.4%,D:5.6%],F:1.7%,M:4.3%,n:2121

        1995    Complete BUSCOs (C)
        1876    Complete and single-copy BUSCOs (S)
        119     Complete and duplicated BUSCOs (D)
        36      Fragmented BUSCOs (F)
        90      Missing BUSCOs (M)
        2121    Total BUSCO groups searched
```


# MAKER2 results

##### Round 1
```
Total:  129665327
Count:  53453
Mean:   2425
Median: 872
Min:    2
Max:    81770
```

##### Round 2
```
Total:  124452492
Count:  33687
Mean:   3694
Median: 2549
Min:    179
Max:    91324
```

# BUSCO results for different isoseq datasets:

##### Agrostemma
`python /home/peter/miniconda3/envs/busco/bin/run_BUSCO.py -i agrostemma.clustered.hq_mod.fasta -o agrostemma -l /home/peter/bioinformatics/busco_db/eudicotyledons_odb10/ -m transcriptome -c 50`

        C:87.1%[S:12.8%,D:74.3%],F:2.2%,M:10.7%,n:2121

        1847    Complete BUSCOs (C)
        271     Complete and single-copy BUSCOs (S)
        1576    Complete and duplicated BUSCOs (D)
        47      Fragmented BUSCOs (F)
        227     Missing BUSCOs (M)
        2121    Total BUSCO groups searched

##### S. conica
`python /home/peter/miniconda3/envs/busco/bin/run_BUSCO.py -i conica.clustered.hq_mod.fasta -o conica -l /home/peter/bioinformatics/busco_db/eudicotyledons_odb10/ -m transcriptome -c 50`

        C:82.7%[S:19.9%,D:62.8%],F:2.3%,M:15.0%,n:2121

        1755    Complete BUSCOs (C)
        423     Complete and single-copy BUSCOs (S)
        1332    Complete and duplicated BUSCOs (D)
        49      Fragmented BUSCOs (F)
        317     Missing BUSCOs (M)
        2121    Total BUSCO groups searched
        
##### S. latifolia
`python /home/peter/miniconda3/envs/busco/bin/run_BUSCO.py -i latifolia.clustered.hq_mod.fasta -o latifolia -l /home/peter/bioinformatics/busco_db/eudicotyledons_odb10/ -m transcriptome -c 50`

        C:87.5%[S:16.4%,D:71.1%],F:2.0%,M:10.5%,n:2121

        1855    Complete BUSCOs (C)
        348     Complete and single-copy BUSCOs (S)
        1507    Complete and duplicated BUSCOs (D)
        43      Fragmented BUSCOs (F)
        223     Missing BUSCOs (M)
        2121    Total BUSCO groups searched
        
##### S. noctiflora
`python /home/peter/miniconda3/envs/busco/bin/run_BUSCO.py -i noctiflora.clustered.hq_mod.fasta -o noctiflora -l /home/peter/bioinformatics/busco_db/eudicotyledons_odb10/ -m transcriptome -c 50`

        C:79.7%[S:28.4%,D:51.3%],F:3.3%,M:17.0%,n:2121

        1691    Complete BUSCOs (C)
        602     Complete and single-copy BUSCOs (S)
        1089    Complete and duplicated BUSCOs (D)
        70      Fragmented BUSCOs (F)
        360     Missing BUSCOs (M)
        2121    Total BUSCO groups searched
        
##### S. vulgaris
`python /home/peter/miniconda3/envs/busco/bin/run_BUSCO.py -i vulgaris.clustered.hq_mod.fasta -o vulgaris -l /home/peter/bioinformatics/busco_db/eudicotyledons_odb10/ -m transcriptome -c 50`

        C:89.6%[S:17.0%,D:72.6%],F:2.1%,M:8.3%,n:2121

        1900    Complete BUSCOs (C)
        361     Complete and single-copy BUSCOs (S)
        1539    Complete and duplicated BUSCOs (D)
        44      Fragmented BUSCOs (F)
        177     Missing BUSCOs (M)
        2121    Total BUSCO groups searched
