## rawdata process
- [x] download.sh
- [x] trim_galore.sh
- [x] mapping.sh,samtool.sh to convert bam to sam
## <s>simulation </s>
- [x] ReferenceKmer.sh to determine real nucleotide freq per 2mer -k 2 
- [ ] bamkmers.sh Extract kmers(2~8mers) of left and right read ends from the real data
- [x] bamlength.sh Extract the length distribution from the BAM
- [ ] simulation.sh
## <s>fragment patterns </s>
- [ ] Frag.py dinucleotide frequency calculation on both sim& real data
- [x] bamlength for sim data 
## WPS apply to BH01~CH01
- [ ] alpha satellite sites
- [ ] heuristic peak calling
- [ ] get in vivo nucleosome occupancy
## Compartment A/B
- [ ] peak to peak spacing
## CTCF & TFs
- [ ] CHO1 WPS at 120-180bp(window 120) 35-80bp(window 16)
- [ ] clustered FIMO prediction
## Nucleosome spacing pattern
- [ ] peak to peak spacing within DHS&TFs
- [ ] diff high expr/low expr
## Reads counts across TSS over multiple gene
- [ ] healthy samples
- [ ] patients
