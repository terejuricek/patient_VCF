# patient_VCF
Complete chromosomal positions for each patient from CRAM files

Doplnenie vcf súborov
- doplnenie chromosomálnych pozící súborov pre daného pacienta
- 39 pacientov EP, LP -> 142 súborov

Postup
- uložíme si zoznam všetkých jedinečných chrom pozicií pre jedného pacienta
- prechádzame cez súbory a doplňame ich na základe missing_positions zoznamu
- doplnenie z cram pre chýbajúcu pozíciu
- manuálne pridanie INFO/GT

for patient in MetaCentrum
    for *vcf vcfs.chrom+pos
        union()
            unique
list of positions(+chrom) per patient

    pac1 ... chr1: [3748, 1234, ...], chr2: ...
    pac2 ...


for vcf in vcfs
    which patient
        get.positions()
            missing_positions()
            list missing_positions()


for missing_positions
    vcf == cram_name

bash bcftools mpileup  -- anotate
    vcf --
    !!!! ručne pridať INFOR/GT

