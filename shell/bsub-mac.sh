### WRITTEN ON MACBOOKPRO

#BSUB -q long
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 20:00
#BSUB -R 'rusage[mem=2048]'
#BSUB -n 64 
#BSUB -R 'span[hosts=1]'

module load R/3.5.1_gcc8.1.0 && module load gcc/8.1.0 && module load pandoc/2.7.2 && Rscript --vanilla -e 'rmarkdown::render("../Rscript/tnseq-timbr.Rmd")'
