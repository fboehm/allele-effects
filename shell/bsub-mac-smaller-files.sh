### WRITTEN ON MACBOOKPRO

#BSUB -q long
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 20:00
#BSUB -R 'rusage[mem=20480]'
#BSUB -n 1
#BSUB -R 'span[hosts=1]'

module load R/3.5.1_gcc8.1.0 && module load gcc/8.1.0 && module load pandoc/2.7.2 && Rscript --vanilla -e 'rmarkdown::render("../Rscript/timbr-small-files.Rmd")'
