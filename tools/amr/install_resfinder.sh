cd $1


git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd resfinder

#BLAST
cd cge
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz

mkdir blast
tar -C blast -zxvf ncbi-blast-2.12.0+-x64-linux.tar.gz
cd blast
cd ncbi-blast-2.12.0+/
mv ./* ../
cd ../
rmdir ncbi-blast-2.12.0+/
cd ../
ln ./blast/bin/blastn


#databases
#go to resfinder dir
cd ../
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder

#indexing databases with KMA
#Go to the database directory
cd db_resfinder
#install, use and remove kma
#python3 INSTALL.py - interactive
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
cd ../
kma/kma_index -i fusidicacid.fsa -o ./fusidicacid
kma/kma_index -i phenicol.fsa -o ./phenicol
kma/kma_index -i glycopeptide.fsa -o ./glycopeptide
kma/kma_index -i trimethoprim.fsa -o ./trimethoprim
kma/kma_index -i oxazolidinone.fsa -o ./oxazolidinone
kma/kma_index -i tetracycline.fsa -o ./tetracycline
kma/kma_index -i quinolone.fsa -o ./quinolone
kma/kma_index -i nitroimidazole.fsa -o ./nitroimidazole
kma/kma_index -i fosfomycin.fsa -o ./fosfomycin
kma/kma_index -i aminoglycoside.fsa -o ./aminoglycoside
kma/kma_index -i macrolide.fsa -o ./macrolide
kma/kma_index -i sulphonamide.fsa -o ./sulphonamide
kma/kma_index -i rifampicin.fsa -o ./rifampicin
kma/kma_index -i colistin.fsa -o ./colistin
kma/kma_index -i beta-lactam.fsa -o ./beta-lactam
kma/kma_index -i disinfectant.fsa -o ./disinfectant
kma/kma_index -i pseudomonicacid.fsa -o ./pseudomonicacid
kma/kma_index -i *.fsa -o ./all

