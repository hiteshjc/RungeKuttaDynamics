#!/bin/csh
# Set MATH alias - takes an arithmetic assignment statement
# as argument, e.g., newvar = var1 + var2
# Separate all items and operators in the expression with blanks
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
set folder=/gpfs/home/hchanglani/Runge_Kuttafaster/
foreach tottime(60.00)
set scratchfolder=/gpfs/research/changlani/MD_Changlani_runs/
set move=conical
set disorder=0.00
#foreach nsamples(500000000)
foreach nsamples(5000)
foreach L (8)
foreach seed(1)
foreach temp(1.80) 

cat >in_kemp_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}.txt << EOFm
seed=${seed}
spin=1.0
L=${L}
temp=${temp}
tmin=${temp}
tmax=${temp}
J1=+3.200
J2=+3.200
J3=+0.019
J4=-0.070
Jnnn=-0.025
gxy=1.0
gz=1.0
deltat=0.02
nstarts=1
disorder=${disorder}
mcmove=${move}
start_config=random
nsamples=${nsamples}
nreplicas=192
omegamin=0.00
omegamax=15.00
omegaspacing=0.05
tottime=${tottime}

EOFm

cat >job_kemp_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples} << EOFm
#!/bin/bash
# Batch Queue Script
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -p changlani_q

cd ${folder}/jobs
source ~/.bashrc
export OMP_NUM_THREADS=24
${folder}/mc ${folder}/jobs/in_kemp_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}.txt > ${scratchfolder}/out_kemp_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}.txt

EOFm

chmod 755 job_kemp_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}
sbatch job_kemp_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}

end
end
end
end
end
