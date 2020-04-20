#!/bin/csh
# Set MATH alias - takes an arithmetic assignment statement
# as argument, e.g., newvar = var1 + var2
# Separate all items and operators in the expression with blanks
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
set folder=/home-2/hchangl1@jhu.edu/data/Runge_Kutta
foreach tottime(60.00)
set scratchfolder=/home-2/hchangl1@jhu.edu/scratch/Spin_1_pyrochlore/Sqw/heisenberg
set move=conical
set disorder=0.00
#foreach nsamples(500000000)
foreach nsamples(50000000)
foreach L (8)
#foreach seed(1)
#foreach seed(2 3 4 5 6 7 8 9 10)
foreach seed(11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
foreach temp(1.80) 

cat >in_heisenberg_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}.txt << EOFm
seed=${seed}
spin=1.0
L=${L}
temp=${temp}
tmin=${temp}
tmax=${temp}
J1=+3.20
J2=+3.20
J3=+0.00
J4=+0.00
Jnnn=+0.00
gxy=1.0
gz=1.0
deltat=0.02
nstarts=192
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

cat >job_heisenberg_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples} << EOFm
#!/bin/bash -l
# Batch Queue Script
#SBATCH --time=48:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --account=olegt

cd ${folder}/jobs
source ~/.bashrc
export OMP_NUM_THREADS=48
${folder}/mc ${folder}/jobs/in_heisenberg_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}.txt > ${scratchfolder}/out_heisenberg_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}.txt

EOFm

chmod 755 job_heisenberg_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}
sleep 5
qsub job_heisenberg_seed_${seed}_disorder_${disorder}_L_${L}_temp_${temp}_tottime_${tottime}_replica_move_${move}_nsamples_${nsamples}

end
end
end
end
end
