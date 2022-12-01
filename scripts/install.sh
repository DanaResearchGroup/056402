# reset .bashrc
cp /home/$USER/.bashrc /home/$USER/.bashrc_backup_1
rm /home/$USER/.bashrc
rm /home/$USER/students_bashrc
wget /home/$USER https://raw.githubusercontent.com/DanaResearchGroup/056402/main/scripts/students_bashrc
cp /home/$USER/students_bashrc /home/$USER/.bashrc
source ~/.bashrc

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

source ~/.bashrc

conda remove --name arc_env --all -y
conda remove --name rmg_env --all -y

cd /home/$USER
mkdir runs
cd runs
mkdir example
cd example
wget . https://raw.githubusercontent.com/DanaResearchGroup/056402/main/RMG%20input%20files/1_minimal/input.py

cd /home/$USER
mkdir Code
cd Code

# Clone repos
git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
git clone https://github.com/ReactionMechanismGenerator/RMG-database.git

cd RMG-Py
conda env create -f environment.yml
source ~/.bashrc
conda activate rmg_env
make
python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator;'

conda deactivate

