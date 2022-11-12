# Create folders
mkdir ~/Code
mkdir ~/Code/runs
mkdir ~/Code/runs/RMG
mkdir ~/Code/runs/RMG/example_1

# Copy RMG minimal example
cd ~/Code/runs/RMG/example_1
cp /home/alon/Code/RMG-Py/examples/rmg/minimal/input.py .

# Edit bashrc
cp ~/.bashrc ~/.bashrc_back1
rm ~/.bashrc
wget https://github.com/DanaResearchGroup/056402/raw/main/Aluf/students_bashrc -O ~/.bashrc

