# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/alon/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/alon/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/alon/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/alon/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<



# RMG-Py
export rmgpy_path='/home/alon/Code/RMG-Py/'
export rmgdb_path='/home/alon/Code/RMG-database/'
export PYTHONPATH=$PYTHONPATH:/home/alon/Code/RMG-Py/

# ARC
export arc_path='/home/alon/Code/ARC/'
export PYTHONPATH=$PYTHONPATH:$arc_path
export PYTHONPATH=$PYTHONPATH:/home/alon/Code/AutoTST/
export PYTHONPATH=$PYTHONPATH:/home/alon/Code/KinBot/
export PYTHONPATH=$PYTHONPATH:/home/alon/Code/TS-GCN/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/alon/miniconda3/envs/sella_env/lib

# T3
export t3_path='/home/alon/Code/T3/'
export PYTHONPATH=$PYTHONPATH:$t3_path

# personalized aliases

alias rc='source ~/.bashrc'
alias rce='nano ~/.bashrc'
alias erc='nano ~/.bashrc'

alias rmge='conda activate rmg_env'
alias arce='conda activate arc_env'
alias tcke='conda activate tck_env'
alias t3e='conda activate t3_env'
alias deact='conda deactivate'

alias rmg='python-jl $rmgpy_path/rmg.py input.py  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'
alias arkane='python-jl $rmgpy_path/Arkane.py input.py  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'
alias arc='python $arc_path/ARC.py input.yml  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'
alias arcrestart='python $arc_path/ARC.py restart.yml  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'
alias restartarc='python $arc_path/ARC.py restart.yml  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'
alias t3='python-jl $t3_path/T3.py input.yml  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'

alias tst='pytest -ra -vv'
alias runs='cd ~/runs'


