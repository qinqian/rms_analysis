__conda_setup="$('/PHShome/qq06/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/PHShome/qq06/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/PHShome/qq06/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/PHShome/qq06/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate deseq

