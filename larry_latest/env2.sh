__conda_setup="$('/data/pinello/SHARED_SOFTWARE/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
   eval "$__conda_setup"
else
   if [ -f "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh" ]; then
       . "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh"
   else
       export PATH="/data/pinello/SHARED_SOFTWARE/anaconda3/bin:$PATH"
   fi  
fi

unset __conda_setup
conda activate /data/pinello/SHARED_SOFTWARE/anaconda3/envs/sc-tutorial
