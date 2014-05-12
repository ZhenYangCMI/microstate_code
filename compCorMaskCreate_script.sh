for i in `ls /home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/FunImg`; do
    3dcalc -a /home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/Masks/${i}_CsfMask_07_91x109x91.nii -b /home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/Masks/${i}_WhiteMask_09_91x109x91.nii -expr 'a+b' -prefix /home2/data/Projects/microstate/DPARSF_preprocessed/data/645/all_10min/preprocessed/Masks/${i}_CsfWhiteMask_91x109x91;
    echo "Image for subject $i has been resampled";
done
