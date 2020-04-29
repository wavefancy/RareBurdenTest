conda create -n burden -yq r-essentials r-base scipy pandas rpy2=2.9.4 docopt
conda activate burden
conda install -yq -c conda-forge r-lmtest tzlocal
conda deactivate
conda clean --all -yq
# conda remove --name burden --all
# conda evn export > conda_burden.yml
# conda env create -n burden -f conda_burden.yml
