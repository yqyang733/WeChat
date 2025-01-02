# AlphaFold结构预测系列工具安装及使用
## AF2-colab安装
```shell
conda create -n colab
conda activate colab
conda install -c conda-forge python=3.10 cudnn==8.2.1.32 cudatoolkit==11.6.0 openmm==7.7.0 pdbfixer -y
conda install -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 mmseqs2=14.7e284 -y
pip install --upgrade pip
cd Colabfold-main
pip install --no-warn-conflicts "colabfold[alphafold]==1.5.3"   tensorflow==2.12.0
pip install https://storage.googleapis.com/jax-releases/cuda11/jaxlib-0.3.25+cuda11.cudnn82-cp310-cp310-manylinux2014_x86_64.whl
pip install jax==0.3.25 chex==0.1.6 biopython==1.79
# Use 'Agg' for non-GUI backend
cd /public/home/jiangyw/miniconda3/envs/colab/lib/python3.10/site-packages/colabfold
sed -i -e "s#from matplotlib import pyplot as plt#import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt#g" plot.py
# modify the default params directory
sed -i -e "s#appdirs.user_cache_dir(__package__ or \"colabfold\")#\"${COLABFOLDDIR}/colabfold\"#g" download.py
# remove cache directory
```