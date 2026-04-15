```shell
# python3.8 in conda env
git clone https://github.com/STOmics/STCellbin.git
mamba create -n STCellbin python=3.8 -y
conda activate STCellbin
mamba install pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 -c pytorch -y
cd STCellbin-main
pip install -r requirements.txt # install
```