# 260310
# python main.py --input_file $input_file --output_file $output_file --layers $layers \
# --setX $setX --ydformat_folder $ydformat_folder

import scanpy 
import anndata 
import pandas 
import numpy 
import scipy

import sys
import os

import argparse

parser = argparse.ArgumentParser(description='Process single-cell data with specified parameters')
parser.add_argument('-i', '--input_file', 
                   type=str,
                   default='/data/work/Dataget/output/GM.h5ad',
                   help='Input h5ad file path [default: %(default)s]')
parser.add_argument('-o', '--output_file',
                   type=str,
                   default='/data/work/Dataget/output/GM',
                   help='Output directory path [default: %(default)s]')
parser.add_argument('--layers',
                   type=str,
                   default='counts,splice,unsplice',
                   help='Comma-separated list of layers to include [default: %(default)s]')
parser.add_argument('--setX',
                   type=str,
                   default='RNA',
                   help='Name of the assay to set as X [default: %(default)s]')
parser.add_argument('--ydformat_folder',
                   type=str,
                   default='/data/work/ydFormat',
                   help='Path to ydFormat folder [default: %(default)s]')
args = parse_args()
input_file=args.input_file
output_file=args.output_file
layers=args.layers.split(',')
setX=args.setX
ydformat_folder=args.ydformat_folder


sys.path.append(ydformat_folder)
import ydAnndata

# 判断文件后缀是否为.h5ad
if input_file.endswith('.h5ad'):
    print('Run .h5ad convert to ydfolder...')
    adata = scanpy.read_h5ad(input_file)
    print(adata)
    ydAnndata.write_ydfolder(
        adata, 
        path=output_file,
        layers=layers,
        verbose=True
    )
else:
    print('Run ydfolder convert to .h5ad')
    adata = ydAnndata.read_ydfolder(
        path=input_file,
        setX=setX,              # 指定主矩阵
        layers=layers,
        verbose=True
    )
    print(adata)
    if '_index' in adata.obs.columns:
        del adata.obs['_index']
    adata.write_h5ad(output_file, compression='gzip')