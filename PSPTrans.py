import torch

from transformers import T5EncoderModel, T5Tokenizer
from transformers import BertModel, BertTokenizer
from transformers import XLNetModel, XLNetTokenizer
from transformers import AlbertModel, AlbertTokenizer

import re
import gc
import os
import pandas as pd
import requests
from tqdm.auto import tqdm
from pandas.core.frame import DataFrame
import numpy as np
import  joblib
import argparse
import sklearn
from pathlib import Path


def run_seq(args):
    seqs = fasta2seqs(args.seqs_fasta)[1]
    protein_names = fasta2seqs(args.seqs_fasta)[0]
    embeddings = get_embeddings(seqs)
    LLPS_propensity = get_LLPS_propensity(embeddings)
    df_llps = DataFrame(data = [protein_names,seqs,LLPS_propensity],index = ['protein_name','seqs','p_llps'] )
    df2file(df_llps,args.output_filename_llps)
    
    
    if args.output_filename_embedding != '422' :
 
        df_embedding1 = DataFrame(data = [protein_names,seqs],index = ['protein_name','seq'] ).T
        df_embedding2 = DataFrame(data = embeddings)
        df_embedding = pd.concat([df_embedding1,df_embedding2],axis = 1)
        df2file(df_embedding.T,args.output_filename_embedding)
    print('All of the embedding is calculated ,OK')
    
def fasta2seqs(filepath):
    f = open(filepath, 'r')
    seqs = []
    names = []
    regions = []
    for line in f.readlines():
        if line[0]  == '>':
            seqs.append(''.join(regions))
            regions = []
            names.append(line.strip())
        elif line[0] != '>':
            regions.append(line.strip())
    seqs.append(''.join(regions))
    return names, seqs[1:]



def get_embeddings(seqs) :

    def embed_dataset(dataset_seqs, shift_left = 0, shift_right = -1):
        inputs_embedding = []     
        for sample in tqdm(dataset_seqs):
            with torch.no_grad():
                ids = tokenizer.batch_encode_plus([sample], add_special_tokens=True, padding=True, is_split_into_words=True, return_tensors="pt")
                embedding = model(input_ids=ids['input_ids'].to(device))[0]
                inputs_embedding.append(embedding[0].detach().cpu().numpy()[shift_left:shift_right])
        return inputs_embedding


    current_work_dir = os.path.dirname(__file__)
    model_name = ''.join([current_work_dir,'/prot_t5_xl_bfd/'])
    #model_name = '/home/lsf/no/prot_t5_xl_bfd/'
    #model_name  = '/prot_t5_xl_bfd/'
    embeddings = [] 
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False )
    model = T5EncoderModel.from_pretrained(model_name)
    gc.collect()
    #print("Number of the parameters is: " + str(int(sum(p.numel() for p in model.parameters())/1000000)) + " Million")
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model = model.eval()
    if torch.cuda.is_available():
        model = model.half() 
    shift_left = 0
    shift_right = -1 
    for i in range(len(seqs)):
        seq = list(seqs[i])
        print('Get the sequence embedding of ',str(i))
        seq_emd = embed_dataset(seq, shift_left, shift_right)
        seq_emd = torch.tensor(seq_emd)
        embs = torch.flatten(seq_emd, start_dim = 1)
        emb = np.zeros(1024)  
        for j in range(len(seq)):
            emb = emb + np.array(embs[j,:])
        emb = emb/len(seq)
        embeddings.append(emb)
    return embeddings

    
def get_LLPS_propensity(embeddings):
    
    current_work_dir = os.path.dirname(__file__)

    model0 = joblib.load(os.path.join(current_work_dir,'model','transformer0.pkl'))
    model1 = joblib.load(os.path.join(current_work_dir,'model','transformer1.pkl'))
    model2 = joblib.load(os.path.join(current_work_dir,'model','transformer2.pkl'))
    model3 = joblib.load(os.path.join(current_work_dir,'model','transformer3.pkl'))
    model4 = joblib.load(os.path.join(current_work_dir,'model','transformer4.pkl'))  

    model_all = [model0,model1,model2,model3,model4]

    results = np.zeros(len(embeddings))
    for i in range(5):
        clf = model_all[i]
        result = clf.predict(embeddings)
        results = results + np.array(result)

    results = results/5

    return results

def df2file(df,file_name):
    path = ''.join([file_name,'.csv'])
    path = Path(path)
    df.T.to_csv(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="run xxx model on fasta files")
    parser.add_argument('--input file name','-i', dest = 'seqs_fasta',help = ('input fasta file name'), type = str )

    parser.add_argument('--output filename for LLPS propensity','-o_llps', dest = 'output_filename_llps', help = ('output csv file name for llps'), type = str)
    parser.add_argument('--output filename for seq_embedding','-o_embedding', default = '422',dest = 'output_filename_embedding', help = ('output csv file name for embedding'), type = str)

    args = parser.parse_args()
    run_seq(args)
    print('OK')