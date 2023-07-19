# How to run the PSPTrans

## Step 0: 
Due to storage space limitations, please download "prot_t5_xl_bfd" [click here]("http://bio-comp.ucas.ac.cn/onlineserver/PSPtrans/file_down_prot_t5_xl_bfd")


## Step 1: 
Make sure you have Python3 in your computer.
How to check and/or install Python3 version - https://realpython.com/installing-python/
Official Python website - https://www.python.org/


## Step2 : 
Set up a virtual enviroment with python=3.9.7  
```conda create -n PSPTrans python=3.9.7```


## Step 3:
Activate the virtual environment  
```conda activate PSPTrans```


## Step 4: 
Download the libraries required to run PSPTrans model in your machine  
```pip3 install -r requirements.txt```

and you need to install the pytorch (CPU or GPU)
Official Python website "https://pytorch.org/get-started/locally/"

for the cpu version  
```pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu```

for the GPU version:
you need to know your computer CUDA version,  the command is  
```nvidia-smi```
and then you can install the corresponding pytorch , for mine(CUDA 11.7)  
```pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu117```


## step 5: 
Run PSPTrans on the example fasta file.
if you just want to know the protensity for LLPS, the command is :   
```python3  llps_lsf.py  -i test.fasta  -o_llps output_file```
if you also want to save the embedding for sequences, the command is :  
```python3  llps_lsf.py  -i test.fasta  -o_llps output_file -o_embedding embedding```

And you can see the result in output_file.csv and embedding.csv(if existed)


## Step 6:
To run PSPTrans on your own fasta files
Replace "test.fasta" with your fasta file in Step 5. You can only run one file at a time.


## Step 7: 
After running, exit or delete your virtual environment.
exit virtual environment  
```deactivate```
delete virtual environment  
```conda env remove -n PSPTrans```


test :
python3  llps_lsf.py  -i lsf_test.fasta  -o_llps lsf0617
