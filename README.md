INTNet
=====
Multi-task multilabel deep neural networks for identification and classification of integrons.

INTNet is designed to identify and classify integron integrases, predict bacterial hosts, and associate ARGs (multi-label) with integrons. This versatile tool supports a range of input types including:
* Long Amino Acid Sequences (Full Length/Contigs)
* Long Nucleotide Sequences
* Short Amino Acid Reads (30-50 aa)
* Short Nucleotide Reads (100-150 nt)

All inputs should be in FASTA format.

**INTNet Components**\
INTNet comprises two specialized models to accommodate different read lengths:
* **INTNet-s**: Optimized for short reads, enhancing prediction accuracy for sequences ranging from 30 to 50 amino acids or 100 to 150 nucleotides.
* **INTNet-l**: Tailored for long sequences, ensuring robust predictions for full-length contigs or long nucleotide sequences.

![alt text](https://github.com/patience111/INTNet/blob/master/pics/INTNet_workflow.jpg)

Installation
------------
clone the program to your local machine\
git clone https://github.com/patience111/INTNet


**1. Setting up environment**


**1.1 Installation with conda**


1.1.1 For **CPU** inference, you could install the program with conda YAML file in the installation directory with the following commands:

```
cd ./installation 
conda env create -f INTNet-CPU.yml -n INTNet-cpu
conda activate INTNet-cpu
```

(This was tested on Ubuntu 16.04, 20.04; Windows 10, macOS(14.1.1))\
 ![alt text](https://github.com/patience111/INTNet/blob/master/pics/test_cpu.jpg)

 1.1.2 For **GPU** inference, you could install the program with conda YAML file in the installation directory with the following commands:</br>
```
cd ./installation
conda env create -f INTNet-GPU.yml -n INTNet-gpu
conda activate INTNet-gpu
```
(This was tested on Ubuntu 16.04, cuda 10.1, Driver Version: 430.64)\
    ![alt text](https://github.com/patience111/INTNet/blob/master/pics/test_gpu_1.jpg)
    ![alt text](https://github.com/patience111/INTNet/blob/master/pics/test_gpu_2.jpg)


**1.2 Or, if you prefer installing dependencies manually, you might find this information useful:**\
The program was tested with the following package version, you can install exactly the same version or other compatible versions.
```
Biopython:  1.79
tensorflow:  2.2.0 
cuda: 10.2 (for GPU using)
cudnn: 7.6.5.32 (for GPU using)
numpy: 1.18.5
scikit-learn: 0.24.1
tqdm: 4.56.0
```
**2. Getting trained models**

```   
cd ./models
bash get-models.sh
```
Quickstart Guide
----------------
***for long sequences***

```
python intnet.py --input input_path_data  --type aa/nt --model argnet-l  --outname output_file_name
```
***for short reads***

```
python intnet.py --input input_path_data --type aa/nt --model argnet-s --outname output_file_name
```

**general options:**

     --input/-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the test file as input </br>
     --type/-t &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;molecular type of your test data (aa for amino acid, nt for nucleotide)</br>
     --model/-m&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the model you assign to make the prediction (argnet-l for long sequences, argnet-s for short reads) </br>
     --outname/-on&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the output file name </br>
