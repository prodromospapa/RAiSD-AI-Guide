#!/bin/bash

# Downloading training raw data
wget https://zenodo.org/records/18841026/files/Mild-bottleneck-training.tar.xz; tar -xvJf Mild-bottleneck-training.tar.xz

# Preparing neutral training data for RAiSD-AI
RAiSD-AI -n TrainingData -I Mild-bottleneck-training/neutral.ms -w 128 -L 100000 -its 50000 -op IMG-GEN -icl neutralTR -bin -typ 1 -frm

# Preparing sweep training data for RAiSD-AI
RAiSD-AI -n TrainingData -I Mild-bottleneck-training/selsweep.ms -w 128 -L 100000 -its 50000 -op IMG-GEN -icl sweepTR -bin -typ 1

# Training a FASTER-NN model with RAiSD-AI
RAiSD-AI -n Model -I RAiSD_Images.TrainingData -op MDL-GEN -e 100 -f

# Downloading test raw data
wget https://zenodo.org/records/18841026/files/Mild-bottleneck-testing.tar.xz; tar -xvJf Mild-bottleneck-testing.tar.xz

# Preparing neutral test data for RAiSD-AI
RAiSD-AI -n TestData -I Mild-bottleneck-testing/neutral.ms -w 128 -L 100000 -its 50000 -op IMG-GEN -icl neutralTE -bin -typ 1 -frm

# Preparing sweep test data for RAiSD-AI
RAiSD-AI -n TestData -I Mild-bottleneck-testing/selsweep.ms -w 128 -L 100000 -its 50000 -op IMG-GEN -icl sweepTE -bin -typ 1

# Classifying test data with RAiSD-AI
RAiSD-AI -n Test -mdl RAiSD_Model.Model -I RAiSD_Images.TestData -op MDL-TST -clp 2 neutralTR=neutralTE sweepTR=sweepTE

# Scanning neutral test data with RAiSD-AI
RAiSD-AI -n neutral -mdl RAiSD_Model.Model -op SWP-SCN -I Mild-bottleneck-testing/neutral.ms -L 100000 -k 0.05 -G 100 -pci 1 1

# Getting threshold from neutral output info file
neutralRunName=neutral
fprThresholdMUVAR=$(grep " muVar " RAiSD_Info.$neutralRunName | grep min | awk -F: '{print $2}' | awk '{print $1}')
fprThresholdMUSFS=$(grep " muSFS" RAiSD_Info.$neutralRunName | grep min | awk -F: '{print $2}' | awk '{print $1}')
fprThresholdMULD=$(grep " muLD " RAiSD_Info.$neutralRunName | grep min | awk -F: '{print $2}' | awk '{print $1}')
fprThresholdMU=$(grep " mu " RAiSD_Info.$neutralRunName | grep min | awk -F: '{print $2}' | awk '{print $1}')
fprThresholdPCL0=$(grep " sweepTR" RAiSD_Info.$neutralRunName | grep min | awk -F: '{print $2}' | awk '{print $1}')
fprThresholdPCL1=$(grep " muvar^sweepTR" RAiSD_Info.$neutralRunName | grep min | awk -F: '{print $2}' | awk '{print $1}')

# Scanning sweep test data with RAiSD-AI
RAiSD-AI -n sweep -mdl RAiSD_Model.Model -op SWP-SCN -I Mild-bottleneck-testing/selsweep.ms -L 100000 -T 50000 -d 2500 -G 100 -pci 1 1 -l 6 var=$fprThresholdMUVAR sfs=$fprThresholdMUSFS ld=$fprThresholdMULD mu=$fprThresholdMU pcl0=$fprThresholdPCL0 pcl1=$fprThresholdPCL1
