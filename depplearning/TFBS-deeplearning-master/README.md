# Deeping learning model for predicting tissue-specific TFBS based on ATAC-seq

Purpose:
Deep learning methods could be applied to predict TFBS in different tissues by combining ATAC-seq  with tsCUT&Tag.In order to solve the problem that it is difficult for tissues to obtain high-quality protoplasts with high transformation efficiency.

Usage:
1. Step 1: preprocess the input sequence.
2. Step 2: train the model of LSTM, TCN and SVM. ① The model training selects the overlapping binding sites of transcription factors and the open chromatin region of ATAC-seq , upstream and downstream 100bp of the peak summit is positive data, and the open chromatin region of ATAC-seq without TFBS is negative data set. 
3. Step 3: Predict each TF. The predicted inputs: The predicted inputs is the open chromatin region sequence of ATAC-seq in csv file (see example ./data/greenleaf_ATAC_knox6predicted_seq.CSV and ./data/tassel_ATAC_tb1predicted_seq.CSV). The prediction model reads a window size of 200bp with a displacement of 100. The outputs： The outputs file (see example ./result/knox6LSTM.csv) contain one columns, the column is the TFBS, 0 is not bound. 
4. Step 4: deal with the outputs of deep learning. The inputs file see example ./result/knox6LSTM.csv. The outputs is the TF binding sites in genome (see example ./result/knox6LSTM_bindingsite.txt).
