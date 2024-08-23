# Deep Learning-based Channel Predictor for RIS-assisted NOMA Communication Systems

## Overview

This repository contains the source code and associated materials for my final graduation project, titled **"Deep Learning-based Channel Predictor for RIS-assisted NOMA Communication Systems"**, presented as part of the requirements for obtaining a degree in Electronics and Computer Engineering at the Polytechnic School of Engineering at Universidade Federal do Rio de Janeiro.

## Abstract

The transition from fifth generation (5G) to 5G Advanced and forthcoming sixth generation (6G) networks necessitates the development of innovative approaches to enhance data rates, minimize latency, and ensure reliable connectivity. This project analyzes the potential of non-orthogonal multiple access (NOMA) and reconfigurable intelligent surfaces (RIS) as pivotal technologies to overcome these challenges. With NOMA increasing spectral efficiency by enabling resource sharing among multiple users and RIS optimizing signal propagation to enhance quality and coverage with lower energy consumption.

To address the complexities of modern networks, artificial intelligence (AI) solutions, particularly deep learning (DL) models such as convolutional neural networks (CNN) and long short term memory (LSTM) networks, are being investigated for their potential in robust feature extraction and sequential data processing in channel estimation. The objective of integrating NOMA and RIS with CNN-LSTM models is to address the high complexity and non-linearity of channel estimation in this communication system in a creative and less complex way.

Our principal contribution is a proposed new DL model, which incorporates additional layers that facilitate greater granularity, enabling more precise magnitude and phase prediction. Moreover, the project employs simulations that approximate real-world scenarios to evaluate the performance of the communication system, emphasizing the bit-error rate (BER) across varying signal-to-noise ratios (SNR) as a figure of merit. Finally, in comparison to the state-of-the-art, our proposed model demonstrates a reduction in average inference time by 17% and a reduction in the number of training parameters by more than 35%.

**Key-words**: NOMA, SIC, RIS, channel estimation, deep learning, CNN, LSTM.

## Project Structure

- **MATLAB/**: Contains the source code for the simulations on MATLAB.
- **Deep Learning Model/**: Jupyter notebooks for exploratory data analysis and model development.

## Getting Started

### Prerequisites

- Python 3.x
- TensorFlow 
- MATLAB

### Installation

Clone this repository:

```bash
git clone https://github.com/eduardo-henriques/Deep-Learning-based-Channel-Predictor-for-RIS-assisted-NOMA-Communication-Systems.git
cd Deep-Learning-based-Channel-Predictor-for-RIS-assisted-NOMA-Communication-Systems
```

## PDF of the Graduation Project

You can download the full text of my graduation project [here]().

## Authors

- **Eduardo Francisco Silva de Lima Henriques** - [LinkedIn](https://www.linkedin.com/in/eduardo-henriques-103b69208/)

## Acknowledgments

I would like to thank my advisors, Prof. Paulo S. R. Diniz and Prof. Rafael da Silva Chaves, for their guidance and support throughout this project. I also express my gratitude to my family and friends for their encouragement during my studies.
