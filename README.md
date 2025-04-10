This repository contains MATLAB code for the research paper:

"Recovering Noisy Measurements Over Inland Water Bodies by Regenerating L1B SAR Altimetry Waveforms Using Segment-Weighted Fully-Focused SAR (swFF-SAR) Processing."

The method improves SAR altimetry waveform quality over small inland water bodies (e.g., lakes and rivers) by reducing the impact of unwanted backscattered signals. It regenerates corrupted L1B waveforms using a segment-weighted FF-SAR approach, leading to more accurate water level time series.

 Requirements
L1A SAR Altimetry Data
The code requires L1A data products, which are large files (~1.8 GB per cycle). These can be downloaded from the Copernicus Data Space:
 https://browser.dataspace.copernicus.eu/

MATLAB
Developed and tested in MATLAB (recommended version: R2021a or later). 

 Included Test Data
For testing purposes, this repository includes data for one case study:

Lake Rathbun

In-situ water level measurements

L2 correction files

LUT files for the SAMOSA model

This allows you to run the full pipeline and test the regeneration algorithm just you need download L1A datasets for this case study

🚀 Running the Code
Clone the repository.

Open the main script (SWFFSAR_Main.m) in MATLAB.

Follow inline comments for input settings and execution flow.

💡 The code is heavily commented to help you to track each step of the processing chain.

📫 Contact
For questions or collaboration, feel free to reach out via the Issues tab or open a pull request.

