# MATLAB Electronic Skin Model

This repository contains the MATLAB code used to simulate and analyze the force response of an electronic skin model. The model consists of 677 masses and 256 nodes, which are arranged in a 16x16 grid, simulating tactile sensing for various textures and hardness levels.

## Files

- **main_simulation.m**: Main script to run the simulation for force sensing.
- **classification.m**: Code for training and testing classification models (e.g., SVM).
- **data/**: Folder containing experimental data and sample datasets used in the simulation.

## Setup Instructions

1. Download or clone this repository.
2. Ensure MATLAB is installed and that you have the necessary toolboxes:
   - Statistics and Machine Learning Toolbox (for classification models).
3. Place the data files in the `data/` directory.
4. Open and run the **main_simulation.m** script to simulate the electronic skin's response.

## How to Run the Code

To run the simulation:
1. Open MATLAB and navigate to the directory containing this repository.
2. Run the `main_simulation.m` script.
3. Modify parameters or data paths as needed to customize the simulation.

## License

This repository is licensed under the MIT License.
