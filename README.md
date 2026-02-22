## **OTFS Channel Estimation: Compressed Sensing Methods**


### **Introduction**

This repository summarizes and implements **Compressed Sensing (CS)-based Orthogonal Time Frequency Space (OTFS) channel estimation methods**. Built in MATLAB, this project provides an end-to-end pipeline for channel estimation under the **integer delay-Doppler taps** assumption.

### **Key Features**

* **Complete MATLAB Pipeline**: Covers OTFS modulation, embedded pilot insertion, doubly dispersive channel transmission.
* **Integer Delay-Doppler Model**: Accurately simulates high-mobility environments under the assumption that path delays and Doppler shifts align with the system's discrete grid resolutions.
* **Compressive Sensing Recovery**: Integrates classic sparse recovery algorithms, specifically the orthogonal matching pursuit (OMP), to efficiently reconstruct channel vectors by exploiting DD domain sparsity.
* **Modular Codebase**: Highly decoupled components allowing researchers to easily test advanced CS variants.

### **Sparse Channel Estimation**


To estimate the channel without corrupting user data, the system utilizes an **Embedded Pilot** surrounded by a **Guard Zone (GZ)** of null symbols.

Due to the massive spectrum overhead and "dimensionality curse" associated with traditional threshold methods in massive MIMO or 6G environments, this pipeline treats estimation as a sparse signal recovery (SSR) problem. The received pilot signal is mapped using a sensing matrix $\mathbf{\Phi}$:

$$\mathbf{y}_p = \mathbf{\Phi} \mathbf{h} + \mathbf{v}.$$

We solve the NP-hard $\ell_0$ norm minimization problem using:

1. **Orthogonal matching pursuit (OMP)**,
2. ...


### 📂 **Repository Structure**

```text
Compressed_Sensing_OTFS_Channel_Estimation/
├── main.m                            # Main entry point for the simulation
├── simulateOTFSChannelEstimation.m   # Channel Estimation pipeline
├── helperOTFSmod.m                   # OTFS modulation
├── Algorithm_OMP.m                   # The OMP algorithm
└── README.md                         # Project documentation
```

### **📚 Citation**

If you find it useful, please consider cite our works:

```code
@ARTICLE{COMML_SL0,
  author={Wang, Xin and Zheng, Can and Hu, Pengjiang and Yang, Junan and Kang, Chung G.},
  journal={IEEE Communications Letters}, 
  title={A Sparsity-Agnostic SL0 Channel Estimation Approach for OTFS Systems}, 
  year={2025},
  volume={29},
  number={5},
  pages={1097-1101}}
```

```code
@INPROCEEDINGS{VTC25fall_ADMM,
  author={Zheng, Can and Wang, Xin and Kang, Chung G.},
  booktitle={Proc. IEEE Vehicular Technology Conference (VTC2025-Fall)}, 
  title={ADMM-Based Delay-Doppler Domain Channel Estimation for OTFS Systems}, 
  year={2025},
  pages={1-5}}
```

```code
@INPROCEEDINGS{APCC_GP,
  author={Zheng, Can and Kang, Chung G.},
  booktitle={Proc. Asia-Pacific Conference on Communications (APCC)}, 
  title={Delay-Doppler Domain Channel Estimation for OTFS System Using the Gradient Projection Method}, 
  year={2025},
  pages={1-4}}
```