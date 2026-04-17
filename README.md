#### **Overview

These are set of analysis for a project that is intersection of Molecular Dynamics Simulation, NMR and Machine Learning

Related studies:

1. Discerning intersecting fusion-activation pathways in the Nipah virus using machine learning

2. Machine learning approaches to evaluate correlation patterns in allosteric signaling: A case study of the PDZ2 domain

3. Water Dynamics at Protein–Protein Interfaces: Molecular Dynamics Study of Virus–Host Receptor Complexes

4. Effect of intrinsic and extrinsic factors on the simulated D-band length of type I collagen

This MATLAB function, `bbTrajPos`, is designed to extract and reorganize the 3D coordinates of the four primary backbone atoms—**Nitrogen ($N$)**, **Alpha Carbon ($C_{\alpha}$)**, **Carbon ($C$)**, and **Oxygen ($O$)**—from a molecular dynamics trajectory. 

The primary goal of this data extraction is to facilitate the calculation of protein backbone dihedral angles, specifically $\phi$ (phi) and $\psi$ (psi), which define the local conformation of the protein chain.

---

### **Function Overview**

| Feature | Description |
| :--- | :--- |
| **Input: `traj`** | The raw trajectory data containing atom coordinates over multiple frames. |
| **Input: `atmGrpRes`** | A mapping matrix where the second column identifies which residue each atom belongs to. |
| **Input: `atmNam`** | A cell array containing the specific names of the atoms (e.g., 'N', 'CA', 'C', 'O'). |
| **Output: `bbTrajPos`** | A 3D matrix of size `[nRes, 12, nFrm]`. |

---

### **Detailed Logic Flow**

1.  **Initialization**: 
    The function determines the number of frames (`nFrm`) and the total number of residues (`nRes`). It initializes `bbTrajPos` as a zero-filled 3D array. The second dimension is **12** because it stores the $x, y, z$ coordinates for four specific atoms ($4 \times 3 = 12$).

2.  **Coordinate Mapping**:
    The function iterates through the atom list to find the backbone atoms for each residue. It maps them into the output array in the following order:
    * **Indices 1–3**: Nitrogen ($N$)
    * **Indices 4–6**: Alpha Carbon ($C_{\alpha}$)
    * **Indices 7–9**: Carbonyl Carbon ($C$)
    * **Indices 10–12**: Carbonyl Oxygen ($O$)

3.  **Search Mechanism**:
    * It looks for the 'N' atom to signal the start of a residue's backbone.
    * It assumes the **CA** atom immediately follows the **N** atom.
    * It enters a nested loop to find the **C** atom within the same residue.
    * It assumes the **O** atom immediately follows the **C** atom.

---

### **Key Use Case**
By structuring the data this way, the output `bbTrajPos` allows for vector-based calculations of the dihedral angles. For example:
* **$\phi$ (Phi)** involves the coordinates of $C_{i-1}$, $N_i$, $C_{\alpha i}$, and $C_i$.
* **$\psi$ (Psi)** involves the coordinates of $N_i$, $C_{\alpha i}$, $C_i$, and $N_{i+1}$.

> **Note on Assumptions**: This code assumes a specific ordering in the input files (where $C_{\alpha}$ always follows $N$, and $O$ always follows $C$). If the topology file uses a different naming convention or atom order, the indexing `a+1` may need to be adjusted.
This Python function, `traj_center_of_mass_res_bb_sc`, calculates the **Center of Mass (CoM)** trajectory for three distinct parts of each protein residue: the whole residue, the backbone, and the side chain.

---

## **Mathematical Context**
The Center of Mass ($R_{CoM}$) for a system of particles is calculated as the weighted average of their positions:

$$R_{CoM} = \frac{\sum_{i=1}^{n} m_i r_i}{\sum_{i=1}^{n} m_i}$$

Where:
* $m_i$ is the atomic mass of atom $i$.
* $r_i$ is the position vector $(x, y, z)$ of atom $i$.

---

## **Inputs and Outputs**

### **Input Parameters**
| Parameter | Type | Description |
| :--- | :--- | :--- |
| `traj` | NumPy Array | The trajectory data in the format `[frames, atoms, coordinates]`. |
| `resId` | Array/List | Residue IDs corresponding to each atom in the system. |
| `atmNm` | Array/List | Atom names (e.g., 'N', 'CA', 'CB') for mass identification. |

### **Returned Arrays**
All outputs are 3D NumPy arrays with dimensions **`[residue, coordinates, frames]`**.
1.  **`CoM_res_traj`**: Movement of the center of mass for the **entire residue**.
2.  **`CoM_bb_traj`**: Movement of the center of mass for the **backbone** ($N, C_{\alpha}, C, O$).
3.  **`CoM_sc_traj`**: Movement of the center of mass for the **side chain** (all atoms except backbone).

---

## **Key Logic & Workflow**

1.  **Mass Assignment**:
    The function defines standard atomic weights for Nitrogen (14.01), Carbon (12.01), Oxygen (16.00), and Sulfur (32.60). It identifies the element by checking the first letter of the atom name (`atmNm[a][0]`).

2.  **Backbone Definition**:
    The backbone is strictly defined as the set of atoms: **N, CA, C, and O**.
    * The total backbone weight is pre-calculated as $N + C + O + C_{\alpha}$ (since both $C$ and $C_{\alpha}$ are Carbon).

3.  **Side Chain Handling**:
    * The side chain is defined as everything that is *not* a backbone atom.
    * **Special Case (Glycine)**: If a residue has only 4 atoms (the backbone), the side chain CoM is set to `NaN` (Not a Number) because Glycine effectively has no side chain beyond a Hydrogen (which is often excluded in these trajectories).

4.  **Dimension Swapping**:
    The function reorders the input `traj` from `[frames, atoms, coords]` to `[atoms, coords, frames]`. This allows the code to perform vectorized summations over all frames simultaneously for each residue, significantly improving calculation speed.

---

##  **Technical Observations**
* **Legacy Code**: Similar to previous scripts, this uses `xrange`, indicating it was written for **Python 2.x**. To run this in Python 3, simply change `xrange` to `range`.
* **Dependency**: The function relies on a helper function `ismember` (likely a custom utility or a port of the MATLAB function) to identify which indices in the atom list match the backbone atom names.
* **Memory Management**: The `del` statements at the end of the loop help clear temporary residue-specific arrays, which is a good practice when dealing with large MD trajectories that consume significant RAM.
This analysis computes the fundamental backbone conformational angles—**$\phi$ (phi)**, **$\psi$ (psi)**, and ****$\omega$ (omega)****—for a protein from a Molecular Dynamics (MD) trajectory. These angles are the primary descriptors of a protein's secondary structure and are used to construct **Ramachandran plots**.

---

## 🧬 Overview of Backbone Dihedrals

The function `diHedCal` iterates through each frame of a simulation and each residue of the protein to calculate the torsion angles formed by the peptide backbone atoms ($N$, $C_\alpha$, and $C$).



### 1. Dihedral Angle Definitions
The script extracts specific coordinate sets to calculate the rotation around backbone bonds:

* **$\phi$ (Phi):** Rotation around the $N — C_\alpha$ bond. It requires the coordinates of $C_{i-1}, N_i, C\alpha_i, C_i$.
* **$\psi$ (Psi):** Rotation around the $C_\alpha — C$ bond. It requires the coordinates of $N_i, C\alpha_i, C_i, N_{i+1}$.
* **$\omega$ (Omega):** Rotation around the peptide bond ($C — N$). It requires the coordinates of $C\alpha_i, C_i, N_{i+1}, C\alpha_{i+1}$.

### 2. Geometric Implementation
The code uses a robust vector-based approach to calculate the torsion angle between two planes:
1.  **Vector Definition:** It defines three vectors ($a, b, c$) connecting the four atoms.
2.  **Cross Products:** It calculates the normals to the planes formed by these vectors.
3.  **Trigonometry:** It uses the dot products and cross products of these normals to determine the angle, ensuring the sign (direction) of the rotation is preserved.
4.  **Conversion:** Angles are converted from radians to **degrees** ($180/\pi$).

---

## 📥 Input Requirements

| Variable | Description | Structure/Assumptions |
| :--- | :--- | :--- |
| `bbTrajPos` | Backbone trajectory positions | 3D array: `(Residue, Atom_Coords, Frame)` |
| `atmGrpRes` | Topological mapping | Maps atoms to specific residue numbers |

**Coordinate Mapping Assumption:**
The code assumes a specific ordering of backbone atoms in the input array:
* **Indices 1–3:** Nitrogen ($N$)
* **Indices 4–6:** Alpha Carbon ($C_\alpha$)
* **Indices 7–9:** Carbonyl Carbon ($C$)

---

## 📤 Analysis Outputs

The function returns three matrices, each of size `(nResidues x nFrames)`:

* **`phi`**: The distribution of $\phi$ angles over time.
* **`psi`**: The distribution of $\psi$ angles over time.
* **`omega`**: The peptide bond dihedrals (usually near **180°** for trans or **0°** for cis configurations).

---

## ⚠️ Boundary Conditions
The analysis automatically skips the **first** and **last** residues of the chain. This is because:
* The first residue lacks a preceding Carbon ($C_{i-1}$) to define $\phi$.
* The last residue lacks a succeeding Nitrogen ($N_{i+1}$) to define $\psi$ and $\omega$.

> **Note:** The `localNewAngle` and `localCalculateTorsionAngle` sub-functions are optimized for numerical stability, handling cases where vectors might be collinear or the angle approaches the $180°$ limit.
# Methyl Group Dipole Autocorrelation Analysis

This notebook performs a specialized analysis of protein dynamics from **Molecular Dynamics (MD) simulations**. It calculates the time-autocorrelation functions of methyl group ($C-CH_3$) dipole moments, which are essential for estimating **NMR $S^2$ order parameters**. These parameters provide a quantitative measure of the spatial restriction and flexibility of side-chain motion.

## 🔬 Analysis Workflow

### 1. Data Initialization and Preprocessing
The script loads molecular trajectory data (atomic coordinates) and topological information (residue names/IDs and atom names). 
* **Input Formats:** Supports both Python (`.npz`) and MATLAB (`.mat`) files via `h5py` and `scipy.io`.
* **Preprocessing:** Ensures the trajectory tensor is correctly shaped as `(frames, atoms, coordinates)`.

### 2. Residue-Specific Methyl Mapping
Because different amino acids have unique side-chain architectures, the code implements logic to map the specific atom pairs defining the $C-CH_3$ bond vector for the following residues:
* **ALA:** $C_\alpha \rightarrow C_\beta$
* **THR:** $C_\beta \rightarrow C_{\gamma2}$
* **ILE:** $C_\beta \rightarrow C_{\gamma2}$ and $C_{\gamma1} \rightarrow C_\delta$
* **VAL:** $C_\beta \rightarrow C_{\gamma1}$ and $C_\beta \rightarrow C_{\gamma2}$
* **LEU:** $C_\gamma \rightarrow C_{\delta1}$ and $C_\gamma \rightarrow C_{\delta2}$

### 3. Dipole Moment Vectorization
For every frame in the simulation, the script calculates the displacement vector for each identified methyl group. These are then converted into **unit vectors** to isolate rotational dynamics from bond-length fluctuations.

### 4. Autocorrelation Calculation ($P_2$ Correlation)
The core of the analysis involves calculating the **Second-order Legendre Polynomial ($P_2$)** of the unit vector's orientation over time:

$$C(\Delta t) = \langle P_2(\vec{\mu}(t) \cdot \vec{\mu}(t+\Delta t)) \rangle = \langle \frac{3 \cos^2 \theta - 1}{2} \rangle$$

Where $\theta$ is the angle between the dipole vector at time $t$ and $t + \Delta t$.

### 5. Performance Optimization
To handle high-frequency trajectory data efficiently, the notebook employs:
* **Vectorization:** Utilizing `numpy` for high-speed array operations.
* **JIT Compilation:** Leveraging `numba` to compile Python loops into optimized machine code, significantly accelerating the calculation of the correlation decay.

---

## 📊 Key Outputs

| File | Description |
| :--- | :--- |
| `dipoles.mat` | Intermediate storage of calculated dipole vectors for all frames and residues. |
| `Cor.mat` | The final autocorrelation functions ($C$) for each methyl group across specified time delays. |
| `nDataPnts` | The count of data points used for each $\Delta t$, used for assessing statistical uncertainty. |

---

## 🛠 Prerequisites

* **Python Libraries:** `numpy`, `scipy`, `h5py`, `numba`.
* **Input Data:** An MD trajectory file (coordinates) and a corresponding topology (residue/atom identifiers).
These MATLAB functions form a coordinated **PDB Trajectory Parser**, specifically designed to handle multi-model PDB files (typically generated by simulation packages like Gromacs). The suite extracts both structural topology and time-resolved atomic coordinates.

---

### 🔍 Code Review & Logic Check

While the logic is generally sound for standard PDB formats, I identified one **critical bug** in the orchestrator function (`PdbTrajRead`) and a few performance considerations:

1.  **The `nargin` Logic Error:** In your current `PdbTrajRead`, if the user does *not* provide `nFrm`, the code calculates the number of frames but **skips** the actual coordinate reading.
    * *Fix:* The call to `PdbTrajReadPos` should happen regardless of whether `nFrm` was provided or calculated.
2.  **Fixed-Width Parsing:** Your use of `fscanf` with specific offsets (e.g., `%*9c%3c`) correctly follows the [PDB fixed-column standard](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html). However, it assumes the input file strictly adheres to these column positions.
3.  **Memory Efficiency:** The code reads the file up to three times (once for topology, once for counting, once for positions). For very large MD trajectories (GBs), this will be slow. It is usually faster to count frames and read topology in a single pass.
4.  **Array Dimensions:** Your trajectory is saved as `(frames, coordinates, atoms)`. This is a $3D$ tensor of size $N_{frames} \times 3 \times N_{atoms}$. This is efficient for looking at a single frame, but if you later want to calculate the autocorrelation of a single atom, you will need to squeeze/permute the dimensions.

---

### 📝 Analysis Description: PDB Trajectory Parser

This suite of functions automates the conversion of raw PDB trajectory data into MATLAB-compatible data structures for biophysical analysis.

#### 1. Topology Extraction (`PdbFirstFrame`)
The parser initializes by scanning the first model of the file. It extracts:
* **Atom Types:** Distinguishes between standard protein residues (`ATOM`) and heteroatoms/ligands (`HETATM`).
* **Naming Conventions:** Captures specific atom names (e.g., `CA`, `CB`) and residue names (e.g., `ALA`, `ILE`) to allow for residue-specific filtering in later analysis.
* **Connectivity Mapping:** Assigns residue IDs to each atom, creating a map of the molecular structure.

#### 2. Frame Synchronization (`PdbNumFrms`)
To ensure memory is allocated efficiently, the tool scans the trajectory for `ENDMDL` tags. This determines the total number of snapshots (frames) available in the simulation, allowing the system to pre-allocate large coordinate tensors.

#### 3. Coordinate Harvesting (`PdbTrajReadPos`)
The core function performs a deep scan of the trajectory, extracting:
* **Temporal Data:** Reads the `TITLE` or header lines to extract the simulation timestamp for each frame.
* **Spatial Data:** Parses the X, Y, and Z coordinates for every atom in every frame, storing them in a $3D$ matrix.

#### 4. Data Structure Output
The final output is a structured set of arrays:
* `traj`: A $3D$ array where `traj(f, :, a)` represents the $(x, y, z)$ coordinates of atom `a` at frame `f`.
* `atmNam` / `resNam`: Cell arrays containing the chemical identity of the parsed system.
* `time`: A vector containing the time evolution of the simulation.
This script, `PdbRead`, is a utility designed to parse molecular dynamics trajectory files (like `.pdb` or `.dcd`) using the **MDAnalysis** library. It extracts structural metadata and atomic coordinates, converting them into a format ready for numerical analysis in **NumPy**.

---

## **Overview**
The code consists of two primary functions that work together to flatten a complex molecular universe into organized arrays. This is particularly useful for preparing data for custom calculations (like the backbone position analysis we discussed earlier).

### **Dependencies**
* **MDAnalysis**: To handle the heavy lifting of reading molecular formats.
* **NumPy**: To store and return data as high-performance arrays.

---

## **Function Details**

### **1. `save_atm_res_seg(trajFileName)`**
This function extracts the "identity" of every atom in the system. 

* **Logic**: It iterates through the atoms, converts the MDAnalysis atom object to a string, and performs manual string slicing/splitting to pull out specific data points.
* **Captured Data**:
    * `atmGrp`: Atom type (e.g., ATOM or HETATM).
    * `atmId`: The unique numerical ID of the atom.
    * `atmNm`: The chemical name (e.g., 'N', 'CA').
    * `resNm`: The residue name (e.g., 'ALA', 'GLY').
    * `resId`: The residue sequence number.
    * `segId`: The chain or segment identifier.

### **2. `save_traj_pos(trajFileName)`**
This function focuses on the "movement" of the atoms over time.

* **Logic**: It initializes a 3D NumPy array of size `[frames, atoms, coordinates]`. It then loops through every frame in the trajectory, updating the `Universe` state and copying the atomic positions into the array.
* **Output**: A high-density array where `traj[0][5][2]` would represent the **Z-coordinate** of the **6th atom** in the **1st frame**.

---

## **Summary of Outputs**

| Array | Data Type | Description |
| :--- | :--- | :--- |
| `atmGrp` | String | Atom/Hetatm classification |
| `atmId` | Integer | Unique Atom ID |
| `atmNm` | String | Atom Name (e.g., Nitrogen) |
| `resNm` | String | Residue Name (e.g., Amino Acid) |
| `resId` | Integer | Residue sequence number |
| `segId` | String | Chain/Segment ID |
| `traj` | Float (3D) | XYZ coordinates across all frames |

---

## **Technical Note**
The script uses some "legacy" Python 2 patterns (like `xrange` and `__future__` imports). If you are running this in a modern **Python 3.x** environment, `xrange` should be replaced with `range` for compatibility. Additionally, the manual string parsing of `atmsAndResStr` is a creative way to bypass deep object inspection, though it relies heavily on the default string formatting of the MDAnalysis version being used.
