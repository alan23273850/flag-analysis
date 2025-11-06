# SAT-based flag qubit circuit verification tool 
## Installation 
Copy this to command line 
```bash
git clone https://github.com/boriswu412/flag-analysis.git
```
Copy this to install the requirement 

```bash
pip install -r requirements.txt
```

## Usage 
In the `config.txt`  there are to rows `qasm_path` and `stab_txt_path`. 

### Qasm file 
#### Qubits
There are five kind of qregs `q` for data qubits, `ancX`, `ancZ` for ancillas which are prepared in X and Z basis, `flagX` ,`flagZ` for flag qubits prepared in X and Z basis. 

#### Gates 
We support three kinds of gates: CNOT, CZ, and NOTNOT(two qubit gate that both qubits aretargets )
Some times we will have hadamard gates on ancillas and flags, instead we can change the basis between X-basis and Z-basis.
We dont write measurent in qasm files. If some qubits are reused in some flag circuits we will treat them as different qubits. 


### Stab.txt
This is for the parity check matrix for the code, first x then z
ps. the order of the rows should match the order of ancillas 



### Run 
To run the program use:
```bash
python flag_analysis.py        
```
#### Output
The whole verfication steps includes four steps:
Step 1 : check if the circuit extracts right syndromes when there are not faults in the circuit.

Step 2:  check every gate in the circuit whether a fault on the gate will cause a high-weight error. If will, store the index of the gate to bad_location  

Step 3: check if there is a fault on one of the gates in bad_location and cause high-weight error, on of the flags should raised.

Step 4: check after a round of flag circuit and a round of circuit without flag the generalised syndrome is unique up to degeneracy.

Each step should be "Success".

For Step 3 and Step 4 is result is fail, there will be a counter example. The example will tell which variable is True and the rest are False.

E.g. for Step 3:
faulty_gate4_z1 = True
faulty_gate4_x0 = True

gate 4 in the flag circuit is a cx gate on qubits [2,7] 
This means when a fault on gate 4 with error X on qubit 2 and Z on qubit 7 will cause a high-weight error but no flag is raised.

E.g. for Step 4:
faulty_gate19_z1_p2 = True
faulty_gate14_z1_p1 = True
faulty_gate19_z0_p2 = True
faulty_gate19_x0_p2 = True
faulty_gate14_x0_p1 = True

gate 19 in the flag circuit is a  cx gate on qubits [8,1] 
gate 14 in the flag circuit is a  cx gate on qubits [7,2]
This means when a fault on gate 19 with error Y on qubit 8 and Z on qubit 1 and another fault on gate 19 with error X on qubit 7 and Z on qubit 2 will have the same generalized syndrome. 
# SAT‑based Flag Qubit Circuit Verification Tool

A tool to verify flag‑based stabilizer extraction circuits using SAT checks.

---

## Installation

Clone the repository and install requirements:

```bash
git clone https://github.com/boriswu412/flag-analysis.git
cd flag-analysis
pip install -r requirements.txt
```

(Optional) create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
```

> Requires Python 3.9+ (tested with 3.10).

---

## Usage

Edit `config.txt` with two rows:
```
qasm_path=path/to/circuit.qasm
stab_txt_path=path/to/stab.txt
```

Then run:
```bash
python flag_analysis.py
```

### QASM file
**Qubits:**
- `q` — data qubits
- `ancX`, `ancZ` — ancilla qubits in X or Z basis
- `flagX`, `flagZ` — flag qubits in X or Z basis

**Gates:**
- Supported: `cx` (CNOT), `cz`, and `notnot` (two‑qubit X‑X gate)
- Hadamard gates on ancillas or flags are replaced by basis changes (X ↔ Z)
- Measurements are omitted from QASM
- Reused qubits in flag circuits are treated as distinct logical qubits

### `stab.txt`
Defines the parity‑check matrix (first X, then Z rows). Row order must match ancilla order.

---

## Running and Output

The verification includes **four steps**:

1. **Syndrome correctness:** verify expected syndromes with no faults.
2. **Fault detection:** find gates that cause high‑weight errors → stored in `bad_location`.
3. **Flag check:** verify faults on bad gates trigger at least one flag.
4. **Syndrome uniqueness:** check combined circuits produce unique syndromes.

Each step should output **“Success.”**

If Step 3 or 4 fails, a counterexample is printed showing which variables are `True`.

Example (Step 3):
```
faulty_gate4_z1 = True
faulty_gate4_x0 = True
```
Gate 4 (`cx [2,7]`): fault with X on 2 and Z on 7 causes a high‑weight error but no flag.

Example (Step 4):
```
faulty_gate19_z1_p2 = True
faulty_gate14_z1_p1 = True
faulty_gate19_z0_p2 = True
faulty_gate19_x0_p2 = True
faulty_gate14_x0_p1 = True
```
Gates 19 (`cx [8,1]`) YZ error  and 14 (`cx [7,2]`) XZ error: these faults yield identical generalized syndromes, violating uniqueness.

