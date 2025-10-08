OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];          // data qubits

qreg ancZ[1];       // ancillas measured in Z (for Z-type stabilizers)
qreg flagX[1];      // flags measured in X (paired with Z-type stabs)





// Define a new 2-qubit gate "notnot"
gate notnot c, t {
    h c;
    cx c, t;
    h c;
}



// =========================
// Stabilizer 1 : 
// ancilla = ancZ[0] (Z basis), flag = flagX[0] (X basis)
// =========================
notnot q[0], ancZ[0];
cx flagX[0], ancZ[0];
cx q[1], ancZ[0];
cx q[2], ancZ[0];
cx flagX[0], ancZ[0];
notnot q[3], ancZ[0];




