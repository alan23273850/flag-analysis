OPENQASM 2.0;
include "qelib1.inc";

qreg q[15];          // data qubits

qreg ancZ[1];       // ancillas measured in Z (for Z-type stabilizers)
qreg flagX[1];      // flags measured in X (paired with Z-type stabs)



// =========================
// Stabilizer 5 (Z syndrome):
// ancilla = ancZ[0] (Z basis), flag = flagX[0] (X basis)
// =========================
cx  q[14], ancZ[0];
cx  flagX[0], ancZ[0];
cx  q[12], ancZ[0];
cx  q[13], ancZ[0];
cx  q[10], ancZ[0];
cx  q[11], ancZ[0];
cx  q[9],  ancZ[0];
cx  q[8],  ancZ[0];
cx  flagX[0], ancZ[0];
cx  q[7],  ancZ[0];



