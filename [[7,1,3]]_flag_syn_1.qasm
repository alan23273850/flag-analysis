OPENQASM 2.0;
include "qelib1.inc";

qreg q[7];          // data qubits
qreg ancZ[1];       // ancillas measured in Z (for Z-type stabilizers)
qreg flagX[1];      // flags measured in X (paired with Z-type stabs)





// =========================
// Stabilizer 1 : ( I I I X X X X )
// ancilla = ancZ[0] (Z basis), flag = flagX[0] (X basis)
// =========================
cx  q[3], ancZ[0];
cx  flagX[0],ancZ[0];
cx  q[4], ancZ[0];
cx  q[5], ancZ[0];
cx  flagX[0],ancZ[0];
cx  q[6], ancZ[0];

