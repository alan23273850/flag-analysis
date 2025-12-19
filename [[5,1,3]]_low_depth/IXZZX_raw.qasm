OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];          // data qubits
qreg ancX[1];       // ancilla

// =========================
// Stabilizer : ( IXZZX )
// ancilla = ancX[0]
// =========================
cx  ancX[0], q[1];
cz  ancX[0], q[2];
cz  ancX[0], q[3];
cx  ancX[0], q[4];

barrier;