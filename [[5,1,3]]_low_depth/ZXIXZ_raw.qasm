OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];          // data qubits
qreg ancX[1];       // ancilla

// =========================
// Stabilizer : ( ZXIXZ )
// ancilla = ancX[0]
// =========================
cz  ancX[0], q[0];
cx  ancX[0], q[1];
cx  ancX[0], q[3];
cz  ancX[0], q[4];

barrier;