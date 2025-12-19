OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];          // data qubits
qreg ancX[1];       // ancilla

// =========================
// Stabilizer : ( IYXXY )
// ancilla = ancX[0]
// =========================
cx  ancX[0], q[2];
cy  ancX[0], q[1];
cy  ancX[0], q[4];
cx  ancX[0], q[3];

barrier;