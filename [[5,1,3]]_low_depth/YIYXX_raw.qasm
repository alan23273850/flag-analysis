OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];          // data qubits
qreg ancX[1];       // ancilla

// =========================
// Stabilizer : ( YIYXX )
// ancilla = ancX[0]
// =========================
cx  ancX[0], q[3];
cy  ancX[0], q[0];
cy  ancX[0], q[2];
cx  ancX[0], q[4];

barrier;