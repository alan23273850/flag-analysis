OPENQASM 2.0;
include "qelib1.inc";

qreg q[5];          // data qubits
qreg ancX[0];       // ancillas measured in X (for X-type stabilizers)
qreg ancZ[4];       // ancillas measured in Z (for Z-type stabilizers)
qreg flagX[4];      // flags measured in X (paired with Z-type stabs)
qreg flagZ[0];      // flags measured in Z (paired with X-type stabs)

creg syn[4];        // syndrome bits (one per stabilizer)
creg flgX[4];       // flag bits for X-basis flags  (used in stabs 4..6)
creg flgZ[0];       // flag bits for Z-basis flags  (used in stabs 1..3)


// Define a new 2-qubit gate "notnot"
gate notnot c, t {
    h c;
    cx c, t;
    h c;
}

barrier q, ancX, ancZ, flagX, flagZ;

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
measure ancZ[0]  -> syn[0];
measure flagX[0] -> flgX[0];

barrier q, ancX, ancZ, flagX, flagZ;

// =========================
// Stabilizer 2 : 
// ancilla = ancZ[1] (Z basis), flag = flagX[1] (X basis)
// =========================
notnot q[1], ancZ[1];
cx flagX[1], ancZ[1];
cx q[2], ancZ[1];
cx q[3], ancZ[1];
cx flagX[1], ancZ[1];
notnot q[4], ancZ[1];
measure ancZ[1]  -> syn[1];
measure flagX[1] -> flgX[1];

barrier q, ancX, ancZ, flagX, flagZ;

// =========================
// Stabilizer 3 : 
// ancilla = ancZ[2] (Z basis), flag = flagX[2] (X basis)
// =========================
notnot q[0], ancZ[2];
cx flagX[2], ancZ[2];
cx q[3], ancZ[2];
cx q[4], ancZ[2];
cx flagX[2], ancZ[2];
notnot q[2], ancZ[2];
measure ancZ[2]  -> syn[2];
measure flagX[2] -> flgX[2];

barrier q, ancX, ancZ, flagX, flagZ;

// =========================
// Stabilizer 4 : 
// ancilla = ancZ[3] (Z basis), flag = flagX[3] (X basis)
// =========================
notnot q[1], ancZ[3];
cx flagX[3], ancZ[3];
cx q[0], ancZ[3];
cx q[4], ancZ[3];
cx flagX[3], ancZ[3];
notnot q[3], ancZ[3];
measure ancZ[3]  -> syn[3];
measure flagX[3] -> flgX[3];

barrier q, ancX, ancZ, flagX, flagZ;
