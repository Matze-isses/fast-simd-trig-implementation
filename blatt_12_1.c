#include "quest.h"
#include "math.h"
#include <stdio.h>


void applyQFT(Qureg qureg) {
  int controls[1];
  int targets[1];

  // Transforms on Q3
  applyHadamard(qureg, 3);


  // Transforms on Q2
  controls[0] = 3;
  targets[0] = 2;
  applyMultiControlledPhaseGadget(qureg, controls, 1, targets, 1, M_PI / 2);
  applyHadamard(qureg, 2);

  // Transforms on Q1
  controls[0] = 3;
  targets[0] = 1;
  applyMultiControlledPhaseGadget(qureg, controls, 1, targets, 1, M_PI / 4);

  controls[0] = 2;
  applyMultiControlledPhaseGadget(qureg, controls, 1, targets, 1, M_PI / 2);
  applyHadamard(qureg, 1);

  // Transforms on Q0
  controls[0] = 3;
  targets[0] = 0;
  applyMultiControlledPhaseGadget(qureg, controls, 1, targets, 1, M_PI / 8);

  controls[0] = 2;
  targets[0] = 0;
  applyMultiControlledPhaseGadget(qureg, controls, 1, targets, 1, M_PI / 4);

  controls[0] = 1;
  applyMultiControlledPhaseGadget(qureg, controls, 1, targets, 1, M_PI / 2);
  applyHadamard(qureg, 0);
}

int main(void) {
    initQuESTEnv();
    Qureg qureg = createForcedQureg(4);

    initClassicalState(qureg, 4);

    applyQFT(qureg);
    reportQureg(qureg);


    printf("Basis Probabilities TASK A:\n");
    printf("|1000>: %f\n", getQuregAmp(qureg, 0));
    printf("|0100>: %f\n", getQuregAmp(qureg, 1));
    printf("|0010>: %f\n", getQuregAmp(qureg, 2));
    printf("|0001>: %f\n", getQuregAmp(qureg, 4));


    // For task b
    initZeroState(qureg);

    // Hadamard auf Qubit 0 → erzeugt (|0> + |1>) / sqrt(2)
    // Resultat: (|0000> + |1000>) / sqrt(2)
    applyHadamard(qureg, 0);

    applyQFT(qureg);

    reportQureg(qureg);
    printf("Basis Probabilities TASK B:\n");
    printf("|1000>: %f\n", getQuregAmp(qureg, 0));
    printf("|0100>: %f\n", getQuregAmp(qureg, 1));
    printf("|0010>: %f\n", getQuregAmp(qureg, 2));
    printf("|0001>: %f\n", getQuregAmp(qureg, 3));

    // qreal prob = calcTotalProb(qureg);
    // reportScalar("Total probability", prob);
    destroyQureg(qureg);
    finalizeQuESTEnv();

    return 0;
}
