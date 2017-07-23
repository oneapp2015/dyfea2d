#include <iostream>
#include "computation.h"

using namespace std;

int main()
{
    // generate mesh
    genMesh1(ND_N, ELM_N, pND, pELM,
        TARGET_L2_START, INTER_L2,
        L1_START, L1_END,
        TARGET_N2, N1
        );
    // initialize mesh data
    initMeshData();

    outputNodeIndexforElem();

    CT = 0.0; // current time
    DT = computeDTFromElementSize();

    while(CT < TT)
    {
        cout << CT << endl;
        CT += DT;
    }
    
    releaseMeshData();
	return 0;
}