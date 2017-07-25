#include <iostream>
#include <iomanip>
#include "computation.h"

using namespace std;

int main()
{
    // generate mesh
    genMesh(ND_N, ELM_N, pND, pELM,
        TARGET_L2_START, INTER_L2, FLYER_L2_END,
        L1_START, L1_END,
        TARGET_N2, FLYER_N2, N1);
    
    /*
    genMesh1(ND_N, ELM_N, pND, pELM,
        INTER_L2, FLYER_L2_END,
        L1_START, L1_END,
        FLYER_N2, N1
        );
    */

    // initialize mesh data
    initMeshData();

    // output node index information
    outputNodeIndexforElem();

    CT = 0.0; // current time
    DT = computeDTFromElementSize();
    int steps = TT/DT;
    int out_step = steps / OUT_INTERVAL;
    int out_index = 0;
    int disp_step = steps / DISP_INTERVAL;
    COUNT = 0;

    // output initial state
    outputState(out_index);

    cout << "let's go!!! " << endl << endl;
    while(CT < TT)
    {
        if(COUNT % disp_step == 0)
        {
            cout << setw(8) << 100.0*CT/TT << " % is completed";
            cout << " ====> CT / TT : " << CT << " / " << TT << endl;
        }
               
        computeNode();
        computeElement();

        CT += DT;
        ++COUNT;

        if(COUNT % out_step == 0)
        {
            ++out_index;
            outputState(out_index);
        }
    }

    outputEnergy();
    
    releaseMeshData();
	return 0;
}