#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "computation.h"

using namespace std;

const double PI = 3.141592653589793;
const double MIN_EPS = 1.0e-12;

Node2D *pND;
int ND_N; // size of pND array
Element2D *pELM;
int ELM_N; // size of pELM array
 
// SOLV_TYPE:
// 0 ==> flyer-target impact
// 1 ==> Taylor bar impact
// 2 ==> blast wave
const int SOLV_TYPE = 0;
// GEO_TYPE:
// 0 ==> plain strain
// 1 ==> cylinder symmetry
const int GEO_TYPE = 0;

// global parameters for SOLV_TYPE ==> 0
// geometry
const double L1_START = -3.0; // cm
const double L1_END = 3.0; // cm
const double TARGET_L2_START = 0.0; // cm
const double INTER_L2 = 0.3; // cm
const double FLYER_L2_END = 0.4; // cm
// mesh generation control
const int TARGET_N2 = 60;
const int FLYER_N2 = 20;
const int N1 = 120;
// impact velocity
const double IMPACT_U = -0.05; // cm/us == 10^4 m/s

// material parameters
double RHO0 = 2.78; // g/cm3
double GAMMA = 2.13;
double S = 1.35;
double C0 = 0.54; // cm/us == 10^4 m/s

double G = 27.0e-2; // 10^2 GPa
double Y0 = 0.7e-2; // 10^2 GPa
double TEN = 2.0e-2; // 10^2 GPa

// time marching parameters
double TT = 4.0; // us
double CT, DT;
double DT_RATIO = 0.8;
int COUNT;

/*
void test()
{
    
    genMesh(
        ND_N, ELM_N, pND, pELM,
        TARGET_L2_START, INTER_L2, FLYER_L2_END,
        L1_START, L1_END,
        TARGET_N2, FLYER_N2, N1
        );
    
    genMesh1(ND_N, ELM_N, pND, pELM,
        TARGET_L2_START, INTER_L2,
        L1_START, L1_END,
        TARGET_N2, N1
        );
}
*/

void releaseMeshData()
{
	delete[] pND;
	delete[] pELM;
}

void outputMeshData()
{
    ofstream outFile("MeshData.dat");

    outFile << scientific;    
    outFile << "*NODE" << endl;
    for(int i=0; i<ND_N; ++i)
    {
        outFile << setw(8) << i+1;
        outFile << setw(20) << pND[i].x0[0];
        outFile << setw(20) << pND[i].x0[1];
        outFile << setw(4) << pND[i].constrain[0];
        outFile << setw(4) << pND[i].constrain[1];
        outFile << endl;        
    }
    outFile << "*ELEMENT" << endl;
    for(int i=0; i<ELM_N; ++i)
    {
        outFile << setw(8) << i+1;
        for(int j=0; j<4; ++j)
        {
            outFile << setw(8) << pELM[i].index[j]+1;
        }
        outFile << endl;
    }
    outFile << "*END" << endl;

    outFile.close();
}

void genMesh(
	int &ND_N, int &ELM_N,
	Node2D* &pND, Element2D* &pELM,
	double target_l2_start, double inter_l2, double flyer_l2_end,
	double l1_start, double l1_end,
	int target_n2, int flyer_n2, int n1
	)
{
    if(GEO_TYPE == 1 && (l1_start < 0.0 || l1_end < 0.0))
    {
        cout << "radius error! radius must be larger than 0. program stopped!" << endl;
        exit(1);            
    }

	ND_N = ( target_n2 + flyer_n2 + 1 )*( n1 + 1 );
	pND = new Node2D [ND_N];

	int i, j, k;
	double target_dx2, flyer_dx2, dx1;
	double target_l2_length, flyer_l2_length, l1_length;

    target_l2_length = fabs(inter_l2 - target_l2_start);
	flyer_l2_length = fabs(flyer_l2_end - inter_l2);	
	l1_length = fabs(l1_end - l1_start);

	target_dx2 = target_l2_length / target_n2;
    flyer_dx2 = flyer_l2_length / flyer_n2;	
	dx1 = l1_length / n1;
	
	k = 0;
	for(i=0; i<target_n2+flyer_n2+1; ++i)
	{
		for(j=0; j<n1+1; ++j)
		{
            // assume all nodes are free
            pND[k].constrain[0] = 0;
            pND[k].constrain[1] = 0;

			pND[k].x0[0] = l1_start + j*dx1;
            pND[k].u[0] = 0.0;
            if(i<=target_n2)
            {
				pND[k].x0[1] = target_l2_start + i*target_dx2;
                pND[k].u[1] = 0.0;
            }
			else
            {
				pND[k].x0[1] = inter_l2 + (i-target_n2)*flyer_dx2;
                pND[k].u[1] = IMPACT_U;
            }
            
            ++k;
		}		
	}

    // add constrain for cylinder symmetry
    if(GEO_TYPE == 1)
    {
        for(i=0; i<ND_N; ++i)
        {
            if(fabs(pND[i].x0[0]) < MIN_EPS)
                pND[i].constrain[0] = 1;
        }
    }

	ELM_N = (target_n2 + flyer_n2)*n1;
	pELM = new Element2D [ELM_N];

	int index1, index2;
	k = 0;
	for(i=0; i<target_n2+flyer_n2; ++i)
	{
		for(j=0; j<n1; ++j)
		{
			index1 = i*(n1+1) + j;
			index2 = (i+1)*(n1+1) + j;
			
			pELM[k].index[0] = index1;
            pELM[k].index[1] = index1 + 1;
            pELM[k].index[2] = index2 + 1;
			pELM[k].index[3] = index2;

			++k;
		}		
	}

    outputMeshData();
}

void genMesh1(
	int &ND_N, int &ELM_N,
	Node2D* &pND, Element2D* &pELM,
	double flyer_l2_start, double flyer_l2_end,
	double l1_start, double l1_end,
	int flyer_n2, int n1
	)
{
    if(GEO_TYPE == 1 && (l1_start < 0.0 || l1_end < 0.0))
    {
        cout << "radius error! radius must be larger than 0. program stopped!" << endl;
        exit(1);            
    }

	ND_N = ( flyer_n2 + 1 )*( n1 + 1 );
	pND = new Node2D [ND_N];

	int i, j, k;
	double flyer_dx2, dx1;
	double flyer_l2_length, l1_length;

	flyer_l2_length = fabs(flyer_l2_end - flyer_l2_start);
	l1_length = fabs(l1_end - l1_start);

	flyer_dx2 = flyer_l2_length / flyer_n2;
	dx1 = l1_length / n1;
	
	k = 0;
	for(i=0; i<flyer_n2+1; ++i)
	{
		for(j=0; j<n1+1; ++j)
		{
            // assume all nodes are free
            pND[k].constrain[0] = 0;
            pND[k].constrain[1] = 0;

			pND[k].x0[0] = l1_start + j*dx1;
            pND[k].x0[1] = flyer_l2_start + i*flyer_dx2; 
            
            pND[k].u[0] = 0.0;
            if(i == 0)
                pND[k].u[1] = 0.0;
            else
                pND[k].u[1] = IMPACT_U;

            ++k;
		}		
	}
    
    // add rigid wall constrain
    for(i=0; i<ND_N; ++i)
    {
        if(fabs(pND[i].x0[1] - flyer_l2_start) < MIN_EPS)
            pND[i].constrain[1] = 1;
    }

    // add constrain for cylinder symmetry
    if(GEO_TYPE == 1)
    {
        for(i=0; i<ND_N; ++i)
        {
            if(fabs(pND[i].x0[0]) < MIN_EPS)
                pND[i].constrain[0] = 1;         
        }
    }

	ELM_N = flyer_n2*n1;
	pELM = new Element2D [ELM_N];

	int index1, index2;
	k = 0;
	for(i=0; i<flyer_n2; ++i)
	{
		for(j=0; j<n1; ++j)
		{
			index1 = i*(n1+1) + j;
			index2 = (i+1)*(n1+1) + j;
			
			pELM[k].index[0] = index1;
            pELM[k].index[1] = index1 + 1;
            pELM[k].index[2] = index2 + 1;
			pELM[k].index[3] = index2;

			++k;
		}		
	}

    outputMeshData();
}

void initMeshData()
{
	int i, j, k;
	for(i=0; i<ND_N; ++i)
	{
		pND[i].x[0] = pND[i].x0[0];
		pND[i].x[1] = pND[i].x0[1];

		pND[i].f_ext[0] = 0.0;
		pND[i].f_ext[1] = 0.0;

		pND[i].f_int[0] = 0.0;
		pND[i].f_int[1] = 0.0;

		pND[i].f_h[0] = 0.0;
		pND[i].f_h[1] = 0.0;

        pND[i].mass = 0.0;
	}

	int index;
	double x[4][2], ct[2];
	double r, area;
	for(i=0; i<ELM_N; ++i)
	{
		pELM[i].rho0 = RHO0;
		pELM[i].rho = RHO0;

		for(j=0; j<4; ++j)
		{
			index = pELM[i].index[j];
			x[j][0] = pND[index].x[0];
			x[j][1] = pND[index].x[1];
		}
		getQuadrangleCenter(ct, x);
		pELM[i].center0[0] = ct[0];
		pELM[i].center0[1] = ct[1];
		pELM[i].center[0] = ct[0];
		pELM[i].center[1] = ct[1];		

        if(GEO_TYPE == 0) // plain strain
			pELM[i].mass = RHO0*getQuadrangleArea(x);
		else // GEO_TYPE == 1, cylinder symmetry
		{
			r = pELM[i].center0[0];
			area = getQuadrangleArea(x);
			pELM[i].mass = RHO0*2.0*PI*r*area;
		}			

		pELM[i].e = 0.0;		
		pELM[i].p = 0.0;		
		pELM[i].q = 0.0;
		
		pELM[i].rot_rate = 0.0;

		for(j=0; j<4; j++)
		{
			pELM[i].sig[j] = 0.0;			
			pELM[i].sd[j] = 0.0;			
			pELM[i].eps_rate[j] = 0.0;
		}

		for(j=0; j<4; j++)
		{
			index = pELM[i].index[j];
			pND[index].mass += pELM[i].mass / 4.0;
		}

		pELM[i].alive = true;
		pELM[i].ten = TEN;
	}
}

double computeDTFromElementSize()
{
    double x[4][2];
    int index;
    double length, length_1;

    // for first element: i = 0
    for(int j=0; j<4; ++j)
	{
		index = pELM[0].index[j];
		x[j][0] = pND[index].x[0];
		x[j][1] = pND[index].x[1];
	}
    length = getQuadrangleDiagLength(x);

    // for rest of the elements: i <= [1, ELM_N-1]
    for(int i=1; i<ELM_N; ++i)
    {
        for(int j=0; j<4; ++j)
		{
			index = pELM[i].index[j];
			x[j][0] = pND[index].x[0];
			x[j][1] = pND[index].x[1];
		}
        length_1 = getQuadrangleDiagLength(x);
        if(length_1 < length)
            length = length_1;
    }

    return DT_RATIO*length / C0;
}

void computeNode()
{
	int i, j;
	double mass;
	double f_ext[2], f_int[2], f_h[2];
	double u_n[2], u_n1[2];
	double x_n[2], x_n1[2];

	for(i=0; i<ND_N; ++i)
	{
		mass = pND[i].mass;
		for(j=0; j<2; ++j)
		{
            if(pND[i].constrain[j] == 0)
            {
                u_n[j] = pND[i].u[j];			
			    f_ext[j] = pND[i].f_ext[j];
			    f_int[j] = pND[i].f_int[j];
			    f_h[j] = pND[i].f_h[j];
			
			    u_n1[j] = u_n[j] + DT*(	f_ext[j] - f_int[j] + f_h[j] ) / mass;
			    pND[i].u[j] = u_n1[j];

			    x_n[j] = pND[i].x[j];			
			
			    x_n1[j] = x_n[j] + DT*u_n1[j];
			    pND[i].x[j] = x_n1[j];
            }
            else
            {
                pND[i].u[j] = 0.0;
            }			
		}
	}

     // clear all node force
    for(i=0; i<ND_N; ++i)
	{
		pND[i].f_ext[0] = 0.0;
		pND[i].f_ext[1] = 0.0;

		pND[i].f_int[0] = 0.0;
		pND[i].f_int[1] = 0.0;		

		pND[i].f_h[0] = 0.0;
		pND[i].f_h[1] = 0.0;
	}
}

void computeElement()
{
	int i;	
	for(i=0; i<ELM_N; ++i)
	{
		if(pELM[i].alive)
		{
			computeElementEpsRotRate(i);
			computeElementSig(i);
		}		
	}	

	for(i=0; i<ELM_N; ++i)
	{
		if(pELM[i].alive)
		{
			computeElementNodeForce(i);
		}		
	}
}

void computeElementEpsRotRate(int elmNum)
{
	int i, j, k;
	double u[4][2];
	double x[4][2];
	double ct[2];
	int index;

	for(i=0; i<4; i++)
	{
		index = pELM[elmNum].index[i];
		for(j=0; j<2; j++)
		{
			u[i][j] = pND[index].u[j];
			x[i][j] = pND[index].x[j];
		}
	}

	getQuadrangleCenter(ct, x);	
	pELM[i].center[0] = ct[0];
	pELM[i].center[1] = ct[1];

	double temp1[2][2];
	double detJ;
	temp1[0][0] = 0.25*( - x[0][0] + x[1][0] + x[2][0] - x[3][0] );
	temp1[0][1] = 0.25*( - x[0][0] - x[1][0] + x[2][0] + x[3][0] );
	temp1[1][0] = 0.25*( - x[0][1] + x[1][1] + x[2][1] - x[3][1] );
	temp1[1][1] = 0.25*( - x[0][1] - x[1][1] + x[2][1] + x[3][1] );
	detJ = temp1[0][0]*temp1[1][1] - temp1[0][1]*temp1[1][0];
    if(fabs(detJ) <= MIN_EPS)
    {
        cout << "detJ calculation error! program stoped!" << endl;
        exit(1);
    }
		
	double temp2[2][2];
	temp2[0][0] = temp1[1][1] / detJ;
	temp2[0][1] = -temp1[0][1] / detJ;
	temp2[1][0] = -temp1[1][0] / detJ;
	temp2[1][1] = temp1[0][0] / detJ;

	double nv[4][2] = 
	{
        {-0.25, -0.25},
        {0.25, -0.25},
        {0.25, 0.25},
        {-0.25, 0.25}		
	};
	double nl[4][2];
	double s[2];
	for(i=0; i<4; ++i)
	{
		computeVector2D_dot_Matrix2x2(s, nv[i], temp2);

		for(j=0; j<2; ++j)
		{
			nl[i][j] = s[j];
			pELM[elmNum].nl[i][j] = s[j];
		}
	}

	double grad_u[2][2];
	double sum;
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			sum = 0.0;
			for(k=0; k<4; k++)
			{
				sum += u[k][i]*nl[k][j]; 
			}
			grad_u[i][j] = sum;
		}
	}

	pELM[elmNum].eps_rate[0] = grad_u[0][0];
	pELM[elmNum].eps_rate[1] = grad_u[1][1];
    pELM[elmNum].eps_rate[2] = 0.5*(grad_u[0][1] + grad_u[1][0]);
    if(GEO_TYPE == 0) // plain strain
        pELM[elmNum].eps_rate[3] = 0.0;
    else // cylinder symmetry
        pELM[elmNum].eps_rate[3] = 0.25
        *(u[0][0] + u[1][0] + u[2][0] + u[3][0]) / ct[0];

	pELM[elmNum].rot_rate = 0.5*(grad_u[0][1] - grad_u[1][0]);	
}

void computeElementSig(int elmNum)
{
    // compute deviate eps rate	
	double theta_rate;
    if(GEO_TYPE == 0) // plain strain
        theta_rate = pELM[elmNum].eps_rate[0] + pELM[elmNum].eps_rate[1];
    else // cylinder symmetry
        theta_rate = pELM[elmNum].eps_rate[0] + pELM[elmNum].eps_rate[1]
        + pELM[elmNum].eps_rate[3];

	double dev_eps_rate[2][2];
	dev_eps_rate[0][0] = pELM[elmNum].eps_rate[0] - theta_rate / 3.0;
	dev_eps_rate[1][1] = pELM[elmNum].eps_rate[1] - theta_rate / 3.0;
	dev_eps_rate[0][1] = pELM[elmNum].eps_rate[2];
	dev_eps_rate[1][0] = pELM[elmNum].eps_rate[2];
	
	double sd_aver[2][2];
	ComputeElementSd(sd_aver, dev_eps_rate, elmNum);

	double v_n, v_n1;
	v_n = 1.0 / pELM[elmNum].rho;		
	v_n1 = v_n + DT*theta_rate*v_n;
	if(v_n1 < 0.0)
	{
		cout << "v_n1 error! Program stopped!" << endl;
		exit(1);
	}
	pELM[elmNum].rho = 1.0 / v_n1;

	double v_aver, dv;	
	v_aver = 0.5*(v_n1 + v_n);
	dv = v_n1 - v_n;
	
	double q;
	q = ComputeElementQ(v_aver, theta_rate, elmNum);

	double e_n, p_n;
	e_n = pELM[elmNum].e;
	p_n = pELM[elmNum].p;
	
	double term1, term2;
	term1 = (0.5*p_n + q)*dv;
	
	term2 = computeMatrix2x2_doubledot_Matrix2x2(sd_aver, dev_eps_rate);
	term2 += ( sd_aver[0][0] + sd_aver[1][1] ) * ( dev_eps_rate[0][0] + dev_eps_rate[1][1] );
	term2 = v_aver * term2 * DT;

	double ea = e_n - term1 + term2;
	double p = ComputeElementP_GruneisenEOS(ea, v_n1, dv);
	double e = ea - 0.5*p*dv;

	pELM[elmNum].p = p;
	pELM[elmNum].e = e;

	pELM[elmNum].sig[0] = -(p + q) + pELM[elmNum].sd[0];
	pELM[elmNum].sig[1] = -(p + q) + pELM[elmNum].sd[1];
	pELM[elmNum].sig[2] = pELM[elmNum].sd[2];
    pELM[elmNum].sig[3] = -(p + q) + pELM[elmNum].sd[3];

	double ten = pELM[elmNum].ten;	
	if(p < 0.0 && pELM[elmNum].sig[1] > ten)
	{
		pELM[elmNum].alive = false;
	}
}

void ComputeElementSd(
    double sd_aver[2][2], double dev_eps_rate[2][2], 
    int elmNum
    )
{
	int i, j;
	double sdj_rate[2][2];
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			sdj_rate[i][j] = 2.0*G*dev_eps_rate[i][j];
		}		
	}

	double sd_n[2][2];
	sd_n[0][0] = pELM[elmNum].sd[0];
	sd_n[1][1] = pELM[elmNum].sd[1];
	sd_n[0][1] = pELM[elmNum].sd[2];
	sd_n[1][0] = pELM[elmNum].sd[2];	
	
	double rot_rate[2][2];
	rot_rate[0][0] = 0.0;
	rot_rate[1][1] = 0.0;
	rot_rate[0][1] = pELM[elmNum].rot_rate;
	rot_rate[1][0] = -pELM[elmNum].rot_rate;	
	
	double s1[2][2], s2[2][2], s3[2][2], s4[2][2];
	computeMatrix2x2_dot_Matrix2x2(s1, sd_n, rot_rate);
	computeMatrix2x2_dot_Matrix2x2(s2, rot_rate, sd_n);
	computeMatrix2x2_sub_Matrix2x2(s3, sdj_rate, s1);
	computeMatrix2x2_plus_Matrix2x2(s4, s3, s2);

	double sd_rate[2][2];
	copyMatrix2x2_to_Matrix2x2(sd_rate, s4);

	double sd_n1[2][2];
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			sd_n1[i][j] = sd_n[i][j] + DT*sd_rate[i][j];
		}
	}

	double dd = computeMatrix2x2_doubledot_Matrix2x2(sd_n1, sd_n1);
	dd += (sd_n1[0][0] + sd_n1[1][1]) * (sd_n1[0][0] + sd_n1[1][1]);

	double sig_eff = sqrt( 1.5*dd );
	if(sig_eff < 0.0)
	{
		cout << "sig_eff error! Program stopped!" << endl;
		exit(1);
	}

	if(sig_eff > Y0)
	{
		double s[2][2];
		double f;

		f = Y0 / sig_eff;
		scaleMatrix2x2(s, f, sd_n1);
		copyMatrix2x2_to_Matrix2x2(sd_n1, s);
	}

	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			sd_aver[i][j] = 0.5*( sd_n[i][j] + sd_n1[i][j] );
		}
	}

	pELM[elmNum].sd[0] = sd_n1[0][0];
	pELM[elmNum].sd[1] = sd_n1[1][1];
	pELM[elmNum].sd[2] = 0.5*( sd_n1[0][1] + sd_n1[1][0] );
    pELM[elmNum].sd[3] = -(sd_n1[0][0] + sd_n1[1][1]);
}

double ComputeElementQ(double v_aver, double theta_rate, int elmNum)
{
	if(theta_rate < 0.0)
	{
		
		double q, q1, q2, c0, rho;
		q1 = 1.5;
		//q2 = 0.06;
		//q2 = 0.5;
		q2 = 0.3;
		c0 = C0;
		rho = 1.0 / v_aver;

		double area_aver;
        if(GEO_TYPE == 0)
            area_aver = v_aver * pELM[elmNum].mass;
        else
            area_aver = v_aver * pELM[elmNum].mass
            / (2.0*PI*pELM[elmNum].center[0]);

		q = q1*q1 * rho * area_aver * theta_rate*theta_rate
            + q2 * rho * sqrt(area_aver) * c0 * fabs(theta_rate);

		pELM[elmNum].q = q;
		return q;
	}
	else
	{
		pELM[elmNum].q = 0.0;
		return 0.0;
	}
}

double ComputeElementP_GruneisenEOS(double ea, double V, double dV)
{
	double gama = GAMMA;
	double rho0 = RHO0;
	double c0 = C0;
	double s = S;

	double V0=1.0/rho0;
	double temp=1.0-V/V0;

	double pH = rho0*c0*c0*temp / ( (1.0-s*temp)*(1.0-s*temp) );
	double eH = 0.5*pH*(V0-V);

	return ( pH + rho0*gama*(ea-eH) ) / ( 1.0 + 0.5*rho0*gama*dV );
}

void computeElementNodeForce(int elmNum)
{
    int i, j;
	double u[4][2];	
	int index;

    for(i=0; i<4; i++)
	{
		index = pELM[elmNum].index[i];
		for(j=0; j<2; j++)
		{
			u[i][j] = pND[index].u[j];
		}
	}

	double rho = pELM[elmNum].rho;
    double qh = 0.1;
	double c0 = C0;
    double gb[4] = {1.0, -1.0, 1.0, -1.0};
	
    double alpha_h, area, volume;
    if(GEO_TYPE == 0) // plain strain
    {
        area = pELM[elmNum].mass / rho;
        alpha_h = qh * rho * sqrt(area) * c0 / 4.0;
    }        
    else // cylinder symmetry
    {
        volume = pELM[elmNum].mass / rho;
        alpha_h = qh * rho * pow(volume, 2.0/3.0) * c0 / 4.0;
    } 

	// get hourglass force
    double fh[4][2];
	double sum;
	for(i=0; i<2; i++)
	{
		sum = 0.0;
		for(j=0; j<4; j++)
		{
			sum += gb[j]*u[j][i];
		}
		
		for(j=0; j<4; j++)
		{
			fh[j][i] = -alpha_h * sum * gb[j];
		}
	}

	double nl[4][2];
	for(i=0; i<4; i++)
	{
		for(j=0; j<2; j++)
		{
			nl[i][j] = pELM[elmNum].nl[i][j];
		}
	}

	double sig[2][2];
    double sig_phai, radius;
	sig[0][0] = pELM[elmNum].sig[0];
	sig[1][1] = pELM[elmNum].sig[1];
	sig[0][1] = pELM[elmNum].sig[2];
	sig[1][0] = pELM[elmNum].sig[2];

	double f[4][2];
	double s[2], v[2];
	for(i=0; i<4; i++)
	{
		v[0] = nl[i][0];
		v[1] = nl[i][1];

		computeVector2D_dot_Matrix2x2(s, v, sig);
        if(GEO_TYPE == 0) // plain strain
        {
            f[i][0] = area * s[0];
		    f[i][1] = area * s[1];
        }
        else // cylinder symmetry
        {
            sig_phai = pELM[elmNum].sig[3];
            radius = pELM[elmNum].center[0];
            f[i][0] = volume * (s[0] + sig_phai*0.25 / radius);
		    f[i][1] = volume * s[1];
        }				
	}

	for(i=0; i<4; i++)
	{
		index = pELM[elmNum].index[i];
		for(j=0; j<2; j++)
		{
			pND[index].f_int[j] += f[i][j];
			pND[index].f_h[j] += fh[i][j];
		}
	}
}

// matrix and vector tools
void computeMatrix2x2_dot_Vector2D(
	double result[2], 
	double mx[2][2], double vec[2]
	)
{
	int i, j;
	double sum;
	for(i=0; i<2; i++)
	{
		sum = 0.0;
		for(j=0; j<2; j++)
		{
			sum += mx[i][j]*vec[j];
		}
		result[i] = sum;
	}
}

void computeVector2D_dot_Matrix2x2(
	double result[2], 
	double vec[2], double mx[2][2]
	)
{
	int i, j;
	double sum;
	for(j=0; j<2; j++)
	{
		sum = 0.0;
		for(i=0; i<2; i++)
		{
			sum += vec[i]*mx[i][j];
		}
		result[j] = sum;
	}
}

void computeMatrix2x2_dot_Matrix2x2(
	double result[2][2], 
	double mx_a[2][2], double mx_b[2][2]
	)
{
	int i, j, k;
	double sum;
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			sum = 0.0;
			for(k=0; k<2; k++)
			{
				sum += mx_a[i][k]*mx_b[k][j];
			}
			result[i][j] = sum;
		}
	}
}

void scaleMatrix2x2(double result[2][2], double s, double mx[2][2])
{
	int i, j;
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			result[i][j] = s*mx[i][j];
		}
	}

}

double computeMatrix2x2_doubledot_Matrix2x2(double mx_a[2][2], double mx_b[2][2])
{
	double mx[2][2];
	mx[0][0] = mx_a[0][0];
	mx[1][1] = mx_a[1][1];
	mx[0][1] = mx_a[1][0];
	mx[1][0] = mx_a[0][1];	
	double s[2][2];
	computeMatrix2x2_dot_Matrix2x2(s, mx, mx_b);
	double sum = s[0][0] + s[1][1];	
	return sum;	
}

void computeMatrix2x2_plus_Matrix2x2(
	double result[2][2], 
	double mx_a[2][2], double mx_b[2][2]
	)
{
	int i, j;	
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			result[i][j] = mx_a[i][j] + mx_b[i][j];
		}
	}

}

void computeMatrix2x2_sub_Matrix2x2(
	double result[2][2], 
	double mx_a[2][2], double mx_b[2][2]
	)
{
	int i, j;	
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			result[i][j] = mx_a[i][j] - mx_b[i][j];
		}
	}
}

void copyMatrix2x2_to_Matrix2x2(double result[2][2], double mx_a[2][2])
{
	int i, j;	
	for(i=0; i<2; i++)
	{
		for(j=0; j<2; j++)
		{
			result[i][j] = mx_a[i][j];
		}
	}	
}

// geometry tools
double getArrayAver(int size, double array[])
{
	double sum = 0.0;
	for(int i=0; i<size; ++i)
	{
		sum += array[i];
	}
	return sum / size;
}

void getQuadrangleCenter(double ct[2], double x[4][2])
{
	double sum;
	for(int j=0; j<2; ++j)
	{
		sum = 0.0;
		for(int i=0; i<4; ++i)
			sum += x[i][j];
		ct[j] = sum / 4.0;
	}
}

double getQuadrangleArea(double x[4][2])
{
	double a0, a1, b0, b1, area;
	double sum = 0.0;
	// area of two triangles
	for(int i=0; i<2; ++i)
	{
		a0 = x[i+1][0] - x[0][0];
	    a1 = x[i+1][1] - x[0][1];
		b0 = x[i+2][0] - x[0][0];
		b1 = x[i+2][1] - x[0][1];
		area = 0.5*(a0*b1 - a1*b0);
		if(area <= MIN_EPS)
		{
			cout << "Triangle area calculation error! program stoped!" << endl;
			exit(1);
		}	
		sum += area;
	}
	return sum;
}

double getQuadrangleDiagLength(double x[4][2])
{
    double t0, t1;
    double len1, len2;
    t0 = (x[2][0] - x[0][0])*(x[2][0] - x[0][0]);
    t1 = (x[2][1] - x[0][1])*(x[2][1] - x[0][1]);
    len1 = sqrt(t0 + t1);

    t0 = (x[3][0] - x[1][0])*(x[3][0] - x[1][0]);
    t1 = (x[3][1] - x[1][1])*(x[3][1] - x[1][1]);
    len2 = sqrt(t0 + t1);

    return len1 < len2 ? len1 : len2;
}

