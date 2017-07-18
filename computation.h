#ifndef COMPUTATION_H
#define COMPUTATION_H

// 2D element
struct Element2D
{
	int index[4];

	double p; 
	double q;
	
	// x-y-z coord index: 0 ==> xx, 1 ==> yy, 2 ==> xy, 3 ==> zz
	// cylinder coord index: 0==> rr, 1 ==> zz, 2 ==> rz, 3 ==> phai-phai
	double sig[4];	
	double sd[4];

	// x-y-z coord index: 0 ==> xx, 1 ==> yy, 2 ==> xy, 3 ==> zz
	// cylinder coord index: 0 ==> rr, 1 ==> zz, 2 ==> rz, 3 ==> phai-phai
	double eps_rate[4];
	double rot_rate; // first row and second col position value

	double rho0;
	double rho;	
	double e;	
	double mass;

	double nl[4][2];

	// x-y-z coord index: 0 ==> x, 1 ==> y
	// cylinder coord index: 0 ==> r, 1 ==> z
	double center0[2];
	double center[2];

	bool alive;

	// tensile strength
	double ten;
};

// 2D node
struct Node2D
{
	double mass;
    
    // movement constrain
    // x-y-z coord index: 0 ==> x, 1 ==> y
	// cylinder coord index: 0 ==> r, 1 ==> z
    // value: 0 ==> free, 1 ==> constrain
    int constrain[2];

	// x-y-z coord index: 0 ==> x, 1 ==> y
	// cylinder coord index: 0 ==> r, 1 ==> z
	double x0[2];
	double x[2];
	double u[2];

	double f_ext[2];
	double f_int[2];
	
	double f_h[2];
};

extern const double PI;
extern const double MIN_EPS;

extern Node2D *pND;
extern int ND_N; // size of pND array
extern Element2D *pELM;
extern int ELM_N; // size of pELM array

extern const int SOLV_TYPE;
extern const int GEO_TYPE;

extern const double L1_START;
extern const double L1_END;
extern const double TARGET_L2_START;
extern const double INTER_L2;
extern const double FLYER_L2_END;
// mesh generation control
extern const int TARGET_N2;
extern const int FLYER_N2;
extern const int N1;

extern const double IMPACT_U;

extern double RHO0;
extern double GAMMA;
extern double S;
extern double C0;

extern double G;
extern double Y0;
extern double TEN;

// time marching parameters
extern double TT;
extern double CT, DT;
extern double DT_RATIO;
extern int COUNT;

//void test();

void genMesh(
	int &ND_N, int &ELM_N,
	Node2D* &pND, Element2D* &pELM,
	double target_l2_start, double inter_l2, double flyer_l2_end,
	double l1_start, double l1_end,
	int target_n2, int flyer_n2, int n1
	);
void genMesh1(
	int &ND_N, int &ELM_N,
	Node2D* &pND, Element2D* &pELM,
	double flyer_l2_start, double flyer_l2_end,
	double l1_start, double l1_end,
	int flyer_n2, int n1
	);
void releaseMeshData();
void outputMeshData();
void initMeshData();

double computeDTFromElementSize();
void computeNode();
void computeElement();
void computeElementEpsRotRate(int elmNum);
void computeElementSig(int elmNum);
void ComputeElementSd(
    double sd_aver[2][2], double dev_eps_rate[2][2], 
    int ElmNum
    );
double ComputeElementQ(double v_aver, double theta_rate, int elmNum);
double ComputeElementP_GruneisenEOS(double EA, double V, double dV);
void computeElementNodeForce(int elmNum);

// matrix and vector tools
void computeMatrix2x_dot_Vector2D(
	double result[2], 
	double mx[2][2], double vec[2]
	);
void computeVector2D_dot_Matrix2x2(
	double result[2], 
	double vec[2], double mx[2][2]
	);
void computeMatrix2x2_dot_Matrix2x2(
	double result[2][2], 
	double mx_a[2][2], double mx_b[2][2]
	);
void scaleMatrix2x2(double result[2][2], double s, double mx[2][2]);
double computeMatrix2x2_doubledot_Matrix2x2(double mx_a[2][2], double mx_b[2][2]);
void computeMatrix2x2_plus_Matrix2x2(
	double result[2][2], 
	double mx_a[2][2], double mx_b[2][2]
	);
void computeMatrix2x2_sub_Matrix2x2(
	double result[2][2], 
	double mx_a[2][2], double mx_b[2][2]
	);
void copyMatrix2x2_to_Matrix2x2(double result[2][2], double mx_a[2][2]);

// geometry tools
double getArrayAver(int size, double array[]);
void getQuadrangleCenter(double ct[2], double x[4][2]);
double getQuadrangleArea(double x[4][2]);
double getQuadrangleDiagLength(double x[4][2]);

#endif