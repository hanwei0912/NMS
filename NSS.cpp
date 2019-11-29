#include<cstdio>
#include<iostream>
#include<cstring>
#include<set>
#include<cmath>
#include<vector>
#include<algorithm>
#include<time.h>
#include "mex.h"

using namespace std;

#define ct(A) cout << #A << ": " << (A) << endl;
#define PI  3.1415926535897932384626433832795
const double eps = 1e-6;


int XDim;    // dimension of X
int FDim;    // dimension of F
int NData;   // number of data
double fes;     // number of the times of evaluation
double *W;   // pointer to weight
double *ideaPoint; // pointer to ideaPoint
double alpha, beta, gamma; // parameters
char *problem;  // name of problem
double *domain; // range of x
double y_obj[4];
double* evaluate(double *ps, int Xdim) {
	// 	double y_obj[4];
	double *x = ps;
	int nx = XDim;
	if (!strcmp("ZDT1", problem)){
		double g = 0;
		for (int n = 1; n<nx; n++)
			g += x[n] * x[n];
		// 		g = 1 + 9*g/(nx-1);

		y_obj[0] = x[0];
		y_obj[1] = (1 - sqrt(g));
		return y_obj;
	}
	if (!strcmp("TEC09_F1", problem)){
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for (j = 1; j < nx; j++)
		{
			yj = x[j] - pow(x[0], (0.5 + 1.5*(j - 1) / (nx - 2.0)));
			yj = yj * yj;
			if (j % 2 == 0)
			{
				sum1 += yj;
				count1++;
			}
			else
			{
				sum2 += yj;
				count2++;
			}
		}
		y_obj[0] = x[0] + 2.0 * sum1 / (double)count1;
		y_obj[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
		return y_obj;
	}
	if (!strcmp("TEC09_F3", problem)){
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1 = sum2 = 0.0;
		count1 = count2 = 0;
		for (j = 1; j < nx; j++)
		{
			if (j % 2 == 0)
			{
				yj = x[j] - 0.8*x[0] * cos(6 * PI*x[0] + (j + 1)*PI / nx);
				yj = yj * yj;
				sum1 += yj;
				count1++;
			}
			else
			{
				yj = x[j] - 0.8*x[0] * sin(6 * PI*x[0] + (j + 1)*PI / nx);
				yj = yj * yj;
				sum2 += yj;
				count2++;
			}
		}
		y_obj[0] = x[0] + 2.0 * sum1 / (double)count1;
		y_obj[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
		return y_obj;
	}
}

double evaluate_sub(double *pf, int Fdim){
	double tmp, max = -9999, sum = 0;
	for (int i = 0; i < Fdim; i++){
// 		sum += (W[i] * obj[i]);
		tmp = W[i] *  (pf[i] - ideaPoint[i]);
		if (tmp>max)
            max = tmp;
	}
	return max;
// 	return sum;
}

struct node {
	int Xdim;
    int Fdim;
	double ps[50];
    double pf[4];
	double value;
	//     double obj[4];
	node() {
		Xdim = 0;
        Fdim=0;
		value = 0.0;
	}
	void update() {
		for (int i = 0; i<Xdim; i++){
			srand((unsigned)time(NULL));
			double r = (rand() % 100) / 100.0;
			if (ps[i]>domain[Xdim + i])
				ps[i] = 2 * domain[Xdim + i] - ps[i] + r*(ps[i] - domain[Xdim + i]);
			if (ps[i]<domain[i])
				ps[i] = 2 * domain[i] - ps[i] - r*(domain[i] - ps[i]);
		}
        double *obj = evaluate(ps, Xdim);
        for (int j=0;j<Fdim;j++){
            pf[j]=obj[j];
        }
		value = evaluate_sub(pf, Fdim);
	}
    void polymutate(){
        double eta_m = 20;
        double rate=0.01;
        double mut_pow = 1.0/(eta_m+1.0); 
        for(int j=0;j<XDim;j++){
            double r=(rand() % 10000) / 10000.0;
            if(r<=rate){
                double y=ps[j];
                double yl=domain[j];
                double yu=domain[XDim+j];
                double delta1=(y-yl)/(yu-yl);
                double delta2=(yu-y)/(yu-yl);

                double rnd=(rand() % 10000) / 10000.0;
                double deltaq;
                if(rnd<=0.5){
                    double xy=1.0-delta1;
                    double val=2.0*rnd+(1.0-2.0*rnd)*pow(xy,eta_m+1.0);
                    deltaq=pow(val,mut_pow)-1.0;
                }
                else{
                    double xy=1.0-delta2;
                    double val=2.0*(1.0-rnd)+2.0*(rnd-0.5)*pow(xy,eta_m+1.0);
                    deltaq=1.0-pow(val,mut_pow);
                }
                y=y+deltaq*(yu-yl);
                if(y>=yl&&y<=yu){
                    ps[j]=y;
                }
            }
        }
    }
	node(double *_ps, int _Xdim, double *_pf, int _Fdim) {
		Xdim = _Xdim;
        Fdim = _Fdim;
		for (int i = 0; i<_Xdim; i++) {
			ps[i] = _ps[i];
		}
        for(int i=0;i<_Fdim;i++){
            pf[i]=_pf[i];
        }
		value=evaluate_sub(pf, Fdim);
	}
	bool operator < (const node &b) const {
		if (fabs(value - b.value) > eps) return value < b.value;
		if (Xdim < b.Xdim) return Xdim < b.Xdim;
		for (int i = 0; i<FDim; i++) {
			if (fabs(pf[i] - b.pf[i]) > eps) return pf[i] < b.pf[i];
		}
	}
	bool operator == (const node &b) const {
		if (Xdim != b.Xdim) return false;
		for (int i = 0; i<Xdim; i++) {
			if (fabs(ps[i] - b.ps[i]) > eps) return false;
		}
		return fabs(value - b.value) < eps;
	}
	node& operator = (const node &b) {
		if (this == &b) return *this;
		Xdim = b.Xdim;
        Fdim = b.Fdim;
		for (int i = 0; i<b.Xdim; i++) {
			ps[i] = b.ps[i];
		}
        for (int i = 0; i<b.Fdim; i++) {
			pf[i] = b.pf[i];
		}
		value = b.value;
		return *this;
	}
	node operator + (const node &b) const {
		node z;
		if (Xdim != b.Xdim) return z;
		z.Xdim = Xdim;
        z.Fdim = Fdim;
		for (int i = 0; i<Xdim; i++) {
			z.ps[i] = ps[i] + b.ps[i];
		}
		z.update();
		return z;
	}
	node operator - (const node &b) const {
		node z;
		if (Xdim != b.Xdim) return z;
		z.Xdim = Xdim;
        z.Fdim = Fdim;
		for (int i = 0; i<Xdim; i++) {
			z.ps[i] = ps[i] - b.ps[i];
		}
		z.update();
		return z;
	}
	node operator * (const double &a) const {
		node z;
		z.Xdim = Xdim;
        z.Fdim = Fdim;
		for (int i = 0; i<Xdim; i++) {
			z.ps[i] = a * ps[i];
		}
		z.update();
		return z;
	}
	void show() {
		ct(Xdim);
		for (int i = 0; i<Xdim; i++) {
			cout << ps[i] << " ";
		}cout << endl;
		ct(value);
		//cout<< #Xdim <<":  "<<Xdim<<endl;
	}
};

set<node> pop;  // store simplex points
set<node> xx;  // store points produce by simplex
struct node sub;  // record the point of subproblem
struct node best;
struct node soso;
struct node worst;
struct node expectation;

void worst_improve(){
    set<node> ::iterator it;
    double tmp[50];
    int ii;
    for (it = pop.begin(), ii = 0; it != pop.end(); it++, ii++){
        if (it == pop.begin()){
            best = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] = best.ps[j];
            }
            /*best.show();
            expectation.show();*/
        }
        else{
            expectation = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] += expectation.ps[j];
            }
            //expectation.show();
        }
        if (ii == NData - 2)
            soso = *it;
        if (ii == NData - 1){
            worst = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] -= worst.ps[j];
            }
        }

    }
    for (int j = 0; j<XDim; j++){
        tmp[j] = tmp[j] / (1.0*NData-1);
        expectation.ps[j] = tmp[j];
    }
    expectation.update();
    xx.insert(expectation);
    fes++;
}

void sub_improve(){
    set<node> ::iterator it;
    double tmp[50];
    int ii;
    for (it = pop.begin(), ii = 0; it != pop.end(); it++, ii++){
        if (it == pop.begin()){
            best = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] = best.ps[j];
            }
            /*best.show();
            expectation.show();*/
        }
        else{
            expectation = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] += expectation.ps[j];
            }
            //expectation.show();
        }
        if (ii == NData - 2)
            soso = *it;
        if (ii == NData - 1){
            worst = sub;
            for (int j = 0; j<XDim; j++){
                tmp[j] -= worst.ps[j];
            }
        }

    }
    for (int j = 0; j<XDim; j++){
        tmp[j] = tmp[j] / (1.0*NData-1);
        expectation.ps[j] = tmp[j];
    }
    expectation.update();
    xx.insert(expectation);
    fes++;
}

void rand_improve(){
    set<node> ::iterator it;
    double tmp[50];
    int ii;
    for (it = pop.begin(), ii = 0; it != pop.end(); it++, ii++){
        if (it == pop.begin()){
            best = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] = best.ps[j];
            }
            /*best.show();
            expectation.show();*/
        }
        else{
            expectation = *it;
            for (int j = 0; j<XDim; j++){
                tmp[j] += expectation.ps[j];
            }
            //expectation.show();
        }
        if (ii == NData - 2)
            soso = *it;
        if (ii == NData - 1){
            worst = sub;
            for (int j = 0; j<XDim; j++){
                tmp[j] -= worst.ps[j];
            }
        }

    }
    for (int j = 0; j<XDim; j++){
        tmp[j] = tmp[j] / (1.0*NData-1);
        expectation.ps[j] = tmp[j];
    }
    expectation.update();
    xx.insert(expectation);
    fes++;
}

void simplex(int flag){
	set<node> ::iterator it;
	if(flag==0){
        worst_improve();
    }
    if(flag==1){
        sub_improve();
    }
    if(flag==2){
        rand_improve();
    }
    struct node x_r;
    x_r = expectation*(1 + alpha) - worst*alpha;
    x_r.polymutate();
    x_r.update();
    xx.insert(x_r);
    fes++;
    if (best.value <= x_r.value && x_r.value <= soso.value){
        it = pop.end();
        pop.erase(--it);
        pop.insert(x_r);
        return ;
    }
    if (x_r.value < best.value){
        struct node x_e;
        x_e = expectation*(1 + alpha*gamma) - worst*(alpha*gamma);
        x_e.polymutate();
        x_e.update();
        xx.insert(x_e);
        fes++;
        if (x_e.value < x_r.value){
            it = pop.end();
            pop.erase(--it);
            pop.insert(x_e);
            return ;
        }
        else{
            it = pop.end();
            pop.erase(--it);
            pop.insert(x_r);
            return ;
        }
    }
    if (soso.value <= x_r.value && x_r.value < worst.value){
        struct node x_co;
        x_co = expectation*(1 + alpha*beta) - worst*(alpha*beta);
        x_co.polymutate();
        x_co.update();
        xx.insert(x_co);
        fes++;
        if (x_co.value <= x_r.value){
            it = pop.end();
            pop.erase(--it);
            pop.insert(x_co);
            return ;
        }
    }
    if (worst.value <= x_r.value){
        struct node x_ci;
        x_ci = expectation*(1 - beta) + worst*beta;
        x_ci.polymutate();
        x_ci.update();
        xx.insert(x_ci);
        fes++;
        if (x_ci.value < worst.value){
            it = pop.end();
            pop.erase(--it);
            pop.insert(x_ci);
            return ;
        }
    }
	return ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[])
{
	int i, j;
    set<node> ::iterator it;
	// Check for proper number of arguments.
	if (nrhs != 9)
	{
		mexErrMsgTxt("Nine input required.");
	}
	else if (nlhs != 5)
	{
		mexErrMsgTxt("Five output arguments");
	}
	// Check for data number
	XDim = mxGetM(prhs[1]);
	NData = mxGetN(prhs[1]);
    
    FDim = mxGetM(prhs[5]);

	problem = mxArrayToString(prhs[6]);
    
	// Copy data.
	double *px = mxGetPr(prhs[1]);
    double *pf = mxGetPr(prhs[8]);
    double *t=mxGetPr(prhs[2]);
    alpha=t[0];
    
    t=mxGetPr(prhs[3]);
    beta=t[0];
    
    t=mxGetPr(prhs[4]);
    gamma=t[0];
    
    W = mxGetPr(prhs[0]);
	ideaPoint = mxGetPr(prhs[5]);
    domain = mxGetPr(prhs[7]);


	for (i = 0; i<NData; i++)
	{
		double *tt, *tf;
        tt = new double[XDim];
        tf = new double[FDim];
		for (j = 0; j<XDim; j++){
            tt[j] = px[i*XDim + j];
        }
        for (j=0;j<FDim;j++){
            tf[j]=pf[i*FDim+j];
        }
		pop.insert(node(tt,XDim,tf,FDim));
        if(i==0){
            sub=node(tt,XDim,tf,FDim);
        }
        delete []tt;
        delete []tf;
	}
    
    fes = 0;
	// Nelder-Mean Simplex
	simplex(0);
    
	// Create matrix for the return arguments.
    double *indx,*indf;
    plhs[0]=mxCreateDoubleMatrix(XDim,NData,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(FDim,NData,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(XDim,fes,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(FDim,fes,mxREAL);
    plhs[4]=mxCreateDoubleMatrix(1,1,mxREAL);
    px = mxGetPr(plhs[0]);
    pf = mxGetPr(plhs[1]);
    indx = mxGetPr(plhs[2]);
    indf = mxGetPr(plhs[3]);
    t = mxGetPr(plhs[4]);
	
	// Copy data to return arguments.
	for (it = pop.begin(), i = 0; it != pop.end(); it++, i++){
        for(j=0;j<XDim;j++) px[i*XDim+j]=(*it).ps[j];
        for(j=0;j<FDim;j++) pf[i*FDim+j]=(*it).pf[j];
        pop.erase(it);
    }
    for(it=xx.begin(),i=0;it!=xx.end();it++,i++){
        for(j=0;j<XDim;j++) indx[i*XDim+j]=(*it).ps[j];
        for(j=0;j<FDim;j++) indf[i*FDim+j]=(*it).pf[j];
        xx.erase(it);
    }
	t[0]=fes;

}