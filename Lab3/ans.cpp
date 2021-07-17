//Aviraj Singh Bevli
//18MA20009
//OR Lab3

//Consider the system of equations to be of the form AX=B
#include<iostream>
#include<vector>
using namespace std;

#define vd vector<double>
#define vvd vector<vector<double>>
#define vi vector<int>
#define pb push_back

void print_2d_mat(vvd tableau){
	for(int i=0;i<tableau.size();i++){
		for(int j=0;j<tableau[i].size();j++) cout<<tableau[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;
}

void printvec(vd vec){
	cout<<"vector is: "<<endl;
	for(int i=0;i<vec.size();i++) cout<<vec[i]<<" ";
	cout<<endl;
}

void take_problem_input(vvd &tableau){
	
	int m = tableau.size()-2;
	int n = tableau[0].size()-m-3;

	vvd A(m, vd(n+m,0));//LHS of all the constraint equations
	vd B(m,0);//RHS of contraints
	vd objective(n,0);//objective function

	cout<<"Enter LHS of all constraint equations(<= type): "<<endl;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cin>>A[i][j];
		}
	}
	cout<<endl<<"A matrix: "<<endl; print_2d_mat(A); cout<<endl;

	cout<<"Enter RHS of the constraint equations: "<<endl;
	for(int i=0;i<m;i++) cin>>B[i];
	//cout<<"B column vector: "<<endl; printvec(B); cout<<endl;
		
	cout<<"Enter objective function: "<<endl;
	for(int i=0;i<n;i++) cin>>objective[i];

	//add slack variables to A matrix
	for(int i=0;i<m;i++){
		for(int j=n;j<n+m;j++) A[i][j]=0;
		A[i][i+n]=1;
	}

	
	//Now we will construct the tableau
	for(int i=0;i<m+2;i++){
		for(int j=0;j<n+m+3;j++) tableau[i][j]=0;
	}
	
	for(int i=0;i<n;i++) tableau[1][2+i] = -1*objective[i];
	for(int i=2;i<m+2;i++) tableau[i][0] = n+i-1;
	for(int i=2;i<n+m+2;i++) tableau[0][i] = i-1;
	for(int i=2;i<m+2;i++){
		for(int j=2;j<n+m+2;j++) tableau[i][j] = A[i-2][j-2];
		tableau[i][n+m+2] = B[i-2];
	}
}

void solve_with_simplex(vvd tableau){
	int m = tableau.size()-2;
	int n = tableau[0].size()-m-3;
	int num_cols = n+m+3;
	int num_rows = m+2;

	double pivot, min_ratio, most_neg;
	/*while(1){
		most_neg=0;
		for(int i=2;i<num_cols-1;i++){

		}
	}*/


	//////////TO BE CONTINUED!!//////////////
	/////////////////////////////////////////
	/////////////////////////////////////////
}

int main(){
	int n,m;
	cout<<"Enter num_vars: "; cin>>n;
	cout<<"Enter num_eqns: "; cin>>m;
	vvd tableau(m+2, vd(n+m+3,0));
	take_problem_input(tableau);//construct the tableau
	cout<<endl<<"tableau: "<<endl; print_2d_mat(tableau);

	//solve_with_simplex(tableau,n,m);
}