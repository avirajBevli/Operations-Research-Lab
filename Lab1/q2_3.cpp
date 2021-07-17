//Aviraj Singh Bevli
//18MA20009
//Solve m equations, n variables(m<n) for basic solutions

#include<iostream>
#include<vector>
#include<algorithm>//To use the function next_permutation directly

using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>
#define vi vector<int>
#define num_iterations 25
#define pb push_back


//FOR DEBUGGING
////////////////////////////
void pritntvec(vi arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

void pritntvecd(vd arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

void printmat(vvd mat){
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++)
			cout<<mat[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;
}
///////////////////////////


//Find dot product of 2 n dimensional vectors
double find_dot_prod(vd arr1, vd arr2){
	int n = arr1.size();
	double ret=0;
	for(int i=0;i<n;i++)
		ret+=(arr1[i]*arr2[i]);
	return ret;
}

//Solve a system of n equations, n variables using Gauss-Siedel method using 0 initialisation
vd solve_gauss_siedel(vvd mat){
	vd temp_arr(mat.size(), 0);
	int n = mat.size();

	for(int i=0;i<num_iterations;i++){
		for(int j=0;j<n;j++){
			double dot = find_dot_prod(temp_arr, mat[j]);
			dot-=(mat[j][j]*temp_arr[j]);
			temp_arr[j] = (mat[j][n] - dot)/mat[j][j];
		}
		//pritntvecd(temp_arr);
	}	

	return temp_arr;
}

//Print the obtained basic solution
void printans(vd vec, vi index_set){
	int n = index_set.size();
	int count=0;
	for(int i=0;i<n;i++){
		if(index_set[i]==0) cout<<"0(non-basic var) ";
		else cout<<vec[count++]<<" ";
	}
	cout<<endl;
}

//calculate factorial of n using DP
int fact(int n, vi &dp){
	if(n==1) return 1;
	if(n==2) return 2;
	if(dp[n]!=1) return dp[n];

	dp[n] = n*fact(n-1,dp);
	return dp[n];
}

//Calculate n combination m
int nCm(int n, int m){
	vi dp(n+1,1);
	dp[2]=2;

	int temp = fact(n,dp);
	int temp2 = fact(m,dp);
	int temp3 = fact(n-m,dp);

	return (temp/(temp2*temp3));
}

//Construct a matrix(from the input set of equations) to be used for the Gauss-Siedel method
vvd form_mat(vi index_set, vvd mat){
	int m = mat.size();
	int n = (mat[0].size())-1;
	vvd ret;
	for(int i=0;i<m;i++){
		vd temp;
		for(int j=0;j<n;j++){
			if(index_set[j]==0) continue;
			temp.pb(mat[i][j]);
		}
		temp.pb(mat[i][n]);
		ret.pb(temp);
	}

	return ret;
}

int main(){
	int n, m;
	cout<<"(m): "; cin>>m;//num_eqns
	cout<<"(n): "; cin>>n;//num_vars
	vector<vector<double>> mat(m , vector<double> (n+1, 0)); 

	for(int i=0;i<m;i++){
		cout<<"Enter equation "<<i<<" coefficients: "<<endl;
		for(int j=0;j<n+1;j++){
			cin>>mat[i][j];
		}
	}
	//printmat(mat);

	vi index_set(n,1);
	for(int i=0;i<n-m;i++) index_set[i]=0;

	cout<<endl<<endl<<"Basic solutions: "<<endl;
	vvd mat_temp;

	int num_cases = nCm(n,m);
	for(int i=0;i<num_cases;i++){
		//cout<<endl<<endl<<"i: "<<i<<endl;
		//pritntvec(index_set);
		mat_temp = form_mat(index_set,mat);
		//printmat(mat_temp);
		
		vd ans = solve_gauss_siedel(mat_temp);
		printans(ans,index_set);

		next_permutation(index_set.begin(), index_set.end());
	}
}