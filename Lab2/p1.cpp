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
#define num_iterations_max 250
#define num_iterations_min 10
#define pb push_back
#define epsilon_threshold 0.00001

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
double absol(double a){
	if(a<0) a*=-1;
	return a;
}

//Find dot product of 2 n dimensional vectors
double find_dot_prod(vd arr1, vd arr2){
	int n = arr1.size();
	double ret=0;
	for(int i=0;i<n;i++)
		ret+=(arr1[i]*arr2[i]);
	return ret;
}

double find_magnitude(vd arr){
	double ret=0;
	for(int i=0;i<arr.size();i++){
		ret+=(arr[i]*arr[i]);
	}
	return ret;
}

//Solve a system of n equations, n variables using Gauss-Siedel method using 0 initialisation
vd solve_gauss_siedel(vvd mat){
	vd temp_arr(mat.size(), 0);
	int n = mat.size();

	double epsilon = 100000;
	double prev = find_magnitude(temp_arr);
	int num_iterns = 0;
	while(epsilon > epsilon_threshold){
		for(int j=0;j<n;j++){
			double dot = find_dot_prod(temp_arr, mat[j]);
			dot-=(mat[j][j]*temp_arr[j]);
			temp_arr[j] = (mat[j][n] - dot)/mat[j][j];
		}
		num_iterns++;
		double curr = find_magnitude(temp_arr);
		epsilon = absol(curr - prev);
		if(num_iterns < num_iterations_min){prev = curr; continue;}
		if(num_iterns > num_iterations_max) break;
		if(epsilon < epsilon_threshold) break;
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

double find_obj_fn_value(vd obj_fn, vd ans, vi index_set){
	int n = index_set.size();
	double ret=0;
	int count=0;
	for(int i=0;i<n;i++){
		if(index_set[i]==0)
			continue;

		if(ans[count]<0){
			cout<<"This is an infeasible solution"<<endl<<endl;
			return -1;
		}
		ret+=(obj_fn[i]*ans[count++]);
	}
	cout<<"objective_fn value: "<<ret<<endl<<endl;
	return ret;
}

int main(){
	int n, m;
	cout<<"(m): "; cin>>m;//num_eqns
	cout<<"(n): "; cin>>n;//num_vars
	vector<vector<double>> mat(m , vector<double> (n+1, 0)); 

	for(int i=0;i<m;i++){
		cout<<"Enter Constraint "<<i<<" coefficients(after inlcuding slack, surplus variables): "<<endl;
		for(int j=0;j<n+1;j++){
			cin>>mat[i][j];
		}
	}
	//printmat(mat);

	cout<<"Enter the objective fucntion to be maximised"<<endl;
	vd objective_fn(n,0);
	for(int i=0;i<n;i++){
		cin>>objective_fn[i];
	}

	vi index_set(n,1);
	for(int i=0;i<n-m;i++) index_set[i]=0;

	cout<<endl<<endl<<"Basic solutions: "<<endl;
	vvd mat_temp;

	double max_obj_fn_value = -1;
	int num_cases = nCm(n,m);
	for(int i=0;i<num_cases;i++){
		//cout<<endl<<endl<<"i: "<<i<<endl;
		//pritntvec(index_set);
		mat_temp = form_mat(index_set,mat);
		//printmat(mat_temp);
		
		vd ans = solve_gauss_siedel(mat_temp);
		printans(ans,index_set);
		
		double obj_fn_value = find_obj_fn_value(objective_fn,ans,index_set);
		if(obj_fn_value > max_obj_fn_value) max_obj_fn_value = obj_fn_value;
		
		next_permutation(index_set.begin(), index_set.end());
	}

	cout<<"The objective_fn's optimum value: "<<max_obj_fn_value<<endl;
}