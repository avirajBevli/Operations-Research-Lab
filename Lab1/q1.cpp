//Aviraj Singh Bevli
//18MA20009
//Contains only gauss seidel method for n equations, n variables

#include<iostream>
#include<vector>
using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>

double find_dot_prod(vd arr1, vd arr2){
	int n = arr1.size();
	double ret=0;
	for(int i=0;i<n;i++)
		ret+=(arr1[i]*arr2[i]);
	return ret;
}

void pritntvecd(vd arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

vd solve_gauss_siedel(vvd mat){
	vd temp_arr(mat.size(), 0); //All variables are initialised to zero
	int n = mat.size();

	int num_iterations=10;
	for(int i=0;i<num_iterations;i++){
		for(int j=0;j<n;j++){
			double dot = find_dot_prod(temp_arr, mat[j]);
			dot-=(mat[j][j]*temp_arr[j]);
			temp_arr[j] = (mat[j][n] - dot)/mat[j][j];
		}
		pritntvecd(temp_arr);
	}	

	return temp_arr;
}

void printans(vd vec){
	int n = vec.size();
	for(int i=0;i<n;i++) cout<<vec[i]<<" ";
	cout<<endl;
}

int main(){
	int n;
	cout<<"(n): "; cin>>n;
	vector<vector<double>> mat(n , vector<double> (n+1, 0)); //n+1 for RHS of all n equations

	for(int i=0;i<n;i++){
		cout<<"Enter equation "<<i<<" coefficients: "<<endl;
		for(int j=0;j<n+1;j++){
			cin>>mat[i][j];
		}
	}

	vd ans = solve_gauss_siedel(mat);
	cout<<"Basic solutions: ";
	printans(ans);
}