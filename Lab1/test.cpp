#include<iostream>
#include<vector>
#include <algorithm> 

using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>
#define vi vector<int>

void pritntvec(vi arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

int fact(int n, vi &dp){
	if(n==1) return 1;
	if(n==2) return 2;
	if(dp[n]!=1) return dp[n];

	dp[n] = n*fact(n-1,dp);
	return dp[n];
}

int nCm(int n, int m){
	vi dp(n+1,1);
	dp[2]=2;

	int temp = fact(n,dp);
	int temp2 = fact(m,dp);
	int temp3 = fact(n-m,dp);

	return (temp/(temp2*temp3));
}

int main(){
	vi arr(5,1);
	arr[0]=0;
	arr[3]=0;
	sort(arr.begin(),arr.end());
	//pritntvec(arr);

	int num_cases = nCm(5,2);
	for(int i=0;i<num_cases;i++){
		pritntvec(arr);
		next_permutation(arr.begin(), arr.end());
	}
}