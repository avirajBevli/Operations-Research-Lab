#include<iostream>
#include<vector>
#include<algorithm>//To use the function next_permutation directly
#include<cmath>//for isnan()

using namespace std;

int main(){
	double d = -9.21;
	int ret = (int)d;
	double left = d - ret;
	cout<<ret<<", "<<left<<endl;
}