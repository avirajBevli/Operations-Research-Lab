//Aviraj Singh Bevli
//18MA20009
//Lab11 submission
//Use the SIMPLEX method to solve a linear optimisation problem generated from a 2 person, zero sum game
#include<bits/stdc++.h>
#include<iostream>
#include<vector>
#include<algorithm>//To use the function next_permutation directly
#include<cmath>//for isnan()

using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>
#define vi vector<int>
#define num_iterations_max 250
#define num_iterations_min 10
#define pb push_back
#define epsilon_threshold 0.00001
#define double_MAXX 1000000
#define zero_threshold 0.0000001
#define counter_overflow_threshold 15

void pritntvec(vi arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

void pritntvecd(vd arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<"  ";
	cout<<endl;
}

void printmat(vvd mat){
	for(int j=1;j<mat[0].size();j++){
		cout.width(13); cout<<"y"<<j;
	}
	cout.width(13); cout<<"RHS";
	cout<<endl;
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			cout.width(13);
			cout<<mat[i][j]<<" ";
		}
		cout<<endl;
	}
}

void print_tableau(vvd tableau){
	cout.width(13); cout<<"C0";
	cout.width(13); cout<<"Basis";
	for(int j=3;j<tableau[0].size();j++){
		cout.width(12); cout<<"x"<<j-2;
	}
	cout.width(13); cout<<"RHS";
	cout<<endl;

	for(int i=0;i<tableau.size();i++){
		for(int j=0;j<tableau[i].size();j++){
			cout.width(13);
			if(j==1){
				cout<<tableau[i][j]-1;
			}
			else{
				cout<<tableau[i][j];
			}
		}
		cout<<endl;
	}
}


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

//Function to be used to decide when to terminate in solve_gauss_siedel
double find_magnitude(vd arr){
	double ret=0;
	for(int i=0;i<arr.size();i++){
		ret+=(arr[i]*arr[i]);
	}
	return ret;
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

//Calculate the value of the objective function for the current values of the variables
double find_obj_fn_value(vd obj_fn, vd ans, vi index_set){
	int n = index_set.size();
	double ret=0;
	int count=0;
	for(int i=0;i<n;i++){
		if(index_set[i]==0)
			continue;

		if(ans[count]<0 || isnan(ans[count])){
			//cout<<"This is an infeasible solution"<<endl<<endl;
			return -1;
		}
		ret+=(obj_fn[i]*ans[count++]);
	}

	if(ret<0) return -1;
	count=0;
	for(int i=0;i<n;i++){
		if(index_set[i]==0) cout<<"0(NB_Var) ";//NB_Var = Non-basic variable
		else cout<<ans[count++]<<" ";
	}

	cout<<",  Objective_fn value: "<<ret<<endl<<endl;
	return ret;
}

//Modify the objective function to take the slack variables into account
void create_maximisation_objective_function(vd &objective_fn, bool is_max, int n){
	for(int i=objective_fn.size();i<n;i++) objective_fn.pb(0);
	if(is_max) return;
	cout<<"objective_fn_size: "<<objective_fn.size()<<endl;
	cout<<"n: "<<n<<endl;
	for(int i=0;i<objective_fn.size();i++) objective_fn[i] = -1*objective_fn[i];
}

//Create the full coefficient matrix which includes the slack variables as well
void create_mat(vi eqn_types, vvd &mat){
	int m = mat.size();
	int n = mat[0].size();

	vd last_column(m,0);
	for(int i=0;i<m;i++){
		last_column[i] = mat[i][n-1];
		mat[i].pop_back();
	}

	//Add the slack variables to the matrix
	for(int i=0;i<m;i++){
		if(eqn_types[i]==1){
			n++;
			for(int j=0;j<m;j++){
				if(j==i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}

		else if(eqn_types[i]==2){
			n++;
			for(int j=0;j<m;j++){
				if(j==i) mat[j].pb(-1);
				else mat[j].pb(0);
			}
		}
	}

	for(int i=0;i<m;i++){
		mat[i].pb(last_column[i]);
	}
}

bool can_simplex_be_used(vi eqn_types){
	for(int i=0;i<eqn_types.size();i++){
		if(eqn_types[i]!=1){cout<<endl<<"SORRY, Simplex method cant be used to solve this system!"<<endl; return 0;}
	}
	cout<<"Using simplex method to solve optimisation question"<<endl;
	return 1;
}

vvd construct_tableau(vvd mat){
	int num_rows = mat.size();
	int num_cols = mat[0].size();
	vector<vector<double>> tableau(num_rows , vector<double> (num_cols+2, 0));//2 extra columns for C0, Basis
	for(int i=0;i<num_rows;i++){
		tableau[i][0] = 0;
		tableau[i][1] = num_cols-num_rows+i+1;//column index of the basis variable inside the tableau
	}
	for(int i=0;i<num_rows;i++){
		for(int j=0;j<num_cols;j++) tableau[i][j+2] = mat[i][j];
	}
	return tableau;
}

vd calculate_deviations(vvd tableau, vd objective_fn){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	vd ret(num_cols-3,0);
	
	vd col(num_rows,0);
	for(int i=0;i<num_rows;i++) col[i] = tableau[i][0];

	for(int col_index = 2; col_index < num_cols-1; col_index++){
		vd temp(num_rows,0);		
		for(int i=0;i<num_rows;i++) temp[i] = tableau[i][col_index];

		ret[col_index-2] = objective_fn[col_index-2] - find_dot_prod(temp, col);
		temp.clear();
	}
	return ret;
}

bool is_termination_reached(vd deviations){
	for(int i=0;i<deviations.size();i++){
		if(deviations[i]>zero_threshold) return 0;
	}
	return 1;
}

int find_max_elem_index(vd arr){
	double max = 0;
	int max_index=-1;
	for(int i=0;i<arr.size();i++){
		if(arr[i]>max){
			max = arr[i];
			max_index = i;
		}
	}
	return max_index;
}

// R1 -> R1 - dR2
void subtract_rows(vvd &tableau, double d, int r1, int r2){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=2;i<num_cols;i++) tableau[r1][i] = tableau[r1][i] - d*tableau[r2][i];
}

// R1 -> R1/d
void divide_row(vvd &tableau, double d, int r){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=2;i<num_cols;i++) tableau[r][i] = tableau[r][i]/d;
}

bool do_multiple_solutions_exist(vd deviations, int num_rows){
	int num_zeros = 0;
	for(int i=0;i<deviations.size();i++){
		if(( deviations[i] > (-1*zero_threshold) ) && ( deviations[i] < zero_threshold )) num_zeros++;
	}
	if(num_zeros > num_rows) return 1;
	return 0;
}

void convert_small_vals_to_zeros(vvd &tableau){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=0;i<num_rows;i++){
		for(int j=2;j<num_cols;j++){
			if( (tableau[i][j] < zero_threshold) && (tableau[i][j] > (-1*zero_threshold)) ){
				tableau[i][j]=0;
			}
		}
	}
}

void solve_with_simplex(vvd tableau, vd objective_fn){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();

	cout<<endl<<endl<<"Beginning Simplex method!!"<<endl; int count=0;
	while(1){
		convert_small_vals_to_zeros(tableau);

		if(count > counter_overflow_threshold){cout<<"Excess iterations, returning!!...."<<endl; return;}
		cout<<endl<<endl<<"Iteration: "<<count<<endl;
		cout<<"The Tableau: "<<endl; print_tableau(tableau);

		vd deviations = calculate_deviations(tableau, objective_fn);
		cout<<"Deviations: "; pritntvecd(deviations);
		bool is_end = is_termination_reached(deviations);
		if(is_end){cout<<"Reached the termination state since all deviations are non positive"<<endl;break;}

		//Identifying the entering variable
		int entering_var_col_index = find_max_elem_index(deviations) + 2;
		if(entering_var_col_index==-1){cout<<"LAPHDA!!"<<endl; break;}

		//Identifying the leaving variable
		int leaving_var_row_index=-1;
		double min_ratio = double_MAXX;
		for(int i=0;i<num_rows;i++){
			if(tableau[i][entering_var_col_index] <= 0) continue;
			double temp = tableau[i][num_cols-1]/tableau[i][entering_var_col_index];
			if(temp < min_ratio){
				min_ratio = temp;
				leaving_var_row_index = i;
			}
		}
		if(leaving_var_row_index==-1){cout<<"Unbounded solution!!!....EXITING..."<<endl; return;}
		cout<<"leaving_var_row_index: "<<leaving_var_row_index<<", entering_var_col_index: "<<entering_var_col_index<<endl;


		//Perform row operations on the tableau to get the tableau for the next iteration ready
		for(int row_index = 0; row_index < num_rows; row_index++){
			if(row_index == leaving_var_row_index) divide_row(tableau, tableau[row_index][entering_var_col_index], row_index);
			else{
				double ratio_temp = tableau[row_index][entering_var_col_index]/tableau[leaving_var_row_index][entering_var_col_index];
				subtract_rows(tableau, ratio_temp, row_index, leaving_var_row_index);
			}
		}

		//Modify 1st and 2nd columns of the tableau
		tableau[leaving_var_row_index][0] = objective_fn[entering_var_col_index-2];
		tableau[leaving_var_row_index][1] = entering_var_col_index;

		deviations.clear();
		count++;
	}

	int num_vars = num_cols - 3;
	vd variable_vals(num_vars,0);
	for(int i=0;i<num_rows;i++)	variable_vals[ tableau[i][1]-2 ] = tableau[i][num_cols-1];

	cout<<endl<<"The Optimisation function is optimised at the point : "; 
	for(int i=0;i<num_rows;i++){
		cout<<"x"<<tableau[i][1]-1<<" = "<<tableau[i][num_cols-1]<<", ";
	}
	cout<<"Rest all xis=0"<<endl;

	double optimal_value = find_dot_prod(objective_fn, variable_vals);	
	cout<<"The value of the objective function at this point is : "<<optimal_value<<endl;


	//CHECKING FOR MULTIPLE SOLUTIONS
	/*vd deviations = calculate_deviations(tableau, objective_fn);
	if(do_multiple_solutions_exist(deviations, tableau.size())){
		cout<<"We are at an optimal point and there are non-basic variables with reduced cost equal to 0, so there are multiple values for the decision variables that allow obtaining the optimal value"<<endl;
		cout<<"Along with the mentioned point, infinite other points exist s.t. the optimisation function is optimised"<<endl;
		return;
	}*/

	cout<<endl<<endl<<"Done solving the obtained LPP by simplex method"<<endl;
	cout<<"The optimal strategies of the 2 players are the following: "<<endl;
	cout<<"For Player B -> ";
}


void modify_payoff(vvd &mat){
	int m = mat.size()-1;
	int n = (mat[0].size())-1;

	//Find and Store row_min, col_max
	for(int i=0;i<m;i++){
		//find row min
		double row_min = DBL_MAX;
		for(int j=0;j<n;j++){
			if(mat[i][j] < row_min){
				row_min = mat[i][j];
			}
		}

		//store row min
		mat[i][n]=row_min;
	}

	for(int j=0;j<n;j++){
		//find col max
		double col_max = DBL_MIN;
		for(int i=0;i<m;i++){
			if(mat[i][j] > col_max){
				col_max = mat[i][j];
			}
		}

		//store col max
		mat[m][j]=col_max;
	}
}

void print_payoff(vvd mat){
	int rows = mat.size();
	int cols = mat[0].size();

	cout.width(14);cout<<"";
	cout.width(13); cout<<"y1";
	for(int j=2;j<cols;j++){
		cout.width(12); cout<<"y"<<j;
	}
	cout.width(13); cout<<"row_min";
	cout<<endl;

	for(int i=0;i<rows;i++){ 
		if(i==rows-1){
			cout.width(14);
			cout<<"col_max";
		}
		else{
			cout.width(13);
			cout<<"x"<<i;
		}
		for(int j=0;j<cols;j++){
			cout.width(13);
			if(i==rows-1 && j==cols-1){
				cout<<" ";
			}
			else{
				cout<<mat[i][j]<<"";
			}
		}
		cout<<endl;
	}
}

vvd develop_pay_off_matrix(int m, int n){
	vvd ret(m+1, vd (n+1,-1));
	cout<<"Enter the payoff matrix ("<<m<<" rows, "<<n<<" cols)"<<endl;
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++) cin>>ret[i][j];
	}
	modify_payoff(ret);
	return ret;
}

double find_k1(vvd payoff_mat){
	int rows = payoff_mat.size();
	int cols = payoff_mat[0].size();
	double ret = DBL_MIN;
	for(int i=0;i<rows-1;i++){
		if(payoff_mat[i][cols-1] > ret){
			ret = payoff_mat[i][cols-1];
		}
	}
	return ret;
}

double find_k2(vvd payoff_mat){
	int rows = payoff_mat.size();
	int cols = payoff_mat[0].size();
	double ret = DBL_MAX;
	for(int j=0;j<cols-1;j++){
		if(payoff_mat[rows-1][j] < ret){
			ret = payoff_mat[rows-1][j];
		}
	}
	return ret;
}


vvd convert_into_LPP(vvd payoff_mat, vd &objective_function){
	int rows = payoff_mat.size();
	int cols = payoff_mat[0].size();
	for(int i=0;i<cols-1;i++) objective_function.pb(1);

	double k1 = find_k1(payoff_mat);
	double k2 = find_k2(payoff_mat);

	double c = max(k1,k2)+2;//randomly pick a c that is bigger than max(k1,k2)

	vvd mat(rows-1, vd(cols, -1));
	for(int i=0;i<rows-1;i++){
		for(int j=0;j<cols;j++){
			if(j==cols-1) mat[i][j]=1;
			else mat[i][j] = payoff_mat[i][j]+c;
		}
	}

	return mat;
}

int main(){
	int n, m;
	cout<<"Choose number of choices for player A (m): "; cin>>m;//num_eqns
	cout<<"Choose number of choices for player B (n): "; cin>>n;//num_vars
	vvd payoff_mat = develop_pay_off_matrix(m,n);

	cout<<"The Payoff tableau: "<<endl;
	print_payoff(payoff_mat);

	vd objective_function;
	vvd mat = convert_into_LPP(payoff_mat, objective_function);
	
	cout<<endl<<endl<<"The coeeficient matrix: "<<endl;
	printmat(mat);

	vvd tableau = construct_tableau(mat);
	cout<<endl<<endl<<"The Simplex tableau: "<<endl;
	print_tableau(tableau);

	solve_with_simplex(tableau, objective_function);
}


/*
3
3
8 4 2 1
2 8 4 1
1 2 8 1
1
1
1
1 1 1 
1
*/