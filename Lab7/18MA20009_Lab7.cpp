//Aviraj Singh Bevli
//18MA20009
//Lab7 submission
//REVISED SIMPLEX METHOD(Can handle any type of linear constraints and linear maximisation/minimisation objective function)

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
#define zero_threshold 0.000001
#define counter_overflow_threshold 15
#define M_value 1000000

///////////////////////////
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
}
///////////////////////////

//Create the full coefficient matrix which includes the slack variables as well
void create_mat(vi eqn_types, vvd &mat, vi &slack_list, vi &surplus_list, vi &artificial_list){
	int m = mat.size();
	int n = mat[0].size();

	vd last_column(m,0);
	for(int i=0;i<m;i++){
		last_column[i] = mat[i][n-1];
		mat[i].pop_back();
	}
	n--;//n is the number of variables till now

	//Add the slack, surplus, artificial variables to the matrix
	for(int i=0;i<m;i++){
		if(eqn_types[i]==1){
			slack_list.pb(n+1); n++;
			for(int j=0;j<m;j++){
				if(j==i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}

		else if(eqn_types[i]==2){
			surplus_list.pb(n+1);
			artificial_list.pb(n+2);
			n=n+2;//add surplus and artificial variable
			for(int j=0;j<m;j++){
				if(j==i){mat[j].pb(-1); mat[j].pb(1);}//subtract surplus and add artificial variable 
				else {mat[j].pb(0); mat[j].pb(0);}
			}
		}

		else if(eqn_types[i]==3){
			artificial_list.pb(n+1); n++;//add artificial variable
			for(int j=0;j<m;j++){
				if(j==i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}
	}

	//Add the RHS of the contraints back again into the matrix
	for(int i=0;i<m;i++){
		mat[i].pb(last_column[i]);
	}
}

//Modify the objective function to take the slack variables into account
void create_maximisation_objective_function(vd &objective_fn, bool is_max, int n, vi slack_list, vi surplus_list, vi artificial_list){
	for(int i=objective_fn.size();i<n;i++) objective_fn.pb(0);
	//cout<<"objective_fn_size: "<<objective_fn.size()<<endl;
	cout<<endl<<endl<<"Total number of variables: "<<n<<endl;
	if(!is_max){
		for(int i=0;i<objective_fn.size();i++) objective_fn[i] = -1*objective_fn[i];
	}
	
	//Use artificial variables to modify the objective function to include the bigM terms
	for(int i=0;i<artificial_list.size();i++)
		objective_fn[artificial_list[i]-1] = -1*M_value;
}

void convert_to_standard_form(vvd &mat, vi &eqn_types){
	int n = mat[0].size() - 1;
	for(int i=0;i<mat.size();i++){
		if(mat[i][n] < 0){
			for(int j=0;j<n+1;j++) mat[i][j]*=-1;
			if(eqn_types[i]==1) eqn_types[i]=2;
			else if(eqn_types[i]==2) eqn_types[i]=1; 
		}
	}
}

vd find_prod_vec_mat(vd vec, vvd mat){
	int n = vec.size();
	vd ret(n);
	for(int i=0;i<n;i++){
		double sum=0;
		for(int j=0;j<n;j++){
			sum+=(vec[j]*mat[j][i]);
		}
		ret[i]=sum;
	}
	return ret;
}

vd calc_simplex_multipliers(vvd tableau, vd objective_fn){
	vd cb_vec;
	for(int i=0;i<tableau.size();i++){
		cb_vec.push_back(objective_fn[tableau[i][0]-1]);
	}

	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	vvd B_inv(tableau.size(), vd (num_cols-2,0));
	for(int i=0;i<num_rows;i++){
		for(int j=1;j<num_cols-1;j++) B_inv[i][j-1] = tableau[i][j];
	}

	cout<<"cb_vec: "; pritntvecd(cb_vec);
	cout<<"B_inv: "<<endl;
	printmat(B_inv);
	
	return find_prod_vec_mat(cb_vec, B_inv);
}

double calculate_cjbar(double cj, vd pi_vec, vd pj){
	double ret=0;
	for(int i=0;i<pj.size();i++){
		ret+=(pi_vec[i]*pj[i]);
	}
	ret*=-1;
	ret+=cj;
	return ret;
}

bool is_artificial_variable(int variable, vi artificial_list){
	for(int i=0;i<artificial_list.size();i++){
		if(artificial_list[i] == variable) return 1;
	}
	return 0;
}

vd find_pi(int index, vvd mat){
	vd ret(mat.size());
	for(int i=0;i<mat.size();i++){
		ret[i] = mat[i][index];
	}
	return ret;
}

vi find_basic_vars_vec(vvd tableau){
	vi ret;
	for(int i=0;i<tableau.size();i++) ret.push_back(tableau[i][0]);
	return ret;
}

bool is_in_basic(int index, vi basic_vars){
	for(int i=0;i<basic_vars.size();i++){
		if(basic_vars[i] == index) return 1;
	}
	return 0;
}

int find_entering_var(vd pi_vec, vd objective_fn, vvd mat, vvd tableau, vi artificial_variable_list){
	int total_vars = objective_fn.size();
	int curr_tableau_row_index=0;
	vector<pair<int,double>> cj_bars;

	for(int i=0;i<total_vars;i++){
		if(is_artificial_variable(i+1, artificial_variable_list)) continue;
		vd pi = find_pi(i, mat);
		if(curr_tableau_row_index < tableau.size()){
			//if((tableau[curr_tableau_row_index][0]-1) == i){
			//	curr_tableau_row_index++;
			//	continue;
			//}
			vi basic_vars = find_basic_vars_vec(tableau);
			if(is_in_basic(i+1,basic_vars)) continue;
			else{
				pair<int, double> temp;
				temp.first = i+1;//wihch variable(one based indexing)
				temp.second = calculate_cjbar(objective_fn[i], pi_vec, pi);
				cj_bars.push_back(temp);
			}
		}
		else{	
			pair<int, double> temp;
			temp.first = i+1;//wihch variable(one based indexing)
			temp.second = calculate_cjbar(objective_fn[i], pi_vec, pi);
			cj_bars.push_back(temp);
		}
	}
	cout<<"CJ_bars are: ";
	for(int i=0;i<cj_bars.size();i++){
		cout<<"("<<cj_bars[i].first<<","<<cj_bars[i].second<<") ";
	} cout<<endl;

	//find the Cj_bar that is the most positive, that is the entering variable
	int max=0;
	int max_var_index=-1;
	for(int i=0;i<cj_bars.size();i++){
		if(cj_bars[i].second > max){
			max = cj_bars[i].second;
			max_var_index = cj_bars[i].first;
		}
	}

	if(max_var_index == -1)
		cout<<endl<<endl<<"Terminating since all CJs are negative!"<<endl;
	return max_var_index;
}

vd find_pivot_column(vvd tableau, vd Pi){
	int n = Pi.size();
	vd ret(n);
	for(int i=0;i<n;i++){
		double sum=0;
		for(int j=0;j<n;j++) sum+=(tableau[i][j+1]*Pi[j]);
		ret[i]=sum;
	}
	return ret;
}

int find_leaving_var_row_index(vvd tableau, vd pivot_column, int ev){
	int num_cols = tableau[0].size();
	double min_ratio = double_MAXX;
	int lv=-1;

	for(int i=0;i<tableau.size();i++){
		if(pivot_column[i] < 0) continue;
		if(((pivot_column[i] > -1*zero_threshold)) && (pivot_column[i] < zero_threshold)) continue;
		double temp = (tableau[i][num_cols-1])/(pivot_column[i]);
		if(temp < min_ratio){
			lv = i;
			min_ratio = temp;
		}
	}
	if(lv==-1) cout<<endl<<endl<<"Couldn't find a positive ratio, hence cant decide leaving variable. Hence the solution of the BVP is unbounded"<<endl;
	return lv;
}

// R1 -> R1 - dR2
void subtract_rows(vvd &tableau, double d, int r1, int r2){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=1;i<num_cols;i++) tableau[r1][i] = tableau[r1][i] - d*tableau[r2][i];
}

// R1 -> R1/d
void divide_row(vvd &tableau, double d, int r){
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=1;i<num_cols;i++) tableau[r][i] = tableau[r][i]/d;
}

void update_tableau(vvd &tableau, int ev, int pivot_row_index, vd pivot_column){
	for(int i=0;i<tableau.size();i++){
		if(i==pivot_row_index){
			continue;
		}
		else{
			double d = pivot_column[i]/pivot_column[pivot_row_index];
			subtract_rows(tableau, d, i, pivot_row_index);
		}
	}
	divide_row(tableau, pivot_column[pivot_row_index], pivot_row_index);
	tableau[pivot_row_index][0] = ev;
}

void solve_with_revised_simplex(vvd mat, vd objective_fn, vi slack_variable_list, vi surplus_variable_list, vi artificial_variable_list, bool is_maximisation){
	//Create Tableau
	int num_rows = mat.size();
	int num_cols = slack_variable_list.size() + artificial_variable_list.size() + 2;
	vvd tableau(num_rows, vd (num_cols, 0));
	int last_col_index = mat[0].size()-1;

	for(int i=0;i<slack_variable_list.size();i++) tableau[i][0] = slack_variable_list[i];
	for(int i=0;i<artificial_variable_list.size();i++) tableau[i+slack_variable_list.size()][0] = artificial_variable_list[i];

	//Initially the B matrix will be just identity 
	for(int i=0;i<num_rows;i++){
		tableau[i][i+1] = 1;
		tableau[i][num_cols-1] = mat[i][last_col_index];
	}

	int it=0;
	while(it < counter_overflow_threshold){
		cout<<endl<<endl<<"Tableau "<<it<<endl; 
		cout<<"//////////////////"<<endl; printmat(tableau); cout<<"//////////////////"<<endl;

		vd pi_vec = calc_simplex_multipliers(tableau, objective_fn);
		cout<<"Simplex Multiplier: "; pritntvecd(pi_vec);

		int ev = find_entering_var(pi_vec, objective_fn, mat, tableau, artificial_variable_list);
		if(ev == -1){
			bool is_there_an_artificial_variable_in_basis=0;
			for(int i=0;i<tableau.size();i++){
				if(is_artificial_variable(tableau[i][0], artificial_variable_list)){
					int num_cols = tableau[0].size(); 
					if((tableau[i][num_cols-1] > zero_threshold) || (tableau[i][num_cols-1] < -1*zero_threshold)) 
						is_there_an_artificial_variable_in_basis = 1;
				}
			}
			if(is_there_an_artificial_variable_in_basis){
				//Means this is an infeasible solution
				cout<<"Since even upon termination a non-zero artificial variable exists, we have an infeasible solution"<<endl;
				return;
			}

			//This tableau is the solution
			cout<<endl<<"We have found the solution. The Optimisation function is optimised at the point : ";
			int num_cols = tableau[0].size(); 
			for(int i=0;i<tableau.size();i++){
				cout<<"x"<<tableau[i][0]<<" = "<<tableau[i][num_cols-1]<<", ";
			}
			cout<<"Rest all xis=0"<<endl;

			///
			double optimal_value = 0;
			for(int i=0;i<tableau.size();i++){
				if(tableau[i][0] <= objective_fn.size()){
					optimal_value+=(tableau[i][num_cols-1]*objective_fn[tableau[i][0]-1]);
				}
			}	
			if(!is_maximisation) optimal_value*=-1;
			cout<<"The value of the objective function at this point is : "<<optimal_value<<endl;
			///

			return;
		}
		cout<<"Entering varable: "<<ev<<endl;
		
		vd Pi = find_pi(ev-1, mat);
		vd pivot_column = find_pivot_column(tableau, Pi);
		cout<<"Pivot column: "; pritntvecd(pivot_column);

		int pivot_row_index = find_leaving_var_row_index(tableau, pivot_column, ev);
		if(pivot_row_index == -1){
			return;
		}
		cout<<"Leaving variable: "<<tableau[pivot_row_index][0]<<endl;

		//Update leaving, entering variable, perform elemenratry row operations on the tableau
		update_tableau(tableau, ev, pivot_row_index, pivot_column);
		it++;
	}
	if(it>=counter_overflow_threshold){
		cout<<"Since even after many iterations, still have not converged, we have an infeasible solution"<<endl;
		return;
	}
}

int main(){
	int n, m;
	cout<<"(m): "; cin>>m;//num_eqns
	cout<<"(n): "; cin>>n;//num_vars
	vector<vector<double>> mat(m , vector<double> (n+1, 0)); 

	for(int i=0;i<m;i++){
		cout<<"Enter Constraint "<<i<<" coefficients: "<<endl;
		for(int j=0;j<n+1;j++){
			cin>>mat[i][j];
		}
	}
	cout<<endl<<"Enter whether the equations are <= type(1), >= type(2), or = type(3): "<<endl;
	vi eqn_types(m,0);
	for(int i=0;i<m;i++){
		cout<<"Equation "<<i<<": "; cin>>eqn_types[i]; 
	}

	//Handle the cases when the RHS of an equation is negaitve
	convert_to_standard_form(mat, eqn_types);

	vi slack_variable_list; vi surplus_variable_list; vi artificial_variable_list;
	create_mat(eqn_types, mat, slack_variable_list, surplus_variable_list, artificial_variable_list);	

	printmat(mat);
	cout<<"slack variables: "; pritntvec(slack_variable_list);
	cout<<"surplus variables: "; pritntvec(surplus_variable_list);
	cout<<"artifiical variables: "; pritntvec(artificial_variable_list);	

	cout<<endl<<"Enter the objective function: "<<endl;
	vd objective_fn(n,0);
	for(int i=0;i<n;i++){
		cin>>objective_fn[i];
	}
	bool is_maximisation;
	cout<<endl<<"Does the objective_fn have to be maximised(1) or minimised(0)? : ";
	cin>>is_maximisation;
	n = mat[0].size();
	create_maximisation_objective_function(objective_fn, is_maximisation, n-1, slack_variable_list, surplus_variable_list, artificial_variable_list);

	cout<<"The coeficient matrix: "<<endl;
	printmat(mat);
	cout<<"The objective function: "<<endl;
	pritntvecd(objective_fn);
	cout<<endl;


	//Solve with Revised Simplex method
	solve_with_revised_simplex(mat, objective_fn, slack_variable_list, surplus_variable_list, artificial_variable_list, is_maximisation);
}


/*
Inputs for questions ->
(A)
3
3
3 5 2 60
4 4 4 72
2 4 5 100
1
1
1
5 10 8
1

(B)
4
2
2 1 1000
1 1 800
1 0 400
0 1 700
1
1
1
1
4 3
1

(C)
2
4
2 2 5 1 12
3 3 4 0 11
1
1
3 3 2 1
1

(D)
3
3
2 3 0 8
0 2 5 10
3 2 4 15
1
1
1
3 5 4
1

(E)
3
4
2 2 5 1 43
3 -3 4 0 11
4 -2 3 -1 25
2
2
2
3 3 2 -1
0

Another sample input
3
3
1 -2 1 11
-4 1 2 3
2 0 -1 -1
1
2
3
-3 1 1
0

*/