//Aviraj Singh Bevli
//18MA20009
//Lab10 solution on Assignment problem(Assigning n number of jobs to m number of opertors such that 
//the total money payabe to operators is minimum
//OR
//the total time used by the operators in minimum

#include<bits/stdc++.h>
using namespace std;
#define pb push_back
#define pi pair<int,int>
#define vpi vector<pair<int,int> >
#define vi vector<int>
#define vvd vector<vector<double> >
#define vb vector<bool>
#define vvb vector<vector<bool> >

//The number of rows = The number of columns in the job-operator assignment matrix after adding dummy rows/columns(if any)
int n_global;

void print2dmat(vvd mat){
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			cout.width(5);
			cout<<mat[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}

//Take inputs of the job-operator assignment table from the user and construct the tableau
vvd construct_talbeau(){
	int num_jobs, num_operators;
	cout<<"Enter the number of jobs: ";
	cin>>num_jobs;
	cout<<"Enter the number of operators: ";
	cin>>num_operators;
	if(num_jobs!=num_operators) 
		cout<<endl<<"This is not a Balanced Assignment problem, will have to make it balanced by introducing dummy variables!"<<endl;
	else
		cout<<endl<<"This is a Balanced Assignment problem"<<endl;
	
	n_global = max(num_jobs, num_operators);
	
	vvd tableau(n_global , vector<double> (n_global, 0));
	cout<<endl<<"Enter job-operator cost table (cell in the ith row and jth column represents the cost of assignment of ith job to the jth operator): "<<endl;
	for(int i=0;i<num_jobs;i++){
		for(int j=0;j<num_operators;j++){
			cin>>tableau[i][j];
		}
	}
	
	return tableau;
}

//Have all zeros in the matrix been covered with lines?
bool are_all_zeros_covered(vvd mat, vb rows_marked, vb cols_marked){
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			if(mat[i][j]==0){
				if(rows_marked[i]==0 && cols_marked[j]==0){
					cout<<"Cell ("<<i<<","<<j<<") identified as a zero not covered with lines"<<endl;
					return 0;
				}
			}
		}
	}
	return 1;
}

//Have all zeros in the matrix been marked?
bool are_all_zeros_marked(vvd mat, vvb bool_mat){
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			if(mat[i][j]==0){
				if(bool_mat[i][j]==0) return 0;
			}
		}
	}
	return 1;
}

//Has optimality condition been reached?
bool is_optimality_reached(vvb bool_mat){
	int count_marked_zeros=0;
	for(int i=0;i<bool_mat.size();i++){
		for(int j=0;j<bool_mat[i].size();j++){
			if(bool_mat[i][j]==1) count_marked_zeros++;
		}
	}
	cout<<"Number of marked cells are "<<count_marked_zeros<<", Number of rows are "<<n_global<<endl;
	if(count_marked_zeros == n_global) return 1;
	return 0;
}

//Print the indices where the bool vector is 1
void printvecb_indices(vb vec){
	for(int i=0;i<vec.size();i++){
		if(vec[i]==1) cout<<i<<", ";
	}
	cout<<endl;
}

//print cells where the bool mat is 1
void printvecvecb_indices(vvb mat){
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			if(mat[i][j]==1) cout<<"("<<i<<","<<j<<"), ";
		}
	}
	cout<<endl;
}

//Find the optimum assignment of jobs to operators (such that the total cost is minimum)
vi find_operator_assignment(vvd mat, vvb cells_marked){
	vi ret_vec(n_global,0);
	for(int i=0;i<n_global;i++){
		for(int j=0;j<n_global;j++){
			if(cells_marked[i][j]){
				ret_vec[i] = j;
			}
		}
	}
	return ret_vec;
}

void do_row_scan(vb &cols_marked, vb &rows_marked, vvb &cells_marked, vvd mat){
	//cout<<"Doing row scanning...";
	for(int i=0;i<n_global;i++){
		int num_zeros = 0;
		for(int j=0;j<n_global;j++){
			if(cols_marked[j]==0){
				if(mat[i][j]==0) 
					num_zeros++;
			}
		}
		if(num_zeros == 1){
			for(int j=0;j<n_global;j++){
				if(cols_marked[j]==0){
					if(mat[i][j]==0) {
						cols_marked[j]=1;
						cells_marked[i][j]=1;
					}
				}
			}
		}
	}
	cout<<"Done Row scanning"<<endl;
}

void do_col_scan(vb &cols_marked, vb &rows_marked, vvb &cells_marked, vvd mat){
	for(int j=0;j<n_global;j++){
		if(cols_marked[j]) continue;

		int num_zeros=0;
		for(int i=0;i<n_global;i++){
			if(rows_marked[i]==0){
				if(mat[i][j]==0)
					num_zeros++;
			}
		}
		if(num_zeros==1){
			for(int i=0;i<n_global;i++){
				if(rows_marked[i]==0){
					if(mat[i][j]==0){
						rows_marked[i]=1;
						cells_marked[i][j]=1;
					}
				}
			}
		}
	}
	cout<<"Done Col scanning"<<endl;
}

void solve_using_Hungerian(vvd tableau){
	vvd mat = tableau;//use it for local computations

	//Phase 1
	cout<<"//////////Phase 1 -> "<<endl;

	//Perform row reductions
	cout<<"Min row elements are: ";
	for(int i=0;i<n_global;i++){
		double row_min = DBL_MAX;
		for(int j=0;j<n_global;j++){
			if(mat[i][j] < row_min){
				row_min = mat[i][j];
			}
		}
		cout<<row_min<<", ";
		for(int j=0;j<n_global;j++){
			mat[i][j]-=row_min;
		}
	}
	cout<<endl<<"Done row reductions, the new mat is -> "<<endl;
	print2dmat(mat);

	//Perform column reductions
	cout<<"Min col elements are: ";
	for(int j=0;j<n_global;j++){
		double col_min = DBL_MAX;
		for(int i=0;i<n_global;i++){
			if(mat[i][j] < col_min){
				col_min = mat[i][j];
			}	
		}
		cout<<col_min<<", ";
		for(int i=0;i<n_global;i++){
			mat[i][j]-=col_min;	
		}
	}
	cout<<endl<<"Done column reductions, the new mat is -> "<<endl;
	print2dmat(mat);

	//Phase 2
	cout<<endl<<"//////////Phase 2 (Optimzation of the problem) -> "<<endl;
	
	vb rows_marked(n_global, 0);
	vb cols_marked(n_global, 0);
	vvb cells_marked(n_global , vector<bool> (n_global, 0));

	int iter_index=0;
	bool does_alternate_optimum_exist = 0;
	while(1){
		iter_index++;
		if(iter_index > 10){
			cout<<"Iterations exceeded, something went wrong !! Breaking!"<<endl;
			break;
		}
		//Initialise rows_marked, cols_marked, cells_marked to all zeros
		for(int i=0;i<n_global;i++){
			rows_marked[i]=0;
			cols_marked[i]=0;
			for(int j=0;j<n_global;j++){
				cells_marked[i][j]=0;
			}
		}

		//(Step 1) -> Draw minimum number of lines to cover all the zeros of mat
		//a) Row scanning
		do_row_scan(cols_marked, rows_marked, cells_marked, mat);

		if(are_all_zeros_covered(mat, rows_marked, cols_marked)){
			cout<<"Since all the zeros are already covered with lines, no need to do column scanning"<<endl;
		}
		else{
			cout<<"All the zeros are not covered with lines, hence Doing column scanning"<<endl;
			//b) Column scanning
			do_col_scan(cols_marked, rows_marked, cells_marked, mat);
		}

		//Check for Alternate Optimum condition
		if(!are_all_zeros_covered(mat, rows_marked, cols_marked)){
			do_row_scan(cols_marked, rows_marked, cells_marked, mat);
			if(are_all_zeros_covered(mat, rows_marked, cols_marked))
				cout<<"Since all the zeros are already covered with lines, no need to do column scanning"<<endl;
			else{
				cout<<"All the zeros are not covered with lines, hence Doing column scanning"<<endl;
				//b) Column scanning
				do_col_scan(cols_marked, rows_marked, cells_marked, mat);
			}
			if(!are_all_zeros_covered(mat, rows_marked, cols_marked)){
				cout<<"Since even after row, column scanning, there are some zeros not covered with lines, means this is the ALTERNATE OPTIMUM CONDITION"<<endl;
				cout<<"We will use random Diagonal selection to break the tie arbitratily"<<endl;
				does_alternate_optimum_exist=1;

				for(int i=0;i<n_global;i++){
					if(rows_marked[i]) continue;
					for(int j=0;j<n_global;j++){
						if(rows_marked[i]) continue;
						if(cols_marked[j]) continue;
						if(mat[i][j]==0 && cells_marked[i][j]==0){
							cols_marked[j]=1;
							cells_marked[i][j]=1;
							break;
						}
					}
				}
			}
		}

		cout<<"Marked rows: "; printvecb_indices(rows_marked);
		cout<<"Marked cols: "; printvecb_indices(cols_marked);
		cout<<"Marked zeros: "; printvecvecb_indices(cells_marked);
		cout<<endl;

		//(Step 2) -> Check for optimality
		if(is_optimality_reached(cells_marked)){
			cout<<"Optimality has been obtained since num_rows is equal to num_cells_marked. Returning optimal solution!"<<endl;
			break;
		}
		else{	
			cout<<"Optimality has not been obtained since num_rows is not equal to num_cells_marked........Continuing iterations"<<endl;
			//break;
		}
		
		//(Step 3) -> Identify the minimum values of the undetected cells
		double min_undetected_val = DBL_MAX;
		pi min_undetected_cell = {-1,-1};
		for(int i=0;i<n_global;i++){
			if(rows_marked[i]) continue;
			for(int j=0;j<n_global;j++){
				if(cols_marked[j]) continue;
				if(mat[i][j] < min_undetected_val){
					min_undetected_val = mat[i][j];
					min_undetected_cell = make_pair(i,j);
				}
			}
		}

		cout<<endl<<"min_undetected_cell : ("<<min_undetected_cell.first<<","<<min_undetected_cell.second<<"), its val: "<<min_undetected_val<<endl;
		
		//(a) Add the minimum undeleted cell value to all the intersection points
		//(b) Subtract the minimum undeleted cell value from all the undeleted cell values
		for(int i=0;i<n_global;i++){
			for(int j=0;j<n_global;j++){
				if(rows_marked[i] && cols_marked[j]){
					mat[i][j]+=min_undetected_val;
				}
				else if((!rows_marked[i]) && (!cols_marked[j])){
					mat[i][j]-=min_undetected_val;
				}
			}
		}
		cout<<endl<<"After modifying the matrix(add minimum value at the intersection points, subtract minimum value from the undeleted cells), the mat now becomes: "<<endl;
		print2dmat(mat);
		//break;
	}

	vi operator_assignment_vec = find_operator_assignment(mat, cells_marked);
	if(does_alternate_optimum_exist)
		cout<<endl<<"Multiple optimal settings exist for this problem, one of which is the following "<<endl;
	else
		cout<<endl<<"Only one optimal setting exists for this problem, which is given by"<<endl;

	cout<<"In the optimal solution, jobs have been allocated to operators as: "<<endl;
	if(operator_assignment_vec.size () != n_global){
		cout<<"SOMETHING IS WRONG BRUHH since all operators have not been assigned jobs"<<endl;
	}

	double total_cost=0;
	for(int i=0;i<n_global;i++){
		int j = operator_assignment_vec[i];
		double cost = tableau[i][j];
		total_cost+=cost;
		cout<<"Job "<<i<<" is assigned to operator "<<j<<" with cost "<<cost;
		if(cost==0) cout<<" {DUMMY}";
		cout<<endl;
	}
	cout<<endl<<"Hence, optimal total cost is "<<total_cost<<endl<<endl;
}

int main(){
	vvd tableau = construct_talbeau();
	cout<<endl<<endl<<"//////////////////////"<<endl<<"The job-operator cost tableau(INPUT) is: "<<endl;
	print2dmat(tableau);

	solve_using_Hungerian(tableau);
}	

/*
Sample question ->
5
5
9 11 14 11 7
6 15 13 13 10
12 13 6 8 8
11 9 10 12 9
7 12 14 10 14

Q1 ->
4
4
8 5 6 2
7 2 1 3
9 7 4 2
9 8 9 6

Q2 ->
5
4
85 70 65 68
93 57 37 97
24 20 25 23
6 1 89 84
10 19 77 38

Q3 ->
4
5
57 33 22 31 98
97 59 22 45 7
56 63 33 84 90
43 97 11 71 96

Q4 ->
5
5
15 86 19 56 54
18 86 19 80 53
1 85 79 35 47
40 86 54 27 99
93 86 78 14 42

*/