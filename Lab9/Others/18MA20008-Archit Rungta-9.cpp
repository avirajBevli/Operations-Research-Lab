#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

// helper function to create 2D arrays of doubles
double **create_2d_doubles(int m, int n)
{
    double **ret = new double *[m]();
    for (int i = 0; i < m; i++)
        ret[i] = new double[n]();
    return ret;
}

// helper function to print 2d arrays
void print_2d(double **ar, int n, int m)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            cout << ar[i][j] << " ";
        cout << "\n";
    }
}

// function to free a 2d array
void free_2d_doubles(double **arr, int m, int n)
{
    for (int i = 0; i < m; i++)
        delete[] arr[i];
    delete[] arr;
}

void print_table(int s, int d, double **cost_table, double **allocation, vector<double> demand, vector<double> supply)
{
    cout << "Printing the cost and current allocation table\n";
    for (int i = 0; i < d; i++)
        cout << setw(9) << "D_" << i + 1;
    cout << setw(10) << "S\n\n";
    for (int i = 0; i < s; i++)
    {
        cout << "S_" << i + 1 << setw(7);;
        for (int j = 0; j < d; j++)
            cout  << cost_table[i][j]<< setw(10);
        cout << "\n"
             << setw(10);
        for (int j = 0; j < d; j++)
            cout << setw(10) << allocation[i][j];
        cout << setw(8) << supply[i] << "\n\n";
    }
    cout << "D" << setw(9);
    for (int i = 0; i < d; i++)
        cout << demand[i]<< setw(10) ;
    cout << "\n";
}

void find_min_index(double **table, vector<int> l1i, vector<int> l2i, int &i, int &j)
{
    int minv = table[l1i[0]][l2i[0]];
    i = l1i[0], j = l2i[0];
    for (int x : l1i)
        for (int y : l2i)
            if (table[x][y] < minv)
                minv = table[x][y], i = x, j = y;
}

double find_cost(int s, int d, double **cost_table, double **allocation) {
    double cost = 0;
    for(int i =0;i<s;i++)
        for(int j =0;j<d;j++)
            cost += cost_table[i][j]*allocation[i][j];
    return cost;
}

double **lcc_method(double **table, vector<double> demand, vector<double> supply, int s, int d)
{
    vector<double> demand_c(demand), supply_c(supply);
    double **ans = create_2d_doubles(s, d);
    vector<int> l1s, l2s;
    for (int i = 0; i < s; i++)
        l1s.push_back(i);
    for (int i = 0; i < d; i++)
        l2s.push_back(i);
    double unmet_demand = 0;
    for (double x : demand)
        unmet_demand += x;
    while (unmet_demand != 0)
    {
        int ss, sd;
        find_min_index(table, l1s, l2s, ss, sd);
        ans[ss][sd] = min(supply[ss], demand[sd]);
        cout << "Minimum index was found to be " << ss+1 << ", " << sd+1<< "\n";
        supply[ss] -= ans[ss][sd];
        demand[sd] -= ans[ss][sd];
        if (supply[ss] == 0)
            l1s.erase(find(l1s.begin(), l1s.end(), ss));
        if (demand[sd] == 0)
            l2s.erase(find(l2s.begin(), l2s.end(), sd));
        unmet_demand -= ans[ss][sd];
        print_table(s,d,table,ans,demand_c,supply_c);
    }
    return ans;
}

double **northwest_corner(double **table, vector<double> demand, vector<double> supply, int s, int d)
{
    vector<double> demand_c(demand), supply_c(supply);
    double **ans = create_2d_doubles(s, d);
    int ss = 0, sd = 0;
    double unmet_demand = 0;
    for (double x : demand)
        unmet_demand += x;
    while (unmet_demand != 0)
    {
        ans[ss][sd] = min(demand[sd], supply[ss]);
        cout << "Current northwest corner is " << ss+1 << ", " << sd+1<< "\n";

        supply[ss] -= ans[ss][sd];
        demand[sd] -= ans[ss][sd];
        unmet_demand -= ans[ss][sd];
        if (supply[ss] == 0)
            ss++;
        if (demand[sd] == 0)
            sd++;
        print_table(s,d,table,ans,demand_c,supply_c);
    }
    return ans;
}

bool get_loop(double **table, double **visited, int &mi, int &mj, int &s, int &d, int i, int j, int sign, vector<tuple<int, int, int>> &path)
{
    if (i < 0 || j < 0 || i >= s || j >= d)
        return false;
    if (mi == i && mj == j)
    {
        if (visited[i][j])
        {
            if (path.size() > 2)
                return true;
            else
                return false;
        }
    }
    else if (visited[i][j])
        return false;
    cout << i << " " << j << "\n";
    visited[i][j] = 1;
    int di =0, dj=0;
    if(path.size()>0)
    {
        path.push_back(make_tuple(i, j, 0));
         di = i - get<0>(path[path.size()-2]), dj = j - get<1>(path[path.size()-2]);
        if (get_loop(table, visited, mi, mj, s, d, i + di, j + dj, sign, path))
            return true;
        path.pop_back();
    }
    if(table[i][j] != 0 || (mi == i && mj == j))
    {
        path.push_back(make_tuple(i, j, sign*-1));
        if (di != 1 && get_loop(table, visited, mi, mj, s, d, i + 1, j, sign * -1, path))
            return true;
        if (dj != 1 && get_loop(table, visited, mi, mj, s, d, i, j + 1, sign * -1, path))
            return true;
        if (di != -1 && get_loop(table, visited, mi, mj, s, d, i - 1, j, sign * -1, path))
            return true;
        if (dj!= -1 && get_loop(table, visited, mi, mj, s, d, i, j - 1, sign * -1, path))
            return true;
        path.pop_back();
    }

    visited[i][j] = 0;
    return false;
}

// return 1 for done, 0 for need more iteration, 2 for degenerate
int modi_iteration(double **table, double **current_soln, vector<double> demand, vector<double> supply, int s, int d)
{
    // cout << "\n";
    // print_2d(current_soln, d, s);
    int allocated = 0;
    for (int i = 0; i < s; i++)
        for (int j = 0; j < d; j++)
            allocated += current_soln[i][j] != 0;
    if (allocated != s + d - 1)
        return 2;
    vector<double> u = vector<double>(s, NAN), v = vector<double>(d, NAN);
    int unsolved = s + d - 1;
    u[0] = 0;
    while (unsolved)
    {
        for (int i = 0; i < s; i++)
            for (int j = 0; j < d; j++)
                if (isnan(v[j]) && !isnan(u[i]) && current_soln[i][j] != 0)
                    v[j] = table[i][j] - u[i], unsolved--;
                else if (!isnan(v[j]) && isnan(u[i]) && current_soln[i][j] != 0)
                    u[i] = table[i][j] - v[j], unsolved--;
    }
    cout << "Calculating u[i] and v[j] assuming u[1] = 0\n";
    for(int i = 0 ; i < s; i++)
        cout << "u[" << i+1 << "] = " <<u[i] << "\n";
    for(int i = 0 ; i < d; i++)
        cout << "v[" << i+1 << "] = " <<v[i] << "\n";
    double maxpenalty = -1;
    int mi, mj;
    cout << "Calculating penalties\n";
    for (int i = 0; i < s; i++)
        for (int j = 0; j < d; j++)
            if (current_soln[i][j] == 0)
            {
                double penalty = u[i] + v[j] - table[i][j];
                cout << "penalty[" << i+1 << "]["<< j+1 << "] = " <<penalty << "\n";
                if (penalty > maxpenalty)
                    maxpenalty = penalty, mi = i, mj = j;
            }
    if (maxpenalty < 0)
    {
        cout << "All negative penalty so we are done\n";
        return 1;
    }
    double **visited = create_2d_doubles(s, d);
    vector<tuple<int, int, int>> path;
    cout << "Finding closed cycle from "<<mi+1 << "," << mj+1 << "\n";
    get_loop(current_soln, visited, mi, mj, s, d, mi, mj, -1, path);
    // cout << "\n" << path.size() << "\n";
    if (path.size() < 4)
        return 2; // currently only tested for path sizes of 4
    double min_route = 1e9;

    cout << "Found path: ";
    for (auto x : path)
    {
        cout << get<0>(x)+1 << "," << get<1>(x)+1;
        if(x == path[path.size()-1]) 
            cout << "\n";
        else
            cout <<  "->";
        if (get<2>(x) == -1)
            min_route = min(min_route, current_soln[get<0>(x)][get<1>(x)]);
    }

    // cout << "min: " << min_route;
    for (auto x : path)
    {
        if (get<2>(x) == -1)
            current_soln[get<0>(x)][get<1>(x)] -= min_route;
        else if(get<2>(x) == 1)
            current_soln[get<0>(x)][get<1>(x)] += min_route;
    }
    return 0;
}

void balance_problem(int& s, int& d, vector<double>& supply, vector<double>& demand, double**& coeff) {
    double ssum=0, dsum=0;
    for(double x : supply)ssum+=x;
    for(double x : demand)dsum += x;
    if(ssum!=dsum) {
        cout << "The input problem is unbalanced. Solving by introducing dummy ";
        int os = s, od =d;
        if(ssum > dsum) {
            cout << "column.\n";
            d++;
            demand.push_back(ssum-dsum);
        }
        else {
            cout << "row.\n";
            s++;
            supply.push_back(dsum-ssum);
        }
        double** c2 = create_2d_doubles(s, d);
        for(int i =0;i<os;i++)
            for(int j =0;j<od;j++)
                c2[i][j] = coeff[i][j];
        free_2d_doubles(coeff, s, d);
        coeff = c2;
    }
}

void transporation_io()
{
    cout << "Enter number of source and demand: ";
    int s, d;
    cin >> s >> d;


    double **coeff = create_2d_doubles(s, d); // store coefficients
    cout << "Enter the cost table: \n";
    vector<double> demand, supply;
    for (int i = 0; i < s; i++)
        for (int j = 0; j < d; j++)
            cin >> coeff[i][j];
    cout << "Enter supply values: ";
    for (int i = 0; i < s; i++)
    {
        int x;
        cin >> x;
        supply.push_back(x);
    }
    cout << "Enter demand values: ";
    for (int i = 0; i < d; i++)
    {
        int x;
        cin >> x;
        demand.push_back(x);
    }
    balance_problem(s, d, supply, demand, coeff);
    cout << "Trying to solve with northwest corner then modified iteration method\n";
    auto res = northwest_corner(coeff, demand, supply, s, d);
    int i = 0;
    int soln = 2;
    print_table(s, d, coeff, res, demand, supply);

    while ((soln = modi_iteration(coeff, res, demand, supply, s, d)) == 0)
    {
        // cout << i++ << "\n";
        print_table(s, d, coeff, res, demand, supply);
    }

    if (soln == 2)
    {
        cout << "Solving failed. Retrying with LCC method and then modified iteration method\n";
        res = lcc_method(coeff, demand, supply, s, d);
        i = 0;
        soln = 0;
        print_table(s, d, coeff, res, demand, supply);

        while ((soln = modi_iteration(coeff, res, demand, supply, s, d)) == 0)
        {
            // cout << i++ << "\n";
            print_table(s, d, coeff, res, demand, supply);
        }
        if (soln == 2)
         cout << "FAILED";
    }
    if(soln!=2)
    {
        cout << "Solution found. ";
        print_table(s, d, coeff, res, demand, supply);
        cout << "Minimum cost is: " << find_cost(s, d, coeff, res) << "\n";
    }


}

int main()
{
    int inp = 1;
    while (inp)
    {
        cout << "Type 0 to quit, 1 to solve Transporation problem using MODI method\n";
        cin >> inp;
        if (inp == 1)
            transporation_io();
        else if (inp != 0)
            cout << "Try again!\n";
    }
}