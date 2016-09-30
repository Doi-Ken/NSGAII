#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<vector>
#include<list>
#include<queue>
#include<stack>
#include<fstream>
#include<iomanip>
#include<ctime>
#include<cmath>
#include"MT.h"
#include"wfg.h"
#include<direct.h>

using namespace std;

#define INT 10000000
#define Q_J 2
#define X_MIN 0
#define X_MAX 1
#define N_BON 10000
#define EPS 0.00000000001

vector< vector<int> > F_i(N_BON);


void polynomial_mutation(double *off_x, int item, double MUT, double MUT_ETA){

	double r = 0.0;
	double shita1 = 0.0;
	double shita2 = 0.0;
	double shitaq = 0.0;
	for (int i = 0; i < item; i++){
		r = nextDoubleIE();
		if (r <= MUT){
			shita1 = (off_x[i] - X_MIN) / (X_MAX - X_MIN);
			shita2 = (X_MAX - off_x[i]) / (X_MAX - X_MIN);

			r = nextDoubleIE();

			if (r <= 0.5){
				shitaq = pow(((2.0 * r) + (1.0 - 2.0 * r) * pow(1.0 - shita1, MUT_ETA + 1.0)), 1.0 / (MUT_ETA + 1.0)) - 1.0;

				off_x[i] += shitaq * (X_MAX - X_MIN);

			}
			else{
				shitaq = 1.0 - pow(2.0 * (1.0 - r) + 2.0 * (r - 0.5) * pow(1.0 - shita2, MUT_ETA + 1.0), 1.0 / (MUT_ETA + 1.0));
				off_x[i] += shitaq * (X_MAX - X_MIN);
			}
			if (off_x[i] < X_MIN){
				off_x[i] = X_MIN;
			}
			if (off_x[i] > X_MAX){
				off_x[i] = X_MAX;
			}

		}
	}
}


void SBXcrossover(double *par_x1, double *par_x2, double **off_x, int off_p_n, int item, double CO, double SBX_ETA){
	int off_p_n1 = off_p_n;
	int off_p_n2 = off_p_n + 1;

	double x1 = 0.0;
	double x2 = 0.0;
	double beta1 = 1.0;
	double beta2 = 1.0;
	double alpha1 = 2.0;
	double alpha2 = 2.0;
	double beta_q_1 = 1.0;
	double beta_q_2 = 1.0;
	double c1, c2;
	double u = nextDoubleIE();

	if (nextDoubleIE() < CO){
		for (int n = 0; n < item; n++){
			if (nextDoubleIE() < 0.5){
				if (abs(par_x2[n] - par_x1[n]) > EPS){
					x1 = 0.0;
					x2 = 0.0;
					beta1 = 1.0;
					beta2 = 1.0;
					alpha1 = 2.0;
					alpha2 = 2.0;
					beta_q_1 = 1.0;
					beta_q_2 = 1.0;

					u = nextDoubleIE();
					if (par_x1[n] > par_x2[n]){
						x2 = par_x1[n];
						x1 = par_x2[n];
					}
					else if (par_x1[n] < par_x2[n]){
						x2 = par_x2[n];
						x1 = par_x1[n];
					}


					beta1 += 2.0 * (x1 - X_MIN) / (x2 - x1);
					beta2 += 2.0 * (X_MAX - x2) / (x2 - x1);

					alpha1 -= pow(beta1, -1.0 * (SBX_ETA + 1.0));
					alpha2 -= pow(beta2, -1.0 * (SBX_ETA + 1.0));

					if (u <= 1.0 / alpha1){
						beta_q_1 = pow(u * alpha1, (1.0 / (SBX_ETA + 1.0)));
					}
					else{
						beta_q_1 = pow(1.0 / (2.0 - u * alpha1), (1.0 / (SBX_ETA + 1.0)));
					}

					if (u <= 1.0 / alpha2){
						beta_q_2 = pow(u * alpha2, (1.0 / (SBX_ETA + 1.0)));
					}
					else{
						beta_q_2 = pow(1.0 / (2.0 - u * alpha2), (1.0 / (SBX_ETA + 1.0)));
					}

					c1 = 0.5 *((x1 + x2) - beta_q_1 * (x2 - x1));
					c2 = 0.5 *((x1 + x2) + beta_q_2 * (x2 - x1));

					if (c1 < X_MIN) c1 = X_MIN;
					if (c1 > X_MAX) c1 = X_MAX;

					if (c2 < X_MIN) c2 = X_MIN;
					if (c2 > X_MAX) c2 = X_MAX;

					//imprement
					if (nextDoubleIE() < 0.5){
						off_x[off_p_n1][n] = c1;
						off_x[off_p_n2][n] = c2;

					}
					else{
						off_x[off_p_n1][n] = c2;
						off_x[off_p_n2][n] = c1;
					}
					//imprement

				}
				else{
					off_x[off_p_n1][n] = par_x1[n];
					off_x[off_p_n2][n] = par_x2[n];
				}
			}
			else{
				off_x[off_p_n1][n] = par_x1[n];
				off_x[off_p_n2][n] = par_x2[n];
			}
		}

	}
	else{
		for (int n = 0; n < item; n++){
			off_x[off_p_n1][n] = par_x1[n];
			off_x[off_p_n2][n] = par_x2[n];
		}

	}


}


bool dominate(int a, int b, double **fit, int ob, double SIGN){

	for (int o = 0; o < ob; o++){
		if (SIGN * fit[a][o] < SIGN * fit[b][o]){
			return false;
		}
	}
	for (int o = 0; o < ob; o++){
		if (SIGN * fit[a][o] > SIGN * fit[b][o]){
			return true;
		}
	}
	return false;
}

void fast_non_dominated_sort(int P[], int *pop_rank, double **fit, int n, int ob, double SIGN){
	vector<int> F;

	F_i.resize(n);

	int **S = new int *[n];
	for (int i = 0; i < n; i++){
		S[i] = new int[n];
		for (int j = 0; j < n; j++){
			S[i][j] = -1;
		}
	}
	int *n_p = new int[n];
	for (int i = 0; i < n; i++){
		n_p[i] = 0;
		pop_rank[i] = 0;
	}

	for (int p = 0; p < n; p++){
		int s = 0;
		for (int q = 0; q < n; q++){
			if (dominate(P[p], P[q], fit, ob, SIGN)){
				S[p][s++] = q;
			}
			else if (dominate(P[q], P[p], fit, ob, SIGN)){
				n_p[p]++;
			}
		}
		if (n_p[p] == 0){
			pop_rank[P[p]] = 1;
			F.push_back(p);
		}
	}
	F_i[0].clear();
	F_i[0] = F;

	int i = 1;
	while (!F.empty()){
		vector<int> Q;

		for (int p = 0; p < F.size(); p++){
			for (int q = 0; S[F[p]][q] != -1; q++){
				n_p[S[F[p]][q]]--;
				if (n_p[S[F[p]][q]] == 0){
					pop_rank[P[S[F[p]][q]]] = i + 1;
					Q.push_back(S[F[p]][q]);
				}
			}
		}
		i++;
		F.clear();
		F = Q;
		F_i[i - 1].clear();
		F_i[i - 1] = F;
		Q.clear();
	}
	F_i.resize(i - 1);
	F.clear();
	for (int i = 0; i < n; i++){
		delete[] S[i];
	}
	delete[] S;
	delete[] n_p;

}

void crowding_distance(vector<int> v, int l, double **fit, double *CD, int ob){

	vector<int> L(l);
	for (int i = 0; i < l; i++){
		L[i] = v[i];
	}
	double temp = 0.0;
	double *temp_array = new double[l];

	for (int i = 0; i < l; i++) CD[L[i]] = 0.0;
	for (int o = 0; o < ob; o++){
		
		//lab_care
		for (int i = 0; i < l; i++){
			L[i] = v[i];
		}
		//lab_care

		for (int i = 0; i < l - 1; i++){
			for (int j = 0; j < l - 1; j++){
				if (fit[L[j + 1]][o] > fit[L[j]][o]){
					swap(L[j + 1], L[j]);
				}
			}
		}


		CD[L[0]] = CD[L[l - 1]] = INT;

		if (abs(fit[L[l - 1]][o] - fit[L[0]][o]) > EPS){
		//if (fit[L[l - 1]][o] != fit[L[0]][o]){
			for (int i = 1; i < l - 1; i++){
				CD[L[i]] += (double)((double)(fit[L[i + 1]][o] - fit[L[i - 1]][o]) / (double)(fit[L[l - 1]][o] - fit[L[0]][o]));

			}
		}
		else{
			if (l > 2){
				for (int i = 1; i < l - 1; i++){
					CD[L[i]] += 0.0;
				}
			}
		}

	}
	L.clear();
	v.clear();

	delete[] temp_array;
}


//continuous optimization
void make_new_pop(double **x, double **off_x, double *CD, int *pop_rank, int pop_n, int item, double MUT, double CO, double MUT_ETA, double SBX_ETA){

	double **temp_off_x = new double*[2];

	temp_off_x[0] = new double[item];
	temp_off_x[1] = new double[item];

	for (int n = 0; n < pop_n; n++){
		int off[2];

		int bin1, bin2;

		bin1 = genrand_int31() % pop_n;
		bin2 = genrand_int31() % pop_n;
		if (pop_rank[bin1] < pop_rank[bin2]){
			off[0] = bin1;
		}
		else if (pop_rank[bin1] > pop_rank[bin2]){
			off[0] = bin2;
		}
		else{
			if (CD[bin1] > CD[bin2]){
				off[0] = bin1;
			}
			else if (CD[bin1] < CD[bin2]){
				off[0] = bin2;
			}
			else{
				if (genrand_int31() % 2){
					off[0] = bin1;
				}
				else{
					off[0] = bin2;
				}
			}
		}

		bin1 = genrand_int31() % pop_n;
		bin2 = genrand_int31() % pop_n;
		if (pop_rank[bin1] < pop_rank[bin2]){
			off[1] = bin1;
		}
		else if (pop_rank[bin1] > pop_rank[bin2]){
			off[1] = bin2;
		}
		else{
			if (CD[bin1] > CD[bin2]){
				off[1] = bin1;
			}
			else if (CD[bin1] < CD[bin2]){
				off[1] = bin2;
			}
			else{
				if (genrand_int31() % 2){
					off[1] = bin1;
				}
				else{
					off[1] = bin2;
				}
			}
		}

		SBXcrossover(x[off[0]], x[off[1]], temp_off_x, 0, item, CO, SBX_ETA);

		if (nextDoubleIE() < 0.5){
			for (int i = 0; i < item; i++){
				off_x[n][i] = temp_off_x[0][i];
			}
		}
		else{
			for (int i = 0; i < item; i++){
				off_x[n][i] = temp_off_x[1][i];
			}
		}
		polynomial_mutation(off_x[n], item, MUT, MUT_ETA);

	}
}

void make_new_pop(int **x, int **off_x, double *CD, int *pop_rank, int pop_n, int item, double MUT, double CO){

	int *temp = new int[item];

	for (int n = 0; n < pop_n; n++){

		int off[2];

		int bin1, bin2;

		bin1 = genrand_int31() % pop_n;
		bin2 = genrand_int31() % pop_n;
		if (pop_rank[bin1] < pop_rank[bin2]){
			off[0] = bin1;
		}
		else if (pop_rank[bin1] > pop_rank[bin2]){
			off[0] = bin2;
		}
		else{
			if (CD[bin1] > CD[bin2]){
				off[0] = bin1;
			}
			else if (CD[bin1] < CD[bin2]){
				off[0] = bin2;
			}
			else{
				if (genrand_int31() % 2){
					off[0] = bin1;
				}
				else{
					off[0] = bin2;
				}
			}
		}

		bin1 = genrand_int31() % pop_n;
		bin2 = genrand_int31() % pop_n;
		if (pop_rank[bin1] < pop_rank[bin2]){
			off[1] = bin1;
		}
		else if (pop_rank[bin1] > pop_rank[bin2]){
			off[1] = bin2;
		}
		else{
			if (CD[bin1] > CD[bin2]){
				off[1] = bin1;
			}
			else if (CD[bin1] < CD[bin2]){
				off[1] = bin2;
			}
			else{
				if (genrand_int31() % 2){
					off[1] = bin1;
				}
				else{
					off[1] = bin2;
				}
			}
		}

		if (nextDoubleIE() <= CO){
			for (int i = 0; i < item; i++){
				temp[i] = genrand_int31() % 2;
			}
		}
		else{
			int k = genrand_int31() % 2;
			for (int i = 0; i < item; i++){
				temp[i] = k;
			}
		}

		for (int i = 0; i < item; i++){
			off_x[n][i] = x[off[temp[i]]][i];
			if (nextDoubleIE() <= MUT){
				off_x[n][i] = (off_x[n][i] + 1) % 2;
			}
		}
	}
	delete[] temp;

}



void sort_cd(int i, double *CD){
	for (int k = 0; k < F_i[i].size() - 1; k++){
		for (int j = 0; j < F_i[i].size() - 1; j++){
			if (CD[F_i[i][j + 1]] > CD[F_i[i][j]]){
				swap(F_i[i][j + 1], F_i[i][j]);
			}
		}
	}
}

int main(int argc, char *argv[]){


	bool combination = false;
	double SIGN = -1.0;

	if (strcmp(argv[1], "-num") != 0){
		cout << "file num " << endl;
		exit(1);
	}

	if (argv[2][0] == '-'){
		cout << "experimental number" << endl;
		exit(1);
	}
	init_genrand(atoi(argv[2]));

	//knapsack
	//knapsack_file_read(); //file_read
	//sorting_profit_per_weight(); //q_j sorting
	//repair_output(); // Step1: repair_output
	//check_input_file_output(x_check); //Step2: knapsack_check
	//for (int o = 0; o < ob; o++){
	//	cout << max_weight[o] << endl;
	//}

	//input
	//ifstream fl("lambda_data_m2_h199_n200.txt");
	int EPN = 6000;

	for (int i = 1; i < argc; i++){
		if (i % 2){
			if (argv[i][0] != '-'){
				cout << "format is disavailable" << endl;
				cout << "-... value -... value -... value ........" << endl;
				exit(1);
			}
		}
	}

	char problem[50] = "wfg1";

	int read_ind = 3;
	if (strcmp(argv[read_ind], "-pro") != 0){
		cout << "default wfg1" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "knapsack") == 0 ||
			strcmp(argv[read_ind], "wfg1") == 0 || strcmp(argv[read_ind], "wfg2") == 0 ||
			strcmp(argv[read_ind], "wfg3") == 0 || strcmp(argv[read_ind], "wfg4") == 0 ||
			strcmp(argv[read_ind], "wfg5") == 0 || strcmp(argv[read_ind], "wfg6") == 0 ||
			strcmp(argv[read_ind], "wfg7") == 0 || strcmp(argv[read_ind], "wfg8") == 0 ||
			strcmp(argv[read_ind], "wfg9") == 0 ||
			strcmp(argv[read_ind], "dtlz1") == 0 || strcmp(argv[read_ind], "dtlz2") == 0 ||
			strcmp(argv[read_ind], "dtlz3") == 0 || strcmp(argv[read_ind], "dtlz4") == 0 ||
			strcmp(argv[read_ind], "dtlz5") == 0 || strcmp(argv[read_ind], "dtlz6") == 0 ||
			strcmp(argv[read_ind], "dtlz7") == 0 || strcmp(argv[read_ind], "dtlz8") == 0 ||
			strcmp(argv[read_ind], "dtlz9") == 0){
			strcpy(problem, argv[read_ind]);
		}
		else{
			cout << "there is no problem." << endl;
			exit(1);
		}

	}


	int ob = 10;

	read_ind++;
	if (strcmp(argv[read_ind], "-obj") != 0){
		cout << "default 10obj" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		ob = atoi(argv[read_ind]);
	}

	int POPSIZE = 220;
	
	read_ind++;
	if (strcmp(argv[read_ind], "-pop") != 0){
		cout << "default pop size 220" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		POPSIZE = atoi(argv[read_ind]);
	}


	int GEN = 10000;
	int VALNUM = 40000;

	read_ind++;
	if (strcmp(argv[read_ind], "-valnum") != 0){
		cout << "default valnum 40000" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		VALNUM = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-gen") != 0){
		cout << "default generation 10000" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		GEN = atoi(argv[read_ind]);
	}


	double SBX_ETA = 15.0;
	double MUT_ETA = 20.0;
	double CO = 0.8;
	double mut = 1.0;

	read_ind++;
	if (strcmp(argv[read_ind], "-sbxeta") != 0){
		cout << "default sbxeta 15.0" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		SBX_ETA = atof(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-muteta") != 0){
		cout << "default muteta 20.0" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		MUT_ETA = atof(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-co") != 0){
		cout << "default crossover probability" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		CO = atof(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-mut") != 0){
		cout << "default mutation probability mut(1.0) / item" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		mut = atof(argv[read_ind]);
	}

	int k_dtlz = 10;
	int k_fac = 2;
	int l_fac = 10;

	read_ind++;
	if (strcmp(argv[read_ind], "-k_fac") != 0){
		cout << "default k_factor 2" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		k_fac = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-l_fac") != 0){
		cout << "default l_factor 10" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		l_fac = atoi(argv[read_ind]);
	}

	read_ind++;
	if (strcmp(argv[read_ind], "-k_dtlz") != 0){
		cout << "default k_dtlz 10" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		k_dtlz = atoi(argv[read_ind]);
	}

	char max_min[10] = "min";
	read_ind++;
	if (strcmp(argv[read_ind], "-maxmin") != 0){
		cout << "default improve_lambda false" << endl;
		read_ind--;
	}
	else{
		read_ind++;
		if (strcmp(argv[read_ind], "max") == 0){
			strcpy(max_min, "max");
		}
	}
	if (strcmp(max_min, "max") == 0){
		SIGN = 1.0;
	}
	else if (strcmp(max_min, "min") == 0){
		SIGN = -1.0;
	}
	else{
		cout << "max_min" << endl;
		exit(1);
		
	}
	int k_val = k_fac;
	int l_val = l_fac;
	int item = k_val + l_val;


	if (strcmp("dtlz1", problem) == 0 || strcmp("dtlz2", problem) == 0 ||
		strcmp("dtlz3", problem) == 0 || strcmp("dtlz4", problem) == 0 ||
		strcmp("dtlz4", problem) == 0 || strcmp("dtlz4", problem) == 0 ||
		strcmp("dtlz5", problem) == 0 || strcmp("dtlz6", problem) == 0 || 
		strcmp("dtlz7", problem) == 0 || strcmp("dtlz8", problem) == 0 ||
		strcmp("dtlz9", problem) == 0){
		item = ob + k_dtlz - 1;
	}
	else if (strcmp("knapsack", problem) == 0){
		combination = true;

		char file_name[50];
		strcpy(file_name, problem);
		strcat(file_name, "_2_500");
		if (ob != 2){
			char ch[5];
			sprintf(ch, "%d", ob);
			strcat(file_name, "to");
			strcat(file_name, ch);
		}
		strcat(file_name, ".txt");
		knapsack_file_read(file_name, item, ob); //file_read
		sorting_profit_per_weight(item, ob); //q_j sorting
		repair_output(item, ob);
		int *cheeck = new int[item];
		for (int i = 0; i < item; i++){
			cheeck[i] = 1;
		}
		check_input_file_output(cheeck, item, ob);

	}

	cout << "number: " << atoi(argv[2]) << endl;
	cout << "problem: " << problem << ", " << ob << "ob" << endl;
	cout << "popsize: " << POPSIZE << endl;
	cout << "valnum: " << VALNUM << ", generation: " << GEN << endl << endl;
	cout << "sbx_eta: " << SBX_ETA << ", mut_eta: " << MUT_ETA << endl;
	cout << "crossover probabilty: " << CO << ", mutation probability: " << mut << " / " << item << endl;
	cout << "k_fac: " << k_fac << ", l_fac: " << l_fac << endl;
	cout << "k_dtlz: " << k_dtlz << endl;
	
	double MUT = mut / item;



	int N = POPSIZE;
	
	double *CD = new double[N];
	double **x = new double *[N];
	int **x_com = new int *[N];
	double **off_x = new double *[N];
	int ** off_x_com = new int *[N];
	int i = 0;
	for (int n = 0; n < N; n++){
		x[n] = new double[item];
		x_com[n] = new int[item];
		off_x[n] = new double[item];
		off_x_com[n] = new int[item];
	}
	if (combination == false){
		for (int n = 0; n < N; n++){
			for (i = 0; i < item; i++){
				x[n][i] = nextDoubleIE();
			}
		}
	}
	else{
		for (int n = 0; n < N; n++){
			for (i = 0; i < item; i++){
				x_com[n][i] = genrand_int31() % 2;
			}
		}
	}
	

	double **fit = new double *[N];
	int *pop_num = new int[N];
	if (combination == false){
		for (int n = 0; n < N; n++){
			fit[n] = new double[ob];
			fitness(x[n], fit[n], problem, ob, k_val, l_val, k_dtlz);
			pop_num[n] = n;
		}
	}
	else{
		for (int n = 0; n < N; n++){
			fit[n] = new double[ob];
			fitness_knap(x_com[n], fit[n], item, ob);
			pop_num[n] = n;
		}
	}

	int *pop_rank = new int[N];

	fast_non_dominated_sort(pop_num, pop_rank, fit, N, ob, SIGN);

	i = 0;
	cout << F_i.size() << endl;
	while (i < F_i.size()){
		crowding_distance(F_i[i], F_i[i].size(), fit, CD, ob);
		i++;
	}

	if (combination == false){
		make_new_pop(x, off_x, CD, pop_rank, N, item, MUT, CO, MUT_ETA, SBX_ETA);
		for (int n = 0; n < N; n++){
			fitness(off_x[n], fit[n], problem, ob, k_val, l_val, k_dtlz);
			pop_num[n] = n;
		}
	}
	else{
		make_new_pop(x_com, off_x_com, CD, pop_rank, N, item, MUT, CO);
		for (int n = 0; n < N; n++){
			fitness_knap(off_x_com[n], fit[n],item,ob);
			pop_num[n] = n;
		}
	}



	fast_non_dominated_sort(pop_num, pop_rank, fit, N, ob, SIGN);

	i = 0;

	while (i < F_i.size()){

		crowding_distance(F_i[i], F_i[i].size(), fit, CD, ob);
		i++;
	}


	char detaildata[100] = "./";
	char datafile[50];
	char outputfile[30] = "graph";
	char bbb[5] = "00";
	sprintf(bbb, "%d", atoi(argv[2]));
	strcat(outputfile, bbb);
	strcpy(datafile, outputfile);
	strcat(outputfile, ".txt");
	//progress_output
	//_mkdir(datafile);

	int val = 0;
	int child_n = N;
	int finflag = 0;
	for (int g = 0; g < GEN; g++){

		//gen_progress
		/*strcpy(detaildata, "./");
		strcat(detaildata, datafile);
		strcat(detaildata, "/");
		strcat(detaildata, "gen");
		char gennumch[8] = "00";
		sprintf(gennumch, "%d", g);
		strcat(detaildata, gennumch);
		strcat(detaildata, ".txt");

		ofstream foutdetail3(detaildata);


		for (int i = 0; i < N; i++){
			for (int o = 0; o < ob; o++){
				foutdetail3 << fit[i][o] << "\t";
			}
			foutdetail3 << endl;
		}

		foutdetail3.close();*/

		double **p_o_fit = new double *[N + child_n];
		double **p_o_x = new double *[N + child_n];
		int **p_o_x_com = new int *[N + child_n];
		int *p_o_pop_num = new int[N + child_n];
		int *p_o_pop_rank = new int[N + child_n];
		double *p_o_CD = new double[N + child_n];
		for (int n = 0; n < N + child_n; n++){
			p_o_x[n] = new double[item];
			p_o_x_com[n] = new int[item];
			p_o_fit[n] = new double[ob];
		}

		int n = 0;
		if (combination == false){
			for (n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					p_o_x[n][i] = x[n][i];
				}
				p_o_pop_num[n] = n;
			}

			for (; n < N + child_n; n++){
				for (int i = 0; i < item; i++){
					p_o_x[n][i] = off_x[n - N][i];
				}
				p_o_pop_num[n] = n;
			}
			for (n = 0; n < N + child_n; n++){
				fitness(p_o_x[n], p_o_fit[n], problem, ob, k_val, l_val, k_dtlz);
			}
		}
		else{
			for (n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					p_o_x_com[n][i] = x_com[n][i];
				}
				p_o_pop_num[n] = n;
			}
			for (; n < N + child_n; n++){
				for (int i = 0; i < item; i++){
					p_o_x_com[n][i] = off_x_com[n - N][i];
				}
				p_o_pop_num[n] = n;
			}
			for (n = 0; n < N + child_n; n++){
				fitness_knap(p_o_x_com[n], p_o_fit[n], item, ob);
			}
		}

		fast_non_dominated_sort(p_o_pop_num, p_o_pop_rank, p_o_fit, N + child_n, ob, SIGN);


		int p = 0;
		int i = 0;
		for (int n = 0; n < N + child_n; n++){
			p_o_CD[n] = 0.0;
		}

		while (p + F_i[i].size() <= N){
			crowding_distance(F_i[i], F_i[i].size(), p_o_fit, p_o_CD, ob);
			for (int j = 0; j < F_i[i].size(); j++){
				if (combination == false){
					for (int it = 0; it < item; it++){
						x[p][it] = p_o_x[F_i[i][j]][it];
					}
				}
				else{
					for (int it = 0; it < item; it++){
						x_com[p][it] = p_o_x_com[F_i[i][j]][it];
					}
				}
				pop_num[p] = p;
				for (int o = 0; o < ob; o++){
					fit[p][o] = p_o_fit[F_i[i][j]][o];
				}
				CD[p] = p_o_CD[F_i[i][j]];
				p++;
			}

			i = i + 1;
		}
		crowding_distance(F_i[i], F_i[i].size(), p_o_fit, p_o_CD, ob);
		sort_cd(i, p_o_CD);
		int j = 0;
		for (; p < N; p++){
			if (combination == false){
				for (int it = 0; it < item; it++){
					x[p][it] = p_o_x[F_i[i][j]][it];
				}
			}
			else{
				for (int it = 0; it < item; it++){
					x_com[p][it] = p_o_x_com[F_i[i][j]][it];
				}
			}
			pop_num[p] = p;
			pop_rank[p] = p_o_pop_rank[F_i[i][j]];
			for (int o = 0; o < ob; o++){
				fit[p][o] = p_o_fit[F_i[i][j]][o];
			}
			CD[p] = p_o_CD[F_i[i][j]];
			j++;
		}

		if (finflag == 0){
			fast_non_dominated_sort(pop_num, pop_rank, fit, child_n, ob, SIGN);
			for (int i = 0; i < F_i.size(); i++){
				crowding_distance(F_i[i], F_i[i].size(), fit, CD, ob);
			}
			if (combination == false){
				make_new_pop(x, off_x, CD, pop_rank, child_n, item, MUT, CO, MUT_ETA, SBX_ETA);
			}
			else{
				make_new_pop(x_com, off_x_com, CD, pop_rank, child_n, item, MUT, CO);
			}
			if (!(g % 200)){
				cout << g << " " << val << endl;
			}
		}
		else{
			cout << g << " " << val << " " << child_n << endl;
			break;
		}

		if (val < VALNUM){
			if (val + N < VALNUM){
				child_n = N;
				val += N;
			}
			else{
				child_n = VALNUM - val;
				val += child_n;
			}
		}
		else{
			finflag = 1;
		}

		for (int n = 0; n < N + child_n; n++){
			delete[] p_o_fit[n];
			delete[] p_o_x[n];
			delete[] p_o_x_com[n];
		}
		delete[] p_o_fit;
		delete[] p_o_x;
		delete[] p_o_pop_num;
		delete[] p_o_pop_rank;
		delete[] p_o_CD;
		delete[] p_o_x_com;

	}


	char filename[30] = "graph";
	char ccc[5] = "00";
	sprintf(ccc, "%d", atoi(argv[2]));
	strcat(filename, ccc);
	strcat(filename, ".txt");


	ofstream fout(filename);
	for (int n = 0; n < N; n++){
		//cout << pop_rank[n] << " ";
		//if (pop_rank[n] == 1){
		for (int o = 0; o < ob - 1; o++){
			//cout << fit[n][o] << " ";
			fout << fit[n][o] << "\t";
		}
		//cout << fit[n][ob - 1];
		fout << fit[n][ob - 1];

		//	}
		//cout << CD[n] << endl;
		fout << endl;
	}
	fout.close();

	//gen_progress
	/*strcpy(detaildata, "./");
	strcat(detaildata, datafile);
	strcat(detaildata, "/");
	strcat(detaildata, "gen");
	char gennumch[8] = "00";
	sprintf(gennumch, "%d", GEN);
	strcat(detaildata, gennumch);
	strcat(detaildata, ".txt");

	ofstream foutdetail3(detaildata);


	for (int i = 0; i < N; i++){
		for (int o = 0; o < ob; o++){
			foutdetail3 << fit[i][o] << "\t";
		}
		foutdetail3 << endl;
	}

	foutdetail3.close();*/


	//cout << val << endl;

	delete[] pop_num;
	delete[] pop_rank;
	delete[] CD;
	for (int n = 0; n < N; n++){
		delete[] x[n];
		delete[] off_x[n];
		delete[] x_com[n];
		delete[] off_x_com[n];
	}
	delete[] x;
	delete[] off_x;
	delete[] x_com;
	delete[] off_x_com;
	for (int n = 0; n < N; n++){
		delete[] fit[n];
	}
	delete[] fit;


	return 0;
}