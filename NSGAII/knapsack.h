#ifndef KSDEF

#define KSDEF
#define Q_J 2

extern double *max_weight; //max_weight
extern int *item_num; // item numbering

extern double **original_weight; //weight[o][i] o: object, i: item
extern double **original_profit; //profit[o][i] o: object, i: item
extern double **original_profit_per_weight; //profit[o][i] / weight[o][i]
extern double *original_q_j; //q_j = max(profit[o][i] / weight[o][i])
extern int *original_x_check; //desicion number
extern int *x_check; //desicion number

extern double **weight; //weight[o][i] o: object, i: item
extern double **profit; //profit[o][i] o: object, i: item
extern double **profit_per_weight; //profit[o][i] / weight[o][i]
extern double *q_j; //q_j = max(profit[o][i] / weight[o][i])

void knapsack_file_read(char *file, int &item, int &ob);
void sorting_profit_per_weight(int item, int ob);
void fitness_knap(int y[], double fit[],int item, int ob);
void repair_output(int item, int ob);
void cheeck_input_file_output(int y[], int item, int ob);

#endif