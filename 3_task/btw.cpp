#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <queue>
#include <stack>
using namespace std;

#define NOTHING -1

struct csr {
    size_t *row_indices;
	int num_rowind;
    int *values;
    size_t *column_indices;
	int num_colind;
};

struct predecessors {
	int* preds;
	int number;
};

void clear(queue<int> &q )
{
   std::queue<int> empty;
   std::swap( q, empty );
}

void clear(stack<int> &s )
{
   std::stack<int> empty;
   std::swap( s, empty );
}

int* get_neighbours(int s, csr mycsr, int& ret_num){
	int* neighbours = (int*) calloc(mycsr.num_rowind, sizeof(neighbours[0]));
	int start, end, j = 0;
	
	start = mycsr.row_indices[s];
	if(s == mycsr.num_rowind - 1) 
		end = mycsr.num_colind - 1; 
	else
		end = mycsr.row_indices[s+1];
	if (end  < 0) end = 0;
	// cout << "s " << s << endl;
	// cout << start << endl;
	// cout << end << endl << endl;

	for(size_t i = start; i < end; i++){
		neighbours[j++] = mycsr.column_indices[i];
	}
	ret_num = end - start;
	if(ret_num < 0) ret_num = 0;
	
	return neighbours;
}

int main(){	
	/* INITIALIZATION */
	
    ifstream f;
	f.open ("csr_graph");
	
	struct csr my_csr;
	size_t val_n, cind_n, rprt_n;
	int *sigma, *d, *V, *neighbours;
	float *delta, *BC;
	int v, w, neig_num;
	queue<int> Q;
	stack<int> S;
	predecessors *P; // predecessors
	
	/* INITIALIZATION */
	
	/* INPUT GRAPH */
	
	f >> val_n; // number of edges
	my_csr.values = (int*) calloc(val_n, sizeof(my_csr.values[0]));
	for(size_t i = 0; i<val_n; i++){
		f >> my_csr.values[i];
	}

	f >> cind_n; // number of col indicies
	my_csr.num_colind = cind_n;
	my_csr.column_indices = (size_t*) calloc(cind_n, sizeof(my_csr.column_indices[0]));
	for(size_t i = 0; i<cind_n; i++){
		f >> my_csr.column_indices[i];
	}
	
	f >> rprt_n; // number of row ptrs
	my_csr.num_rowind = rprt_n;
	my_csr.row_indices = (size_t*) calloc(rprt_n, sizeof(my_csr.row_indices[0]));	
	for(size_t i = 0; i<rprt_n; i++){
		f >> my_csr.row_indices[i];
	}
	
	/* INPUT GRAPH */
	
	/* BTW ALGORITHM */
	
	BC = (float*) calloc(rprt_n, sizeof(BC[0]));
	for(size_t i = 0; i<rprt_n; i++){
		BC[i] = 0.0;
	}
	
	P = (predecessors*) calloc(rprt_n, sizeof(P[0]));
	for(int i = 0; i<rprt_n; i++){
		P[i].preds = (int*) calloc(rprt_n, sizeof(P[0].preds[0]));
		P[i].number = 0;
	}
	
	sigma = (int*) calloc(rprt_n, sizeof(sigma[0]));
	d = (int*) calloc(rprt_n, sizeof(d[0]));
	delta = (float*) calloc(rprt_n, sizeof(delta[0]));
	
	for(size_t s = 0; s < rprt_n; s++){
		// empty stack
		clear(S);
		// empty predecessors list
		for(size_t i = 0; i < rprt_n; i++){
			for(size_t j = 0; j < rprt_n; j++){
				P[i].preds[j] = NOTHING;
			}
			P[i].number = 0;
		}

		// empty sigma and d arrays
		for(size_t i = 0; i < rprt_n; i++){
			sigma[i] = 0;
			d[i] = -1;
		}
		sigma[s] = 1;
		d[s] = 0;
		clear(Q); // empty queue
		Q.push(s); // enqueue s -> Q
		
		while(!Q.empty()){
			v = Q.front();
			Q.pop();
			S.push(v);
			neighbours = get_neighbours(v, my_csr, neig_num);
			for(size_t i = 0; i < neig_num; i++){ 			// w = neighbours[i]
				w = neighbours[i];
				if(d[w] < 0){
					Q.push(w);
					d[w] = d[v] + 1;
					if(d[w] == d[v] + 1){
						sigma[w] += sigma[v];
						P[w].preds[P[w].number++] = v;	// append v -> P[w];
					}
				}
			}
		}
		
		for(size_t i = 0; i < rprt_n; i++){
			delta[i] = 0;
		}
		
		while(!S.empty()){
			w = S.top();
			S.pop();
			// cout << "w " << w << " P[w].number " << P[w].number << endl;
			for(size_t i = 0; i < P[w].number; i++){
				v = P[w].preds[i];
				// cout << "v " << v << " delta[v]: " << delta[v] << " sigma[v] " << sigma[v] << " sigma[w] " << sigma[w] << " sigma[v]/sigma[w]: " << sigma[v]/sigma[w] << endl;
				delta[v] += (float)sigma[v]/(float)sigma[w] * (1 + delta[w]);
			}
			if(w != s){
				BC[w] += delta[w];
			}
		}
	}
	
	/* BTW ALGORITHM */
	
	/* OUTPUT RESULTS */
	
	FILE* out;
	out = fopen ("btwcheck_cpp","w");
	
	for(size_t i = 0; i < rprt_n - 1; i++){
		BC[i] /= 2.0;
		printf("%.1f ", BC[i]);
		fprintf(out, "%.2f ", BC[i]);
		// cout << BC[i] << " ";
		// out << BC[i] << " ";
	}
	cout << endl;
	
	/* OUTPUT RESULTS */
	
	/* CLEANING */

	free(d);
	free(P);
	free(BC);
	free(sigma);
	free(delta);
	free(my_csr.values);
	free(my_csr.column_indices);
	free(my_csr.row_indices);
    f.close();
	fclose(out); 
	
	/* CLEANING */
	
	return 0;
}