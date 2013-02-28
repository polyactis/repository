// SLiM: a forward population genetic simulation with selection and linkage
//
// version 1.2 (September 5th, 2012)
//
// Copyright (C) 2012  Philipp Messer
//
// compile by:
//
// g++ -fast ./slim.cpp -lgsl -lgslcblas -o slim // iMac
// g++ -O3   ./slim.cpp -lgsl -lgslcblas -o slim // BioX2
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version (http://www.gnu.org/licenses/).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>


using namespace std;

const gsl_rng *rng; 


class event
{
  // type of events:
  //
  // t P i n [j]:  add subpopulation i of size n [drawn from j]
  // t N i n:      set size of subpopulation i to n
  // t M i j x:    set fraction x of subpopulation i that originates as migrants from j
  //
  // t S i n:      output sample of n randomly drawn genomes from subpopulation i
  // t F:          output list of all mutations that have become fixed so far
  // t A [file]:   output state of entire population [into file]
  // t T m:        follow trajectory of mutation m (specified by mutation type) from generation t on

public:
  
  char   t;         // event type
  vector<string> s; // vector of strings with parameters of event
  int np;           // number of parameters

  event(char T, vector<string> S)
  {
    t = T;
    s = S;
    np = s.size();

    string options = "PNMSFAT";
    if(options.find(t)==string::npos) 
      { 
	cerr << "ERROR (initialize): invalid event type \"" << t;
	for(int i=0; i<np; i++) { cerr << " " << s[i]; }
	cerr << "\"" << endl;
	exit(1); 
      }
  }  
};



class mutation
{
public:

  int   t; // mutation type identifier
  int   x; // position
  float s; // selection coefficient
  float h; // dominance coefficient
  int   n; // prevalence in population

  mutation(void) { n = 0; }

  mutation(int T, int X, float S, float H) 
  { 
    t = T;
    x = X;
    s = S;
    h = H;
    n = 0;
  }

  mutation(string line)
  {
    string sub;
    istringstream iss(line); 
    iss >> sub;  
    iss >> sub; sub.erase(0,1); t = atoi(sub.c_str());
    iss >> sub; x = atoi(sub.c_str())-1;
    iss >> sub; s = atof(sub.c_str());
    iss >> sub; h = atof(sub.c_str());
    n=0;
  }

  void print() { cout << "m" << t << " " << x+1 << " " << s << " " << h << " "<< n; }

  void print_non() { cout << "m" << t << " " << x+1 << " " << s << " " << h; }

  void print(ofstream& outfile) { outfile << "m" << t << " " << x+1 << " " << s << " " << h << " "<< n; }
};


class introduced_mutation : public mutation
{
public:

  int i;   // subpopulation into which mutation is to be introduced
  int nAA; // number of homozygotes
  int nAa; // number of heterozygotes

  introduced_mutation(int T, int X, float S, float H, int I, int NAA, int NAa)
  {
    t = T;
    x = X;
    s = S;
    h = H;
    i = I;
    nAA = NAA;
    nAa = NAa;
    n = 0;
  }
};


class mutation_type
{
  // a mutation type is specified by the DFE and the dominance coefficient
  //
  // DFE options: f: fixed (s) 
  //              e: exponential (mean s)
  //              g: gamma distribution (mean s,shape)
  //
  // examples: synonymous, nonsynonymous, adaptive, etc.

public:
    
  double h;         // dominance coefficient of mutations of this type
  char   d;         // DFE (f: fixed, g: gamma, e: exponential)
  vector<double> p; // DFE parameters

  mutation_type(double H, char D, vector<double> P)
  {
    h = H;
    d = D;
    p = P;

    string s = "fge";

    if(s.find(d)==string::npos)  { cerr << "ERROR (initialize): invalid mutation type parameters" << endl; exit(1); }
    if(p.size()==0)              { cerr << "ERROR (initialize): invalid mutation type parameters" << endl; exit(1); }
  }

  float draw_s()
  {
    switch(d)
      {
      case 'f': return p[0]; 
      case 'g': return gsl_ran_gamma(rng,p[1],p[0]/p[1]);
      case 'e': return gsl_ran_exponential(rng,p[0]);
      default: exit(1);
      }
  }
};




class genomic_element
{
  // a genomic element has a genomic element type identifier (i), start (s) and end (e) position

public:

  int i,s,e;

  genomic_element(int I, int S, int E) { i = I; s = S; e = E; }
};


class genomic_element_type
{
  // a genomic element type is specified by a vector of the mutation type identifiers off all 
  // mutation types than can occur in such elements and a vector of their relative fractions.
  // examples: exon, intron, utr, intergenic, etc.

private:

  gsl_ran_discrete_t* LT;

public:

  vector<int>    m; // mutation types identifiers in this element
  vector<double> g; // relative fractions of each mutation type

  genomic_element_type(vector<int> M, vector<double> G)
  {
    m = M;
    g = G;  

    if(m.size() != g.size()) { exit(1); }
    double A[m.size()]; for(int i=0; i<m.size(); i++) { A[i] = g[i]; }
    LT = gsl_ran_discrete_preproc(G.size(),A);
  }

  int draw_mutation_type() { return m[gsl_ran_discrete(rng,LT)]; }
};



class partial_sweep
{
 public:

  int t;
  int x;
  float p;

  partial_sweep(int T, int X, float P)
  {
    t = T; x = X; p = P;
  }
};


class chromosome : public vector<genomic_element>
{
  // the chromosome is a vector of genomic elements (type, start, end)

private:

  gsl_ran_discrete_t* LT_M; // mutation
  gsl_ran_discrete_t* LT_R; // recombination

public:

  map<int,mutation_type>        mutation_types;
  map<int,genomic_element_type> genomic_element_types;
  vector<int>                   rec_x;
  vector<double>                rec_r;

  double M; // overall mutation rate
  double R; // overall recombination rate

  void initialize_rng()
  {
    if(size() == 0)       { cerr << "ERROR (initialize): empty chromosome" << endl; exit(1); }
    if(rec_r.size() == 0) { cerr << "ERROR (initialize): recombination rate not specified" << endl; exit(1); }
    if(!(M>=0))           { cerr << "ERROR (initialize): invalid mutation rate" << endl; exit(1); }

    for(int i=0; i<size(); i++)
      {
	if(genomic_element_types.count(operator[](i).i)==0) 
	  { 
	    cerr << "ERROR (initialize): genomic element type " << operator[](i).i << " not defined" << endl; exit(1); 
	  }
      }

    for (map<int,genomic_element_type>::iterator it = genomic_element_types.begin(); it!=genomic_element_types.end(); it++)
      {
	for(int j=0; j<it->second.m.size(); j++)
	  {
	    if(mutation_types.count(it->second.m[j]) == 0)
	      {
	        cerr << "ERROR (initialize): mutation type " << it->second.m[j] << " not defined" << endl; exit(1); 
	      }
	  }
      }  

    double A[size()]; int L = 0;
    for(int i=0; i<size(); i++) 
      { 
	int l = operator[](i).e - operator[](i).s + 1.0; 
	A[i] = (double)l; L += l;
      }
    LT_M = gsl_ran_discrete_preproc(size(),A); M = M*(double)L;

    double B[rec_r.size()];
    B[0] = rec_r[0]*(double)rec_x[0]; R += B[0];
    for(int i=1; i<rec_r.size(); i++) { B[i] = rec_r[i]*(double)(rec_x[i]-rec_x[i-1]); R+= B[i];}
    LT_R = gsl_ran_discrete_preproc(rec_r.size(),B);
  }


  int draw_n_mut() { return gsl_ran_poisson(rng,M); }
 
  int draw_n_rec() { return gsl_ran_poisson(rng,R); }
    

  mutation draw_new_mut()
  {
    int g = gsl_ran_discrete(rng,LT_M); // genomic element
    genomic_element_type ge_type = genomic_element_types.find(operator[](g).i)->second; // genomic element type

    int mut_type_id = ge_type.draw_mutation_type(); // mutation type id
    mutation_type mut_type = mutation_types.find(mut_type_id)->second; // mutation type

    int   x = operator[](g).s + gsl_rng_uniform_int(rng,operator[](g).e - operator[](g).s + 1); // position    
    float s = mut_type.draw_s(); // selection coefficient

    return mutation(mut_type_id,x,s,mut_type.h);
  }


  vector<int> draw_rec_breakpoints(int n)
  {
    vector<int> r;

    for(int i=0; i<n; i++)
      {
	int x = 0;
	int j = gsl_ran_discrete(rng,LT_R);

	if(j==0) { x += gsl_rng_uniform_int(rng,rec_x[j]); }
	else     { x  = rec_x[j-1] + gsl_rng_uniform_int(rng,rec_x[j]-rec_x[j-1]); }

	r.push_back(x);
      }
    return r;
  }
};



class genome : public multimap<int,int> 
{ 
  // keys:   mutation positions
  // values: mutation ids
};




class subpopulation
{
  // a subpopulation is described by the vector G of 2N genomes
  // individual i is constituted by the two genomes 2*i and 2*i+1

private:
  
  gsl_ran_discrete_t* LT;   

public:

  int N;

  vector<genome> G_parent;
  vector<genome> G_child;

  map<int,double> m; // m[i]: fraction made up of migrants from subpopulation i per generation
 

  subpopulation(int n)
  {
    N = n;
    G_parent.resize(2*N); G_child.resize(2*N);
    double A[N]; for(int i=0; i<N; i++) { A[i] = 1.0; }
    LT = gsl_ran_discrete_preproc(N,A);
  }


  int draw_individual() { return gsl_ran_discrete(rng,LT); }


  void update_fitness(vector<mutation>& M)
  {
    // calculate fitnesses in parent population and create new lookup table
    
    gsl_ran_discrete_free(LT);
    double A[(int)(G_parent.size()/2)]; for(int i=0; i<(int)(G_parent.size()/2); i++) { A[i] = W(2*i,2*i+1,M); }
    LT = gsl_ran_discrete_preproc((int)(G_parent.size()/2),A);
  }

  double W(int i, int j, vector<mutation>& M)
  {
    // calculate the fitness of the individual constituted by genomes i and j in the parent population
    //
    // algorithm:
    //
    // make a local copy of genome j
    // iterate through all mut in genome i
    // for each mut in i, check whether mut is also present in j
    // if not present: add fitness 2hs, then move to next mut in i
    // if present: add fitness 2s, remove mut from genome j, move to next mut in i
    // after traversing all mut in i, iterate through all remaining mut in j and add 2hs 

    genome g = G_parent[j]; // make local copy of genome j
    double w = 1.0;

    multimap<int,int>::iterator it1;
    multimap<int,int>::iterator it2;

    // iterate through genome i

    for(it1 = G_parent[i].begin(); it1 != G_parent[i].end(); it1++ )
      {
	if(g.count(it1->first) == 0) // heterozygous
	  { 
	    w = w * (1.0 + 2.0 * M[it1->second].h * M[it1->second].s); 
	  }
	else                         // homozygous
	  {
	    // iterate through all mutations in genome j with same key

	    pair<multimap<int,int>::iterator,multimap<int,int>::iterator> range;
	    range = g.equal_range(it1->first); it2 = range.first;

	    while(it2 != range.second)
	      {
		if(it1->second == it2->second) 
		  { 
		    w = w * (1.0 + 2.0 *M[it1->second].s); 
		    g.erase(it2);
		    it2 = range.second;
		  }
		else{ it2++; }
	      }
	  }
      }

    // iterate through remaining mutations in genome j (all heterozygous)

    for(it2 = g.begin(); it2 != g.end(); it2++) { w = w * (1.0 + 2.0*M[it2->second].h*M[it2->second].s); }
    
    if(w<0) { w = 0.0; }
    return w;
  }
};





class population : public map<int,subpopulation>
{
  // the population is a map of subpopulations

private:

  map<int,int> M_map; // map between mutation ids in parent (key) and child (value) population


public: 

  vector<mutation> M_parent;
  vector<mutation> M_child;
  vector<mutation> Substitutions;

  map<int,subpopulation>::iterator it;

  vector<string> parameters;

  void add_subpopulation(int i, unsigned int N) 
  { 
    // add new empty subpopulation i of size N

    if(count(i)!=0) { cerr << "ERROR (add subpopulation): subpopulation p"<< i << " already exists" << endl; exit(1); }
    if(N<1)         { cerr << "ERROR (add subpopulation): subpopulation p"<< i << " empty" << endl; exit(1); }

    insert(pair<int,subpopulation>(i,subpopulation(N))); 
  }


  void add_subpopulation(int i, int j, unsigned int N) 
  { 
    // add new subpopulation i of size N individuals drawn from source subpopulation j

    if(count(i)!=0) { cerr << "ERROR (add subpopulation): subpopulation p"<< i << " already exists" << endl; exit(1); }
    if(count(j)==0) { cerr << "ERROR (add subpopulation): source subpopulation p"<< j << " does not exists" << endl; exit(1); }
    if(N<1)         { cerr << "ERROR (add subpopulation): subpopulation p"<< i << " empty" << endl; exit(1); }

    insert(pair<int,subpopulation>(i,subpopulation(N))); 

    for(int p=0; p<find(i)->second.N; p++)
      {
	// draw individual from subpopulation j and assign to be a parent in i  

	int m = find(j)->second.draw_individual();

	for(multimap<int,int>::iterator it = find(j)->second.G_parent[2*m].begin(); it != find(j)->second.G_parent[2*m].end(); it++)
	  {
	    find(i)->second.G_parent[2*p].insert(pair<int,int>(it->first,it->second));
	    M_parent[it->second].n++;
	  }
	for(multimap<int,int>::iterator it = find(j)->second.G_parent[2*m+1].begin(); it != find(j)->second.G_parent[2*m+1].end(); it++)
	  {
	    find(i)->second.G_parent[2*p+1].insert(pair<int,int>(it->first,it->second));
	    M_parent[it->second].n++;
	  }
      }
  }


  void change_size(int i, unsigned int N) 
  {
    if(count(i)==0) { cerr << "ERROR (change size): no subpopulation p"<< i << endl; exit(1); }

    if(N==0) // remove subpopulation i 
      {
	erase(i); 
	for(it = begin(); it != end(); it++) { it->second.m.erase(i); } 
      }
    else { find(i)->second.N = N; find(i)->second.G_child.resize(2*N); }
  }


  void set_migration(int i, int j, double m) 
  { 
    // set fraction m of i that originates as migrants from j per generation  

    if(count(i)==0) { cerr << "ERROR (set migration): no subpopulation p"<< i << endl; exit(1); }
    if(count(j)==0) { cerr << "ERROR (set migration): no subpopulation p"<< j << endl; exit(1); }

    if(find(i)->second.m.count(j) !=0 ) { find(i)->second.m.erase(j); }

    find(i)->second.m.insert(pair<int,double>(j,m)); 
  }


  void execute_event(event& E, int g, vector<int>& FM)
  {
    char type = E.t;

    if(type == 'P') // add subpopulation
      { 
	if(E.np==2) // empty subpopulation
	  { 
	    string sub = E.s[0]; sub.erase(0,1);

	    int i = atoi(sub.c_str());
	    int n = (int)atof(E.s[1].c_str());
	    add_subpopulation(i,n);
	  }
	      
	if(E.np==3) // drawn from source population
	  {
	    string sub1 = E.s[0]; sub1.erase(0,1);
	    string sub2 = E.s[2]; sub2.erase(0,1);

	    int i = atoi(sub1.c_str());
	    int j = atoi(sub2.c_str());
	    int n = (int)atof(E.s[1].c_str());
	    add_subpopulation(i,j,n);
	  } 
      }
	  
    if(type == 'N') // change subpopulation size
      { 
	string sub = E.s[0]; sub.erase(0,1);

	int i = atoi(sub.c_str());
	int n = (int)atof(E.s[1].c_str());

	change_size(i,n);
      }
	  
    if(type == 'M') // change migration rate
      {
	string sub1 = E.s[0]; sub1.erase(0,1);
	string sub2 = E.s[1]; sub2.erase(0,1);

	int    i = atoi(sub1.c_str());
	int    j = atoi(sub2.c_str());
	double m = atof(E.s[2].c_str());
	set_migration(i,j,m); 
      }

    if(type == 'A') // output state of entire population
      {
	if(E.s.size()==0) 
	  {
	    cout << "#OUT: " << g << " A" << endl;
	    cout << "Populations:" << endl;
	    for(it = begin(); it != end(); it++) {  cout << "p" << it->first << " " << it->second.N << endl; }
	    print_all(); 
	  }	
	if(E.s.size()==1)
	  {
	    ofstream outfile;
	    outfile.open (E.s[0].c_str());

	    for(int i=0; i<parameters.size(); i++) { outfile << parameters[i] << endl; }

	    if(outfile.is_open()) 
	      { 
		outfile << "#OUT: " << g << " A " << E.s[0].c_str() << endl;
		outfile << "Populations: " << endl;
		for(it = begin(); it != end(); it++) {  outfile << "p" << it->first << " " << it->second.N << endl; }
		print_all(outfile);
		outfile.close(); 
	      }
	    else { cerr << "ERROR (output): could not open "<< E.s[0].c_str() << endl; exit(1); }
	  }
      }

    if(type == 'S') // output subpopulation sample
      {
	string sub = E.s[0]; sub.erase(0,1);

	int    i = atoi(sub.c_str());
	int    n = atoi(E.s[1].c_str());   
	cout << "#OUT: " << g << " S p" << i << " " << n << endl;
	print_sample(i,n);
      }

    if(type == 'F') // output list of fixed mutations
      {
	cout << "#OUT: " << g << " F " << endl;
	cout << "Mutations:" << endl;
	for(int i=0; i<Substitutions.size(); i++) { cout << i+1 << " "; Substitutions[i].print(); cout << endl; }
      }

    if(type == 'T') // track mutation-types
      {
	string sub = E.s[0]; sub.erase(0,1);
	FM.push_back(atoi(sub.c_str()));
      }
  }


  void introduce_mutation(introduced_mutation M) 
  {
    // introduce user-defined mutation
 
    if(count(M.i)==0) { cerr << "ERROR (predetermined mutation): subpopulation "<< M.i << " does not exists" << endl; exit(1); }
    if(find(M.i)->second.G_child.size()/2 < M.nAA + M.nAa) 
      { 
	cerr << "ERROR (predetermined mutation): not enough individuals in subpopulation "<< M.i << endl; exit(1); 
      }

    M_child.push_back(mutation(M.t,M.x,M.s,M.h));
    int id = M_child.size()-1;

    // introduce homozygotes

    for(int j=0; j<M.nAA; j++) 
      { 
	find(M.i)->second.G_child[2*j].insert(pair<int,int>(M.x,id));
	find(M.i)->second.G_child[2*j+1].insert(pair<int,int>(M.x,id));
	M_child[id].n += 2;
      }

    // introduce heterozygotes

    for(int j=M.nAA; j<M.nAA+M.nAa; j++) 
      { 
	find(M.i)->second.G_child[2*j].insert(pair<int,int>(M.x,id));
	M_child[id].n += 1;
      }
  }

  
  void track_mutations(int g, vector<int>& TM, vector<partial_sweep>& PS)
  {
    // output trajectories of followed mutations and set s=0 for partial sweeps 

    int N = 0; for(it = begin(); it != end(); it++) { N += it->second.N; }

    for(int i=0; i<M_child.size(); i++)
      {
	// follow mutations

	for(int j=0; j<TM.size(); j++)
	  {
	    if(M_child[i].t == TM[j]) 
	      { 
		cout << "#OUT: " << g << " T "; M_child[i].print(); cout << endl;
	      }
	  }
	
	// set s=0 once partial sweeps have reached their target frequency

	for(int j=0; j<PS.size(); j++)
	  {
	    if(M_child[i].x == PS[j].x && M_child[i].t == PS[j].t) 
	      { 
		if(((float)M_child[i].n)/(2*N) >= PS[j].p)
		  {
		    M_child[i].s = 0;
		    PS.erase(PS.begin()+j); j--;
		  }
	      }
	  }
      }
  }


  void evolve_subpopulation(int i, chromosome& chr)
  {
    int g1,g2,p1,p2,n_mut_1,n_mut_2;


    // create map of shuffled children ids

    int child_map[find(i)->second.N];          
    for(int j = 0; j < find(i)->second.N; j++) { child_map[j] = j; }
    gsl_ran_shuffle(rng,child_map,find(i)->second.N,sizeof(int));


    int c = 0; // counter over all N children (will get mapped to child_map[c])

    // migration, loop over all source populations
    
    map<int,double>::iterator it;

    for(map<int,double>::iterator it = find(i)->second.m.begin(); it != find(i)->second.m.end(); it++)
      {
	int n_migrants = (int)(it->second * find(i)->second.N + 0.5);
        
	for(int m=0; m<n_migrants; m++) 
	  {
	    if(c>=find(i)->second.N) { cerr << "ERROR (evolve subpopulation): too many migrants in subpopulation "<< i << endl; exit(1); } 

	    g1 = 2*child_map[c];   // child genome 1
	    g2 = 2*child_map[c]+1; // child genome 2

	    // draw parents in source population

	    p1 = gsl_rng_uniform_int(rng,find(it->first)->second.G_parent.size()/2);
	    p2 = gsl_rng_uniform_int(rng,find(it->first)->second.G_parent.size()/2);

	    // recombine and mutate

	    recombine(i,g1,it->first,2*p1,2*p1+1,chr.draw_rec_breakpoints(chr.draw_n_rec()));
	    recombine(i,g2,it->first,2*p2,2*p2+1,chr.draw_rec_breakpoints(chr.draw_n_rec()));

	    n_mut_1 = chr.draw_n_mut();
	    n_mut_2 = chr.draw_n_mut();

	    for(int x=0; x<n_mut_1; x++) { add_new_mut(i,g1,chr.draw_new_mut()); }
	    for(int x=0; x<n_mut_2; x++) { add_new_mut(i,g2,chr.draw_new_mut()); }

	    c++;
	  }
      }
	    
    // remainder

    while(c<find(i)->second.N) 
      {
	g1 = 2*child_map[c];   // child genome 1
	g2 = 2*child_map[c]+1; // child genome 2

	p1 = find(i)->second.draw_individual(); // parent 1
	p2 = find(i)->second.draw_individual(); // parent 2

	recombine(i,g1,i,2*p1,2*p1+1,chr.draw_rec_breakpoints(chr.draw_n_rec()));
	recombine(i,g2,i,2*p2,2*p2+1,chr.draw_rec_breakpoints(chr.draw_n_rec()));

	n_mut_1 = chr.draw_n_mut();
	n_mut_2 = chr.draw_n_mut();

	for(int x=0; x<n_mut_1; x++) { add_new_mut(i,g1,chr.draw_new_mut()); }
	for(int x=0; x<n_mut_2; x++) { add_new_mut(i,g2,chr.draw_new_mut()); }

	c++;
      }
  }


  void recombine(int i, int c, int j, int p1, int p2, vector<int> r)
  {
    // child genome c in subpopulation i is assigned outcome of recombination at breakpoints r 
    // between parent genomes p1 and p2 from subpopulation j
    // 
    // example r = (r1,r2)
    // 
    // mutations (      x < r1) assigned from p1
    // mutations (r1 <= x < r2) assigned from p2
    // mutations (r2 <= x     ) assigned from p1
    //
    // p1 and p2 are swapped in half of the cases to assure random assortement

    if(gsl_rng_uniform_int(rng,2)==0) { int swap = p1; p1 = p2; p2 = swap; } // swap p1 and p2

    find(i)->second.G_child[c].clear();

    if(r.size()>0)
      {
	sort(r.begin(),r.end());
	r.erase(unique(r.begin(),r.end()),r.end());

	int p = p1; int p_temp = p2; int swap;
	multimap<int,int>::iterator it = find(j)->second.G_parent[p].begin();

	for(int x=0; x<r.size(); x++) // assign segments from p1 and p2
	  {
	    while(it != find(j)->second.G_parent[p].lower_bound(r[x])) { add_parent_mut(i,c,it->second); it++; }
	    swap = p; p = p_temp; p_temp = swap;
	    it = find(j)->second.G_parent[p].lower_bound(r[x]);
	    if(x+1 == r.size()) { while(it != find(j)->second.G_parent[p].end()) { add_parent_mut(i,c,it->second); it++; } }
	  }
      }
    else // copy p1 into c 
      {
	for(multimap<int,int>::iterator it = find(j)->second.G_parent[p1].begin(); it != find(j)->second.G_parent[p1].end(); it++)
	  {
	    add_parent_mut(i,c,it->second);
	  }
      }
   }


  void swap_generations(int g)
  {
    // find and remove fixed mutations from the children in all subpopulations

    vector<int> S = find_fixed(g); remove(S);

    // make children the new parents

    for(it = begin(); it != end(); it++) 
      { 
	it->second.G_parent.swap(it->second.G_child);
	it->second.G_child.resize(2*it->second.N); 
      }
    
    M_parent.swap(M_child);

    M_child.clear();
    M_map.clear();

    // update parent fitness values in all subpopulations

    for(it = begin(); it != end(); it++) { it->second.update_fitness(M_parent); }
  }


  vector<int> find_fixed(int g)
  {
    // find mutations that are fixed in all child subpopulations and return vector with their ids

    vector<int> S;
    int N = 0;
    for(it = begin(); it != end(); it++) { N += it->second.N; }
    for(int i=0; i<M_child.size(); i++) 
      { 
	if(M_child[i].n==2*N) 
	  { 
	    S.push_back(i); 
	    M_child[i].n=g; // the variable n now stores the fixation generation
	    Substitutions.push_back(M_child[i]); 
	  } 
      }
    return S;
  }


  void remove(vector<int>& S)
  {
    // remove mutations in S from the child genomes in all subpopulations

    for(it = begin(); it != end(); it++) // subpopulations
      {
	for(int j=0; j<2*it->second.N; j++) // child genomes
	  {
	    for(int k=0; k<S.size(); k++) // mutations
	      {
		multimap<int,int>::iterator it2;
		pair<multimap<int,int>::iterator,multimap<int,int>::iterator> range;
		range = it->second.G_child[j].equal_range(M_child[S[k]].x); it2=range.first;

		while(it2!=range.second)
		  {
		    if(it2->second == S[k]) { it->second.G_child[j].erase(it2); it2 = range.second; }
		    else { it2++; }
		  }
	      }
	  }
      }
  }


  void add_parent_mut(int i, int j, int id_parent)
  {
    // add parent mutation id_parent to child genome j in subpopulation i

    int id;

    if(M_map.count(id_parent)>0) { id = M_map[id_parent]; } // mutation not yet mapped
    else
      {
	M_child.push_back(M_parent[id_parent]);    // add mutation to M_child vector
	id =  M_child.size()-1;                    // mut_id is position of mutation in M_child
	M_child[id].n  = 0;                        // reset mutation prevalence
	M_map.insert(pair<int,int>(id_parent,id)); // add parent and child ids to map
      }

    find(i)->second.G_child[j].insert(pair<int,int>(M_child[id].x,id)); // add mutation to genome
    M_child[id].n++; // increase mutation prevalence
  }


  void add_new_mut(int i, int j, mutation m)
  {
    // add new mutation to child genome j in subpopulation i

    M_child.push_back(m);
    int id = M_child.size()-1;
    find(i)->second.G_child[j].insert(pair<int,int>(m.x,id));
    M_child[id].n = 1;
  }


  void print_all()
  {
    // print all mutations and all genomes 

    cout << "Mutations:"  << endl;

    for(int i=0; i<M_child.size(); i++) { cout << i+1 << " "; M_child[i].print(); cout << endl; }

    cout << "Genomes:" << endl;

    for(it = begin(); it != end(); it++) // subpopulations
      {
	for(int i=0; i<2*it->second.N; i++) // child genomes
	  {
	    cout << "p" << it->first << ":" << i+1 << " "; print_genome(it->first,i);
	  }
      }
  }


  void print_all(ofstream& outfile)
  {
    // print all mutations and all genomes to file 

    outfile << "Mutations:"  << endl;

    for(int i=0; i<M_child.size(); i++) { outfile << i+1 << " "; M_child[i].print(outfile); outfile << endl; }

    outfile << "Genomes:" << endl;

    for(it = begin(); it != end(); it++) // subpopulations
      {
	for(int i=0; i<2*it->second.N; i++) // child genomes
	  {
	    outfile << "p" << it->first << ":" << i+1 << " "; print_genome(it->first,i,outfile);
	  }
      }
  }


  void print_sample(int i, int n)
  {
    if(count(i)==0) { cerr << "ERROR (output): subpopulation p"<< i << " does not exists" << endl; exit(1); }

    vector<int> sample; 
    
    map<int,int> m;
    map<int,int>::iterator im;
    
    multimap<int,int>::iterator ig;

    for(int k=0; k<n; k++) 
      { 
	int j = gsl_rng_uniform_int(rng,find(i)->second.G_child.size());
	sample.push_back(j);
	ig = find(i)->second.G_child[j].begin();
	while(ig != find(i)->second.G_child[j].end()) 
	  { 
	    if(m.count(ig->second)==0) { m.insert(pair<int,int>(ig->second,1)); }
	    else                       { m[ig->second]++; }
	    ig++; 
	  }
      }

    cout << "Mutations:"  << endl;

    im = m.begin(); 
    while(im != m.end()) 
      { 
	cout << im->first+1 << " ";
	M_child[im->first].print_non(); 
	cout << " " << im->second << endl; 
	im++; 
      }

    cout << "Genomes:"  << endl;
    
    for(int k=0; k<n; k++) 
      {
	cout << "p" << i << ":" << sample[k]+1 << " "; print_genome(i, sample[k]);
      }
  }


  void print_genome(int i, int j)
  {
    // print child genome j from subpopulation i

    if(count(i)==0) { cerr << "ERROR (output): subpopulation "<< i << " does not exists" << endl; exit(1); }

    multimap<int,int>::iterator it = find(i)->second.G_child[j].begin();
   
    while(it != find(i)->second.G_child[j].end()) 
      {
	cout << it->second+1; it++;
	if(it != find(i)->second.G_child[j].end()) { cout << " "; }
      }

    cout << endl;
   }

  
  void print_genome(int i, int j, ofstream& outfile)
  {
    // print child genome j from subpopulation i to file

    if(count(i)==0) { cerr << "ERROR (output): subpopulation "<< i << " does not exists" << endl; exit(1); }

    multimap<int,int>::iterator it = find(i)->second.G_child[j].begin();
   
    while(it != find(i)->second.G_child[j].end()) 
      {
	outfile << it->second+1; it++;
	if(it != find(i)->second.G_child[j].end()) { outfile << " "; }
      }

    outfile << endl;
   }
};


void get_line(ifstream& infile,string& line)
{
  getline(infile,line);
  if(line.find("/")!= string::npos) { line.erase(line.find("/")); } // remove all after "/"
  line.erase(0,line.find_first_not_of(' ') ); // remove leading whitespaces
  line.erase(line.find_last_not_of(' ') + 1); // remove trailing whitespaces
};


void input_error(int type, string line)
{
  cerr << endl;

  if(type==0) // parameter file
    {
      cerr << "ERROR (parameter file): could not open: " << line << endl << endl;
    }

  else if(type==1) // mutation rate
    {
      cerr << "ERROR (parameter file): invalid mutation rate: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#MUTATION RATE" << endl;
      cerr << "<u>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#MUTATION RATE" << endl;
      cerr << "1.5e-8" << endl << endl;
    }

  else if(type==2) // mutation type
    {
      cerr << "ERROR (parameter file): invalid mutation type: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#MUTATION TYPES" << endl;
      cerr << "<mutation-type-id> <h> <DFE-type> [DFE parameters]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#MUTATION TYPES" << endl;
      cerr << "m1 0.2 g -0.05 0.2" << endl;
      cerr << "m2 0.0 f 0.0" << endl;
      cerr << "m3 0.5 e 0.01" << endl << endl;
    }

  else if(type==3) // genomic element type
    {
      cerr << "ERROR (parameter file): invalid genomic element type: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#GENOMIC ELEMENT TYPES" << endl;
      cerr << "<element-type-id> <mut-type> <x> [<mut-type> <x>...]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#GENOMIC ELEMENT TYPES" << endl;
      cerr << "g1 m3 0.8 m2 0.01 m1 0.19" << endl << endl;
    }

  else if(type==4) // chromosome organization
    {
      cerr << "ERROR (parameter file): invalid chromosome organization: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#CHROMOSOME ORGANIZATION" << endl;
      cerr << "<element-type> <start> <end>" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#CHROMOSOME ORGANIZATION" << endl;
      cerr << "g1 1000 1999" << endl << endl;
    }

  else if(type==5) // recombination rate
    {
      cerr << "ERROR (parameter file): invalid recombination rate: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#RECOMBINATION RATE" << endl;
      cerr << "<interval-end> <r>" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#RECOMBINATION RATE" << endl;
      cerr << "10000 1e-8" << endl;
      cerr << "20000 4.5e-8" << endl << endl;
    }
  
  else if(type==6) // generations
    {
      cerr << "ERROR (parameter file): invalid generations: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#GENERATIONS" << endl;
      cerr << "<t>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#GENERATIONS" << endl;
      cerr << "10000" << endl << endl;
    }

  else if(type==7) // demography and structure
    {
      cerr << "ERROR (parameter file): invalid demography and structure: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#DEMOGRAPHY AND STRUCTURE" << endl;
      cerr << "<time> <event-type> [event parameters]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "DEMOGRAPHY AND STRUCTURE" << endl;
      cerr << "1 P p1 1000" << endl;
      cerr << "1000 P p2 100 p1" << endl;
      cerr << "2000 N p1 1e4" << endl;
      cerr << "2000 M p2 p1 0.01" << endl << endl;
    }

  else if(type==8) // output
    {
      cerr << "ERROR (parameter file): invalid output: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#OUTPUT" << endl;
      cerr << "<time> <output-type> [output parameters]" << endl;
      cerr << "..." << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "OUTPUT" << endl;
      cerr << "2000 A outfile" << endl;
      cerr << "1000 S p1 10" << endl;
      cerr << "2000 F" << endl;
      cerr << "1 T m3" << endl << endl;

    }

   else if(type==9) // initialization
    {
      cerr << "ERROR (parameter file): invalid initialization: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#INITIALIZATION" << endl;
      cerr << "<filename>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#INITIALIZATION" << endl;
      cerr << "outfile" << endl << endl;
    }

   else if(type==10) // seed
    {
      cerr << "ERROR (parameter file): invalid seed: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#SEED" << endl;
      cerr << "<seed>" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#SEED" << endl;
      cerr << "141235" << endl << endl;
    }

   else if(type==11) // predetermined mutation
    {
      cerr << "ERROR (parameter file): invalid predetermined mutations: " << line << endl << endl;
      cerr << "Required syntax:" << endl << endl;
      cerr << "#PREDETERMINED MUTATIONS" << endl;
      cerr << "<time> <mut-type> <x> <s> <h> <pop> <nAA> <nAa> [P <f>]" << endl << endl;
      cerr << "Example:" << endl << endl;
      cerr << "#PREDETERMINED MUTATIONS" << endl;
      cerr << "5000 m7 45000 0.05 0.5 p1 0 1" << endl;
      cerr << "8000 m9 2e5 0.02 1.0 p1 0 10 P 0.5" << endl << endl;
    }

   else if(type==12) // unknown parameter
     {
       cerr << "ERROR (parameter file): unknown parameter: " << line << endl << endl;
     }

   else if(type==13) // no population defined
     {
       cerr << "ERROR (parameter file): no population to simulate:" << endl << endl;
     }

  exit(1);
};


void check_input_file(char* file)
{
  int mutation_types = 0;
  int mutation_rate  = 0;
  int genomic_element_types = 0;
  int chromosome_organization = 0;
  int recombination_rate = 0;
  int generations = 0;
  int population = 0;

  ifstream infile (file);
  if (!infile.is_open()) { input_error(0,string(file)); }

  string line; string sub;

  get_line(infile,line);

  while(!infile.eof())
    {
      if(line.find('#') != string::npos) 
	{
	  if(line.find("MUTATION RATE") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  if(line.find_first_not_of("1234567890.e-") != string::npos ) { input_error(1,line); }
		  else { mutation_rate++; }
		}
	      get_line(infile,line); } }


	  else if(line.find("MUTATION TYPES") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.compare(0,1,"m") != 0) { good = 0; } sub.erase(0,1);

		  if(sub.find_first_not_of("1234567890") != string::npos )   { good = 0; } // id
		  if(iss.eof()) { good = 0; } iss >> sub;
 		  if(sub.find_first_not_of("1234567890.-") != string::npos ) { good = 0; } // h
		  if(iss.eof()) { good = 0; } iss >> sub; 
		  if(sub.find_first_not_of("fge") != string::npos )          { good = 0; } // DFE-type
		  
		  if(sub.compare("f")==0 || sub.compare("e")==0) // one parameter
		    { 
		      if(iss.eof()) { good = 0; } iss >> sub;
		      if(sub.find_first_not_of("1234567890.-") != string::npos ) { good = 0; }
		      if(!iss.eof()) { good = 0; }
		    }
		  if(sub.compare("g")==0) // two parameters
		    {
		      if(iss.eof()) { good = 0; } iss >> sub;
		      if(sub.find_first_not_of("1234567890.-") != string::npos ) { good = 0; }
		      if(iss.eof()) { good = 0; } iss >> sub;
		      if(sub.find_first_not_of("1234567890.-") != string::npos ) { good = 0; }
		      if(!iss.eof()) { good = 0; }
		    }
		   
		  if(good==0) { input_error(2,line); }
		  else { mutation_types++; }
		}
	      get_line(infile,line); } }


	  else if(line.find("GENOMIC ELEMENT TYPES") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.compare(0,1,"g") != 0) { good = 0; } sub.erase(0,1);

		  if(sub.find_first_not_of("1234567890") != string::npos )   { good = 0; } // id
		  if(iss.eof()) { good = 0; }

 		  while(!iss.eof())
		    {
		      iss >> sub;
		      if(sub.compare(0,1,"m") != 0) { good = 0; } sub.erase(0,1); // mutation type id
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; } 
		      if(iss.eof()) { good = 0; } iss >> sub; 
		      if(sub.find_first_not_of("1234567890e.") != string::npos ) { good = 0; } // fraction
		    }

		  if(good==0) { input_error(3,line); }
		  else { genomic_element_types++; }
		}
	      get_line(infile,line); } }


	  else if(line.find("CHROMOSOME ORGANIZATION") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.compare(0,1,"g") != 0) { good = 0; } sub.erase(0,1);

		  if(sub.find_first_not_of("1234567890") != string::npos )   { good = 0; } // id
		  if(iss.eof()) { good = 0; } iss >> sub;
		  if(sub.find_first_not_of("1234567890e") != string::npos )   { good = 0; } // start
		  if(iss.eof()) { good = 0; } iss >> sub;
		  if(sub.find_first_not_of("1234567890e") != string::npos )   { good = 0; } // end
		  if(!iss.eof()) { good = 0; }
 		  

		  if(good==0) { input_error(4,line); }
		  else { chromosome_organization++; }
		}
	      get_line(infile,line); } }

	  
	  else if(line.find("RECOMBINATION RATE") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.find_first_not_of("1234567890e") != string::npos )   { good = 0; } // end
		  if(iss.eof()) { good = 0; } iss >> sub;
		  if(sub.find_first_not_of("1234567890e.-") != string::npos )   { good = 0; } // rate
		  if(!iss.eof()) { good = 0; }
 		  

		  if(good==0) { input_error(5,line); }
		  else { recombination_rate++; }
		}
	      get_line(infile,line); } }


	  else if(line.find("GENERATIONS") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.find_first_not_of("1234567890e") != string::npos )   { good = 0; } // T
		  if(!iss.eof()) { good = 0; } iss >> sub;

		  if(good==0) { input_error(6,line); }
		  else { generations++; }
		}
	      get_line(infile,line); } }


	  else if(line.find("DEMOGRAPHY AND STRUCTURE") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.find_first_not_of("1234567890e") != string::npos )   { good = 0; } // t
		  if(iss.eof()) { good = 0; } iss >> sub;
		  if(sub.find_first_not_of("PMN") != string::npos )  { good = 0; } // event type

		  if(sub.compare("P")==0) // two or three positive integers
		    { 
		      if(iss.eof()) { good = 0; } iss >> sub; // p1
		      if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }

		      if(iss.eof()) { good = 0; } iss >> sub; // N
		      if(sub.find_first_not_of("1234567890e") != string::npos ) { good = 0; }
		      
		      if(!iss.eof()) // p2
			{ 
			  iss >> sub;
			  if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
			  if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
			  if(!iss.eof()) { good = 0; }
			}

		      population++;
		    }

		  if(sub.compare("N")==0) // two positive integers
		    { 
		      if(iss.eof()) { good = 0; } iss >> sub; // p
		      if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		      if(iss.eof()) { good = 0; } iss >> sub; // N
		      if(sub.find_first_not_of("1234567890e") != string::npos ) { good = 0; }		      
		      if(!iss.eof()) { good = 0; }
		    }

		  if(sub.compare("M")==0) // two positive integers and a double
		    { 
		      if(iss.eof()) { good = 0; } iss >> sub; // p
		      if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		      if(iss.eof()) { good = 0; } iss >> sub; // p
		      if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		      if(iss.eof()) { good = 0; } iss >> sub; // M
		      if(sub.find_first_not_of("1234567890.-e") != string::npos ) { good = 0; }
		      if(!iss.eof()) { good = 0; }
		    }
		  
		  if(good==0) { input_error(7,line); }
		}
	      get_line(infile,line); } }


	  else if(line.find("OUTPUT") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  
		  if(sub.find_first_not_of("1234567890e") != string::npos )   { good = 0; } // t
		  if(iss.eof()) { good = 0; } iss >> sub;
		  if(sub.find_first_not_of("ASFT") != string::npos )  { good = 0; } // event type

		  if(sub.compare("A")==0) // no parameter of filename
		    { 
		      if(!iss.eof()) { iss >> sub; if(!iss.eof()) { good = 0; } }
		    }

		  if(sub.compare("S")==0) // two positive integers
		    { 
		      if(iss.eof()) { good = 0; } iss >> sub; // p1
		      if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		      if(iss.eof()) { good = 0; } iss >> sub; // p2
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }		      
		      if(!iss.eof()) { good = 0; }
		    }

		   if(sub.compare("F")==0) // no parameter
		    { 
		      if(!iss.eof()) { good = 0; }
		    }

		   if(sub.compare("T")==0) // follow 
		    { 
		      if(iss.eof()) { good = 0; } iss >> sub; // p1
		      if(sub.compare(0,1,"m") != 0) { good = 0; } sub.erase(0,1);
		      if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		    }

		   if(good==0) { input_error(8,line); }
		}
	      get_line(infile,line); } }


	  else if(line.find("INITIALIZATION") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  if(!iss.eof()) { good = 0; }

		  if(good==0) { input_error(9,line); }

		  population++;
		}
	      get_line(infile,line); } }


	   else if(line.find("SEED") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub;
		  if(sub.find_first_not_of("1234567890-") != string::npos ) { good = 0; }
		  if(!iss.eof()) { good = 0; }

		  if(good==0) { input_error(10,line); }
		}
	      get_line(infile,line); } }


	  else if(line.find("PREDETERMINED MUTATIONS") != string::npos) { get_line(infile,line);
	    while(line.find('#') == string::npos && !infile.eof()) {
	      if(line.length()>0)
		{
		  int good = 1;

		  istringstream iss(line); iss >> sub; // time
		  if(sub.find_first_not_of("1234567890e") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // id
		  if(sub.compare(0,1,"m") != 0) { good = 0; } sub.erase(0,1);
		  if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // x
		  if(sub.find_first_not_of("1234567890e") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // s 
		  if(sub.find_first_not_of("1234567890.-e") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // h
		  if(sub.find_first_not_of("1234567890.-e") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // sub
		  if(sub.compare(0,1,"p") != 0) { good = 0; } sub.erase(0,1);
		  if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // nAA
		  if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		  if(iss.eof()) { good = 0; } iss >> sub; // nAa
		  if(sub.find_first_not_of("1234567890") != string::npos ) { good = 0; }
		  
		  if(!iss.eof())
		    {
		      iss >> sub;
		      if(sub.find_first_not_of("P") != string::npos ) { good = 0; }
		      if(iss.eof()) { good = 0; } iss >> sub; // freq
		      if(sub.find_first_not_of("1234567890.-e") != string::npos ) { good = 0; }
		    }

		  if(!iss.eof()) { good = 0; }

		  if(good==0) { input_error(11,line); }
		}
	      get_line(infile,line); } }


	  else { input_error(12,line); }
	}
      else { get_line(infile,line); }
    }

  if(mutation_rate!=1) { input_error(1,string()); }
  if(mutation_types<1) { input_error(2,string()); }
  if(genomic_element_types<1) { input_error(3,string()); }
  if(chromosome_organization<1) { input_error(4,string()); }
  if(recombination_rate<1) { input_error(5,string()); }
  if(generations<1) { input_error(6,string()); }
  if(population<1) { input_error(13,string()); }
};



void initialize_from_file(population& P, const char* file)
{
  // initialize population from file

  string line; 
  string sub; 

  ifstream infile (file);

  if(!infile.is_open()) { cerr << "ERROR (initialize): could not open initialization file" << endl; exit(1); }

  get_line(infile,line);

  while(line.find("Populations") == string::npos && !infile.eof()) { get_line(infile,line); } 

  get_line(infile,line);

  while(line.find("Mutations") == string::npos && !infile.eof())
    { 
      istringstream iss(line); iss >> sub; sub.erase(0,1);  
      int i = atoi(sub.c_str()); iss >> sub;  
      int n = atoi(sub.c_str());
      P.add_subpopulation(i,n);
      get_line(infile,line);      
    }

  get_line(infile,line);

  while(line.find("Genomes") == string::npos && !infile.eof()) 
    {     
      mutation m(line);
      P.M_parent.push_back(m);
      get_line(infile,line); 
    }

  get_line(infile,line);

  while(!infile.eof())
    {
      istringstream iss(line); iss >> sub; sub.erase(0,1);
      int pos = sub.find_first_of(":"); 
      int p = atoi(sub.substr(0,pos+1).c_str()); sub.erase(0,pos+1);  
      int i = atoi(sub.c_str());

      while(iss >> sub) 
	{ 
	  int m_id = atoi(sub.c_str())-1;
	  int m_x  = P.M_parent[m_id].x;

	  P.find(p)->second.G_parent[i-1].insert(pair<int,int>(m_x,m_id));
	  P.M_parent[m_id].n++; 
	}
      get_line(infile,line);
    }

  for(P.it=P.begin(); P.it!=P.end(); P.it++) { P.it->second.update_fitness(P.M_parent); }
};


void initialize(population& P, char* file, chromosome& chr, int &T, multimap<int,event>& E, multimap<int,event>& O, multimap<int,introduced_mutation>& IM, vector<partial_sweep>& PS, vector<string>& parameters)
{
  string line; 
  string sub; 

  ifstream infile (file);

  long pid=getpid();
  time_t *tp,t; tp=&t; time(tp); t+=pid;
  int seed=t;


  get_line(infile,line);

  while(!infile.eof())
    {
      if(line.find('#') != string::npos) 
	{
	  if(line.find("MUTATION RATE") != string::npos) { get_line(infile,line);
	    parameters.push_back("#MUTATION RATE");
	    while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);
		      istringstream iss(line); iss >> sub;  chr.M = atof(sub.c_str());
		    }
		  get_line(infile,line); } }

	  else if(line.find("MUTATION TYPES") != string::npos) { get_line(infile,line);
	    parameters.push_back("#MUTATION TYPES");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: i h t p1 p2 ... (identifier, dominance coefficient DFE type, DFE parameters) 

		      int i; float h; char t; vector<double> p; istringstream iss(line);
		      iss >> sub; sub.erase(0,1); i = atoi(sub.c_str());

		      if(chr.mutation_types.count(i)>0) 
			{  
			  cerr << "ERROR (initialize): mutation type "<< i << " already defined" << endl; exit(1);
			}

		      iss >> sub; h = atof(sub.c_str());
		      iss >> sub; t = sub.at(0);
		      while(iss >> sub) { p.push_back(atof(sub.c_str())); }
		      chr.mutation_types.insert(pair<int,mutation_type>(i,mutation_type(h,t,p)));
		    } get_line(infile,line); } }

	  else if(line.find("GENOMIC ELEMENT TYPES") != string::npos) { get_line(infile,line);
	    parameters.push_back("#GENOMIC ELEMENT TYPES");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: i m1 g1 [m2 g2 ...] (identifier, mut type class, fraction)
		      
		      int i; vector<int> m; vector<double> g; istringstream iss(line);
		      iss >> sub; sub.erase(0,1); i = atoi(sub.c_str());
		      while(iss >> sub) 
			{ 
			  sub.erase(0,1);
			  m.push_back(atoi(sub.c_str())); iss >> sub;
			  g.push_back(atof(sub.c_str()));
			}

		      if(chr.genomic_element_types.count(i)>0) 
			{  
			  cerr << "ERROR (initialize): genomic element type "<< i << " already defined" << endl; exit(1);
			}

		      chr.genomic_element_types.insert(pair<int,genomic_element_type>(i,genomic_element_type(m,g)));
		    } get_line(infile,line); } }
	  
	  else if(line.find("CHROMOSOME ORGANIZATION") != string::npos) { get_line(infile,line);
	    parameters.push_back("#CHROMOSOME ORGANIZATION");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: i s e (genomic element type identifier, start, end)
		      
		      int i,s,e; istringstream iss(line);
		      iss >> sub; sub.erase(0,1); i = atoi(sub.c_str());
		      iss >> sub; s = (int)atof(sub.c_str())-1;
		      iss >> sub; e = (int)atof(sub.c_str())-1;
		      chr.push_back(genomic_element(i,s,e));
		    } get_line(infile,line); } }

	  else if(line.find("RECOMBINATION RATE") != string::npos) { get_line(infile,line);
	    parameters.push_back("#RECOMBINATION RATE");
	       while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: x r (interval end, rec rate in events per bp)

		      int x; double r; istringstream iss(line);
		      iss >> sub; x = (int)atof(sub.c_str())-1;
		      iss >> sub; r = atof(sub.c_str());
		      chr.rec_x.push_back(x); chr.rec_r.push_back(r);
		    } get_line(infile,line); } }

	  else if(line.find("GENERATIONS") != string::npos) { get_line(infile,line);
	    parameters.push_back("#GENERATIONS");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);
		      istringstream iss(line); iss >> sub;  T = (int)atof(sub.c_str());
		    }
		  get_line(infile,line); } }

	  else if(line.find("DEMOGRAPHY AND STRUCTURE") != string::npos) { get_line(infile,line);
	    parameters.push_back("#DEMOGRAPHY AND STRUCTURE");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: t event_type event_parameters

		      int t; char c; vector<string> s; istringstream iss(line); 
		      iss >> sub; t = (int)atof(sub.c_str());
		      iss >> sub; c = sub.at(0);

		      while(iss >> sub) { s.push_back(sub.c_str()); }
		      E.insert(pair<int,event>(t,event(c,s)));
		    }
		  get_line(infile,line); } }

	  else if(line.find("INITIALIZATION") != string::npos) { get_line(infile,line);
	    parameters.push_back("#INITIALIZATION");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      istringstream iss(line); iss >> sub;  
		      initialize_from_file(P,sub.c_str());
		    }
		  get_line(infile,line); } }

	  else if(line.find("OUTPUT") != string::npos) { get_line(infile,line); 
	    parameters.push_back("#OUTPUT");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: t event_type event_paramaters

		      int t; char c; vector<string> s; istringstream iss(line); 
		      iss >> sub; t = (int)atof(sub.c_str());
		      iss >> sub; c = sub.at(0);

		      while(iss >> sub) { s.push_back(sub.c_str()); }
		      O.insert(pair<int,event>(t,event(c,s)));
		    }
		  get_line(infile,line); } }

	  else if(line.find("PREDETERMINED MUTATIONS") != string::npos) { get_line(infile,line);
	    parameters.push_back("#PREDETERMINED MUTATIONS");
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      parameters.push_back(line);

		      // FORMAT: t type x s h i nAA nAa [P <freq>]

		      istringstream iss(line); 

		      iss >> sub; int time = (int)atof(sub.c_str());
		      iss >> sub; sub.erase(0,1); int t = atoi(sub.c_str());
		      iss >> sub; int x    = (int)atof(sub.c_str())-1;
		      iss >> sub; float s  = atof(sub.c_str());
		      iss >> sub; float h  = atof(sub.c_str());
		      iss >> sub; sub.erase(0,1); int i = atoi(sub.c_str());
		      iss >> sub; int nAA  = (int)atof(sub.c_str());
		      iss >> sub; int nAa  = (int)atof(sub.c_str());

		      introduced_mutation M(t,x,s,h,i,nAA,nAa);

		      IM.insert(pair<int,introduced_mutation>(time,M));
		      
		      while(iss >> sub) 
			{ 
			  if(sub.find('P')!=string::npos) 
			    {
			      iss >> sub; float p = atof(sub.c_str());
			      PS.push_back(partial_sweep(t,x,p));
			    }
			}
		    }
		  get_line(infile,line); } }

	  else if(line.find("SEED") != string::npos) { get_line(infile,line);
	      while(line.find('#') == string::npos && !infile.eof()) {
		  if(line.length()>0)
		    {
		      istringstream iss(line); iss >> sub;  seed = atoi(sub.c_str());
		    }
		  get_line(infile,line); } }

	  
	}
      else { get_line(infile,line); }
    }

  // initialize chromosome

  chr.initialize_rng();

  // initialize rng

  rng=gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng,(long)seed); 

  parameters.push_back("#SEED");
  stringstream ss; ss << seed;
  parameters.push_back(ss.str());

  // parameter output

  for(int i=0; i<P.parameters.size(); i++) { cout << parameters[i] << endl; }
}







int main(int argc,char *argv[])
{
  // initialize simulation parameters

  if(argc<=1) { cerr << "usage: slim <parameter file>" << endl; exit(1); } 

  char* input_file = argv[1];
  check_input_file(input_file);

  int T; chromosome chr;

  population P; map<int,subpopulation>::iterator itP;

  P.parameters.push_back("#INPUT PARAMETER FILE");
  P.parameters.push_back(input_file);

  // demographic and structure events

  multimap<int,event> E; 
  multimap<int,event>::iterator itE;


  // output events (time, output)
  
  multimap<int,event> O; 
  multimap<int,event>::iterator itO;


  // user-defined mutations that will be introduced (time,mutation)

  multimap<int,introduced_mutation> IM; 
  multimap<int,introduced_mutation>::iterator itIM;


  // tracked mutation-types

  vector<int> TM; 


  // mutations undergoing partial sweeps

  vector<partial_sweep> PS;

  initialize(P,input_file,chr,T,E,O,IM,PS,P.parameters); 
 
  // evolve over t generations

  for(int g=1; g<=T; g++)
    { 
      // execute demographic and substructure events in this generation 

      pair<multimap<int,event>::iterator,multimap<int,event>::iterator> rangeE = E.equal_range(g);
      for(itE = rangeE.first; itE != rangeE.second; itE++) { P.execute_event(itE->second,g,TM); }
   
      // evolve all subpopulations

      for(itP = P.begin(); itP != P.end(); itP++) { P.evolve_subpopulation(itP->first,chr); }     
            
      // introduce user-defined mutations

      pair<multimap<int,introduced_mutation>::iterator,multimap<int,introduced_mutation>::iterator> rangeIM = IM.equal_range(g);
      for(itIM = rangeIM.first; itIM != rangeIM.second; itIM++) { P.introduce_mutation(itIM->second); }

      // execute output events

      pair<multimap<int,event>::iterator,multimap<int,event>::iterator> rangeO = O.equal_range(g);
      for(itO = rangeO.first; itO != rangeO.second; itO++) { P.execute_event(itO->second,g,TM); }

      // track particular mutation-types and set s=0 for partial sweeps when completed
      
      if(TM.size()>0 || PS.size()>0) { P.track_mutations(g,TM,PS); }

      // swap generations

      P.swap_generations(g);   
    }
}
