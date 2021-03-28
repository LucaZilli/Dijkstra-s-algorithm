#include <limits>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
//#include <omp.h> //Include it if you want to get the time of the algorithm
#include <cassert>
  
using namespace std;

// Number of nodes in the graph
#define V 6

// A function that read a matrix from file and return it as a std vector
vector < vector<double> > ReadRows( const char* filename , int N,int M);

// A function which needs to find the node with minimum distance value in V-S
int minDistance(double dist[], bool sptSet[]);
  
// Function to print shortest path from source to j using SP array
void printPath(int SP[], int j);

// A recursive function to print the distance array
void printSolution(double dist[], int n,
                   int SP[]);
  
// Function that implements Dijkstra's algorithm
// The graph is represented using the adjacency matrix
void dijkstra(double** graph, int src);
  




// main program 
int main()
{

    vector < vector <double> > graph;
    //Uploading an adjacent matrix
    graph=ReadRows("Graph6",V,V);
    
    int rows=V,cols=V;
    double** G = new double*[rows];
    for (int i = 0; i < rows; ++i)
        G[i] = new double[cols];
    
    for(int j=0;j<V;j++)
        for(int i=0;i<V;i++)
            G[j][i]=graph[j][i];
    
    dijkstra(G, 0);
  
    return 0;
}







vector < vector<double> > ReadRows( const char* filename , int N,int M) {

    vector < vector<double> > v;
    vector<double> x_y;
    ifstream fin (filename);
    assert(fin && "file non esiste");
    if( !fin )
    {
        cout << "file non esistente" << endl;
        exit(0);
    }
    else
    {
        for( int k=0 ; k<N ; k++) {
            x_y.clear();
            for( int j=0 ; j<M ; j++){
                double valore=0;
                fin >> valore;
                x_y.push_back(valore);
               }
                v.push_back(x_y);
        
                //assert(!(fin.eof()) && "lettura file terminata");
        
                if( fin.eof() ) {
                    //cout << "ho finito i dati" << endl;
                    //exit(1);
                }
        }
           
    }
    return v;
}





void printPath(int SP[], int j)
{
      
    // Base Case : If j is source
    if (SP[j] == - 1)
        return;
  
    printPath(SP, SP[j]);
  
    printf("%d ", j);
    
}



void printSolution(double dist[], int n,
                   int SP[])
{
    int src = 0;
    printf("Vertex\t\t Distance\t\tPath");
    for (int i = 1; i < V; i++)
    {
        printf("\n%d -> %d \t\t %f\t\t%d ",
                      src, i, dist[i], src);
        printPath(SP, i);
    }
    cout<<endl;
}
    
int minDistance(double dist[], bool S[])
{
    // Initialize min value
    double min = std::numeric_limits<double>::max();
    int min_index;
  
    for (int v = 0; v < V; v++)
        if (S[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
  
    return min_index;
}

void dijkstra(double** graph, int src)
{
    //double start = omp_get_wtime();
    double dist[V]; // dist[i] will contain the shortest distance from src to i
    bool S[V]; // sptSet[i] will be true if vertex i is included in S, else will be false
  
  
    // ShortestPath array to store the shortest path
    int ShortestPath[V];
    
    //Initialize SP[0]=-1
    ShortestPath[0] = -1;
    //Include src in S
    S[src] = true;
    // Distance of source vertex from itself is always 0
    dist[src] = 0.;
    
    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 1; i < V; i++){
        S[i] = false;
        dist[i] = std::numeric_limits<double>::max();
    }
    

    for (int v = 0; v < V; v++)

        // Update dist[v]
        if (!S[v] && (graph[src][v]>0.) && dist[src] != std::numeric_limits<double>::max()
            && dist[src] + graph[src][v] < dist[v]){
        
            ShortestPath[v] = src;
            dist[v] = dist[src] + graph[src][v];
        }
  
    // Find shortest path for all vertices
    for (int count = 1; count < V - 1; count++) {
        // Pick the Argimn(D[v]) where v is contained in V-S (S[v]=False)
        int u = minDistance(dist, S);
  
        // Add the node in S
        S[u] = true;
  
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int v = 0; v < V; v++)
  
            // Update dist[v] only if is not in S and total weight of path from src to v through u is smaller than current value of dist[v]
            if (!S[v] && (graph[u][v]>0.) && dist[u] != std::numeric_limits<double>::max()
                && dist[u] + graph[u][v] < dist[v]){
                
                ShortestPath[v] = u;
                dist[v] = dist[u] + graph[u][v];
            }
    }
  
    //double end = omp_get_wtime();
    cout<<endl<<"time= "<<end-start<<endl;
    
    printSolution(dist, V, ShortestPath);
}



