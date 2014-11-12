#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>

using namespace std;

typedef double EdgeWeight;
typedef int Node;

struct Vertex {
	Vertex(): community(-1), weighted_degree(0){}
	int community;
	EdgeWeight weighted_degree;
};

typedef vector<Vertex> VertexList;
typedef vector<pair<Node, EdgeWeight> > EdgeList;
typedef vector<EdgeList> GRA;

class Graph{
	private:
		GRA graph;
		VertexList vl;
		int vsize;
		int esize;
		EdgeWeight total_weight;
	public:
		Graph();
		Graph(int num, EdgeWeight total_weight);
		Graph(istream& in);
		~Graph();
		void readGraph(istream&);
		void writeGraph();
		void addVertex(Vertex);
		void addEdge(Node, Node, const EdgeWeight&);
		const int& getNumVertices() const;
		const int& getNumEdges() const;
		const EdgeWeight& getTotalWeight() const;
		VertexList& getVertices();
		int getDegree(Node src);
		EdgeList& getNeighbors(Node src);
		void setCommunity(Node src, int community);
		int getCommunity(Node src);
		EdgeWeight getWeightedDegree(Node src);
		void setWeightedDegree(Node src, EdgeWeight weight);
		void resize(int size);
};

inline
void Graph::addVertex(Vertex v){
	vl.push_back(v);
	vsize++;
}

inline
void Graph::addEdge(Node sid, Node tid, const EdgeWeight& weight) {
	graph[sid].push_back(make_pair(tid, weight));
	esize++;
}

inline
const int& Graph::getNumVertices() const{
	return vsize;
}

inline
const int& Graph::getNumEdges() const{
	return esize;
}

inline
const EdgeWeight& Graph::getTotalWeight() const{
	return total_weight;
}

inline
VertexList& Graph::getVertices(){
	return vl;
}

inline
int Graph::getDegree(Node src) {
	return graph[src].size();
}

inline
EdgeList& Graph::getNeighbors(Node src){
	return graph[src];
}

inline
void Graph::setCommunity(Node src, int comm){
	vl[src].community = comm;
}

inline
int Graph::getCommunity(Node src){
	return vl[src].community;
}

inline
EdgeWeight Graph::getWeightedDegree(Node src){
	return vl[src].weighted_degree;
}

inline
void Graph::setWeightedDegree(Node src, EdgeWeight weight){
	vl[src].weighted_degree = weight;
}

#endif
