#include "Graph.h"

Graph::Graph(){
	vl = VertexList(0);
	vsize = 0;
	esize = 0;
	EdgeList el;
	graph = GRA(vsize, el);
	total_weight = 0;
}

Graph::Graph(istream& in) {
	readGraph(in);
}

Graph::~Graph() {}

void Graph::readGraph(istream& in){
	string buf;
	vsize = 0;
	int id;
	while(getline(in, buf)){
		if(buf.find("#") != string::npos)
			continue;
		istringstream(buf.substr(0, buf.find("	"))) >> id;
		if(vsize < id)
			vsize = id;
		buf.erase(0, buf.find("	")+1);
		
		istringstream(buf.substr(0, buf.find("	"))) >> id;
		if(vsize < id)
			vsize = id;
	}
	id = 0;
	EdgeList el;
	vl = VertexList(++vsize);	
	graph = GRA(vsize, el);
	
	//build edges
	in.clear();
	in.seekg(0, ios::beg);
	Node sid, tid;
	EdgeWeight weight = 1;
	while(getline(in, buf)){
		if(buf.find("#") != string::npos)
			continue;
		istringstream(buf.substr(0, buf.find("	"))) >> sid;
		buf.erase(0, buf.find("	")+1);
		
		istringstream(buf.substr(0, buf.find("	"))) >> tid;
		buf.erase(0, buf.find("	")+1);

		if(buf.find("	") != string::npos)
			istringstream(buf.substr(0, buf.find("	"))) >> weight;
		addEdge(sid, tid, weight);
		if(tid != sid){
			esize--;
			addEdge(tid, sid, weight);
		}
		vl[sid].weighted_degree += weight;
		vl[tid].weighted_degree += weight;
		total_weight += weight;
	}
	while(id < vsize){
		setCommunity(id, id);
		id++;
	}
}

void Graph::writeGraph(){
	cout<<"number of vertices: "<<getNumVertices()<<"  number of edges: "<<getNumEdges()<<endl;
	for(Node i = 0; i < getNumVertices(); i++){
		for(EdgeList::iterator eit = getNeighbors(i).begin(); eit != getNeighbors(i).end(); eit++){
			cout<<i<<" "<<(*eit).first<<" "<<(*eit).second<<endl;
		}
	}
}

void Graph::resize(int num){
	graph.clear();
	vl.resize(num);
	EdgeList el;
	graph = GRA(vsize, el);
	vsize = num;
	esize = 0;
}
