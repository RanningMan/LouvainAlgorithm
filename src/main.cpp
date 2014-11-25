#include "Graph.h"
#include "hash_map.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sys/time.h>

void writeCommunity(vector<int>& ntoc){
	int i = 0;
	for(vector<int>::iterator it = ntoc.begin(); it != ntoc.end(); it++,i++){
		cout<<i<<"  "<<*it<<endl;
	}
}

void writeCommunity(vector<vector<int> >& community){
	int i = 0;
	for(vector<vector<int> >::iterator it = community.begin(); it != community.end(); it++, i++){
	int max_community_size=-1;
	int max_community_id=-1;
	for(vector<vector<int> >::iterator it = community.begin(); it != community.end(); it++, i++){
		if(int(it->size())>max_community_size){
			max_community_size=int(it->size());
			max_community_id=i;
		}
		cout<<i<<": ";
		for(vector<int>::iterator iit = (*it).begin(); iit != (*it).end(); iit++){
			cout<<*iit<<" ";
		}
		cout<<endl;
	}
	cout<<"Community "<<max_community_id<<" has the biggest size: "<<max_community_size<<endl;
}

EdgeWeight getEdgeWeightBetweenVertices(Graph& g, Node i, Node j){ //O(m)
	for(EdgeList::iterator eit = g.getNeighbors(i).begin(); eit != g.getNeighbors(i).end(); eit++){
		if((*eit).first == j)
			return (*eit).second;
	}
	return 0;
}

EdgeWeight getEdgeWeightBetweenCommunity(Graph& g, Node i, int comm){ //O(m)
	EdgeWeight ew = 0;
	for(EdgeList::iterator eit = g.getNeighbors(i).begin(); eit != g.getNeighbors(i).end(); eit++){
		if(i == (*eit).first)continue;
		if(g.getCommunity((*eit).first) == comm)
			ew += (*eit).second;
	}
	return ew;
}

void remove(Graph& g, Node i, vector<EdgeWeight>& in, vector<EdgeWeight>& tot){ //O(m)
	int comm = g.getCommunity(i);
	in[comm] -= (getEdgeWeightBetweenVertices(g,i,i) + 2 * getEdgeWeightBetweenCommunity(g,i,comm));
	tot[comm] -= g.getWeightedDegree(i);
	g.setCommunity(i, -1);
}

void add(Graph& g, Node i, int comm, vector<EdgeWeight>& in, vector<EdgeWeight>& tot){ //O(m)
	in[comm] += (getEdgeWeightBetweenVertices(g,i,i) + 2 * getEdgeWeightBetweenCommunity(g,i,comm));
	tot[comm] += g.getWeightedDegree(i);
	g.setCommunity(i, comm);
}

double modularity(Graph& g, vector<EdgeWeight>& in, vector<EdgeWeight>& tot){ //O(n)
	double sum = 0;
	double m = g.getTotalWeight();
	for (int i = 0; i < g.getNumVertices(); i++){
		if(tot[i] > 0){
			sum += in[i] / (2 * m) - tot[i] * tot[i] / (4 * m * m);
		}
	}
	return sum;
}

double modularityGain(Graph& g, Node i, int comm, EdgeWeight comm_weight, vector<EdgeWeight>& in, vector<EdgeWeight>& tot){ //O(1)
	double m = g.getTotalWeight();
	double n2c = comm_weight;
	double wd = g.getWeightedDegree(i);
	return (n2c - tot[comm] * wd / (2 * m)) / m;
}

//the first phase, stop when some condition reached
double detectCommunity(Graph& g, vector<EdgeWeight>& in, vector<EdgeWeight>& tot, float EPSILON, int PASS_NUM){
	int max, pass = 0;
	int pass = 0;
	double gain = 0;
	double maxgain = 0;
	double mod = 1;
	double pre_mod = 0;
	int pre_comm;
	bool changed = false;
	vector<int> visited;
	vector<EdgeWeight> comm_weight(g.getNumVertices(), 0);
	for (int i = 0; i < g.getNumVertices(); i++){
		in[i]  = getEdgeWeightBetweenVertices(g,i, i);
		tot[i] = g.getWeightedDegree(i);
	}
	do {
		changed=false;
		pre_mod = modularity(g, in, tot);
		pass++;
		for(Node i = 0; i < g.getNumVertices(); i++){
			pre_comm = g.getCommunity(i);
			for(EdgeList::iterator eit = g.getNeighbors(i).begin(); eit != g.getNeighbors(i).end(); eit++){
				comm_weight[g.getCommunity((*eit).first)] += (*eit).second;
				visited.push_back(g.getCommunity((*eit).first));
			}
			remove(g,i,in,tot);
			max = i;
			for(EdgeList::iterator eit = g.getNeighbors(i).begin(); eit != g.getNeighbors(i).end(); eit++){
				if(i == (*eit).first)continue;
				gain = modularityGain(g, i, g.getCommunity((*eit).first), comm_weight[g.getCommunity((*eit).first)], in, tot);
				if(gain > maxgain){
					maxgain = gain;
					max = (*eit).first;
				}
			}
			if(max == i){
				add(g,i,pre_comm,in,tot);
			}
			else{
				add(g,i,g.getCommunity(max),in,tot);
			remove(g,i,in,tot);
			for(EdgeList::iterator eit = g.getNeighbors(i).begin(); eit != g.getNeighbors(i).end(); eit++){
				if(i == (*eit).first)continue;
				comm_weight[g.getCommunity((*eit).first)] += (*eit).second;
				visited.push_back(g.getCommunity((*eit).first));
			}
			int max = pre_comm;
			for(vector<int>::iterator cmit=visited.begin(); cmit!=visited.end(); cmit++){
				gain = modularityGain(g, i, *cmit, comm_weight[*cmit], in, tot);
				if(gain > maxgain){
					maxgain = gain;
					max = *cmit;
				}
			}
			
			add(g,i,max,in,tot);
			
			if(max != pre_comm){
				changed = true;
			}
			gain = 0;
			maxgain = 0;
			for(vector<int>::iterator it = visited.begin(); it != visited.end(); it++){
				comm_weight[*it] = 0;
			}
			visited.clear();
		}
		mod = modularity(g, in, tot);
		cout<<"mod is "<<mod<<", pre_mod is "<<pre_mod<<", pass is "<<pass<<", mod implement is "<<(mod - pre_mod)<<endl;
	} while(changed &&((mod - pre_mod) > EPSILON) && pass < PASS_NUM); 
	return mod;
}

//the second phase
int constructNewGraph(Graph& g, vector<EdgeWeight>& in, vector<EdgeWeight>& tot, vector<int>& ntoc){
	int num_comm = 0;	
	int cur_comm, v_comm;
	vector<int> ctoc(g.getNumVertices(), 0);	
	int num_vertices = g.getNumVertices();	
	for(Node i = 0; i < g.getNumVertices(); i++){
		if(tot[i] > 0)
			ctoc[i] = num_comm++;
	}
	for(vector<int>::iterator it = ntoc.begin(); it != ntoc.end(); it++){
		*it = ctoc[*it];
	}	

	EdgeList commEdge[num_comm];	
	vector<hash_map<int, EdgeWeight> > newGraphEdgeSet (num_comm, hash_map<int, EdgeWeight>());
	for(Node i = 0; i < g.getNumVertices(); i++){	
		v_comm = g.getCommunity(i);
		for(EdgeList::iterator eit = g.getNeighbors(i).begin(); eit != g.getNeighbors(i).end(); eit++){
			cur_comm = g.getCommunity((*eit).first);
			if(cur_comm != v_comm){
				commEdge[ctoc[v_comm]].push_back(make_pair(ctoc[cur_comm], (*eit).second));
				if(newGraphEdgeSet[ctoc[cur_comm]].find(ctoc[v_comm])==newGraphEdgeSet[ctoc[cur_comm]].end()){
					newGraphEdgeSet[ctoc[cur_comm]][ctoc[v_comm]]=(*eit).second;
				}else{
					newGraphEdgeSet[ctoc[cur_comm]][ctoc[v_comm]]+=(*eit).second;
				}
			}
		}		
	}	
	
	g.resize(num_comm);
	for(Node i = 0; i < num_vertices; i++){
		if(tot[i] > 0){
			g.setWeightedDegree(ctoc[i], tot[i]);
			g.setCommunity(ctoc[i], ctoc[i]);
			g.addEdge(ctoc[i], ctoc[i], in[i]);
		}
	}	
	for(int i = 0; i < num_comm; i++){
		for(EdgeList::iterator it = commEdge[i].begin(); it != commEdge[i].end(); it++){
			g.addEdge(i, (*it).first, (*it).second);
		}
	}	
			
			//Here is an attempt to add no self edge
			//g.setWeightedDegree(ctoc[i], tot[i]-in[i]);
			//g.setCommunity(ctoc[i], ctoc[i]);
		}
	}	
	
	for(vector<hash_map<int, EdgeWeight> >::iterator it=newGraphEdgeSet.begin(); it!=newGraphEdgeSet.end(); it++){
		for(hash_map<int, EdgeWeight>::iterator jt=it->begin(); jt!=it->end(); jt++){
			g.addEdge(distance(newGraphEdgeSet.begin(), it), jt->first, jt->second);
		}
	}
	g.writeGraph();
	return num_comm;
}

void LouvainAlgorithm(Graph& g, vector<vector<int> >& community, float EPSILON, int PASS_NUM){
	vector<EdgeWeight> in, tot;
	double mod = 0;
	double pre_mod = 0;
	int pass = 0;
	int com = 0;
	vector<int> ntoc(g.getNumVertices(), 0);
	int i = 0;
	for(vector<int>::iterator it = ntoc.begin(); it != ntoc.end(); it++,i++){
		*it = i;
	}
	while(true){
		in.resize(g.getNumVertices());
		tot.resize(g.getNumVertices());
		pass++;
		cout<<"pass"<<pass<<":"<<endl;
		pre_mod = mod;
		cout<<"first phase"<<endl;
		mod = detectCommunity(g, in, tot, EPSILON, PASS_NUM);
		if((mod - pre_mod) < EPSILON || pass > PASS_NUM)
			break;			
		for(vector<int>::iterator it = ntoc.begin(); it != ntoc.end(); it++){
			*it = g.getCommunity(*it);
		}
		cout<<"second phase"<<endl;
		com = constructNewGraph(g, in, tot, ntoc);
		cout<<"num of community is "<<com<<endl;
		in.clear();
		tot.clear();
	}
	community.resize(com);
	for(int i = 0; i < (int)ntoc.size(); i++){
		community[ntoc[i]].push_back(i);
	}
	writeCommunity(community);	
	cout<<"mod is "<<mod<<endl;
	cout<<"total_pass"<<pass<<endl;
}

void usage(){
}

int main(int argc, char **argv){
	string input_graphfile;
	int input_para_counter = 1;
	float EPSILON = 0.000001;
	int	PASS_NUM = 100;
	
	while (input_para_counter < argc) {
		if (strcmp("-h", argv[input_para_counter]) == 0) {
			usage();
			return 1;
		}
		
		if(strcmp("-e", argv[input_para_counter]) == 0) {
			input_para_counter++;
			EPSILON = atof(argv[input_para_counter++]);
		}
		else if(strcmp("-p", argv[input_para_counter]) == 0){
			input_para_counter++;
			PASS_NUM = atoi(argv[input_para_counter++]);
		}
		else{
			input_graphfile = argv[input_para_counter++];
		}
	}		
	
	ifstream infile_graph(input_graphfile.c_str());
	if (!infile_graph) {
		cout << "Error: Cannot open " << input_graphfile<<endl;
		return -1;
	}
	struct timeval after_time, before_time;
	float graph_time = 0;
	float community_time = 0;
	gettimeofday(&before_time, NULL); 
	Graph g(infile_graph);
	gettimeofday(&after_time, NULL);
	graph_time = (after_time.tv_sec - before_time.tv_sec)*1000.0+(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
	vector<vector<int> > community;
	gettimeofday(&before_time, NULL); 
	LouvainAlgorithm(g, community, EPSILON, PASS_NUM);
	gettimeofday(&after_time, NULL);
	community_time = (after_time.tv_sec - before_time.tv_sec)*1000.0+(after_time.tv_usec - before_time.tv_usec)*1.0/1000.0;
	cout<<"graph construction time: "<<graph_time<<endl;
	cout<<"community finding time: "<<community_time<<endl;
}
